#include "contact.h"

// Definition of static member variables
size_t contact_patch_t::num_nodes_in_patch;
double contact_patches_t::bs;
size_t contact_patches_t::n;

void contact_patch_t::update_nodes(const mesh_t &mesh, const node_t &nodes, // NOLINT(*-make-member-function-const)
                                   const corner_t &corner)
{
    // Constructing the points, vel_points, and internal_force arrays
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < num_nodes_in_patch; j++)
        {
            const size_t node_gid = this->nodes_gid(j);

            points(i, j) = nodes.coords(0, node_gid, i);
            vel_points(i, j) = nodes.vel(0, node_gid, i);

            // looping over the corners
            for (size_t corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++)
            {
                size_t corner_gid = mesh.corners_in_node(node_gid, corner_lid);
                internal_force(i, j) += corner.force(corner_gid, i);
            }

        }  // end local node loop
    }  // end dimension loop
}

void contact_patches_t::initialize(const mesh_t &mesh, const CArrayKokkos<size_t> &bdy_contact_patches,
                                   const node_t &nodes)
{
    // Contact is only supported in 3D
    if (mesh.num_dims != 3)
    {
        std::cerr << "Error: contact is only supported in 3D" << std::endl;
        exit(1);
    }

    // Set up the patches that will be checked for contact
    patches_gid = bdy_contact_patches;
    num_contact_patches = patches_gid.size();

    // Set up array of contact_patch_t
    contact_patches = CArrayKokkos<contact_patch_t>(num_contact_patches);

    CArrayKokkos<size_t> nodes_in_patch(num_contact_patches, mesh.num_nodes_in_patch);
    FOR_ALL_CLASS(i, 0, num_contact_patches,
                  j, 0, mesh.num_nodes_in_patch, {
                      nodes_in_patch(i, j) = mesh.nodes_in_patch(patches_gid(i), j);
                  });  // end parallel for
    Kokkos::fence();

    for (int i = 0; i < num_contact_patches; i++)
    {
        contact_patch_t &contact_patch = contact_patches(i);
        contact_patch.gid = patches_gid(i);

        // Make contact_patch.nodes_gid equal to the row of nodes_in_patch(i)
        // This line is what is limiting the parallelism
        contact_patch.nodes_gid = CArrayKokkos<size_t>(mesh.num_nodes_in_patch);
        for (size_t j = 0; j < mesh.num_nodes_in_patch; j++)
        {
            contact_patch.nodes_gid(j) = nodes_in_patch(i, j);
        }
    }  // end for

    // Setting up the iso-parametric coordinates for all patch objects
    if (mesh.num_nodes_in_patch == 4)
    {
        contact_patch_t::num_nodes_in_patch = 4;
        double xi_temp[4] = {-1.0, 1.0, 1.0, -1.0};
        double eta_temp[4] = {-1.0, -1.0, 1.0, 1.0};

        // Allocating memory for the points and vel_points arrays
        for (int i = 0; i < num_contact_patches; i++)
        {
            contact_patch_t &contact_patch = contact_patches(i);
            contact_patch.points = CArrayKokkos<double>(3, 4);
            contact_patch.vel_points = CArrayKokkos<double>(3, 4);
            contact_patch.internal_force = CArrayKokkos<double>(3, 4);

            contact_patch.xi = CArrayKokkos<double>(4);
            contact_patch.eta = CArrayKokkos<double>(4);

            for (int j = 0; j < 4; j++)
            {
                contact_patch.xi(j) = xi_temp[j];
                contact_patch.eta(j) = eta_temp[j];
            }
        }

    } else
    {
        std::cerr << "Error: higher order elements are not yet tested for contact" << std::endl;
        exit(1);
    }  // end if

    // Determine the bucket size. This is defined as 1.001*min_node_distance
    CArrayKokkos<double> node_distances(num_contact_patches, contact_patch_t::num_nodes_in_patch);
    FOR_ALL_CLASS(i, 0, num_contact_patches,
                  j, 0, contact_patch_t::num_nodes_in_patch, {
                      const contact_patch_t &contact_patch = contact_patches(i);
                      if (j < contact_patch_t::num_nodes_in_patch - 1)
                      {
                          const size_t n1 = contact_patch.nodes_gid(j); // current node
                          const size_t n2 = contact_patch.nodes_gid(j + 1); // next node

                          double sum_sq = 0.0;

                          for (int k = 0; k < 3; k++)
                          {
                              sum_sq += pow(nodes.coords(0, n1, k) - nodes.coords(0, n2, k), 2);
                          }

                          node_distances(i, j) = sqrt(sum_sq);
                      } else
                      {
                          const size_t n1 = contact_patch.nodes_gid(0); // current node
                          const size_t n2 = contact_patch.nodes_gid(
                              contact_patch_t::num_nodes_in_patch - 1); // next node

                          double sum_sq = 0.0;

                          for (int k = 0; k < 3; k++)
                          {
                              sum_sq += pow(nodes.coords(0, n1, k) - nodes.coords(0, n2, k), 2);
                          }

                          node_distances(i, j) = sqrt(sum_sq);
                      }
                  });
    Kokkos::fence();

    double result = 0.0;
    double local_min = 1.0e10;
    REDUCE_MIN(i, 0, num_contact_patches,
               j, 0, contact_patch_t::num_nodes_in_patch, local_min, {
                   if (local_min > node_distances(i, j))
                   {
                       local_min = node_distances(i, j);
                   }
               }, result);

    contact_patches_t::bs = 1.001*result;

    // Find the total number of nodes (this is should always be less than or equal to mesh.num_bdy_nodes)
    size_t local_max_index = 0;
    size_t max_index = 0;
    REDUCE_MAX(i, 0, num_contact_patches,
               j, 0, contact_patch_t::num_nodes_in_patch, local_max_index, {
                   if (local_max_index < contact_patches(i).nodes_gid(j))
                   {
                       local_max_index = contact_patches(i).nodes_gid(j);
                   }
               }, max_index);
    Kokkos::fence();

    CArrayKokkos<size_t> node_count(max_index + 1);
    for (int i = 0; i < num_contact_patches; i++)
    {
        contact_patch_t &contact_patch = contact_patches(i);
        FOR_ALL(j, 0, contact_patch_t::num_nodes_in_patch, {
            size_t node_gid = contact_patch.nodes_gid(j);
            if (node_count(node_gid) == 0)
            {
                node_count(node_gid) = 1;
                contact_patches_t::n += 1;
            }
        });
    }

    // Construct nodes_gid
    nodes_gid = CArrayKokkos<size_t>(contact_patches_t::n);
    size_t node_lid = 0;
    for (int i = 0; i < num_contact_patches; i++)
    {
        contact_patch_t &contact_patch = contact_patches(i);
        for (int j = 0; j < contact_patch_t::num_nodes_in_patch; j++)
        {
            size_t node_gid = contact_patch.nodes_gid(j);
            if (node_count(node_gid) == 1)
            {
                node_count(node_gid) = 2;
                nodes_gid(node_lid) = node_gid;
                node_lid += 1;
            }
        }
    }
}


void contact_patches_t::sort(const mesh_t &mesh, const node_t &nodes, const corner_t &corner)
{
    // Update the points and vel_points arrays for each contact patch
    FOR_ALL_CLASS(i, 0, num_contact_patches, {
        contact_patches(i).update_nodes(mesh, nodes, corner);
    });  // end parallel for
    Kokkos::fence();

    // Grouping all the coordinates, velocities, and accelerations
    CArrayKokkos<double> points(3, contact_patch_t::num_nodes_in_patch*num_contact_patches);
    CArrayKokkos<double> velocities(3, contact_patch_t::num_nodes_in_patch*num_contact_patches);
    CArrayKokkos<double> accelerations(3, contact_patch_t::num_nodes_in_patch*num_contact_patches);

    FOR_ALL(i, 0, num_contact_patches, {
        const contact_patch_t &contact_patch = contact_patches(i);
        for (int j = 0; j < contact_patch_t::num_nodes_in_patch; j++)
        {
            const size_t node_gid = contact_patch.nodes_gid(j);
            const double &mass = nodes.mass(node_gid);
            for (int k = 0; k < 3; k++)
            {
                points(k, i*contact_patch_t::num_nodes_in_patch + j) = contact_patch.points(k, j);
                velocities(k, i*contact_patch_t::num_nodes_in_patch + j) = fabs(contact_patch.vel_points(k, j));
                accelerations(k, i*contact_patch_t::num_nodes_in_patch + j) = fabs(
                    contact_patch.internal_force(k, j)/mass);
            }
        }
    });  // end parallel for
    Kokkos::fence();

    // Find the max and min for each dimension
    double local_x_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_x_max, {
        if (local_x_max < points(0, i))
        {
            local_x_max = points(0, i);
        }
    }, x_max);
    double local_y_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_y_max, {
        if (local_y_max < points(1, i))
        {
            local_y_max = points(1, i);
        }
    }, y_max);
    double local_z_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_z_max, {
        if (local_z_max < points(2, i))
        {
            local_z_max = points(2, i);
        }
    }, z_max);
    double local_x_min = 0.0;
    REDUCE_MIN_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_x_min, {
        if (local_x_min > points(0, i))
        {
            local_x_min = points(0, i);
        }
    }, x_min);
    double local_y_min = 0.0;
    REDUCE_MIN_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_y_min, {
        if (local_y_min > points(1, i))
        {
            local_y_min = points(1, i);
        }
    }, y_min);
    double local_z_min = 0.0;
    REDUCE_MIN_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_z_min, {
        if (local_z_min > points(2, i))
        {
            local_z_min = points(2, i);
        }
    }, z_min);
    double local_vx_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_vx_max, {
        if (local_vx_max < velocities(0, i))
        {
            local_vx_max = velocities(0, i);
        }
    }, vx_max);
    double local_vy_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_vy_max, {
        if (local_vy_max < velocities(1, i))
        {
            local_vy_max = velocities(1, i);
        }
    }, vy_max);
    double local_vz_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_vz_max, {
        if (local_vz_max < velocities(2, i))
        {
            local_vz_max = velocities(2, i);
        }
    }, vz_max);
    double local_ax_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_ax_max, {
        if (local_ax_max < accelerations(0, i))
        {
            local_ax_max = accelerations(0, i);
        }
    }, ax_max);
    double local_ay_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_ay_max, {
        if (local_ay_max < accelerations(1, i))
        {
            local_ay_max = accelerations(1, i);
        }
    }, ay_max);
    double local_az_max = 0.0;
    REDUCE_MAX_CLASS(i, 0, contact_patch_t::num_nodes_in_patch*num_contact_patches, local_az_max, {
        if (local_az_max < accelerations(2, i))
        {
            local_az_max = accelerations(2, i);
        }
    }, az_max);
    Kokkos::fence();

    // If the max velocity is zero, then we want to set it to a small value. The max velocity and acceleration are used
    // for creating a capture box around the contact patch. We want there to be at least some thickness to the box.
    double *vel_max[3] = {&vx_max, &vy_max, &vz_max};
    for (auto &i: vel_max)
    {
        if (*i == 0.0)
        {
            *i = 1.0e-3; // Set to a small value
        }
    }

    // Define Sx, Sy, and Sz
    Sx = floor((x_max - x_min)/bs) + 1; // NOLINT(*-narrowing-conversions)
    Sy = floor((y_max - y_min)/bs) + 1; // NOLINT(*-narrowing-conversions)
    Sz = floor((z_max - z_min)/bs) + 1; // NOLINT(*-narrowing-conversions)

    // Initializing the nbox, lbox, nsort, and npoint arrays
    size_t nb = Sx*Sy*Sz;  // total number of buckets
    nbox = CArrayKokkos<size_t>(nb);
    lbox = CArrayKokkos<size_t>(contact_patches_t::n);
    nsort = CArrayKokkos<size_t>(contact_patches_t::n);
    npoint = CArrayKokkos<size_t>(nb);
    CArrayKokkos<size_t> nsort_lid(contact_patches_t::n);

    // Find the bucket id for each node by constructing lbox
    FOR_ALL_CLASS(i, 0, contact_patches_t::n, {
        size_t node_gid = nodes_gid(i);
        double x = nodes.coords(0, node_gid, 0);
        double y = nodes.coords(0, node_gid, 1);
        double z = nodes.coords(0, node_gid, 2);

        size_t Si_x = floor((x - x_min)/bs);
        size_t Si_y = floor((y - y_min)/bs);
        size_t Si_z = floor((z - z_min)/bs);

        lbox(i) = Si_z*Sx*Sy + Si_y*Sx + Si_x;
        nbox(lbox(i)) += 1;  // increment nbox
    });
    Kokkos::fence();

    // Calculate the pointer for each bucket into a sorted list of nodes
    for (size_t i = 1; i < nb; i++) {
        npoint(i) = npoint(i - 1) + nbox(i - 1);
    }

    // Zero nbox
    FOR_ALL_CLASS(i, 0, nb, {
        nbox(i) = 0;
    });
    Kokkos::fence();

    // Sort the slave nodes according to their bucket id into nsort
    for (int i = 0; i < contact_patches_t::n; i ++) {
        nsort_lid(nbox(lbox(i)) + npoint(lbox(i))) = i;
        nbox(lbox(i)) += 1;
    }

    // Change nsort to reflect the global node id's
    FOR_ALL_CLASS(i, 0, contact_patches_t::n, {
        nsort(i) = nodes_gid(nsort_lid(i));
    });
    Kokkos::fence();
}
