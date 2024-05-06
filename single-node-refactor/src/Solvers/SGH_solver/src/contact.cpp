#include "contact.h"

// Definition of static member variables
size_t contact_patch_t::num_nodes_in_patch;

void contact_patch_t::update_nodes(const node_t &nodes) // NOLINT(*-make-member-function-const)
{
    // Constructing the points and vel_points arrays
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < num_nodes_in_patch; j++)
        {
            points(i, j) = nodes.coords(0, this->nodes_gid(j), i);
            vel_points(i, j) = nodes.vel(0, this->nodes_gid(j), i);
        }
    }
}

void contact_patches_t::initialize(const mesh_t &mesh, const CArrayKokkos<size_t> &bdy_contact_patches)
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
}


void contact_patches_t::sort(const node_t &nodes)
{
    // Update the points and vel_points arrays for each contact patch
    FOR_ALL_CLASS(i, 0, num_contact_patches, {
        contact_patches(i).update_nodes(nodes);
    });  // end parallel for
    Kokkos::fence();
}
