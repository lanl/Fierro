#include "contact.h"

void contact_patch_t::initialize()
{

}

void contact_patches_t::initialize(const mesh_t &mesh, const CArrayKokkos<size_t> &bdy_contact_patches)
{
    if (mesh.num_dims != 3)
    {
        std::cerr << "Error: contact is only supported in 3D" << std::endl;
        exit(1);
    }

    patches_gid = bdy_contact_patches;
    num_contact_patches = patches_gid.size();

    contact_patches = CArrayKokkos<contact_patch_t>(num_contact_patches);

    CArrayKokkos<size_t> nodes_in_patch(num_contact_patches, mesh.num_nodes_in_patch);
    FOR_ALL_CLASS(i, 0, num_contact_patches,
                  j, 0, mesh.num_nodes_in_patch, {
                      nodes_in_patch(i, j) = mesh.nodes_in_patch(patches_gid(i), j);
                  });
    Kokkos::fence();

    for (int i = 0; i < num_contact_patches; i++)
    {
        contact_patch_t &contact_patch = contact_patches(i);
        contact_patch.gid = patches_gid(i);

        // Make contact_patch.nodes_gid equal to the row of nodes_in_patch(i)
        contact_patch.nodes_gid = CArrayKokkos<size_t>(mesh.num_nodes_in_patch);
        for (size_t j = 0; j < mesh.num_nodes_in_patch; j++)
        {
            contact_patch.nodes_gid(j) = nodes_in_patch(i, j);
        }
    }
    Kokkos::fence();
}
