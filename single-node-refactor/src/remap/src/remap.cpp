/**********************************************************************************************
� 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and
to permit others to do so.
This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or other materials
provided with the distribution.
3.  Neither the name of the copyright holder nor the names of its contributors may be used
to endorse or promote products derived from this software without specific prior
written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************************/

#include "remap.h"
#include "matar.h"
#include "mesh.h"
#include "material.h"
#include "state.h"

// NOTES:
// Option1: purely cell centric
//   loop over all patches in mesh
//   get elem_gid1 and elem_gid2 for a patch
//   loop over mats in these elems
//   use reverse map to find mat_storage_lid_1 and mat_storage_lid_2 for given mat_id
//   advect all vars for that material (requires atomics)

// Option2: material centric
//    loop over all mat_elems
//    loop over the patches
//    get neighbor elem_gid
//    find which bin my mat_id is in
//    advect all vars for that material (requires atomics)

// Option3: material centric patch (best option)
//    loop over all mat_patches
//    get neighboring elem_gids
//    find which bin the mat_id is in
//    advect all vars for that material (requires atomics)

// ----
//   saving fluxes to patches and doing a second loop over elements and 
//   patches in elem is thread safe



void advect_mat_var(const Mesh_t& mesh,
                    const DCArrayKokkos<double>& node_coords,     // coords after mesh smoothing
                    const DCArrayKokkos<double>& node_coords_n0,  // coords_n0 are prior to remap
                    const DRaggedRightArrayKokkos<double>& MaterialPoints_den, 
                    const DRaggedRightArrayKokkos<double>& MaterialPoints_sie, 
                    const DRaggedRightArrayKokkos<double>& MaterialPoints_stress, 
                    const DRaggedRightArrayKokkos<double>& MaterialPoints_volfrac,     // a remapped var
                    const DRaggedRightArrayKokkos<double>& MaterialPoints_geo_volfrac, // interface recon
                    const DRaggedRightArrayKokkos<size_t>& MaterialToMeshMaps_elem,
                    const DRaggedRightArrayKokkos<size_t>& MeshToMaterialMaps_mat_id,
                    const DRaggedRightArrayKokkos<size_t>& MeshToMaterialMaps_mat_storage_lid,
                    const DCArrayKokkos<size_t>& MeshToMaterialMaps_num_mats_in_elem,
                    const DCArrayKokkos<size_t>& MaterialToMeshMaps_patches,
                    const size_t num_mat_patches,
                    const size_t mat_id,
                    const double fuzz,
                    const double small)
    {

        // Loop over the patches on the mesh
        FOR_ALL(mat_patch_lid, 0, num_mat_patches, {

            // get patches gid
            size_t patch_gid = MaterialToMeshMaps_patches(mat_id, mat_patch_lid); 

            // get elem gids
            const size_t elem_gid_0 = mesh.elems_in_patch(patch_gid, 0); 
            const size_t elem_gid_1 = mesh.elems_in_patch(patch_gid, 1);
            
            // the material elem storage index
            size_t mat_elem_lid_0;
            size_t mat_elem_lid_1;

            // loop over mat_lid in elem_gid_0
            for (size_t mat_lid=0; mat_lid<MeshToMaterialMaps_num_mats_in_elem(elem_gid_0); mat_lid++){
                if (MeshToMaterialMaps_mat_id(elem_gid_0, mat_lid)==mat_id){
                    mat_elem_lid_0 = MeshToMaterialMaps_mat_storage_lid(elem_gid_0, mat_lid);
                } // end if
            } // end for mat_lid

            // loop over mat_lid in elem_gid_1
            for (size_t mat_lid=0; mat_lid<MeshToMaterialMaps_num_mats_in_elem(elem_gid_1); mat_lid++){
                if (MeshToMaterialMaps_mat_id(elem_gid_1, mat_lid)==mat_id){
                    mat_elem_lid_1 = MeshToMaterialMaps_mat_storage_lid(elem_gid_1, mat_lid);
                } // end if
            } // end for mat_lid


            // --- calculate surfaced area ---
            double avg_coords[3];
            avg_coords[0] = 0.0;
            avg_coords[1] = 0.0;
            avg_coords[2] = 0.0;
            for (size_t node_lid=0; node_lid<mesh.num_nodes_in_patch; node_lid++){
                
                // get the node id
                size_t node_gid = mesh.nodes_in_patch(patch_gid, node_lid);


                for (size_t dim=0; dim<mesh.num_dims; dim++){
                    avg_coords[dim] += node_coords(node_gid, dim);
                } // end for dim

            } // end for

            for (size_t dim=0; dim<mesh.num_dims; dim++){
                avg_coords[dim] /= (double)mesh.num_nodes_in_patch;
            } // end for dim


            // allocate the corner surface normals and forces 
            double corn_patch_area_normal_1D[3];
            ViewCArrayKokkos <double> corn_patch_area_normal(&corn_patch_area_normal_1D[0],3);
            
            // initialize corner patch vector to zero
            for (size_t dim=0; dim<mesh.num_dims; dim++){
                corn_patch_area_normal(dim) = 0.0;
            } // end for dim

            double vec_a[3];
            double vec_b[3];
            for (size_t node_lid=0; node_lid<mesh.num_nodes_in_patch; node_lid++){
                
                // get the node ids for the triangle
                size_t node_gid_0 = mesh.nodes_in_patch(patch_gid, node_lid);
                size_t node_gid_1;
                if(node_lid<mesh.num_nodes_in_patch-1){
                    node_gid_1 = mesh.nodes_in_patch(patch_gid, node_lid+1);
                } else {
                    node_gid_1 = mesh.nodes_in_patch(patch_gid, 0);
                } // end if
        
                for (size_t dim=0; dim<mesh.num_dims; dim++){
                    vec_a[dim] = node_coords(node_gid_0, dim) - avg_coords[dim];
                    vec_b[dim] = node_coords(node_gid_1, dim) - avg_coords[dim];
                } // end for dim

                // calculating the cross product of the 2 vectors, 1/2 multiply is later
                corn_patch_area_normal(0) += (vec_a[1]*vec_b[2] - vec_a[2]*vec_b[1]);
                corn_patch_area_normal(1) += (vec_a[2]*vec_b[0] - vec_a[0]*vec_b[2]);
                corn_patch_area_normal(2) += (vec_a[0]*vec_b[1] - vec_a[1]*vec_b[0]);

            } // end for node_lid in patch

            // each node gets 1/4 of the surface area
            for (size_t dim=0; dim<mesh.num_dims; dim++){
                corn_patch_area_normal(dim) /= 8;  //  comes from 1/2 * 1/4
            } // end for dim

        });

    } // end advect function

