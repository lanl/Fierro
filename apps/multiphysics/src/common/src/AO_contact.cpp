/**********************************************************************************************
© 2020. Triad National Security, LLC. All rights reserved.
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

#include "AO_contact.hpp"

void AO_contact_initialize(DCArrayKokkos <double>& bdy_node_coords,
                           const size_t num_bdy_nodes,
                           DCArrayKokkos <size_t>& num_nodes_in_bounding_boxes,
                           const size_t num_bdy_surfs,
                           CArrayKokkos <double>& bounding_boxes)
{
    // allocating member variables that are statically sized
    bdy_node_coords = DCArrayKokkos <double> (num_bdy_nodes, 3, "bdy_node_coords");
    num_nodes_in_bounding_boxes = DCArrayKokkos <size_t> (num_bdy_surfs, "num_nodes_in_bounding_boxes");
    bounding_boxes = CArrayKokkos <double> (num_bdy_surfs, 6, "bounding_boxes");

    

    return;
} // end AO_contact_initialize

void get_bounding_boxes(
                            )
{
    
}

void AO_contact_sort(DCArrayKokkos <double>& bdy_node_coords,
                     const size_t num_bdy_nodes,
                     const MPICArrayKokkos <double>& node_coords,
                     const CArrayKokkos <size_t>& bdy_nodes,
                     swage::PointCloud_t& bdy_node_point_cloud,
                     const size_t num_bins)
{
    // updating bdy_node_coords
    FOR_ALL(i, 0, (int)num_bdy_nodes,
            j, 0, 3, 
            {
                bdy_node_coords(i,j) = node_coords(bdy_nodes(i), j);
            });
    Kokkos::fence();

    // getting domain bounds
    bdy_node_point_cloud.get_bounds_point_cloud(bdy_node_point_cloud.xmin, bdy_node_point_cloud.ymin, bdy_node_point_cloud.zmin,
                                                bdy_node_point_cloud.xmax, bdy_node_point_cloud.ymax, bdy_node_point_cloud.zmax,
                                                bdy_node_coords);

    // building underlying bin mesh
    const double bin_mesh_sizing_factor_x = (bdy_node_point_cloud.xmax - bdy_node_point_cloud.xmin) / num_bins * 0.1; // setting bin mesh to extend 10 percent of a bin
    const double bin_mesh_sizing_factor_y = (bdy_node_point_cloud.ymax - bdy_node_point_cloud.ymin) / num_bins * 0.1; // setting bin mesh to extend 10 percent of a bin
    const double bin_mesh_sizing_factor_z = (bdy_node_point_cloud.zmax - bdy_node_point_cloud.zmin) / num_bins * 0.1; // setting bin mesh to extend 10 percent of a bin
    bdy_node_point_cloud.build_bin_mesh(bdy_node_point_cloud.xmin - bin_mesh_sizing_factor_x,
                                        bdy_node_point_cloud.ymin - bin_mesh_sizing_factor_y,
                                        bdy_node_point_cloud.zmin - bin_mesh_sizing_factor_z,
                                        bdy_node_point_cloud.xmax + bin_mesh_sizing_factor_x,
                                        bdy_node_point_cloud.ymax + bin_mesh_sizing_factor_y,
                                        bdy_node_point_cloud.zmax + bin_mesh_sizing_factor_z,
                                        num_bins, num_bins, num_bins);

    // getting a bounding box for each boundary surface


    return;
} // end AO_contact_sort