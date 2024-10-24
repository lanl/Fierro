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

#ifndef NO_STRENGTH_H
#define NO_STRENGTH_H



// -----------------------------------------------------------------------------
// This is the user material model function for the stress tensor
// ------------------------------------------------------------------------------
namespace NoStrengthModel {

    static void init_strength_state_vars(
        const DCArrayKokkos <double> &MaterialPoints_eos_state_vars,
        const DCArrayKokkos <double> &MaterialPoints_strength_state_vars,
        const RaggedRightArrayKokkos <double> &eos_global_vars,
        const RaggedRightArrayKokkos <double> &strength_global_vars,
        const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const size_t num_material_points,
        const size_t mat_id)
    {

        // walk over all elements that have this material
        FOR_ALL(mat_points_lid, 0, num_material_points, {
            
            // get elem gid
            size_t elem_gid = MaterialToMeshMaps_elem(mat_points_lid); // might be used with some models

            // first index is matpt, second index is the number of vars
            size_t num_strength_state_vars = MaterialPoints_strength_state_vars.dims(1); 
            
            for(size_t var=0; var<num_strength_state_vars; var++){
                MaterialPoints_strength_state_vars(mat_points_lid,var) = 0.0;
            } // end for

        });  // end parallel for

    }  // end of init_strength_state_vars


    KOKKOS_FUNCTION
    static void calc_stress(
        const DCArrayKokkos<double>  &vel_grad,
        const DCArrayKokkos <double> &node_coords,
        const DCArrayKokkos <double> &node_vel,
        const DCArrayKokkos<size_t>  &nodes_in_elem,
        const DCArrayKokkos<double>  &MaterialPoints_pres,
        const DCArrayKokkos<double>  &MaterialPoints_stress,
        const DCArrayKokkos<double>  &MaterialPoints_sspd,
        const DCArrayKokkos <double> &MaterialPoints_eos_state_vars,
        const DCArrayKokkos <double> &MaterialPoints_strength_state_vars,
        const double MaterialPoints_den,
        const double MaterialPoints_sie,
        const DCArrayKokkos<double>& MaterialPoints_shear_modulii,
        const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const RaggedRightArrayKokkos <double> &eos_global_vars,
        const RaggedRightArrayKokkos <double> &strength_global_vars,
        const double vol,
        const double dt,
        const double rk_alpha,
        const double time,
        const size_t cycle,
        const size_t MaterialPoints_lid,
        const size_t mat_id,
        const size_t gauss_gid,
        const size_t elem_gid)
    {

        // -----------------------------------------------------------------------------
        // Required variables are here
        // ------------------------------------------------------------------------------

        // -----------------------------------------------------------------------------
        // The user must coding goes here
        // ------------------------------------------------------------------------------

        return;
    } // end of user mat

    
    static void destroy(
        const DCArrayKokkos <double> &MaterialPoints_eos_state_vars,
        const DCArrayKokkos <double> &MaterialPoints_strength_state_vars,
        const RaggedRightArrayKokkos <double> &eos_global_vars,
        const RaggedRightArrayKokkos <double> &strength_global_vars,
        const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const size_t num_material_points,
        const size_t mat_ids)
    {

    } // end destory

} // end namespace




#endif // end Header Guard