

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

#ifndef NO_DISSIPATION_H
#define NO_DISSIPATION_H

namespace NoDissipationModel {

 
    KOKKOS_FUNCTION
    static void calc_dissipation (const ViewCArrayKokkos<size_t> elem_node_gids,
                                  const RaggedRightArrayKokkos <double>& dissipation_global_vars,
                                  const DCArrayKokkos<double>& GaussPoints_vel_grad,
                                  const DCArrayKokkos<bool>&   MaterialPoints_eroded,
                                  const DCArrayKokkos<double>& node_vel,
                                  const DCArrayKokkos<double>& MaterialPoints_den,
                                  const DCArrayKokkos<double>& MaterialPoints_sspd,
                                  const ViewCArrayKokkos<double>& disp_corner_forces,
                                  const ViewCArrayKokkos<double>& area_normal,
                                  const RaggedRightArrayKokkos<size_t>& elems_in_elem,
                                  const CArrayKokkos<size_t>& num_elems_in_elem,
                                  const double vol,
                                  const double fuzz,
                                  const double small,
                                  const double elem_gid,
                                  const size_t mat_point_lid,
                                  const size_t mat_id)
                                  {

                                    return;
                                  }
    
    } // end namespace


#endif // end Header Guard