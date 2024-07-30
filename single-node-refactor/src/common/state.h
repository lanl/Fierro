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
#ifndef STATE_H
#define STATE_H

#include "matar.h"

using namespace mtr;

/////////////////////////////////////////////////////////////////////////////
///
/// \struct node_t
///
/// \brief Stores state information associated with a node
///
/////////////////////////////////////////////////////////////////////////////
struct node_t
{
    DCArrayKokkos<double> coords; ///< Nodal coordinates
    MPIArrayKokkos<double> vel;  ///< Nodal velocity
    DCArrayKokkos<double> mass; ///< Nodal mass

    // initialization method (num_rk_storage_bins, num_nodes, num_dims)
    void initialize(size_t num_rk, size_t num_nodes, size_t num_dims)
    {
        this->coords = DCArrayKokkos<double>(num_rk, num_nodes, num_dims, "node_coordinates");
        this->vel    = MPIArrayKokkos<double>(num_rk, num_nodes, num_dims, "node_velocity");
        this->mass   = DCArrayKokkos<double>(num_nodes, "node_mass");
    }; // end method
}; // end node_t

/////////////////////////////////////////////////////////////////////////////
///
/// \struct Zone_t
///
/// \brief Stores state information associated with zone index space
///
/////////////////////////////////////////////////////////////////////////////
struct Zone_t
{
    DCArrayKokkos<double> sie; ///< coefficients for the sie polynomial field

    // initialization method (num_rk_storage_bins, num_nodes, num_dims)
    void initialize(size_t num_rk, size_t num_nodes, size_t num_dims)
    {
        this->sie = DCArrayKokkos<double>(num_rk, num_nodes, num_dims, "zone_sie");
    }; // end method

}; // end zone_t



/////////////////////////////////////////////////////////////////////////////
///
/// \struct GaussPoint_t
///
/// \brief Stores state information associated with the Gauss point
///
/////////////////////////////////////////////////////////////////////////////
struct GaussPoint_t
{
    DCArrayKokkos<double> vol;  ///< GAussPoint volume
    DCArrayKokkos<double> div;  ///< GAussPoint divergence of velocity

    DCArrayKokkos<size_t> mat_id; ///< MaterialPoint material index

    DCArrayKokkos<bool> eroded; ///< MaterialPoint eroded or not

    // initialization method (num_rk_storage_bins, num_cells, num_dims)
    void initialize(size_t num_rk, size_t num_elems, size_t num_dims)
    {
        this->vol    = DCArrayKokkos<double>(num_elems, "gauss_point_volume");
        this->div    = DCArrayKokkos<double>(num_elems, "gauss_point_div");
        this->mat_id = DCArrayKokkos<size_t>(num_elems, "gauss_point_mat_id");
        this->eroded = DCArrayKokkos<bool>(num_elems, "gauss_point_eroded");
    }; // end method

};  // end GuassPoint_t


/////////////////////////////////////////////////////////////////////////////
///
/// \struct MaterialPoint_t
///
/// \brief Stores state information associated with the material
///
/////////////////////////////////////////////////////////////////////////////
struct MaterialPoint_t
{
    DCArrayKokkos<double> den;    ///< MaterialPoint density
    DCArrayKokkos<double> pres;   ///< MaterialPoint pressure
    DCArrayKokkos<double> stress; ///< MaterialPoint stress
    DCArrayKokkos<double> sspd;   ///< MaterialPoint sound speed
    DCArrayKokkos<double> sie;    ///< MaterialPoint specific internal energy
    DCArrayKokkos<double> mass;   ///< MaterialPoint mass
    DCArrayKokkos<double> statev; ///< MaterialPoint state variable

    // initialization method (num_rk_storage_bins, num_cells, num_dims)
    void initialize(size_t num_rk, size_t num_elems, size_t num_dims)
    {
        this->den    = DCArrayKokkos<double>(num_elems, "material_point_density");
        this->pres   = DCArrayKokkos<double>(num_elems, "material_point_pressure");
        this->stress = DCArrayKokkos<double>(num_rk, num_elems, num_dims, num_dims, "material_point_stress");
        this->sspd   = DCArrayKokkos<double>(num_elems, "material_point_sspd");
        this->sie    = DCArrayKokkos<double>(num_rk, num_elems, "material_point_sie");  // only used with DG
        this->mass   = DCArrayKokkos<double>(num_elems, "material_point_mass");
    }; // end method

}; // end MaterialPoint

/////////////////////////////////////////////////////////////////////////////
///
/// \struct corner_t
///
/// \brief Stores state information associated with a corner
///
/////////////////////////////////////////////////////////////////////////////
struct corner_t
{
    DCArrayKokkos<double> force; ///< Corner force
    DCArrayKokkos<double> mass; ///< Corner mass

    // initialization method (num_corners, num_dims)
    void initialize(size_t num_corners, size_t num_dims)
    {
        this->force = DCArrayKokkos<double>(num_corners, num_dims, "corner_force");
        this->mass  = DCArrayKokkos<double>(num_corners, "corner_mass");
    }; // end method
}; // end corner_t

#endif
