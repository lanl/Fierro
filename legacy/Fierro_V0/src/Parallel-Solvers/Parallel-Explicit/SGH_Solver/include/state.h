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
#ifndef STATE_H
#define STATE_H

#include "matar.h"

using namespace mtr;

// node_state
struct node_t
{
    // Position
    DCArrayKokkos<double> coords;

    // velocity
    DCArrayKokkos<double> vel;

    // mass at nodes
    DCArrayKokkos<double> mass;

    // Includes Ghost Positions
    // DCArrayKokkos <double> all_coords;

    // Includes Ghost velocities
    // DCArrayKokkos <double> all_vel;

    // Includes Ghost masses
    // DCArrayKokkos <double> all_mass;

    // initialization method (num_rk_storage_bins, num_nodes, num_dims)
    void initialize(size_t num_rk, size_t num_nodes, size_t num_dims)
    {
        this->coords = DCArrayKokkos<double>(num_rk, num_nodes, num_dims);
        this->vel    = DCArrayKokkos<double>(num_rk, num_nodes, num_dims);
        this->mass   = DCArrayKokkos<double>(num_nodes);
    }; // end method
}; // end node_t

// elem_state
struct elem_t
{
    CArray<double> den; ///< element density
    CArray<double> pres; ///< element pressure
    CArray<double> stress; ///< element stress
    CArray<double> sspd; ///< element sound speed
    CArray<double> sie; ///< specific internal energy
    CArray<double> vol; ///< element volume
    CArray<double> div; ///< divergence of velocity
    CArray<double> mass; ///< element mass
    CArray<size_t> mat_id; ///< element material id

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn initialize
    ///
    /// \brief Initialize state data associated with the element
    ///
    /// This function initializes storage for density, pressure, stress,
    /// sound speed, specific internal energy, volume, divergence of velocity,
    /// mass, and a material ID for each element in the mesh.
    ///
    /// \param Number of time integration bins required for the Runge Kutta scheme
    /// \param Number of elements in the mesh
    /// \param Number of spatial dimnesions
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize(size_t num_rk, size_t num_elems, size_t num_dims)
    {
        this->den    = CArray<double>(num_elems);
        this->pres   = CArray<double>(num_elems);
        this->stress = CArray<double>(num_rk, num_elems, num_dims, num_dims);
        this->sspd   = CArray<double>(num_elems);
        this->sie    = CArray<double>(num_rk, num_elems);
        this->vol    = CArray<double>(num_elems);
        this->div    = CArray<double>(num_elems);
        this->mass   = CArray<double>(num_elems);
        this->mat_id = CArray<size_t>(num_elems);
    }; // end method
}; // end elem_t

// corner_state
struct corner_t
{
    CArray<double> force; ///< Force acting on a corner
    CArray<double> mass; ///< Partitioned mass of the corner

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn initialize
    ///
    /// \brief Initialize state data associated with the corners
    ///
    /// This function initializes storage for force and mass associated with each
    /// corner
    ///
    /// \param Number of corners in the mesh
    /// \param Number of spatial dimnesions
    ///
    /////////////////////////////////////////////////////////////////////////////
    void initialize(size_t num_corners, size_t num_dims)
    {
        this->force = CArray<double>(num_corners, num_dims);
        this->mass  = CArray<double>(num_corners);
    }; // end method
}; // end corner_t

#endif
