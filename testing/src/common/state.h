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
    CArray<double> coords;

    // velocity
    CArray<double> vel;

    // mass at nodes
    CArray<double> mass;

    // initialization method (num_rk_storage_bins, num_nodes, num_dims)
    void initialize(size_t num_rk, size_t num_nodes, size_t num_dims)
    {
        this->coords = CArray<double>(num_rk, num_nodes, num_dims);
        this->vel    = CArray<double>(num_rk, num_nodes, num_dims);
        this->mass   = CArray<double>(num_nodes);
    }; // end method
}; // end node_t

// elem_state
struct elem_t
{
    // den
    CArray<double> den;

    // pres
    CArray<double> pres;

    // stress
    CArray<double> stress;

    // sspd
    CArray<double> sspd;

    // sie
    CArray<double> sie;

    // vol
    CArray<double> vol;

    // divergence of velocity
    CArray<double> div;

    // mass of elem
    CArray<double> mass;

    // mat ids
    CArray<size_t> mat_id;

    // state variables
    CArray<double> statev;

    // initialization method (num_rk_storage_bins, num_cells, num_dims)
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
    // force
    CArray<double> force;

    // mass of corner
    CArray<double> mass;

    // initialization method (num_corners, num_dims)
    void initialize(size_t num_corners, size_t num_dims)
    {
        this->force = CArray<double>(num_corners, num_dims);
        this->mass  = CArray<double>(num_corners);
    }; // end method
}; // end corner_t

namespace model
{
// strength model types
enum strength_tag
{
    none = 0,
    hypo = 1,         // hypoelastic plastic model
    hyper = 2,        // hyperelastic plastic model
};
} // end of namespace

namespace model_init
{
// strength model setup
enum strength_setup_tag
{
    input = 0,
    user_init = 1,
};
} // end of namespace

// material model parameters
struct material_t
{
    // statev(0) = gamma
    // statev(1) = minimum sound speed
    // statev(2) = specific heat c_v
    // statev(3) = ref temperature
    // statev(4) = ref density
    // statev(5) = ref specific internal energy

    // eos fcn pointer
    void (*eos_model) (const DViewCArrayKokkos<double>& elem_pres,
                       const DViewCArrayKokkos<double>& elem_stress,
                       const size_t elem_gid,
                       const size_t mat_id,
                       const DViewCArrayKokkos<double>& elem_state_vars,
                       const DViewCArrayKokkos<double>& elem_sspd,
                       const double den,
                       const double sie) = NULL;

    // strength fcn pointer
    void (*strength_model) (const DViewCArrayKokkos<double>& elem_pres,
                            const DViewCArrayKokkos<double>& elem_stress,
                            const size_t elem_gid,
                            const size_t mat_id,
                            const DViewCArrayKokkos<double>& elem_state_vars,
                            const DViewCArrayKokkos<double>& elem_sspd,
                            const double den,
                            const double sie,
                            const ViewCArrayKokkos<double>&  vel_grad,
                            const ViewCArrayKokkos<size_t>&  elem_node_gids,
                            const DViewCArrayKokkos<double>& node_coords,
                            const DViewCArrayKokkos<double>& node_vel,
                            const double vol,
                            const double dt,
                            const double alpha) = NULL;

    // hypo or hyper elastic plastic model
    model::strength_tag strength_type;

    // setup the strength model via the input file for via a user_setup
    model_init::strength_setup_tag strength_setup = model_init::input;

    size_t num_state_vars;

    double q1;    // acoustic coefficient in Riemann solver for compresion
    double q1ex;  // acoustic coefficient in Riemann solver for expansion
    double q2;    // linear coefficient in Riemann solver for compression
    double q2ex;  // linear coefficient in Riemann solver for expansion
}; // end material_t

#endif
