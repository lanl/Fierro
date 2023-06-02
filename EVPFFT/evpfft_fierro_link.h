#pragma once

#include "matar.h"

using namespace mtr;

void init_evpfft(const DCArrayKokkos <double> &file_state_vars,
                 const size_t num_state_vars,
                 const size_t mat_id,
                 const size_t num_elems);

void destroy_evpfft(const DCArrayKokkos <double> &file_state_vars,
                    const size_t num_state_vars,
                    const size_t mat_id,
                    const size_t num_elems);

void evpfft_strength_model(const DViewCArrayKokkos <double> &elem_pres,
                           const DViewCArrayKokkos <double> &elem_stress,
                           const size_t elem_gid,
                           const size_t mat_id,
                           const DViewCArrayKokkos <double> &elem_state_vars,
                           const DViewCArrayKokkos <double> &elem_sspd,
                           const double den,
                           const double sie,
                           const ViewCArrayKokkos <double> &vel_grad,
                           const ViewCArrayKokkos <size_t> &elem_node_gids,
                           const DViewCArrayKokkos <double> &node_coords,
                           const DViewCArrayKokkos <double> &node_vel,
                           const double vol,
                           const double dt,
                           const double rk_alpha,
                           const size_t cycle);

