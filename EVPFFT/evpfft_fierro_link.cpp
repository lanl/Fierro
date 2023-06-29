#include "evpfft_fierro_link.h"
#include "evpfft.h"
#include <vector>

#if BUILD_EVPFFT_FIERRO

// each element will have it own evpfft model
std::vector<EVPFFT*> elem_evpfft;

void init_evpfft(const DCArrayKokkos <double> &file_state_vars,
                 const size_t num_state_vars,
                 const size_t mat_id,
                 const size_t num_elems)
{
    if (mat_id > 0) {
        throw std::runtime_error("evpfft-fierro-link is not implemented for multiple materials yet.");
    }

    real_t stress_scale = 1.0; // 1.0e-5; // used to convert MPa to MegaBar
    real_t time_scale = 1.0; // 1.0e+6; // used to convert second to microsecond

    CommandLineArgs cmd;
    cmd.nn = {8,8,8};
    cmd.input_filename = "evpfft.in";
    cmd.micro_filetype = 0;
    cmd.check_cmd_args();

    elem_evpfft = std::vector<EVPFFT*>(num_elems);

    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++)
    {
      elem_evpfft[elem_gid] = new EVPFFT(cmd,stress_scale,time_scale);
    }

}

void destroy_evpfft(const DCArrayKokkos <double> &file_state_vars,
                    const size_t num_state_vars,
                    const size_t mat_id,
                    const size_t num_elems)
{
    for (auto evpfft_ptr : elem_evpfft) {
        delete evpfft_ptr;
    }
}

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
                           const size_t cycle)
{
    real_t dt_rk = dt; // since using rk_num_stages = 1

    // Note EVPFFT uses F-layout while Fierro uses C-layout
    FArray <double> Fvel_grad(3,3);
    FArray <double> Fstress(3,3);
    // Transpose vel_grad
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        Fvel_grad(i,j) = 0.0; //vel_grad(i,j);
        Fstress(i,j) = elem_stress.host(1,elem_gid,i,j);
      }
    }

Fvel_grad(0,0) = 4000;
Fvel_grad(0,2) = -2000;
Fvel_grad(1,1) = 6000;
Fvel_grad(2,2) = -10000;

    elem_evpfft[elem_gid]->solve(Fvel_grad.pointer(), Fstress.pointer(), dt_rk, cycle, elem_gid);

    // Transpose stress. Not needed, stress is symmetric. But why not.
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        elem_stress.host(1,elem_gid,i,j) = Fstress(i,j);
      }
    }

    // write into elem_state_vars for output
    // evm,evmp,dvm,dvmp,svm
    elem_state_vars.host(elem_gid,0) = elem_evpfft[elem_gid]->evm;
    elem_state_vars.host(elem_gid,1) = elem_evpfft[elem_gid]->evmp;
    elem_state_vars.host(elem_gid,2) = elem_evpfft[elem_gid]->dvm;
    elem_state_vars.host(elem_gid,3) = elem_evpfft[elem_gid]->dvmp;
    elem_state_vars.host(elem_gid,4) = elem_evpfft[elem_gid]->svm;
}

#endif
