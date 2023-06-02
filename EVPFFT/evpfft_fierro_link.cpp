#include "evpfft_fierro_link.h"
#include "evpfft.h"
#include <vector>

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
    real_t dt_rk = dt*0.5;
    elem_evpfft[elem_gid]->solve(&vel_grad(0,0), &elem_stress.host(1,elem_gid,0,0), dt_rk, cycle, elem_gid);

    // write into elem_state_vars for output
    // evm,evmp,dvm,dvmp,svm
    elem_state_vars.host(elem_gid,0) = elem_evpfft[elem_gid]->evm;
    elem_state_vars.host(elem_gid,1) = elem_evpfft[elem_gid]->evmp;
    elem_state_vars.host(elem_gid,2) = elem_evpfft[elem_gid]->dvm;
    elem_state_vars.host(elem_gid,3) = elem_evpfft[elem_gid]->dvmp;
    elem_state_vars.host(elem_gid,4) = elem_evpfft[elem_gid]->svm;
}
