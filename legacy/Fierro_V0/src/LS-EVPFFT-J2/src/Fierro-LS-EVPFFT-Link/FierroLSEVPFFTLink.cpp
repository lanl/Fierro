#include "FierroLSEVPFFTLink.h"
#include <string>
#include "evpfft.h"


MPI_Comm evpfft_mpi_comm = MPI_COMM_NULL;
// to hold all evpfft in each element
std::vector<EVPFFT*> elem_evpfft;

namespace FierroLSEVPFFTLink
{
  void init_strength_state_vars(
    const DCArrayKokkos <material_t> &material,
    const DViewCArrayKokkos <size_t> &elem_mat_id,
    const DCArrayKokkos <double> &eos_state_vars,
    const DCArrayKokkos <double> &strength_state_vars,
    const DCArrayKokkos <double> &eos_global_vars,
    const DCArrayKokkos <double> &strength_global_vars,
    const DCArrayKokkos <double> &elem_user_output_vars,
    const size_t num_elems)
    {
      printf("Executing FierroLSEVPFFTLink::init_strength_state_vars ...\n");

      // First, lets create a new communicator with each rank having its own communicator containing only itself.
      if (evpfft_mpi_comm == MPI_COMM_NULL) {
        int global_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
        MPI_Comm_split(MPI_COMM_WORLD, global_rank, global_rank, &evpfft_mpi_comm);
      }
      //...

      // assign correct size to elem_evpfft, all element are nullptr during initialization
      elem_evpfft = std::vector<EVPFFT*> (num_elems, nullptr);

      for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {

        size_t mat_id = elem_mat_id.host(elem_gid);

        // only fill the elem_evpfft of elements that use this model
        if (material.host(mat_id).strength_model == STRENGTH_MODEL::ls_evpfft) {
          
          // evpfft only runs on host so check to see
          if (material.host(mat_id).strength_run_location == RUN_LOCATION::device) {
            throw std::runtime_error("EVPFFT only runs on Host");
          }

          // Input files for for multiple materials should be names as evpfft1.in, evpfft2.in, etc.
          std::string filename = "evpfft" + std::to_string(mat_id+1) + ".in";

          real_t stress_scale = 1.0; // 1.0e-5; // used to convert MPa to MegaBar
          real_t time_scale = 1.0; // 1.0e+6; // used to convert second to microsecond

          CommandLineArgs cmd;
          cmd.input_filename = filename; //"evpfft.in";
          cmd.micro_filetype = 0;
          cmd.check_cmd_args();

          // create EVPFFT model in element that used evpfft
          elem_evpfft[elem_gid] = new EVPFFT(evpfft_mpi_comm,
                                             cmd,
                                             stress_scale,
                                             time_scale);

        } // end if (material.host(mat_id).strength_model...

      } // end for (size_t elem_gid = 0... 

      return;     
    }

  void calc_stress (
    const DViewCArrayKokkos <double> &elem_pres,
    const DViewCArrayKokkos <double> &elem_stress,
    const size_t elem_gid,
    const size_t mat_id,
    const DCArrayKokkos <double> &eos_state_vars,
    const DCArrayKokkos <double> &strength_state_vars,
    const DCArrayKokkos <double> &eos_global_vars,
    const DCArrayKokkos <double> &strength_global_vars,
    const DCArrayKokkos <double> &elem_user_output_vars,
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
    const size_t cycle,
    const size_t rk_level,
    const double time)
    {
      if (elem_evpfft[elem_gid] == nullptr)
      {
        throw std::runtime_error("LSEVPFFT not initialized in this element");
      }

      real_t dt_rk = dt; // since using rk_num_stages = 1

      // Note EVPFFT uses F-layout while Fierro uses C-layout
      FArray <double> Fvel_grad(3,3);
      FArray <double> Fstress(3,3);
      // Transpose vel_grad
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          Fvel_grad(i,j) = vel_grad(i,j);
          Fstress(i,j) = elem_stress.host(rk_level,elem_gid,i,j);
        }
      }

      double udotAccTh = strength_global_vars.host(mat_id,0); // Linear Aprox. Threshold
      elem_evpfft[elem_gid]->solve(Fvel_grad.pointer(), Fstress.pointer(), dt_rk, cycle, elem_gid, udotAccTh);

      // Transpose stress. Not needed, stress is symmetric. But why not.
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          elem_stress.host(rk_level,elem_gid,i,j) = Fstress(i,j);
        }
      }

      // write into elem_state_vars for output
      // evm,evmp,dvm,dvmp,svm
      elem_user_output_vars.host(elem_gid,0) = elem_evpfft[elem_gid]->evm;
      elem_user_output_vars.host(elem_gid,1) = elem_evpfft[elem_gid]->evmp;
      elem_user_output_vars.host(elem_gid,2) = elem_evpfft[elem_gid]->dvm;
      elem_user_output_vars.host(elem_gid,3) = elem_evpfft[elem_gid]->dvmp;
      elem_user_output_vars.host(elem_gid,4) = elem_evpfft[elem_gid]->svm;

      return;
    }

  void destroy(
    const DCArrayKokkos <material_t> &material,
    const DViewCArrayKokkos <size_t> &elem_mat_id,
    const DCArrayKokkos <double> &eos_state_vars,
    const DCArrayKokkos <double> &strength_state_vars,
    const DCArrayKokkos <double> &eos_global_vars,
    const DCArrayKokkos <double> &strength_global_vars,
    const DCArrayKokkos <double> &elem_user_output_vars,
    const size_t num_elems)
    {
      printf("Executing FierroEVPFFTLink::destroy ...\n");

      for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
        delete elem_evpfft[elem_gid];
      }

      return;
    }
    
} // end namespace FierroLSEVPFFTLink

