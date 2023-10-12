#include "UserStrengthModel.h"
#include <string>
#include "evpfft.h"

MPI_Comm evpfft_mpi_comm = MPI_COMM_NULL;

KOKKOS_FUNCTION
UserStrengthModel::UserStrengthModel(
  const DCArrayKokkos <material_t> &material,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t mat_id,
  const size_t elem_gid)
{
    // First, lets create a new communicator with each rank having its own communicator containing only itself.
    if (evpfft_mpi_comm == MPI_COMM_NULL) {
      int global_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
      MPI_Comm_split(MPI_COMM_WORLD, global_rank, global_rank, &evpfft_mpi_comm);
    }
    //...

    // Input files for for multiple materials should be names as evpfft1.in, evpfft2.in, etc.
    std::string filename = "evpfft" + std::to_string(mat_id+1) + ".in";

    // get dimensions
    int N1 = global_vars.host(mat_id,0);
    int N2 = global_vars.host(mat_id,1);
    int N3 = global_vars.host(mat_id,2);

    real_t stress_scale = 1.0; // 1.0e-5; // used to convert MPa to MegaBar
    real_t time_scale = 1.0; // 1.0e+6; // used to convert second to microsecond

    CommandLineArgs cmd;
    //cmd.nn = {N1,N2,N3};
    cmd.input_filename = filename; //"evpfft.in";
    cmd.micro_filetype = 0;
    cmd.check_cmd_args();

    evpfft_ptr = new EVPFFT(evpfft_mpi_comm, cmd,stress_scale,time_scale);
}

KOKKOS_FUNCTION
UserStrengthModel::~UserStrengthModel()
{
    delete evpfft_ptr;
}

KOKKOS_FUNCTION
int UserStrengthModel::calc_stress(
  const DViewCArrayKokkos <double> &elem_pres,
  const DViewCArrayKokkos <double> &elem_stress,
  const size_t elem_gid,
  const size_t mat_id,
  const DCArrayKokkos <double> &global_vars,
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
  const size_t rk_level)
{
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

    double udotAccTh = global_vars.host(mat_id,3); // Linear Aprox. Threshold
    evpfft_ptr->solve(Fvel_grad.pointer(), Fstress.pointer(), dt_rk, cycle, elem_gid, udotAccTh);

    // Transpose stress. Not needed, stress is symmetric. But why not.
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        elem_stress.host(rk_level,elem_gid,i,j) = Fstress(i,j);
      }
    }

    // write into elem_state_vars for output
    // evm,evmp,dvm,dvmp,svm
    elem_user_output_vars.host(elem_gid,0) = evpfft_ptr->evm;
    elem_user_output_vars.host(elem_gid,1) = evpfft_ptr->evmp;
    elem_user_output_vars.host(elem_gid,2) = evpfft_ptr->dvm;
    elem_user_output_vars.host(elem_gid,3) = evpfft_ptr->dvmp;
    elem_user_output_vars.host(elem_gid,4) = evpfft_ptr->svm;

    return 0;
}

