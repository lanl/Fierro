#include "UserMatModel.h"
#include "command_line_args.h" 
#include "evpfft.h"


void init_user_mat_model(std::shared_ptr<UserMatModel>* elem_user_mat_model,
                         const size_t* elem_mat_id,
                         const size_t num_elems)
{
  /*
      Initialize user_mat_model in this function
      Different user_mat_model per material can be initialized here using elem_mat_id(elem_gid)
  */

  real_t stress_scale = 1.0; // 1.0e-5; // used to convert MPa to MBar
  real_t time_scale = 1.0; // 1.0e+6; // used to convert second to microsecond

  CommandLineArgs cmd;
  cmd.nn = {8,8,8};
  cmd.input_filename = "fft.in";
  cmd.micro_filetype = 0;
  cmd.check_cmd_args();
 
  for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++)
  {
    elem_user_mat_model[elem_gid] = std::make_shared<EVPFFT>(cmd, stress_scale, time_scale);
    //elem_user_mat_model[elem_gid] = std::shared_ptr<UserMatModel>(new EVPFFT(cmd, stress_scale));
  }
}

void destroy_user_mat_model()
{}

void solve_user_mat_model()
{}
