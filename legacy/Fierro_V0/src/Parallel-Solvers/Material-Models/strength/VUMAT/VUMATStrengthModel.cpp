#include "VUMATStrengthModel.h"
#include "VUMAT.h"
#include <vector>
#include <memory>


// to hold all vumat in each element
std::vector<VUMAT*> elem_vumat;

namespace VUMATStrengthModel
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
    /*
    In this function, initialize the strength_state_vars according to the location this model needs it.
    DO NOT strength_state_vars.update_host or strength_state_vars.update_device.
    This is because different elements might have different models, run location, and 
    initialization location leading to overwriting of data.
    */

    // assign correct size to elem_vumat, all element are nullptr during initialization
    elem_vumat = std::vector<VUMAT*> (num_elems, nullptr);

    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {

      size_t mat_id = elem_mat_id.host(elem_gid);

      // only fill the elem_vumat of elements that use this model
      if (material.host(mat_id).strength_model == STRENGTH_MODEL::vumat) {
        
        // currently vumat only runs on host so check to see
        if (material.host(mat_id).strength_run_location == RUN_LOCATION::device) {
          throw std::runtime_error("VUMAT only runs on Host");
        }

        // create VUMAT model in element that used vumat
        elem_vumat[elem_gid] = new VUMAT(material,
                                         eos_state_vars,
                                         strength_state_vars,
                                         eos_global_vars,
                                         strength_global_vars,
                                         elem_user_output_vars,
                                         mat_id,
                                         elem_gid);

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
    if (elem_vumat[elem_gid] == nullptr)
    {
      throw std::runtime_error("VUMAT not initialized in this element");
    }

    elem_vumat[elem_gid]->calc_stress(elem_pres,
                                      elem_stress,
                                      elem_gid,
                                      mat_id,
                                      eos_state_vars,
                                      strength_state_vars,
                                      eos_global_vars,
                                      strength_global_vars,
                                      elem_user_output_vars,
                                      elem_sspd,
                                      den,
                                      sie,
                                      vel_grad,
                                      elem_node_gids,
                                      node_coords,
                                      node_vel,
                                      vol,
                                      dt,
                                      rk_alpha,
                                      cycle,
                                      rk_level,
                                      time);

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
    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
      delete elem_vumat[elem_gid];
    }
  } 
     
} // end namespace VUMATStrengthModel