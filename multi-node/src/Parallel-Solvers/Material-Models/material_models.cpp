#include "material_models.h"
#include "IdealGasEOSModel.h"
#include "IdealGasStrengthModel.h"
#include "UserEOSModel.h"
#include "UserStrengthModel.h"

/* EOSParent */
KOKKOS_FUNCTION
EOSParent::EOSParent() {}

KOKKOS_FUNCTION
EOSParent::~EOSParent() {}

/* sterngth_parent */
KOKKOS_FUNCTION
StrengthParent::StrengthParent() {}

KOKKOS_FUNCTION
StrengthParent::~StrengthParent() {}

template <typename P, typename C>
P* allocate_memory(RUN_LOCATION run_loc)
{
  P* ptr = nullptr;
  if (run_loc == RUN_LOCATION::device) {
    ptr = (C*)Kokkos::kokkos_malloc<Kokkos::DefaultExecutionSpace::memory_space>(sizeof(C));
  }
  else if (run_loc == RUN_LOCATION::host) {
    ptr = (C*)Kokkos::kokkos_malloc<Kokkos::HostSpace>(sizeof(C));
  }
  else {
    throw std::runtime_error("run_loc should be device or host");
  }
  return ptr;
}

template <typename P>
void free_memory(P* model_ptr, RUN_LOCATION run_loc)
{
  static_assert(std::is_same<EOSParent, P>::value or std::is_same<StrengthParent, P>::value, 
                "Invalid type. Type must be EOSParent or StrengthParent");

  if (model_ptr == nullptr)
    return;

  if (run_loc == RUN_LOCATION::device) {
    Kokkos::kokkos_free<Kokkos::DefaultExecutionSpace::memory_space>(model_ptr);
    model_ptr = nullptr;
  }
  else if (run_loc == RUN_LOCATION::host) {
    Kokkos::kokkos_free<Kokkos::HostSpace>(model_ptr);
    model_ptr = nullptr;
  }
  else {
    throw std::runtime_error("run_loc should be device or host");
  }
  
}

void init_strength_model(
  DCArrayKokkos <strength_t> &elem_strength,
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t num_elems)
{
  /*
  Material model strength should be initialized here.
  */

  // Allocate memory on GPU or CPU depending on model run_location
  for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {
    size_t mat_id = elem_mat_id.host(elem_gid);
    auto run_loc = material.host(mat_id).strength_run_location;
    auto strength_model = material.host(mat_id).strength_model;
    switch (strength_model) {
      case STRENGTH_MODEL::ideal_gas:
        elem_strength.host(elem_gid).model = allocate_memory <StrengthParent, IdealGasStrengthModel> (run_loc);
        break;
      case STRENGTH_MODEL::user_strength_model:
        elem_strength.host(elem_gid).model = allocate_memory <StrengthParent, UserStrengthModel> (run_loc);
        break;
      default:
        break;
    } // end switch
  }
  elem_strength.update_device();

  // Create model on GPU using `placement new` for models that run on device
  FOR_ALL(elem_gid, 0, num_elems, {
    size_t mat_id = elem_mat_id(elem_gid);
    auto run_loc = material(mat_id).strength_run_location;
    auto strength_model = material(mat_id).strength_model;
    if (run_loc == RUN_LOCATION::device) {
      switch (strength_model) {
        case STRENGTH_MODEL::ideal_gas:
          new ((IdealGasStrengthModel*)elem_strength(elem_gid).model) IdealGasStrengthModel(material, global_vars, elem_user_output_vars, mat_id, elem_gid);
          break;
        case STRENGTH_MODEL::user_strength_model:
          new ((UserStrengthModel*)elem_strength(elem_gid).model) UserStrengthModel(material, global_vars, elem_user_output_vars, mat_id, elem_gid);
          break;
        default:
          break;
      } // end switch
    } // end if
  });
  Kokkos::fence();


  // Create model on CPU using `placement new` for models that run on host
  for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {
    size_t mat_id = elem_mat_id.host(elem_gid);
    auto run_loc = material.host(mat_id).strength_run_location;
    auto strength_model = material.host(mat_id).strength_model;
    if (run_loc == RUN_LOCATION::host) {
      switch (strength_model) {
        case STRENGTH_MODEL::ideal_gas:
          new ((IdealGasStrengthModel*)elem_strength.host(elem_gid).model) IdealGasStrengthModel(material, global_vars, elem_user_output_vars, mat_id, elem_gid);
          break;
        case STRENGTH_MODEL::user_strength_model:
          new ((UserStrengthModel*)elem_strength.host(elem_gid).model) UserStrengthModel(material, global_vars, elem_user_output_vars, mat_id, elem_gid);
          break;
        default:
          break;
      } // end switch
    } // end if
  }

} // end init_user_strength_model


void init_eos_model(
  DCArrayKokkos <eos_t> &elem_eos,
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t num_elems)
{
  /*
  Material EOS model should be initialized here.
  */

  // Allocate memory on GPU or CPU depending on model run_location
  for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {
    size_t mat_id = elem_mat_id.host(elem_gid);
    auto run_loc = material.host(mat_id).eos_run_location;
    auto eos_model = material.host(mat_id).eos_model;
    switch (eos_model) {
      case EOS_MODEL::ideal_gas:
        elem_eos.host(elem_gid).model = allocate_memory <EOSParent, IdealGasEOSModel> (run_loc);
        break;
      case EOS_MODEL::user_eos_model:
        elem_eos.host(elem_gid).model = allocate_memory <EOSParent, UserEOSModel> (run_loc);
        break;
      default:
        break;
    } // end switch
  }
  elem_eos.update_device();

  // Create model on GPU using `placement new` for models that run on device
  FOR_ALL(elem_gid, 0, num_elems, {
    size_t mat_id = elem_mat_id(elem_gid);
    auto run_loc = material(mat_id).eos_run_location;
    auto eos_model = material(mat_id).eos_model;
    if (run_loc == RUN_LOCATION::device) {
      switch (eos_model) {
        case EOS_MODEL::ideal_gas:
          new ((IdealGasEOSModel*)elem_eos(elem_gid).model) IdealGasEOSModel(material, global_vars, elem_user_output_vars, mat_id, elem_gid);
          break;
        case EOS_MODEL::user_eos_model:
          new ((UserEOSModel*)elem_eos(elem_gid).model) UserEOSModel(material, global_vars, elem_user_output_vars, mat_id, elem_gid);
          break;
        default:
          break;
      } // end switch
    } // end if
  });
  Kokkos::fence();


  // Create model on CPU using `placement new` for models that run on host
  for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {
    size_t mat_id = elem_mat_id.host(elem_gid);
    auto run_loc = material.host(mat_id).eos_run_location;
    auto eos_model = material.host(mat_id).eos_model;
    if (run_loc == RUN_LOCATION::host) {
      switch (eos_model) {
        case EOS_MODEL::ideal_gas:
          new ((IdealGasEOSModel*)elem_eos.host(elem_gid).model) IdealGasEOSModel(material, global_vars, elem_user_output_vars, mat_id, elem_gid);
          break;
        case EOS_MODEL::user_eos_model:
          new ((UserEOSModel*)elem_eos.host(elem_gid).model) UserEOSModel(material,  global_vars, elem_user_output_vars, mat_id, elem_gid);
          break;
        default:
          break;
      } // end switch
    } // end if 
  }

} // end init_user_eos_model

void destroy_strength_model(
  DCArrayKokkos <strength_t> &elem_strength,
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> elem_user_output_vars,
  const size_t num_elems)
{
  /*
    All memory cleanup related to the user material model should be done in this fuction.
    Fierro calls `destroy_strength_model` at the end of a simulation.
  */

  // Destroy GPU model
  FOR_ALL(elem_gid, 0, num_elems, {
    size_t mat_id = elem_mat_id(elem_gid);
    auto run_loc = material(mat_id).strength_run_location;
    if (run_loc == RUN_LOCATION::device and elem_strength(elem_gid).model != nullptr) {  
        elem_strength(elem_gid).model->~StrengthParent();
    }
  });
  Kokkos::fence();

  // Destroy CPU model
  for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
    size_t mat_id = elem_mat_id.host(elem_gid);
    auto run_loc = material.host(mat_id).strength_run_location;
    if (run_loc == RUN_LOCATION::host and elem_strength.host(elem_gid).model != nullptr) {  
        elem_strength.host(elem_gid).model->~StrengthParent();
    }
  }

  // Free CPU memory
  for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
    size_t mat_id = elem_mat_id.host(elem_gid);
    auto run_loc = material.host(mat_id).strength_run_location;
    if (elem_strength.host(elem_gid).model != nullptr) {
      free_memory(elem_strength.host(elem_gid).model, run_loc);
    }
  }

  return;
} // end destroy_strength_model


void destroy_eos_model(
  DCArrayKokkos <eos_t> &elem_eos,
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t num_elems)
{
  /*
    All memory cleanup related to the user material model should be done in this fuction.
    Fierro calls `destroy_eos_model()` at the end of a simulation.
  */

  // Destroy GPU model
  FOR_ALL(elem_gid, 0, num_elems, {
    size_t mat_id = elem_mat_id(elem_gid);
    auto run_loc = material(mat_id).eos_run_location;
    if (run_loc == RUN_LOCATION::device and elem_eos(elem_gid).model != nullptr) {  
        elem_eos(elem_gid).model->~EOSParent();
    }
  });
  Kokkos::fence();

  // Destroy CPU model
  for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
    size_t mat_id = elem_mat_id.host(elem_gid);
    auto run_loc = material.host(mat_id).eos_run_location;
    if (run_loc == RUN_LOCATION::host and elem_eos.host(elem_gid).model != nullptr) {  
        elem_eos.host(elem_gid).model->~EOSParent();
    }
  }

  // Free CPU memory
  for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
    size_t mat_id = elem_mat_id.host(elem_gid);
    auto run_loc = material.host(mat_id).eos_run_location;
    if (elem_eos.host(elem_gid).model != nullptr) {
      free_memory(elem_eos.host(elem_gid).model, run_loc);
    }
  }

  return;
} // end destroy_eos_model

