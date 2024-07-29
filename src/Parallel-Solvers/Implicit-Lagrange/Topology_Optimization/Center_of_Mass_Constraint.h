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
 
#ifndef CENTER_OF_MASS_CONSTRAINT_TOPOPT_H
#define CENTER_OF_MASS_CONSTRAINT_TOPOPT_H

#include "matar.h"
#include "elements.h"
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"

#include "ROL_Types.hpp"
#include <ROL_TpetraMultiVector.hpp>
#include "ROL_Constraint.hpp"
#include "ROL_Elementwise_Reduce.hpp"
#include "FEA_Module_Inertial.h"

class CenterOfMassConstraint_TopOpt : public ROL::Constraint<real_t> {
  
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<>::node_type Node;
  typedef Tpetra::Map<LO, GO, Node> Map;
  typedef Tpetra::MultiVector<real_t, LO, GO, Node> MV;
  typedef ROL::Vector<real_t> V;
  typedef ROL::TpetraMultiVector<real_t,LO,GO,Node> ROL_MV;
  
  using traits = Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>;
  using array_layout    = typename traits::array_layout;
  using execution_space = typename traits::execution_space;
  using device_type     = typename traits::device_type;
  using memory_traits   = typename traits::memory_traits;
  using global_size_t = Tpetra::global_size_t;

  typedef Kokkos::View<real_t*, Kokkos::LayoutRight, device_type, memory_traits> values_array;
  typedef Kokkos::View<GO*, array_layout, device_type, memory_traits> global_indices_array;
  typedef Kokkos::View<LO*, array_layout, device_type, memory_traits> indices_array;
  
  //typedef Kokkos::DualView<real_t**, Kokkos::LayoutLeft, device_type>::t_dev vec_array;
  typedef MV::dual_view_type::t_dev vec_array;
  typedef MV::dual_view_type::t_host host_vec_array;
  typedef Kokkos::View<const real_t**, array_layout, HostSpace, memory_traits> const_host_vec_array;
  typedef MV::dual_view_type dual_vec_array;

private:

  FEA_Module_Inertial *FEM_;
  ROL::Ptr<ROL_MV> ROL_Element_Masses;
  ROL::Ptr<ROL_MV> ROL_Element_Moments;
  ROL::Ptr<ROL_MV> ROL_Gradients, ROL_Mass_Gradients;
  Teuchos::RCP<MV> constraint_gradients_distributed;
  Teuchos::RCP<MV> mass_gradients_distributed;
  real_t initial_center_of_mass;
  bool inequality_flag_;
  real_t constraint_value_, normalization_value;
  int constraint_component_;

  ROL::Ptr<const MV> getVector( const V& x ) {
    return dynamic_cast<const ROL_MV&>(x).getVector();
  }

  ROL::Ptr<MV> getVector( V& x ) {
    return dynamic_cast<ROL_MV&>(x).getVector();
  }

public:
  bool nodal_density_flag_;
  int last_comm_step, current_step, last_solve_step;
  std::string my_fea_module = "Inertial";

  CenterOfMassConstraint_TopOpt(FEA_Module *FEM, bool nodal_density_flag, int constraint_component, real_t constraint_value=0, bool inequality_flag=true) 
  {
    FEM_ = dynamic_cast<FEA_Module_Inertial*>(FEM);
    nodal_density_flag_ = nodal_density_flag;
    last_comm_step = last_solve_step = -1;
    current_step = 0;
    inequality_flag_ = inequality_flag;
    constraint_value_ = constraint_value;
    constraint_component_ = constraint_component;
    normalization_value = constraint_value_;
    if(constraint_value_==0) normalization_value = 1;
    int num_dim = FEM_->simparam->num_dims;
    ROL_Element_Masses = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Masses);
    if(constraint_component_ == 0)
    ROL_Element_Moments = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_x);
    if(constraint_component_ == 1)
    ROL_Element_Moments = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_y);
    if(constraint_component_ == 2)
    ROL_Element_Moments = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_z);

    const_host_vec_array design_densities = FEM_->design_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    
    //sum per element results across all MPI ranks
    ROL::Elementwise::ReductionSum<real_t> sumreduc;
    real_t initial_moment;
    real_t initial_mass;

    if(FEM_->mass_init) { initial_mass = FEM_->mass; }
    else{
      FEM_->compute_element_masses(design_densities,true);
      //sum per element results across all MPI ranks
      ROL::Elementwise::ReductionSum<real_t> sumreduc;
      FEM_->mass = initial_mass = ROL_Element_Masses->reduce(sumreduc);
      FEM_->mass_init = true;
    }
    
    if(FEM_->com_init[constraint_component_]) { initial_center_of_mass = FEM_->center_of_mass[constraint_component_]; }
    else{
      FEM_->compute_element_moments(design_densities,true, constraint_component_, false);

      //sum per element results across all MPI ranks
      ROL::Elementwise::ReductionSum<real_t> sumreduc;
      initial_moment = ROL_Element_Moments->reduce(sumreduc);
      FEM_->com_init[constraint_component_] = true;
      FEM_->center_of_mass[constraint_component_] = initial_center_of_mass = initial_moment/initial_mass;
    }

    //debug print
    if(FEM_->myrank==0){
      if(constraint_component_ == 0)
      std::cout << "INITIAL COM X: " << initial_center_of_mass << std::endl;
      if(constraint_component_ == 1)
      std::cout << "INITIAL COM Y: " << initial_center_of_mass << std::endl;
      if(constraint_component_ == 2)
      std::cout << "INITIAL COM Z: " << initial_center_of_mass << std::endl;
    }
  
    if(FEM_->mass_gradients_distributed.is_null()) 
      FEM_->mass_gradients_distributed = Teuchos::rcp(new MV(FEM_->map, 1));
    mass_gradients_distributed = FEM_->mass_gradients_distributed; 

    if(FEM_->center_of_mass_gradients_distributed.is_null()) 
      FEM_->center_of_mass_gradients_distributed = Teuchos::rcp(new MV(FEM_->map, 1));
    constraint_gradients_distributed = FEM_->center_of_mass_gradients_distributed;
  }

  /* --------------------------------------------------------------------------------------
   Update solver state variables to synchronize with the current design variable vector, z
  ----------------------------------------------------------------------------------------- */

  void update(const ROL::Vector<real_t> &z, ROL::UpdateType type, int iter = -1 ) {
    current_step++;
  }

  /* --------------------------------------------------------------------------------------
   Update constraint value (c) with the current design variable vector, z
  ----------------------------------------------------------------------------------------- */
  
  void value(ROL::Vector<real_t> &c, const ROL::Vector<real_t> &z, real_t &tol ) override {
    //std::cout << "Started constraint value on task " <<FEM_->myrank <<std::endl;
    ROL::Ptr<const MV> zp = getVector(z);
    ROL::Ptr<std::vector<real_t>> cp = dynamic_cast<ROL::StdVector<real_t>&>(c).getVector();
    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

    //communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    
    FEM_->compute_element_moments(design_densities,false,constraint_component_);
    FEM_->com_update[constraint_component_] = current_step;
    
    //sum per element results across all MPI ranks
    ROL::Elementwise::ReductionSum<real_t> sumreduc;
    real_t current_moment = ROL_Element_Moments->reduce(sumreduc);
    real_t current_com;

    //compute mass
    real_t current_mass;
    if(FEM_->mass_update == current_step&&0) { current_mass = FEM_->mass; }
    else{
      FEM_->compute_element_masses(design_densities,false);
      //sum per element results across all MPI ranks
      ROL::Elementwise::ReductionSum<real_t> sumreduc;
      FEM_->mass = current_mass = ROL_Element_Masses->reduce(sumreduc);
      FEM_->mass_update = current_step;
    }

    current_com = FEM_->center_of_mass[constraint_component_] = current_moment/current_mass;
    //debug print
    if(FEM_->myrank==0){
      if(constraint_component_ == 0)
      std::cout << "CURRENT COM X " << current_com << std::endl;
      if(constraint_component_ == 1)
      std::cout << "CURRENT COM Y " << current_com << std::endl;
      if(constraint_component_ == 2)
      std::cout << "CURRENT COM Z " << current_com << std::endl;
    }
    
    if(inequality_flag_){
      (*cp)[0] = current_com/normalization_value;
    }
    else{
      (*cp)[0] = (current_com - constraint_value_)/normalization_value;
    }

    //std::cout << "Ended constraint value on task " <<FEM_->myrank <<std::endl;
  }

  /* ----------------------------------------------------------------------------------------
    Update adjoint of the constraint jacobian (ajv), involves gradient vector
    and constraint adjoint (v) for the current design variable vector, z
  ------------------------------------------------------------------------------------------*/

  void applyAdjointJacobian(ROL::Vector<real_t> &ajv, const ROL::Vector<real_t> &v, const ROL::Vector<real_t> &x, real_t &tol) override {
    //std::cout << "Started constraint adjoint grad on task " <<FEM_->myrank << std::endl;
     //get Tpetra multivector pointer from the ROL vector
    ROL::Ptr<const MV> zp = getVector(x);
    ROL::Ptr<const std::vector<real_t>> vp = dynamic_cast<const ROL::StdVector<real_t>&>(v).getVector();
    ROL::Ptr<MV> ajvp = getVector(ajv);
    
    //ROL::Ptr<ROL_MV> ROL_Element_Volumes;

    //get local view of the data
    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    //host_vec_array constraint_gradients = constraint_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    host_vec_array constraint_gradients = ajvp->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    host_vec_array mass_gradients = mass_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    //host_vec_array dual_constraint_vector = vp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

    //communicate ghosts
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    
    int rnum_elem = FEM_->rnum_elem;

    //compute mass
    real_t current_mass;
    if(FEM_->mass_update == current_step&&0) { current_mass = FEM_->mass; }
    else{
      FEM_->compute_element_masses(design_densities,false);
      //sum per element results across all MPI ranks
      ROL::Elementwise::ReductionSum<real_t> sumreduc;
      FEM_->mass = current_mass = ROL_Element_Masses->reduce(sumreduc);
      FEM_->mass_update = current_step;
    }

    //compute center of mass
    real_t current_com;
    if(FEM_->com_update[constraint_component_] == current_step) current_com = FEM_->center_of_mass[constraint_component_];
    else{
    FEM_->compute_element_moments(design_densities,false,constraint_component_);
    FEM_->com_update[constraint_component_] = current_step;
    
    //sum per element results across all MPI ranks
    ROL::Elementwise::ReductionSum<real_t> sumreduc;
    real_t current_moment = ROL_Element_Moments->reduce(sumreduc);
    FEM_->center_of_mass[constraint_component_] = current_com = current_moment/current_mass;
    }

    FEM_->compute_moment_gradients(design_densities, constraint_gradients, constraint_component_);
    FEM_->compute_nodal_gradients(design_densities, mass_gradients);

      //debug print of gradient
      //std::ostream &out = std::cout;
      //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
      //if(FEM_->myrank==0)
      //*fos << "Gradient data :" << std::endl;
      //ajvp->describe(*fos,Teuchos::VERB_EXTREME);
      //*fos << std::endl;
      //std::fflush(stdout);
    for(int i = 0; i < FEM_->nlocal_nodes; i++){
      constraint_gradients(i,0) /= current_mass;
      constraint_gradients(i,0) -= mass_gradients(i,0)*current_com/current_mass;
      constraint_gradients(i,0) *= (*vp)[0]/normalization_value;
    }
    
    
    //std::cout << "Ended constraint adjoint grad on task " <<FEM_->myrank  << std::endl;
    //debug print
    //std::cout << "Constraint Gradient value " << std::endl;
  }

  /* ----------------------------------------------------------------------------------------
    Update the constraint jacobian (jv), involves gradient vector
    and design vector differential (v) for the current design variable vector, x
  ------------------------------------------------------------------------------------------*/

  void applyJacobian(ROL::Vector<real_t> &jv, const ROL::Vector<real_t> &v, const ROL::Vector<real_t> &x, real_t &tol) override {
    //std::cout << "Started constraint grad on task " <<FEM_->myrank  << std::endl;
    //get Tpetra multivector pointer from the ROL vector
    ROL::Ptr<const MV> zp = getVector(x);
    ROL::Ptr<std::vector<real_t>> jvp = dynamic_cast<ROL::StdVector<real_t>&>(jv).getVector();
    
    //ROL::Ptr<ROL_MV> ROL_Element_Volumes;

    //get local view of the data
    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    host_vec_array constraint_gradients = constraint_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    host_vec_array mass_gradients = mass_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

    //communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
    //communicate ghosts
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    
    int rnum_elem = FEM_->rnum_elem;

    //compute mass
    real_t current_mass;
    if(FEM_->mass_update == current_step&&0) { current_mass = FEM_->mass; }
    else{
      FEM_->compute_element_masses(design_densities,false);
      //sum per element results across all MPI ranks
      ROL::Elementwise::ReductionSum<real_t> sumreduc;
      FEM_->mass = current_mass = ROL_Element_Masses->reduce(sumreduc);
      FEM_->mass_update = current_step;
    }

    //compute center of mass
    real_t current_com;
    if(FEM_->com_update[constraint_component_] == current_step) current_com = FEM_->center_of_mass[constraint_component_];
    else{
    FEM_->compute_element_moments(design_densities,false,constraint_component_);
    FEM_->com_update[constraint_component_] = current_step;
    
    //sum per element results across all MPI ranks
    ROL::Elementwise::ReductionSum<real_t> sumreduc;
    real_t current_moment = ROL_Element_Moments->reduce(sumreduc);
    FEM_->center_of_mass[constraint_component_] = current_com = current_moment/current_mass;
    }

    FEM_->compute_moment_gradients(design_densities, constraint_gradients, constraint_component_);
    FEM_->compute_nodal_gradients(design_densities, mass_gradients);

    for(int i = 0; i < FEM_->nlocal_nodes; i++){
      constraint_gradients(i,0) /= current_mass;
      constraint_gradients(i,0) -= mass_gradients(i,0)*current_com/current_mass;
    }

    ROL_Gradients = ROL::makePtr<ROL_MV>(constraint_gradients_distributed);
    real_t gradient_dot_v = ROL_Gradients->dot(v);
    //debug print
    //std::cout << "Constraint Gradient value " << gradient_dot_v << std::endl;

    (*jvp)[0] = gradient_dot_v/normalization_value;
    //std::cout << "Ended constraint grad on task " <<FEM_->myrank  << std::endl;
  }
  
  /* ----------------------------------------------------------------------------------------
    Update adjoint of the constraint Hessian vector product (ahuv), product of hessian and
    differential design vector (v), and constraint adjoint (u) for the 
    current design variable vector, z
  ------------------------------------------------------------------------------------------*/

  void applyAdjointHessian(ROL::Vector<real_t> &ahuv, const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &v, const ROL::Vector<real_t> &z, real_t &tol) {
    
    //get Tpetra multivector pointer from the ROL vector
    ROL::Ptr<const MV> zp = getVector(z);
    ROL::Ptr<const std::vector<real_t>> up = dynamic_cast<const ROL::StdVector<real_t>&>(u).getVector();
    ROL::Ptr<const MV> vp = getVector(v);
    
    //ROL::Ptr<ROL_MV> ROL_Element_Volumes;

    //get local view of the data
    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    const_host_vec_array v_view = vp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    host_vec_array moment_gradients = constraint_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    host_vec_array mass_gradients = mass_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

    //communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
    //communicate ghosts
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    
    int rnum_elem = FEM_->rnum_elem;

    //compute mass
    real_t current_mass;
    if(FEM_->mass_update == current_step&&0) { current_mass = FEM_->mass; }
    else{
      FEM_->compute_element_masses(design_densities,false);
      //sum per element results across all MPI ranks
      ROL::Elementwise::ReductionSum<real_t> sumreduc;
      FEM_->mass = current_mass = ROL_Element_Masses->reduce(sumreduc);
      FEM_->mass_update = current_step;
    }

    //compute center of mass
    real_t current_com;
    if(FEM_->com_update[constraint_component_] == current_step) current_com = FEM_->center_of_mass[constraint_component_];
    else{
    FEM_->compute_element_moments(design_densities,false,constraint_component_);
    FEM_->com_update[constraint_component_] = current_step;
    
    //sum per element results across all MPI ranks
    ROL::Elementwise::ReductionSum<real_t> sumreduc;
    real_t current_moment = ROL_Element_Moments->reduce(sumreduc);
    FEM_->center_of_mass[constraint_component_] = current_com = current_moment/current_mass;
    }

    // Unwrap hv
    ROL::Ptr<MV> ahuvp = getVector(ahuv);
  
    host_vec_array constraint_adjoint_hess = ahuvp->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    real_t matrix_product;
    FEM_->compute_moment_gradients(design_densities, moment_gradients, constraint_component_);
    FEM_->compute_nodal_gradients(design_densities, mass_gradients);
    
    ROL_Gradients = ROL::makePtr<ROL_MV>(constraint_gradients_distributed);
    ROL_Mass_Gradients = ROL::makePtr<ROL_MV>(mass_gradients_distributed);
    real_t dot_product_moment = ROL_Gradients->dot(v);
    real_t dot_product_mass = ROL_Mass_Gradients->dot(v);
    for(int i = 0; i < FEM_->nlocal_nodes; i++){
        constraint_adjoint_hess(i,0) = -dot_product_moment*mass_gradients(i,0)*(*up)[0]/current_mass/current_mass-
                      moment_gradients(i,0)*(*up)[0]*dot_product_mass/current_mass/current_mass+
                      2*current_com*dot_product_mass*mass_gradients(i,0)*(*up)[0]/current_mass/current_mass;
        constraint_adjoint_hess(i,0) /= normalization_value;
    }
    
  }
};

#endif // end header guard
