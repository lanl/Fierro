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
 
#ifndef MOMENT_OF_INERTIA_CONSTRAINT_TOPOPT_H
#define MOMENT_OF_INERTIA_CONSTRAINT_TOPOPT_H

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

class MomentOfInertiaConstraint_TopOpt : public ROL::Constraint<real_t> {
  
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
  ROL::Ptr<ROL_MV> ROL_Element_Moments1;
  ROL::Ptr<ROL_MV> ROL_Element_Moments2;
  ROL::Ptr<ROL_MV> ROL_Element_Moments_of_Inertia;
  ROL::Ptr<ROL_MV> ROL_Element_Moments_of_Inertia_xx, ROL_Element_Moments_of_Inertia_yy, ROL_Element_Moments_of_Inertia_zz;
  ROL::Ptr<ROL_MV> ROL_Gradients;
  Teuchos::RCP<MV> constraint_gradients_distributed;
  //Teuchos::RCP<MV> center_of_mass_gradients_distributed;
  //Teuchos::RCP<MV> mass_gradients_distributed;
  real_t initial_moment_of_inertia, initial_Mxx, initial_Myy, initial_Mzz;
  bool inequality_flag_, use_initial_coords_;
  real_t constraint_value_;
  real_t current_mass;
  real_t current_center_of_mass[3];
  int inertia_component_;
  int com1, com2;

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

  MomentOfInertiaConstraint_TopOpt(FEA_Module *FEM, bool nodal_density_flag, int inertia_component, real_t constraint_value = 0, bool inequality_flag=true, bool use_initial_coords=false) 
  {
    
    FEM_ = dynamic_cast<FEA_Module_Inertial*>(FEM);
    if(inertia_component > 5) { std::cout << "SPECIFIED MOMENT OF INERTIA CONSTRAINT COMPONENT CANNOT EXCEED 5" << std::endl; }
    nodal_density_flag_ = nodal_density_flag;
    use_initial_coords_ = use_initial_coords;
    last_comm_step = last_solve_step = -1;
    current_step = 0;
    inequality_flag_ = inequality_flag;
    constraint_value_ = constraint_value;
    inertia_component_ = inertia_component;
    int num_dim = FEM_->simparam->num_dims;
    ROL_Element_Masses = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Masses);

    if(inertia_component_ == 0)
    ROL_Element_Moments_of_Inertia = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_of_Inertia_xx);
    if(inertia_component_ == 1)
    ROL_Element_Moments_of_Inertia = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_of_Inertia_yy);
    if(inertia_component_ == 2)
    ROL_Element_Moments_of_Inertia = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_of_Inertia_zz);
    if(inertia_component_ == 3)
    ROL_Element_Moments_of_Inertia = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_of_Inertia_xy);
    if(inertia_component_ == 4)
    ROL_Element_Moments_of_Inertia = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_of_Inertia_xz);
    if(inertia_component_ == 5)
    ROL_Element_Moments_of_Inertia = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_of_Inertia_yz);

    ROL_Element_Moments_of_Inertia_xx = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_of_Inertia_xx);
    ROL_Element_Moments_of_Inertia_yy = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_of_Inertia_yy);
    ROL_Element_Moments_of_Inertia_zz = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_of_Inertia_zz);

    const_host_vec_array design_densities = FEM_->design_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    
    //compute initial mass
    real_t initial_mass;

    if(FEM_->mass_init) { initial_mass = FEM_->mass; }
    else{
      FEM_->compute_element_masses(design_densities,true,use_initial_coords_);
      //sum per element results across all MPI ranks
      ROL::Elementwise::ReductionSum<real_t> sumreduc;
      FEM_->mass = initial_mass = ROL_Element_Masses->reduce(sumreduc);
      FEM_->mass_init = true;
    }

    //compute initial center of mass
    real_t initial_center_of_mass[3];
    real_t initial_moment;
    
    //compute the components of the center of mass that are required
    if(inertia_component_==0) {
      com1 = 1;
      com2 = 2;
    }
    else if(inertia_component_==1) {
      com1 = 0;
      com2 = 2;
    }
    else if(inertia_component_==2) {
      com1 = 0;
      com2 = 1;
    }
    else if(inertia_component_==3) {
      com1 = 0;
      com2 = 1;
    }
    else if(inertia_component_==4) {
      com1 = 0;
      com2 = 2;
    }
    else if(inertia_component_==5) {
      com1 = 1;
      com2 = 2;
    }
    
    
    if(com1 == 0)
    ROL_Element_Moments1 = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_x);
    if(com1 == 1)
    ROL_Element_Moments1 = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_y);
    if(com1 == 2)
    ROL_Element_Moments1 = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_z);

    if(FEM_->com_init[com1]) { initial_center_of_mass[com1] = FEM_->center_of_mass[com1]; }
    else{
      FEM_->compute_element_moments(design_densities,true, com1, use_initial_coords_);

      //sum per element results across all MPI ranks
      ROL::Elementwise::ReductionSum<real_t> sumreduc;
      initial_moment = ROL_Element_Moments1->reduce(sumreduc);
      FEM_->com_init[com1] = true;
      FEM_->center_of_mass[com1] = initial_center_of_mass[com1] = initial_moment/initial_mass;
    }
    
    if(com2 == 0)
    ROL_Element_Moments2 = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_x);
    if(com2 == 1)
    ROL_Element_Moments2 = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_y);
    if(com2 == 2)
    ROL_Element_Moments2 = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Moments_z);

    if(FEM_->com_init[com2]) { initial_center_of_mass[com2] = FEM_->center_of_mass[com2]; }
    else{
      FEM_->compute_element_moments(design_densities,true, com2, use_initial_coords_);

      //sum per element results across all MPI ranks
      ROL::Elementwise::ReductionSum<real_t> sumreduc;
      initial_moment = ROL_Element_Moments2->reduce(sumreduc);
      FEM_->com_init[com2] = true;
      FEM_->center_of_mass[com2] = initial_center_of_mass[com2] = initial_moment/initial_mass;
    }
    
    
    FEM_->compute_element_moments_of_inertia(design_densities,true, inertia_component_, use_initial_coords_);
    
    //sum per element results across all MPI ranks
    ROL::Elementwise::ReductionSum<real_t> sumreduc;
    initial_moment_of_inertia = ROL_Element_Moments_of_Inertia->reduce(sumreduc);

    if(inertia_component_ == 3){
      FEM_->compute_element_moments_of_inertia(design_densities,true, 0, use_initial_coords_);
      initial_Mxx = ROL_Element_Moments_of_Inertia_xx->reduce(sumreduc);
    }
    if(inertia_component_ == 4){
      FEM_->compute_element_moments_of_inertia(design_densities,true, 0, use_initial_coords_);
      initial_Mxx = ROL_Element_Moments_of_Inertia_xx->reduce(sumreduc);
    }
    if(inertia_component_ == 5){
      FEM_->compute_element_moments_of_inertia(design_densities,true, 2, use_initial_coords_);
      initial_Mzz = ROL_Element_Moments_of_Inertia_zz->reduce(sumreduc);
    }
    
    //prevent numerical issues
    if(inertia_component_ < 3 && constraint_value_==0) { std::cout << "ERROR PRINCIPAL MOMENTS OF INERTIA CANNOT BE CONSTRAINED TO ZERO " << std::endl; }

    //debug print
    if(FEM_->myrank==0){
      if(inertia_component_ == 0)
        std::cout << "INITIAL MOMENT OF INERTIA XX: " << initial_moment_of_inertia << std::endl;
      if(inertia_component_ == 1)
        std::cout << "INITIAL MOMENT OF INERTIA YY: " << initial_moment_of_inertia << std::endl;
      if(inertia_component_ == 2)
        std::cout << "INITIAL MOMENT OF INERTIA ZZ: " << initial_moment_of_inertia << std::endl;
      if(inertia_component_ == 3){
        std::cout << "INITIAL MOMENT OF INERTIA XY: " << initial_moment_of_inertia << std::endl;
        std::cout << "INITIAL MOMENT OF INERTIA XX: " << initial_Mxx << std::endl;
      }
      if(inertia_component_ == 4){
        std::cout << "INITIAL MOMENT OF INERTIA XZ: " << initial_moment_of_inertia << std::endl;
        std::cout << "INITIAL MOMENT OF INERTIA XX: " << initial_Mxx << std::endl;
      }
      if(inertia_component_ == 5){
        std::cout << "INITIAL MOMENT OF INERTIA YZ: " << initial_moment_of_inertia << std::endl;
        std::cout << "INITIAL MOMENT OF INERTIA ZZ: " << initial_Mzz << std::endl;
      }
    }
    constraint_gradients_distributed = Teuchos::rcp(new MV(FEM_->map, 1));
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
    
    real_t current_mass;
    real_t current_center_of_mass[3];
    update_com_and_mass(design_densities, current_mass, current_center_of_mass);
    
    FEM_->compute_element_moments_of_inertia(design_densities,false,inertia_component_, use_initial_coords_);
    
    //sum per element results across all MPI ranks
    ROL::Elementwise::ReductionSum<real_t> sumreduc;
    real_t current_moment_of_inertia = ROL_Element_Moments_of_Inertia->reduce(sumreduc);
    //debug print
    if(FEM_->myrank==0){
      if(inertia_component_ == 0)
        std::cout << "MOMENT OF INERTIA XX RATIO: " << current_moment_of_inertia/initial_moment_of_inertia << std::endl;
      if(inertia_component_ == 1)
        std::cout << "MOMENT OF INERTIA YY RATIO: " << current_moment_of_inertia/initial_moment_of_inertia << std::endl;
      if(inertia_component_ == 2)
        std::cout << "MOMENT OF INERTIA ZZ RATIO: " << current_moment_of_inertia/initial_moment_of_inertia << std::endl;
      if(inertia_component_ == 3)
        std::cout << "MOMENT OF INERTIA XY RATIO (NORMALIZED BY M_XX): " << current_moment_of_inertia/initial_Mxx << std::endl;
      if(inertia_component_ == 4)
        std::cout << "MOMENT OF INERTIA XZ RATIO (NORMALIZED BY M_XX): " << current_moment_of_inertia/initial_Mxx << std::endl;
      if(inertia_component_ == 5)
        std::cout << "MOMENT OF INERTIA YZ RATIO (NORMALIZED BY M_ZZ): " << current_moment_of_inertia/initial_Mzz << std::endl;
    }
    
    if(inequality_flag_){
      if(inertia_component_ < 3){
        (*cp)[0] = current_moment_of_inertia/initial_moment_of_inertia;
      }
      else if(inertia_component_ == 3||inertia_component_== 4){
        (*cp)[0] = current_moment_of_inertia/initial_Mxx;
      }
      else{
        (*cp)[0] = current_moment_of_inertia/initial_Mzz;
      }
    }
    else{
      if(inertia_component_ < 3){
        (*cp)[0] = current_moment_of_inertia/initial_moment_of_inertia - constraint_value_;
      }
      else if(inertia_component_ == 3||inertia_component_== 4){
        (*cp)[0] = current_moment_of_inertia/initial_Mxx - constraint_value_;
      }
      else{
        (*cp)[0] = current_moment_of_inertia/initial_Mzz - constraint_value_;
      }
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
    //host_vec_array mass_gradients = mass_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    //host_vec_array center_of_mass_gradients = center_of_mass_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    //host_vec_array dual_constraint_vector = vp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

    //communicate ghosts
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    
    int rnum_elem = FEM_->rnum_elem;
    update_com_and_mass(design_densities, current_mass, current_center_of_mass);

    /*
    //compute mass and com derivates needed by chain rule
    FEM_->compute_nodal_gradients(design_densities, mass_gradients);

      //debug print of gradient
      //std::ostream &out = std::cout;
      //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
      //if(FEM_->myrank==0)
      //*fos << "Gradient data :" << std::endl;
      //ajvp->describe(*fos,Teuchos::VERB_EXTREME);
      //*fos << std::endl;
      //std::fflush(stdout);
    
    int com1, com2;
    
    //compute the components of the center of mass that are required
    if(inertia_component_==0) {
      com1 = 1;
      com2 = 2;
    }
    else if(inertia_component_==1) {
      com1 = 0;
      com2 = 2;
    }
    else if(inertia_component_==2) {
      com1 = 0;
      com2 = 1;
    }
    else if(inertia_component_==3) {
      com1 = 0;
      com2 = 1;
    }
    else if(inertia_component_==4) {
      com1 = 0;
      com2 = 2;
    }
    else if(inertia_component_==5) {
      com1 = 1;
      com2 = 2;
    }
    
    FEM_->compute_moment_gradients(design_densities, center_of_mass_gradients, com1);
    FEM_->compute_moment_gradients(design_densities, center_of_mass_gradients, com2);

    for(int i = 0; i < FEM_->nlocal_nodes; i++){
      center_of_mass_gradients(i,com1) /= current_mass;
      center_of_mass_gradients(i,com1) -= mass_gradients(i)*current_center_of_mass[com1]/current_mass;
    }

    for(int i = 0; i < FEM_->nlocal_nodes; i++){
      center_of_mass_gradients(i,com2) /= current_mass;
      center_of_mass_gradients(i,com2) -= mass_gradients(i)*current_center_of_mass[com2]/current_mass;
    }
    */
    
    FEM_->compute_moment_of_inertia_gradients(design_densities, constraint_gradients, inertia_component_, use_initial_coords_);
      //debug print of gradient
      //std::ostream &out = std::cout;
      //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
      //if(FEM_->myrank==0)
      //*fos << "Gradient data :" << std::endl;
      //ajvp->describe(*fos,Teuchos::VERB_EXTREME);
      //*fos << std::endl;
      //std::fflush(stdout);

    //for(int i = 0; i < FEM_->nlocal_nodes; i++){
      //constraint_gradients(i,0) -= 2*current_center_of_mass[com1]*current_mass*center_of_mass_gradients(i, com1);
      //constraint_gradients(i,0) -= 2*current_center_of_mass[com2]*current_mass*center_of_mass_gradients(i, com2);

      //constraint_gradients(i,0) += 2*current_center_of_mass[com1]*current_mass*center_of_mass_gradients(i, com1);
      //constraint_gradients(i,0) += 2*current_center_of_mass[com2]*current_mass*center_of_mass_gradients(i, com2);
    //}

    
      if(inertia_component_ < 3){
        for(int i = 0; i < FEM_->nlocal_nodes; i++)
          constraint_gradients(i,0) *= (*vp)[0]/initial_moment_of_inertia;
      }
      else if(inertia_component_== 3||inertia_component_== 4){
        for(int i = 0; i < FEM_->nlocal_nodes; i++)
          constraint_gradients(i,0) *= (*vp)[0]/initial_Mxx;
      }
      else{
        for(int i = 0; i < FEM_->nlocal_nodes; i++)
          constraint_gradients(i,0) *= (*vp)[0]/initial_Mzz;
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
    //host_vec_array mass_gradients = mass_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    //host_vec_array center_of_mass_gradients = center_of_mass_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

    //communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
    //communicate ghosts
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    
    int rnum_elem = FEM_->rnum_elem;
    real_t current_mass;
    real_t current_center_of_mass[3];
    update_com_and_mass(design_densities, current_mass, current_center_of_mass);
    
    FEM_->compute_moment_of_inertia_gradients(design_densities, constraint_gradients, inertia_component_, use_initial_coords_);
    
      if(inertia_component_ < 3){
        for(int i = 0; i < FEM_->nlocal_nodes; i++)
          constraint_gradients(i,0) /= initial_moment_of_inertia;
      }
      else if (inertia_component_== 3||inertia_component_== 4){
        for(int i = 0; i < FEM_->nlocal_nodes; i++)
          constraint_gradients(i,0) /= initial_Mxx;
      }
      else{
        for(int i = 0; i < FEM_->nlocal_nodes; i++)
          constraint_gradients(i,0) /= initial_Mzz;
      }
    

    ROL_Gradients = ROL::makePtr<ROL_MV>(constraint_gradients_distributed);
    real_t gradient_dot_v = ROL_Gradients->dot(v);
    //debug print
    //std::cout << "Constraint Gradient value " << gradient_dot_v << std::endl;

    (*jvp)[0] = gradient_dot_v;
    //std::cout << "Ended constraint grad on task " <<FEM_->myrank  << std::endl;
  }

  void update_com_and_mass(const_host_vec_array design_densities, real_t &current_mass, real_t *current_center_of_mass){

    //compute mass
    if(FEM_->mass_update == current_step&&0) { current_mass = FEM_->mass; }
    else{
      FEM_->compute_element_masses(design_densities,false,use_initial_coords_);
      //sum per element results across all MPI ranks
      ROL::Elementwise::ReductionSum<real_t> sumreduc;
      FEM_->mass = current_mass = ROL_Element_Masses->reduce(sumreduc);
      FEM_->mass_update = current_step;
    }

    //compute initial center of mass
    real_t current_moment;
    
    if(FEM_->com_update[com1] == current_step&&0) { current_center_of_mass[com1] = FEM_->center_of_mass[com1]; }
    else{
      FEM_->compute_element_moments(design_densities,false, com1,use_initial_coords_);
      //sum per element results across all MPI ranks
      ROL::Elementwise::ReductionSum<real_t> sumreduc;
      current_moment = ROL_Element_Moments1->reduce(sumreduc);
      FEM_->com_update[com1] = current_step;
      FEM_->center_of_mass[com1] = current_center_of_mass[com1] = current_moment/current_mass;
    }

    if(FEM_->com_update[com2] = current_step&&0) { current_center_of_mass[com2] = FEM_->center_of_mass[com2]; }
    else{
      FEM_->compute_element_moments(design_densities,false, com2,use_initial_coords_);
      //sum per element results across all MPI ranks
      ROL::Elementwise::ReductionSum<real_t> sumreduc;
      current_moment = ROL_Element_Moments2->reduce(sumreduc);
      FEM_->com_update[com2] = current_step;
      FEM_->center_of_mass[com2] = current_center_of_mass[com2] = current_moment/current_mass;
    }
  }
  
  /* ----------------------------------------------------------------------------------------
    Update adjoint of the constraint Hessian vector product (ahuv), product of hessian and
    differential design vector (v), and constraint adjoint (u) for the 
    current design variable vector, z
  ------------------------------------------------------------------------------------------*/

  void applyAdjointHessian(ROL::Vector<real_t> &ahuv, const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &v, const ROL::Vector<real_t> &z, real_t &tol) {
    
    // Unwrap hv
    ROL::Ptr<MV> ahuvp = getVector(ahuv);

    ahuvp->putScalar(0);
    
  }
  /*

  void hessVec_21( ROL::Vector<real_t> &hv, const ROL::Vector<real_t> &v, 
                   const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
                     
    // Unwrap g
    ROL::Ptr<MV> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const MV> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const MV> up = getVector(u);
    ROL::Ptr<const MV> zp = getVector(z);
 
    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      std::MV<real_t> U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      std::MV<real_t> V;
      V.assign(vp->begin(),vp->end());
      FEM_->set_boundary_conditions(V);
      FEM_->apply_adjoint_jacobian(*hvp,U,*zp,V);
      for (size_t i=0; i<hvp->size(); i++) {
        (*hvp)[i] *= 2.0;
      }
    }
    
  }

  void hessVec_22( ROL::Vector<real_t> &hv, const ROL::Vector<real_t> &v, 
                   const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
                     
    ROL::Ptr<MV> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const MV> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const MV> up = getVector(u);
    ROL::Ptr<const MV> zp = getVector(z);
    
    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      MV U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      MV V;
      V.assign(vp->begin(),vp->end());
      FEM_->set_boundary_conditions(V);
      FEM_->apply_adjoint_jacobian(*hvp,U,*zp,*vp,U);
    }
    
  }
  */
};

#endif // end header guard
