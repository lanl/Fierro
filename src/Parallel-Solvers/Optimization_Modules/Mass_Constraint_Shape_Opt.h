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
 
#ifndef MASS_CONSTRAINT_SHAPEOPT_H
#define MASS_CONSTRAINT_SHAPEOPT_H

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

class MassConstraint_ShapeOpt : public ROL::Constraint<real_t> {
  
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
  ROL::Ptr<ROL_MV> ROL_Gradients;
  Teuchos::RCP<MV> constraint_gradients_distributed;
  real_t initial_mass;
  real_t current_mass;
  bool inequality_flag_, use_initial_coords_;
  real_t constraint_value_;

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

  MassConstraint_ShapeOpt(FEA_Module *FEM, real_t constraint_value=0, bool inequality_flag=true, bool use_initial_coords=false) 
  {
    FEM_ = dynamic_cast<FEA_Module_Inertial*>(FEM);
    use_initial_coords_ = use_initial_coords;
    last_comm_step = last_solve_step = -1;
    current_step = 0;
    inequality_flag_ = inequality_flag;
    constraint_value_ = constraint_value;
    ROL_Element_Masses = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Masses);
    const_host_vec_array design_coords = FEM_->design_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    
    FEM_->compute_element_masses(design_coords,true,use_initial_coords_);
    FEM_->mass_init = true;
    
    //sum per element results across all MPI ranks
    ROL::Elementwise::ReductionSum<real_t> sumreduc;
    FEM_->mass = initial_mass = ROL_Element_Masses->reduce(sumreduc);
    //debug print
    if(FEM_->myrank==0)
      std::cout << "INITIAL MASS: " << initial_mass << std::endl;
    if(FEM_->mass_gradients_distributed.is_null())
      FEM_->mass_gradients_distributed = Teuchos::rcp(new MV(FEM_->map, 1));
    constraint_gradients_distributed = FEM_->mass_gradients_distributed;
  }

  /* --------------------------------------------------------------------------------------
   Update solver state variables to synchronize with the current design variable vector, z
  ----------------------------------------------------------------------------------------- */

  void update(const ROL::Vector<real_t> &z, ROL::UpdateType type, int iter = -1 ) {
    current_step++;

    ROL::Ptr<const MV> zp = getVector(z);
    const_host_vec_array design_coords = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

    if (type == ROL::UpdateType::Initial)  {

      //initial design density data was already communicated for ghost nodes in init_design()
    }
    else if (type == ROL::UpdateType::Accept) {
      // z was accepted as the new iterate
    }
    else if (type == ROL::UpdateType::Revert) {
      // z has been rejected as the new iterate
      // Revert to cached value
      FEM_->comm_variables(zp);
    }
    else if (type == ROL::UpdateType::Trial) {
      // This is a new value of x
      FEM_->comm_variables(zp);
    }
    else { // ROL::UpdateType::Temp
      // This is a new value of x used for,
      // e.g., finite-difference checks
      FEM_->comm_variables(zp);
    }
  }

  /* --------------------------------------------------------------------------------------
   Update constraint value (c) with the current design variable vector, z
  ----------------------------------------------------------------------------------------- */

  void value(ROL::Vector<real_t> &c, const ROL::Vector<real_t> &z, real_t &tol ) override {
    //std::cout << "Started constraint value on task " <<FEM_->myrank <<std::endl;
    ROL::Ptr<const MV> zp = getVector(z);
    ROL::Ptr<std::vector<real_t>> cp = dynamic_cast<ROL::StdVector<real_t>&>(c).getVector();
    const_host_vec_array design_coords = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

    //communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
    /*
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    */
    FEM_->compute_element_masses(design_coords,false,use_initial_coords_);
    
    //sum per element results across all MPI ranks
    ROL::Elementwise::ReductionSum<real_t> sumreduc;
    current_mass = ROL_Element_Masses->reduce(sumreduc);
    FEM_->mass = current_mass;
    FEM_->mass_update = current_step;
    
    //debug print
    if(FEM_->myrank==0)
    std::cout << "SYSTEM MASS RATIO: " << current_mass/initial_mass << std::endl;
    
    if(inequality_flag_)
      (*cp)[0] = current_mass/initial_mass;
    else
      (*cp)[0] = current_mass/initial_mass - constraint_value_;

    //std::cout << "Ended constraint value on task " <<FEM_->myrank <<std::endl;
  }

  /* ----------------------------------------------------------------------------------------
    Update adjoint of the constraint jacobian (ajv), involves gradient vector
    and constraint adjoint (v) for the current design variable vector, z
  ------------------------------------------------------------------------------------------*/

  void applyAdjointJacobian(ROL::Vector<real_t> &ajv, const ROL::Vector<real_t> &v, const ROL::Vector<real_t> &z, real_t &tol) override {
    //std::cout << "Started constraint adjoint grad on task " <<FEM_->myrank << std::endl;
     //get Tpetra multivector pointer from the ROL vector
    ROL::Ptr<const MV> zp = getVector(z);
    ROL::Ptr<const std::vector<real_t>> vp = dynamic_cast<const ROL::StdVector<real_t>&>(v).getVector();
    ROL::Ptr<MV> ajvp = getVector(ajv);
    int num_dim = FEM_->num_dim;
    
    //ROL::Ptr<ROL_MV> ROL_Element_Volumes;

    //get local view of the data
    const_host_vec_array design_coords = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    //host_vec_array constraint_gradients = constraint_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    host_vec_array constraint_gradients = ajvp->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    //host_vec_array dual_constraint_vector = vp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

    //communicate ghosts
    /*
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    */
    int rnum_elem = FEM_->rnum_elem;

    FEM_->compute_nodal_gradients(design_coords, constraint_gradients, use_initial_coords_);
      //debug print of gradient
      //std::ostream &out = std::cout;
      //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
      //if(FEM_->myrank==0)
      //*fos << "Gradient data :" << std::endl;
      //ajvp->describe(*fos,Teuchos::VERB_EXTREME);
      //*fos << std::endl;
      //std::fflush(stdout);
    for(int i = 0; i < FEM_->nlocal_nodes; i++){
      constraint_gradients(i,0) *= (*vp)[0]/initial_mass;
      constraint_gradients(i,1) *= (*vp)[0]/initial_mass;
      if(num_dim==3){
        constraint_gradients(i,2) *= (*vp)[0]/initial_mass;
      }
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
    int num_dim = FEM_->num_dim;
    
    //ROL::Ptr<ROL_MV> ROL_Element_Volumes;

    //get local view of the data
    const_host_vec_array design_coords = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    host_vec_array constraint_gradients = constraint_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

    //communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
    //communicate ghosts
    /*
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    */

    FEM_->compute_nodal_gradients(design_coords, constraint_gradients);
    for(int i = 0; i < FEM_->nlocal_nodes; i++){
      constraint_gradients(i,0) /= initial_mass;
      constraint_gradients(i,1) /= initial_mass;
      if(num_dim==3){
        constraint_gradients(i,2) /= initial_mass;
      }
    }
    

    ROL_Gradients = ROL::makePtr<ROL_MV>(constraint_gradients_distributed);
    real_t gradient_dot_v = ROL_Gradients->dot(v);
    //debug print
    //std::cout << "Constraint Gradient value " << gradient_dot_v << std::endl;

    (*jvp)[0] = gradient_dot_v;
    //std::cout << "Ended constraint grad on task " <<FEM_->myrank  << std::endl;
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

};

#endif // end header guard
