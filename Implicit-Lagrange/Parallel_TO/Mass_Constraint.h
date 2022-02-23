#ifndef MASS_CONSTRAINT_TOPOPT_H
#define MASS_CONSTRAINT_TOPOPT_H

#include "matar.h"
#include "element_types/elements.h"
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Kokkos_View.hpp>
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"

#include "ROL_Types.hpp"
#include <ROL_TpetraMultiVector.hpp>
#include "ROL_Constraint.hpp"
#include "ROL_Elementwise_Reduce.hpp"
#include "Parallel_Nonlinear_Solver.h"

class MassConstraint_TopOpt : public ROL::Constraint<real_t> {
  
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

  Parallel_Nonlinear_Solver *FEM_;
  ROL::Ptr<ROL_MV> ROL_Element_Masses;
  ROL::Ptr<ROL_MV> ROL_Gradients;
  Teuchos::RCP<MV> constraint_gradients_distributed;
  real_t initial_mass;
  bool inequality_flag_;
  real_t constraint_value_;

  ROL::Ptr<const MV> getVector( const V& x ) {
    return dynamic_cast<const ROL_MV&>(x).getVector();
  }

  ROL::Ptr<MV> getVector( V& x ) {
    return dynamic_cast<ROL_MV&>(x).getVector();
  }

public:
  bool nodal_density_flag_;
  size_t last_comm_step, current_step, last_solve_step;

  MassConstraint_TopOpt(Parallel_Nonlinear_Solver *FEM, bool nodal_density_flag, bool inequality_flag=true, real_t constraint_value=0) 
    : FEM_(FEM) {

    nodal_density_flag_ = nodal_density_flag;
    last_comm_step = last_solve_step = -1;
    current_step = 0;
    inequality_flag_ = inequality_flag;
    constraint_value_ = constraint_value;
    ROL_Element_Masses = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Masses);
    const_host_vec_array design_densities = FEM_->design_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    
    FEM_->compute_element_masses(design_densities,true);
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

  void update(const ROL::Vector<real_t> &z, ROL::UpdateType type, int iter = -1 ) {
    current_step++;

    ROL::Ptr<const MV> zp = getVector(z);
    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

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

  void value(ROL::Vector<real_t> &c, const ROL::Vector<real_t> &z, real_t &tol ) override {
    //std::cout << "Started constraint value on task " <<FEM_->myrank <<std::endl;
    ROL::Ptr<const MV> zp = getVector(z);
    ROL::Ptr<std::vector<real_t>> cp = dynamic_cast<ROL::StdVector<real_t>&>(c).getVector();
    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

    //communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
    /*
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    */
    FEM_->compute_element_masses(design_densities,false);
    
    //sum per element results across all MPI ranks
    ROL::Elementwise::ReductionSum<real_t> sumreduc;
    real_t current_mass = ROL_Element_Masses->reduce(sumreduc);
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
    //host_vec_array dual_constraint_vector = vp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

    //communicate ghosts
    /*
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    */
    int rnum_elem = FEM_->rnum_elem;

    if(nodal_density_flag_){
      FEM_->compute_nodal_gradients(design_densities, constraint_gradients);
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
      }
    }
    else{
      //update per element volumes
      FEM_->compute_element_volumes();
      //ROL_Element_Volumes = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Volumes);
      //local view of element volumes
      const_host_vec_array element_volumes = FEM_->Global_Element_Volumes->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
      for(int ig = 0; ig < rnum_elem; ig++)
        constraint_gradients(ig,0) = element_volumes(ig,0)*(*vp)[0]/initial_mass;
    }
    
    //std::cout << "Ended constraint adjoint grad on task " <<FEM_->myrank  << std::endl;
    //debug print
    //std::cout << "Constraint Gradient value " << std::endl;
  }
  
  void applyJacobian(ROL::Vector<real_t> &jv, const ROL::Vector<real_t> &v, const ROL::Vector<real_t> &x, real_t &tol) override {
    //std::cout << "Started constraint grad on task " <<FEM_->myrank  << std::endl;
     //get Tpetra multivector pointer from the ROL vector
    ROL::Ptr<const MV> zp = getVector(x);
    ROL::Ptr<std::vector<real_t>> jvp = dynamic_cast<ROL::StdVector<real_t>&>(jv).getVector();
    
    //ROL::Ptr<ROL_MV> ROL_Element_Volumes;

    //get local view of the data
    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    host_vec_array constraint_gradients = constraint_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

    //communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
    //communicate ghosts
    /*
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    */
    int rnum_elem = FEM_->rnum_elem;

    if(nodal_density_flag_){
      FEM_->compute_nodal_gradients(design_densities, constraint_gradients);
      for(int i = 0; i < FEM_->nlocal_nodes; i++){
        constraint_gradients(i,0) /= initial_mass;
      }
    }
    else{
      //update per element volumes
      FEM_->compute_element_volumes();
      //ROL_Element_Volumes = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Volumes);
      //local view of element volumes
      const_host_vec_array element_volumes = FEM_->Global_Element_Volumes->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
      for(int ig = 0; ig < rnum_elem; ig++)
        constraint_gradients(ig,0) = element_volumes(ig,0)/initial_mass;
    }

    ROL_Gradients = ROL::makePtr<ROL_MV>(constraint_gradients_distributed);
    real_t gradient_dot_v = ROL_Gradients->dot(v);
    //debug print
    //std::cout << "Constraint Gradient value " << gradient_dot_v << std::endl;

    (*jvp)[0] = gradient_dot_v;
    //std::cout << "Ended constraint grad on task " <<FEM_->myrank  << std::endl;
  }
  
  
  void hessVec( ROL::Vector<real_t> &hv, const ROL::Vector<real_t> &v, const ROL::Vector<real_t> &z, real_t &tol ) {
    
    // Unwrap hv
    ROL::Ptr<MV> hvp = getVector(hv);

    hvp->putScalar(0);
    
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
