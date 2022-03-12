#ifndef HEAT_CAPACITY_POTENTIAL_MINIMIZE_TOPOPT_H
#define HEAT_CAPACITY_POTENTIAL_MINIMIZE_TOPOPT_H

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
#include <Kokkos_View.hpp>
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"

#include "ROL_Types.hpp"
#include <ROL_TpetraMultiVector.hpp>
#include "ROL_Objective.hpp"
#include "ROL_Elementwise_Reduce.hpp"
#include "FEA_Module_Heat_Conduction.h"

class HeatCapacityPotentialMinimize_TopOpt : public ROL::Objective<real_t> {
  
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

  FEA_Module_Heat_Conduction *FEM_;
  ROL::Ptr<ROL_MV> ROL_RHS;
  ROL::Ptr<ROL_MV> ROL_Temperatures;
  ROL::Ptr<ROL_MV> ROL_Gradients;
  Teuchos::RCP<MV> constraint_gradients_distributed;
  Teuchos::RCP<MV> all_node_temperatures_distributed_temp;

  bool useLC_; // Use linear form of compliance.  Otherwise use quadratic form.

  ROL::Ptr<const MV> getVector( const V& x ) {
    return dynamic_cast<const ROL_MV&>(x).getVector();
  }

  ROL::Ptr<MV> getVector( V& x ) {
    return dynamic_cast<ROL_MV&>(x).getVector();
  }

public:
  bool nodal_density_flag_;
  size_t last_comm_step, last_solve_step, current_step;

  HeatCapacityPotentialMinimize_TopOpt(FEA_Module_Heat_Conduction *FEM, bool nodal_density_flag) 
    : FEM_(FEM), useLC_(true) {
      nodal_density_flag_ = nodal_density_flag;
      last_comm_step = last_solve_step = -1;
      current_step = 0;
      
      //deep copy solve data into the cache variable
      FEM_->all_cached_node_temperatures_distributed = Teuchos::rcp(new MV(*(FEM_->all_node_temperatures_distributed), Teuchos::Copy));
      all_node_temperatures_distributed_temp = FEM_->all_node_temperatures_distributed;

      constraint_gradients_distributed = Teuchos::rcp(new MV(FEM_->map, 1));

      ROL_RHS = ROL::makePtr<ROL_MV>(FEM_->Global_Nodal_RHS);
      ROL_Temperatures = ROL::makePtr<ROL_MV>(FEM_->node_temperatures_distributed);

      real_t current_heat_capacity_potential = ROL_Temperatures->dot(*ROL_RHS);
      std::cout.precision(10);
      if(FEM_->myrank==0)
      std::cout << "INITIAL HEAT CAPACITY POTENTIAL " << current_heat_capacity_potential << std::endl;
  }

  void update(const ROL::Vector<real_t> &z, ROL::UpdateType type, int iter = -1 ) {
    //debug
    std::ostream &out = std::cout;
    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

    current_step++;
    ROL::Ptr<const MV> zp = getVector(z);
    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

    if (type == ROL::UpdateType::Initial)  {
      // This is the first call to update
      //first linear solve was done in FEA class run function already

      //initial design density data was already communicated for ghost nodes in init_design()
    }
    else if (type == ROL::UpdateType::Accept) {
      // u_ was set to u=S(x) during a trial update
      // and has been accepted as the new iterate
      /*assign temp pointer to the cache multivector (not the cache pointer) storage for a swap of the multivectors;
        this just avoids deep copy */
      all_node_temperatures_distributed_temp = FEM_->all_cached_node_temperatures_distributed;
      // Cache the accepted value
      FEM_->all_cached_node_temperatures_distributed = FEM_->all_node_temperatures_distributed;
    }
    else if (type == ROL::UpdateType::Revert) {
      // u_ was set to u=S(x) during a trial update
      // and has been rejected as the new iterate
      // Revert to cached value
      FEM_->comm_variables(zp);
      FEM_->all_node_temperatures_distributed = FEM_->all_cached_node_temperatures_distributed;
    }
    else if (type == ROL::UpdateType::Trial) {
      // This is a new value of x
      FEM_->all_node_temperatures_distributed = all_node_temperatures_distributed_temp;
      //communicate density variables for ghosts
      FEM_->comm_variables(zp);
      //update deformation variables
      FEM_->update_linear_solve(zp);
      if(FEM_->myrank==0)
      *fos << "called Trial" << std::endl;
    }
    else { // ROL::UpdateType::Temp
      // This is a new value of x used for,
      // e.g., finite-difference checks
      if(FEM_->myrank==0)
      *fos << "called Temp" << std::endl;
      FEM_->all_node_temperatures_distributed = all_node_temperatures_distributed_temp;
      FEM_->comm_variables(zp);
      FEM_->update_linear_solve(zp);
    }

    //decide to output current optimization state
    if(current_step%FEM_->simparam->optimization_output_freq==0)
      FEM_->Solver_Pointer_->tecplot_writer();
  }

  real_t value(const ROL::Vector<real_t> &z, real_t &tol) {
    //std::cout << "Started obj value on task " <<FEM_->myrank  << std::endl;
    ROL::Ptr<const MV> zp = getVector(z);
    real_t c = 0.0;

    //debug print
    //std::ostream &out = std::cout;
    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    //if(FEM_->myrank==0)
    //*fos << "Value function z:" << std::endl;
    //zp->describe(*fos,Teuchos::VERB_EXTREME);
    //*fos << std::endl;
    //std::fflush(stdout);

    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    //communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
    /*
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    
    if(last_solve_step!=current_step){
      //std::cout << "UPDATED TEMPERATURES" << std::endl;
      FEM_->update_linear_solve(zp);
      last_solve_step = current_step;
    }
    */
    //debug print of displacements
    //std::ostream &out = std::cout;
    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    //if(FEM_->myrank==0)
    //*fos << "Displacement data :" << std::endl;
    //FEM_->node_temperatures_distributed->describe(*fos,Teuchos::VERB_EXTREME);
    //*fos << std::endl;
    //std::fflush(stdout);
    
    ROL_RHS = ROL::makePtr<ROL_MV>(FEM_->Global_Nodal_RHS);
    ROL_Temperatures = ROL::makePtr<ROL_MV>(FEM_->node_temperatures_distributed);

    real_t current_heat_capacity_potential = ROL_Temperatures->dot(*ROL_RHS);
    std::cout.precision(10);
    if(FEM_->myrank==0)
    std::cout << "CURRENT HEAT CAPACITY POTENTIAL " << current_heat_capacity_potential << std::endl;

    //std::cout << "Ended obj value on task " <<FEM_->myrank  << std::endl;
    return current_heat_capacity_potential;
  }

  //void gradient_1( ROL::Vector<real_t> &g, const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
    //g.zero();
  //}
  
  void gradient( ROL::Vector<real_t> &g, const ROL::Vector<real_t> &z, real_t &tol ) {
    //std::cout << "Started obj gradient on task " <<FEM_->myrank  << std::endl;
    //get Tpetra multivector pointer from the ROL vector
    ROL::Ptr<const MV> zp = getVector(z);
    ROL::Ptr<MV> gp = getVector(g);

    //communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
    //FEM_->gradient_print_sync=1;
    /*
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    
    if(last_solve_step!=current_step){
      FEM_->update_linear_solve(zp);
      last_solve_step = current_step;
    }
    */
    //FEM_->gradient_print_sync=0;
    //get local view of the data
    host_vec_array objective_gradients = gp->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

    if(nodal_density_flag_){
      FEM_->compute_adjoint_gradients(design_densities, objective_gradients);
      //debug print of gradient
      //std::ostream &out = std::cout;
      //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
      //if(FEM_->myrank==0)
      //*fos << "Gradient data :" << std::endl;
      //gp->describe(*fos,Teuchos::VERB_EXTREME);
      //*fos << std::endl;
      //std::fflush(stdout);
      //for(int i = 0; i < FEM_->nlocal_nodes; i++){
        //objective_gradients(i,0) *= -1;
      //}
    }
    else{
      
      //update per element volumes
      //FEM_->compute_element_volumes();
      //ROL_Element_Volumes = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Volumes);
      //local view of element volumes
      //const_host_vec_array element_volumes = FEM_->Global_Element_Volumes->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
      //for(int ig = 0; ig < rnum_elem; ig++)
        //objective_gradients(ig,0) = element_volumes(ig,0);
        
    }
    //std::cout << "Objective Gradient called"<< std::endl;
    //debug print of design variables
    //std::ostream &out = std::cout;
    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    //if(FEM_->myrank==0)
    //*fos << "Gradient data :" << std::endl;
    //gp->describe(*fos,Teuchos::VERB_EXTREME);
    //*fos << std::endl;
    //std::fflush(stdout);
    //std::cout << "ended obj gradient on task " <<FEM_->myrank  << std::endl;
  }
  
  
  void hessVec( ROL::Vector<real_t> &hv, const ROL::Vector<real_t> &v, const ROL::Vector<real_t> &z, real_t &tol ) {
    //debug
    std::ostream &out = std::cout;
    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    // Unwrap hv
    ROL::Ptr<MV> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const MV> vp = getVector(v);
    ROL::Ptr<const MV> zp = getVector(z);
    
    host_vec_array objective_hessvec = hvp->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    const_host_vec_array direction_vector = vp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

    FEM_->compute_adjoint_hessian_vec(design_densities, objective_hessvec, vp);
    //if(FEM_->myrank==0)
    //std::cout << "hessvec" << std::endl;
    //vp->describe(*fos,Teuchos::VERB_EXTREME);
    //hvp->describe(*fos,Teuchos::VERB_EXTREME);
    if(FEM_->myrank==0)
    *fos << "Called Hessianvec" << std::endl;
    FEM_->hessvec_count++;
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
