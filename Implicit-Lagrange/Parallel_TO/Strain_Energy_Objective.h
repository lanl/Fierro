#ifndef STRAIN_ENERGY_OBJECTIVE_TOPOPT_H
#define STRAIN_ENERGY_OBJECTIVE_TOPOPT_H

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
#include "ROL_Objective.hpp"
#include "ROL_Elementwise_Reduce.hpp"
#include "Parallel_Nonlinear_Solver.h"

class StrainEnergyObjective_TopOpt : public ROL::Objective<real_t> {
  
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
  ROL::Ptr<ROL_MV> ROL_Force;
  ROL::Ptr<ROL_MV> ROL_Displacements;
  ROL::Ptr<ROL_MV> ROL_Gradients;
  Teuchos::RCP<MV> constraint_gradients_distributed;
  real_t target_strain_energy_;

  bool useLC_; // Use linear form of compliance.  Otherwise use quadratic form.

  ROL::Ptr<const MV> getVector( const V& x ) {
    return dynamic_cast<const ROL_MV&>(x).getVector();
  }

  ROL::Ptr<MV> getVector( V& x ) {
    return dynamic_cast<ROL_MV&>(x).getVector();
  }

public:
  bool nodal_density_flag_;
  size_t last_comm_step, current_step, last_solve_step;

  StrainEnergyObjective_TopOpt(Parallel_Nonlinear_Solver *FEM, bool nodal_density_flag, real_t target_strain_energy) 
    : FEM_(FEM), useLC_(true) {
      nodal_density_flag_ = nodal_density_flag;
      target_strain_energy_ = target_strain_energy;
      last_comm_step = last_solve_step = -1;
      current_step = 0;
      constraint_gradients_distributed = Teuchos::rcp(new MV(FEM_->map, 1));
  }

  void update(const ROL::Vector<real_t> &z, ROL::UpdateType type, int iter = -1 ) {
    current_step++;
  }

  real_t value(const ROL::Vector<real_t> &z, real_t &tol) {
    ROL::Ptr<const MV> zp = getVector(z);
    real_t c = 0.0;

    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    
    //communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    if(last_solve_step!=current_step){
      FEM_->update_linear_solve(zp);
      last_solve_step = current_step;
    }

    ROL_Force = ROL::makePtr<ROL_MV>(FEM_->Global_Nodal_Forces);
    ROL_Displacements = ROL::makePtr<ROL_MV>(FEM_->node_displacements_distributed);

    real_t current_strain_energy = ROL_Displacements->dot(*ROL_Force);
    std::cout.precision(10);
    if(FEM_->myrank==0)
    std::cout << "CURRENT STRAIN ENERGY " << current_strain_energy << std::endl;
    return (current_strain_energy - target_strain_energy_)*(current_strain_energy - target_strain_energy_)/2;
  }

  //void gradient_1( ROL::Vector<real_t> &g, const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
    //g.zero();
  //}
  
  void gradient( ROL::Vector<real_t> &g, const ROL::Vector<real_t> &z, real_t &tol ) {
    //get Tpetra multivector pointer from the ROL vector
    ROL::Ptr<const MV> zp = getVector(z);
    ROL::Ptr<MV> gp = getVector(g);

    //communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
    FEM_->gradient_print_sync=1;
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    if(last_solve_step!=current_step){
      FEM_->update_linear_solve(zp);
      last_solve_step = current_step;
    }
    
    FEM_->gradient_print_sync=0;
    //get local view of the data
    host_vec_array objective_gradients = gp->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

    int rnum_elem = FEM_->rnum_elem;
    real_t current_strain_energy = ROL_Displacements->dot(*ROL_Force);

    if(nodal_density_flag_){
      FEM_->compute_adjoint_gradients(design_densities, objective_gradients);
      //debug print of gradient
      //std::ostream &out = std::cout;
      //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
      //if(FEM_->myrank==0)
      //*fos << "Gradient data :" << std::endl;
      //FEM_->all_node_displacements_distributed->describe(*fos,Teuchos::VERB_EXTREME);
      //*fos << std::endl;
      //std::fflush(stdout);
      //multiply by difference in energy to get final gradient values
      for(int i = 0; i < FEM_->nlocal_nodes; i++){
        objective_gradients(i,0) *= (current_strain_energy - target_strain_energy_);
      }
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
  }
  
  /*
  void hessVec_12( ROL::Vector<real_t> &hv, const ROL::Vector<real_t> &v, 
                   const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
    
    // Unwrap hv
    ROL::Ptr<MV> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const MV> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const MV> up = getVector(u);
    ROL::Ptr<const MV> zp = getVector(z);

    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      MV KU(up->size(),0.0);
      MV U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      FEM_->apply_jacobian(KU,U,*zp,*vp);
      for (size_t i=0; i<up->size(); i++) {
        (*hvp)[i] = 2.0*KU[i];
      }
    }
    
  }

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
