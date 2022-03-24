#ifndef HEAT_CAPACITY_POTENTIAL_CONSTRAINT_TOPOPT_H
#define HEAT_CAPACITY_POTENTIAL_CONSTRAINT_TOPOPT_H

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
#include "ROL_Constraint.hpp"
#include "ROL_Elementwise_Reduce.hpp"
#include "FEA_Module_Heat_Conduction.h"

class HeatCapacityPotentialConstraint_TopOpt : public ROL::Constraint<real_t> {
  
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
  ROL::Ptr<ROL_MV> ROL_Heat;
  ROL::Ptr<ROL_MV> ROL_Temperatures;
  ROL::Ptr<ROL_MV> ROL_Gradients;
  Teuchos::RCP<MV> all_node_temperatures_distributed_temp;
  Teuchos::RCP<MV> constraint_gradients_distributed;
  real_t initial_heat_capacity_potential_;
  bool inequality_flag_;
  real_t constraint_value_;

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
  std::string my_fea_module = "Heat_Conduction";

  HeatCapacityPotentialConstraint_TopOpt(FEA_Module *FEM, bool nodal_density_flag, bool inequality_flag=true, real_t constraint_value=0) 
    : useLC_(true) {
      FEM_ = dynamic_cast<FEA_Module_Heat_Conduction*>(FEM);
      nodal_density_flag_ = nodal_density_flag;
      last_comm_step = last_solve_step = -1;
      current_step = 0;
      inequality_flag_ = inequality_flag;
      constraint_value_ = constraint_value;
      constraint_gradients_distributed = Teuchos::rcp(new MV(FEM_->map, 1));

      //deep copy solve data into the cache variable
      FEM_->all_cached_node_temperatures_distributed = Teuchos::rcp(new MV(*(FEM_->all_node_temperatures_distributed), Teuchos::Copy));
      all_node_temperatures_distributed_temp = FEM_->all_node_temperatures_distributed;

      ROL_Heat = ROL::makePtr<ROL_MV>(FEM_->Global_Nodal_Heat);
      ROL_Temperatures = ROL::makePtr<ROL_MV>(FEM_->node_temperatures_distributed);

      initial_heat_capacity_potential_ = ROL_Temperatures->dot(*ROL_Heat)/2;
      std::cout.precision(10);
      if(FEM_->myrank==0)
        std::cout << "INITIAL HEAT CAPACITY POTENTIAL " << initial_heat_capacity_potential_ << std::endl;
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
  }

  void value(ROL::Vector<real_t> &c, const ROL::Vector<real_t> &z, real_t &tol ) {
    ROL::Ptr<const MV> zp = getVector(z);
    ROL::Ptr<std::vector<real_t>> cp = dynamic_cast<ROL::StdVector<real_t>&>(c).getVector();

    real_t current_heat_capacity_potential = ROL_Temperatures->dot(*ROL_Heat)/2;
    
    if(FEM_->myrank==0)
      std::cout << "CURRENT HEAT CAPACITY POTENTIAL RATIO " << current_heat_capacity_potential/initial_heat_capacity_potential_ << std::endl;

    if(inequality_flag_)
      (*cp)[0] = current_heat_capacity_potential/initial_heat_capacity_potential_;
    else
      (*cp)[0] = current_heat_capacity_potential/initial_heat_capacity_potential_ - constraint_value_;
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
    int nlocal_nodes = FEM_->nlocal_nodes;

    if(nodal_density_flag_){
      FEM_->compute_adjoint_gradients(design_densities, constraint_gradients);
      //debug print of gradient
      //std::ostream &out = std::cout;
      //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
      //if(FEM_->myrank==0)
      //*fos << "Gradient data :" << std::endl;
      //ajvp->describe(*fos,Teuchos::VERB_EXTREME);
      //*fos << std::endl;
      //std::fflush(stdout);
      for(int i = 0; i < nlocal_nodes; i++){
        constraint_gradients(i,0) *= (*vp)[0]/initial_heat_capacity_potential_;
      }
    }
    else{
      FEM_->compute_adjoint_gradients(design_densities, constraint_gradients);
      //debug print of gradient
      //std::ostream &out = std::cout;
      //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
      //if(FEM_->myrank==0)
      //*fos << "Gradient data :" << std::endl;
      //ajvp->describe(*fos,Teuchos::VERB_EXTREME);
      //*fos << std::endl;
      //std::fflush(stdout);
      for(int i = 0; i < rnum_elem; i++){
        constraint_gradients(i,0) *= (*vp)[0]/initial_heat_capacity_potential_;
      }
    }
    //std::cout << "Ended constraint adjoint grad on task " <<FEM_->myrank  << std::endl;
    //debug print
    //std::cout << "Constraint Gradient value " << std::endl;
  }
  
  void applyJacobian(ROL::Vector<real_t> &jv, const ROL::Vector<real_t> &v, const ROL::Vector<real_t> &x, real_t &tol) {
    //get Tpetra multivector pointer from the ROL vector
    ROL::Ptr<const MV> zp = getVector(x);
    ROL::Ptr<std::vector<real_t>> jvp = dynamic_cast<ROL::StdVector<real_t>&>(jv).getVector();

    //get local view of the data
    host_vec_array constraint_gradients = constraint_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

    int rnum_elem = FEM_->rnum_elem;
    int nlocal_nodes = FEM_->nlocal_nodes;

    FEM_->compute_adjoint_gradients(design_densities, constraint_gradients);
    if(nodal_density_flag_){
      for(int i = 0; i < nlocal_nodes; i++){
        constraint_gradients(i,0) /= initial_heat_capacity_potential_;
      }
    }
    else{
      for(int i = 0; i < rnum_elem; i++){
        constraint_gradients(i,0) /= initial_heat_capacity_potential_;
      }
    }
    ROL_Gradients = ROL::makePtr<ROL_MV>(constraint_gradients_distributed);
    real_t gradient_dot_v = ROL_Gradients->dot(v);
    //debug print
    //std::cout << "Constraint Gradient value " << gradient_dot_v << std::endl;

    (*jvp)[0] = gradient_dot_v;
  }

  void applyAdjointHessian(ROL::Vector<real_t> &ahuv, const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &v, const ROL::Vector<real_t> &z, Real &tol) {
    // Unwrap hv
    ROL::Ptr<MV> ahuvp = getVector(ahuv);

    // Unwrap v
    ROL::Ptr<const MV> vp = getVector(v);
    ROL::Ptr<const std::vector<real_t>> up = dynamic_cast<const ROL::StdVector<real_t>&>(u).getVector();
    ROL::Ptr<const MV> zp = getVector(z);
    
    host_vec_array constraint_hessvec = ahuvp->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    const_host_vec_array direction_vector = vp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    
    int rnum_elem = FEM_->rnum_elem;
    int nlocal_nodes = FEM_->nlocal_nodes;

    FEM_->compute_adjoint_hessian_vec(design_densities, constraint_hessvec, vp);
    for(int i = 0; i < nlocal_nodes; i++){
      constraint_hessvec(i,0) *= (*up)[0]/initial_heat_capacity_potential_;
    }
    //if(FEM_->myrank==0)
    //std::cout << "hessvec" << std::endl;
    //vp->describe(*fos,Teuchos::VERB_EXTREME);
    //hvp->describe(*fos,Teuchos::VERB_EXTREME);
    if(FEM_->myrank==0)
    *(FEM_->fos) << "Called Heat Capacity Potential Constraint Hessianvec" << std::endl;
    FEM_->hessvec_count++;
  }
};

#endif // end header guard
