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
 
#ifndef KINETIC_ENERGY_MINIMIZE_TOPOPT_H
#define KINETIC_ENERGY_MINIMIZE_TOPOPT_H

#include "matar.h"
#include "elements.h"
#include <string>
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
#include "ROL_Objective.hpp"
#include "ROL_Elementwise_Reduce.hpp"
#include "FEA_Module_SGH.h"
#include "Simulation_Parameters_Dynamic_Optimization.h"

class KineticEnergyMinimize_TopOpt : public ROL::Objective<real_t> {
  
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
  typedef Kokkos::View<const real_t**, array_layout, device_type, memory_traits> const_vec_array;
  typedef MV::dual_view_type dual_vec_array;

private:

  FEA_Module_SGH *FEM_;
  ROL::Ptr<ROL_MV> ROL_Force;
  ROL::Ptr<ROL_MV> ROL_Velocities;
  ROL::Ptr<ROL_MV> ROL_Gradients;
  Teuchos::RCP<MV> all_node_velocities_distributed_temp;

  bool useLC_; // Use linear form of energy.  Otherwise use quadratic form.

  ROL::Ptr<const MV> getVector( const V& x ) {
    return dynamic_cast<const ROL_MV&>(x).getVector();
  }

  ROL::Ptr<MV> getVector( V& x ) {
    return dynamic_cast<ROL_MV&>(x).getVector();
  }

public:
  bool nodal_density_flag_, time_accumulation;
  int last_comm_step, last_solve_step, current_step;
  std::string my_fea_module = "SGH";
  real_t objective_accumulation;

  KineticEnergyMinimize_TopOpt(FEA_Module *FEM, bool nodal_density_flag) 
    : useLC_(true) {
      FEM_ = dynamic_cast<FEA_Module_SGH*>(FEM);
      nodal_density_flag_ = nodal_density_flag;
      last_comm_step = last_solve_step = -1;
      current_step = 0;
      time_accumulation = true;
      objective_accumulation = 0;
      
      //deep copy solve data into the cache variable
      FEM_->all_cached_node_velocities_distributed = Teuchos::rcp(new MV(*(FEM_->all_node_velocities_distributed), Teuchos::Copy));
      all_node_velocities_distributed_temp = FEM_->all_node_velocities_distributed;

      //ROL_Force = ROL::makePtr<ROL_MV>(FEM_->Global_Nodal_Forces);
      ROL_Velocities = ROL::makePtr<ROL_MV>(FEM_->node_velocities_distributed);

      //real_t current_kinetic_energy = ROL_Velocities->dot(*ROL_Force)/2;
      //std::cout.precision(10);
      //if(FEM_->myrank==0)
      //std::cout << "INITIAL KINETIC ENERGY " << current_kinetic_energy << std::endl;
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
      all_node_velocities_distributed_temp = FEM_->all_cached_node_velocities_distributed;
      // Cache the accepted value
      FEM_->all_cached_node_velocities_distributed = FEM_->all_node_velocities_distributed;
    }
    else if (type == ROL::UpdateType::Revert) {
      // u_ was set to u=S(x) during a trial update
      // and has been rejected as the new iterate
      // Revert to cached value
      FEM_->comm_variables(zp);
      FEM_->all_node_velocities_distributed = FEM_->all_cached_node_velocities_distributed;
    }
    else if (type == ROL::UpdateType::Trial) {
      // This is a new value of x
      FEM_->all_node_velocities_distributed = all_node_velocities_distributed_temp;
      //communicate density variables for ghosts
      FEM_->comm_variables(zp);
      //update deformation variables
      FEM_->update_forward_solve(zp);
      if(FEM_->myrank==0)
      *fos << "called Trial" << std::endl;
    }
    else { // ROL::UpdateType::Temp
      // This is a new value of x used for,
      // e.g., finite-difference checks
      if(FEM_->myrank==0)
      *fos << "called Temp" << std::endl;
      FEM_->all_node_velocities_distributed = all_node_velocities_distributed_temp;
      FEM_->comm_variables(zp);
      FEM_->update_forward_solve(zp);
    }

    //decide to output current optimization state
    if(current_step%FEM_->simparam_dynamic_opt.optimization_output_freq==0)
      FEM_->Explicit_Solver_Pointer_->parallel_tecplot_writer();
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
      //std::cout << "UPDATED velocities" << std::endl;
      FEM_->update_linear_solve(zp);
      last_solve_step = current_step;
    }
    */
    //debug print of velocities
    //std::ostream &out = std::cout;
    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    //if(FEM_->myrank==0)
    //*fos << "Displacement data :" << std::endl;
    //FEM_->node_velocities_distributed->describe(*fos,Teuchos::VERB_EXTREME);
    //*fos << std::endl;
    //std::fflush(stdout);
    
    //ROL_Force = ROL::makePtr<ROL_MV>(FEM_->Global_Nodal_Forces);
    ROL_Velocities = ROL::makePtr<ROL_MV>(FEM_->node_velocities_distributed);

    std::cout.precision(10);
    if(FEM_->myrank==0)
    std::cout << "CURRENT TIME INTEGRAL OF KINETIC ENERGY " << objective_accumulation << std::endl;

    //std::cout << "Ended obj value on task " <<FEM_->myrank  << std::endl;
    return objective_accumulation;
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
    //FEM_->gradient_print_sync=0;
    //get local view of the data
    vec_array objective_gradients = gp->getLocalView<device_type> (Tpetra::Access::ReadWrite);
    const_vec_array design_densities = zp->getLocalView<device_type> (Tpetra::Access::ReadOnly);

    FEM_->compute_topology_optimization_gradient_full(design_densities, objective_gradients);
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
  
  /*
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
    *fos << "Called Strain Energy Hessianvec" << std::endl;
    FEM_->hessvec_count++;
  }
  */
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
