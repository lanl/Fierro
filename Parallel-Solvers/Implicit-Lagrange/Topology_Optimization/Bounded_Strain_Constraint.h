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
 
#ifndef BOUNDED_STRAIN_CONSTRAINT_TOPOPT_H
#define BOUNDED_STRAIN_CONSTRAINT_TOPOPT_H

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
#include "FEA_Module_Elasticity.h"

class BoundedStrainConstraint_TopOpt : public ROL::Constraint<real_t> {
  
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

  FEA_Module_Elasticity *FEM_;
  ROL::Ptr<MV> Element_Masses;
  ROL::Ptr<ROL_MV> ROL_Element_Masses;
  real_t maximum_strain_;

  bool useLC_; // Use linear form of compliance.  Otherwise use quadratic form.

  ROL::Ptr<const MV> getVector( const V& x ) {
    return dynamic_cast<const ROL_MV&>(x).getVector();
  }

  ROL::Ptr<MV> getVector( V& x ) {
    return dynamic_cast<ROL_MV&>(x).getVector();
  }

public:
  bool nodal_density_flag_;
  int last_comm_step, current_step, last_solve_step;
  std::string my_fea_module = "Elasticity";

  BoundedStrainConstraint_TopOpt(FEA_Module *FEM, bool nodal_density_flag, real_t maximum_strain) 
    : useLC_(true) {
      FEM_ = dynamic_cast<FEA_Module_Elasticity*>(FEM);
      nodal_density_flag_ = nodal_density_flag;
      Element_Masses = ROL::makePtr<MV>(FEM_->element_map,1,true);
      maximum_strain_ = maximum_strain;
      last_comm_step = last_solve_step = -1;
      current_step = 0;
  }

  void update(const ROL::Vector<real_t> &z, ROL::UpdateType type, int iter = -1 ) {
    current_step++;
  }

  void value(ROL::Vector<real_t> &c, const ROL::Vector<real_t> &z, real_t &tol ) {
    ROL::Ptr<const MV> zp = getVector(z);
    //ROL::Ptr<MV> cp = getVector(c);
    ROL::Ptr<std::vector<real_t>> cp = dynamic_cast<ROL::StdVector<real_t>&>(c).getVector();
    int num_dim = FEM_->simparam->num_dim;
    int strain_count;
    if(num_dim==3) strain_count = 6;
    else strain_count = 3;
    

    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    const_host_vec_array design_displacements = FEM_->node_displacements_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    //host_vec_array constraint_view = cp->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

    //communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
    if(last_comm_step!=current_step){
      FEM_->comm_variables(zp);
      last_comm_step = current_step;
    }
    if(last_solve_step!=current_step){
      FEM_->update_linear_solve(zp, current_step);
      last_solve_step = current_step;
    }

    FEM_->compute_nodal_strains();

    //get local view of strains
    const_host_vec_array local_nodal_strains = FEM_->node_strains_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    int nlocal_nodes = FEM_->nlocal_nodes;

    //debug print of strain variables
    std::ostream &out = std::cout;
    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    if(FEM_->myrank==0)
    *fos << "Global Strain data :" << std::endl;
    FEM_->node_strains_distributed->describe(*fos,Teuchos::VERB_EXTREME);
    *fos << std::endl;
    std::fflush(stdout);
    
    //find maximum absolute value of strain component in the model
    real_t maximum_strain = 0;
    real_t global_maximum_strain = 0;
    //loop through all stress components evaluated at node positions (assumes Lagrange interpolation)
    for(int inode = 0; inode < nlocal_nodes; inode++)
     for(int istrain = 0; istrain < strain_count; istrain++)
       if(fabs(local_nodal_strains(inode,istrain)) > maximum_strain)
         maximum_strain = fabs(local_nodal_strains(inode,istrain));

    //std::cout << "STRAIN VALUE " << maximum_strain << std::endl;
    //collective comm to find max
    MPI_Allreduce(&maximum_strain, &global_maximum_strain, 1, MPI_DOUBLE, MPI_MAX, FEM_->world);
    //debug print
    //std::cout << "GLOBAL STRAIN VALUE " << global_maximum_strain << std::endl;
    (*cp)[0] = maximum_strain_ - global_maximum_strain;
    //debug print
    std::cout << "CONSTRAINT VALUE " << maximum_strain_ - global_maximum_strain << std::endl;
  }
  /*
  void gradient_1( ROL::Vector<real_t> &g, const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
    g.zero();
  }

  void gradient_2( ROL::Vector<real_t> &g, const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
    //get Tpetra multivector pointer from the ROL vector
    ROL::Ptr<const MV> up = getVector(u);
    ROL::Ptr<const MV> zp = getVector(z);
    ROL::Ptr<MV> gp = getVector(g);
    
    ROL::Ptr<ROL_MV> ROL_Element_Volumes;

    //communicate ghosts here again? check if variables were updated since last communication

    //get local view of the data
    host_vec_array objective_gradients = gp->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    const_host_vec_array design_densities = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    const_host_vec_array design_displacement = up->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

    int rnum_elem = FEM_->rnum_elem;

    if(nodal_density_flag_){
      FEM_->compute_nodal_gradients(design_densities, objective_gradients);
    }
    else{
      //update per element volumes
      FEM_->compute_element_volumes();
      //ROL_Element_Volumes = ROL::makePtr<ROL_MV>(FEM_->Global_Element_Volumes);
      //local view of element volumes
      const_host_vec_array element_volumes = FEM_->Global_Element_Volumes->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
      for(int ig = 0; ig < rnum_elem; ig++)
        objective_gradients(ig,0) = element_volumes(ig,0);
    }
    
  }
  
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
