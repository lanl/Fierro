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
 
#ifndef AREA_NORMALS_SHAPEOPT_H
#define AREA_NORMALS_SHAPEOPT_H

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
#include "FEA_Module_Dynamic_Elasticity.h"
#include "Explicit_Solver.h"

class AreaNormals : public ROL::Objective<real_t> {
  
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

  Explicit_Solver *Explicit_Solver_Pointer_;
  FEA_Module_SGH *FEM_SGH_;
  FEA_Module_Dynamic_Elasticity *FEM_Dynamic_Elasticity_;
  ROL::Ptr<ROL_MV> ROL_Force;
  ROL::Ptr<ROL_MV> ROL_Velocities;
  ROL::Ptr<ROL_MV> ROL_Gradients;

  ROL::Ptr<const MV> getVector( const V& x ) {
    return dynamic_cast<const ROL_MV&>(x).getVector();
  }

  ROL::Ptr<MV> getVector( V& x ) {
    return dynamic_cast<ROL_MV&>(x).getVector();
  }

public:
  bool nodal_density_flag_;
  int last_comm_step, last_solve_step, current_step;
  size_t nvalid_modules;
  std::vector<FEA_MODULE_TYPE> valid_fea_modules; //modules that may interface with this objective function
  FEA_MODULE_TYPE set_module_type;
  int set_node_, set_dim_;
  //std::string my_fea_module = "SGH";

  AreaNormals(Explicit_Solver *Explicit_Solver_Pointer, int set_node, int set_dim, bool nodal_density_flag) 
   {
      Explicit_Solver_Pointer_ = Explicit_Solver_Pointer;
      
      valid_fea_modules.push_back(FEA_MODULE_TYPE::SGH);
      valid_fea_modules.push_back(FEA_MODULE_TYPE::Dynamic_Elasticity);
      nvalid_modules = valid_fea_modules.size();

      const Simulation_Parameters& simparam = Explicit_Solver_Pointer_->simparam;
      for (const auto& fea_module : Explicit_Solver_Pointer_->fea_modules) {
        for(int ivalid = 0; ivalid < nvalid_modules; ivalid++){
          if(fea_module->Module_Type==FEA_MODULE_TYPE::SGH){
            FEM_SGH_ = dynamic_cast<FEA_Module_SGH*>(fea_module);
            set_module_type = FEA_MODULE_TYPE::SGH;
          }
          if(fea_module->Module_Type==FEA_MODULE_TYPE::Dynamic_Elasticity){
            FEM_Dynamic_Elasticity_ = dynamic_cast<FEA_Module_Dynamic_Elasticity*>(fea_module);
            set_module_type = FEA_MODULE_TYPE::Dynamic_Elasticity;
          }
        }
      }
      nodal_density_flag_ = nodal_density_flag;
      last_comm_step = last_solve_step = -1;
      current_step = 0;
      set_node_ = set_node;
      set_dim_ = set_dim;

      //real_t current_kinetic_energy = ROL_Velocities->dot(*ROL_Force)/2;
      //std::cout.precision(10);
      //if(FEM_->myrank==0)
      //std::cout << "INITIAL KINETIC ENERGY " << current_kinetic_energy << std::endl;
  }

  void update(const ROL::Vector<real_t> &z, ROL::UpdateType type, int iter = -1 ) {
  }

  real_t value(const ROL::Vector<real_t> &z, real_t &tol) {
    
    double area_normal_array[24];
    size_t nodes_in_elem_array[8];
    double node_coords_array[2*8*3];
    ViewCArrayKokkos <size_t> elem_node_gids(nodes_in_elem_array, 8);
    for(int init = 0; init < 8; init++){
      elem_node_gids(init) = init;
    }
    ViewCArrayKokkos <double> area_normal(area_normal_array, 8, 3);
    const size_t rk_level = FEM_SGH_->simparam->dynamic_options.rk_num_bins - 1;
    //std::cout << "Started obj value on task " <<FEM_->myrank  << std::endl;
    ROL::Ptr<const MV> zp = getVector(z);

    const_host_vec_array design_coordinates = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    DViewCArrayKokkos<double> node_coords(node_coords_array,2,8,3);
    for(int inode = 0; inode < 8; inode++){
      size_t local_node_id = FEM_SGH_->nodes_in_elem(0,inode);
      node_coords(1,inode,0) = design_coordinates(local_node_id,0);
      node_coords(1,inode,1) = design_coordinates(local_node_id,1);
      node_coords(1,inode,2) = design_coordinates(local_node_id,2);
    }

    // std::cout.precision(10);
    // if(Explicit_Solver_Pointer_->myrank==0)
    // std::cout << "CURRENT TIME INTEGRAL OF KINETIC ENERGY " << objective_accumulation << std::endl;
    FEM_SGH_->get_bmatrix(area_normal,
                      0,
                      node_coords,
                      elem_node_gids,
                      1);

    real_t objective_value = -area_normal(set_node_, set_dim_);
    //std::cout << "Ended obj value on task " <<FEM_->myrank  << std::endl;
    return objective_value;
  }

  //void gradient_1( ROL::Vector<real_t> &g, const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
    //g.zero();
  //}
  
  void gradient( ROL::Vector<real_t> &g, const ROL::Vector<real_t> &z, real_t &tol ) {
    //std::cout << "Started obj gradient on task " <<FEM_->myrank  << std::endl;
    //get Tpetra multivector pointer from the ROL vector
    double area_normal_gradients_array[8*8*3*3];
    size_t nodes_in_elem[8];
    double node_coords_array[2*8*3];
    ViewCArrayKokkos <size_t> elem_node_gids(nodes_in_elem, 8);
    for(int init = 0; init < 8; init++){
      elem_node_gids(init) = init;
    }
    
    ViewCArrayKokkos <double> area_normal_gradients(area_normal_gradients_array, 8, 3, 8, 3);
    const size_t rk_level = FEM_SGH_->simparam->dynamic_options.rk_num_bins - 1;
    ROL::Ptr<const MV> zp = getVector(z);
    ROL::Ptr<MV> gp = getVector(g);
    const_host_vec_array design_coordinates = zp->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    host_vec_array design_gradients = gp->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    DViewCArrayKokkos<double> node_coords(node_coords_array,2,8,3);
    for(int inode = 0; inode < 8; inode++){
      size_t local_node_id = FEM_SGH_->nodes_in_elem(0,inode);
      node_coords(1,inode,0) = design_coordinates(local_node_id,0);
      node_coords(1,inode,1) = design_coordinates(local_node_id,1);
      node_coords(1,inode,2) = design_coordinates(local_node_id,2);
    }

    for (int igradient = 0; igradient < FEM_SGH_->nlocal_nodes+FEM_SGH_->nghost_nodes; igradient++){
      design_gradients(igradient,0) = 0;
      design_gradients(igradient,1) = 0;
      design_gradients(igradient,2) = 0;
    }


    FEM_SGH_->get_bmatrix_gradients(area_normal_gradients,
                0,
                node_coords,
                elem_node_gids,
                1);

    for(int inode = 0; inode < 8; inode++){
      size_t local_node_id = FEM_SGH_->nodes_in_elem(0,inode);
      for(int idim = 0; idim < 3; idim++)
        design_gradients(local_node_id,idim) = -area_normal_gradients(set_node_,set_dim_,inode,idim);
    }
  }
  
};

#endif // end header guard
