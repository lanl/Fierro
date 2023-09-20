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
 
#ifndef MULTI_OBJECTIVE_TOPOPT_H
#define MULTI_OBJECTIVE_TOPOPT_H

#include "matar.h"
#include <string>
#include <vector>
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

class MultiObjective_TopOpt : public ROL::Objective<real_t> {
  
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
  
  Teuchos::RCP<MV> multi_derivative;
  ROL::Ptr<ROL_MV> ROL_Multi_Derivative;
  std::vector<ROL::Ptr<ROL::Objective<real_t>>> Multi_Objective_Terms_;
  std::vector<real_t> Multi_Objective_Weights_;
  int nobjectives;
  bool derivative_allocated;

  ROL::Ptr<const MV> getVector( const V& x ) {
    return dynamic_cast<const ROL_MV&>(x).getVector();
  }

  ROL::Ptr<MV> getVector( V& x ) {
    return dynamic_cast<ROL_MV&>(x).getVector();
  }

public:

  MultiObjective_TopOpt(std::vector<ROL::Ptr<ROL::Objective<real_t>>> Multi_Objective_Terms, std::vector<real_t> Multi_Objective_Weights){  
    Multi_Objective_Terms_ = Multi_Objective_Terms;
    Multi_Objective_Weights_ = Multi_Objective_Weights;
    nobjectives = Multi_Objective_Terms_.size();
    derivative_allocated = false;
  }

  void update(const ROL::Vector<real_t> &z, ROL::UpdateType type, int iter = -1 ) {
    for(int iobjective = 0; iobjective < nobjectives; iobjective++){
      Multi_Objective_Terms_[iobjective]->update(z, type, iter);
    }
  }

  real_t value(const ROL::Vector<real_t> &z, real_t &tol) {
    //std::cout << "Started obj value on task " <<FEM_->myrank  << std::endl;
    ROL::Ptr<const MV> zp = getVector(z);
    real_t c = 0.0;
    for(int iobjective = 0; iobjective < nobjectives; iobjective++){
      c += Multi_Objective_Weights_[iobjective]*Multi_Objective_Terms_[iobjective]->value(z, tol);
    }
    //std::cout << "Ended obj value on task " <<FEM_->myrank  << std::endl;
    return c;
  }

  //void gradient_1( ROL::Vector<real_t> &g, const ROL::Vector<real_t> &u, const ROL::Vector<real_t> &z, real_t &tol ) {
    //g.zero();
  //}
  
  void gradient( ROL::Vector<real_t> &g, const ROL::Vector<real_t> &z, real_t &tol ) {
    //get Tpetra multivector pointer from the ROL vector
    ROL::Ptr<const MV> zp = getVector(z);
    ROL::Ptr<MV> gp = getVector(g);

    if(!derivative_allocated){
      //derivative storage
      multi_derivative = Teuchos::rcp(new MV(gp->getMap(), 1));
      derivative_allocated = true;
    }
    multi_derivative->putScalar(0);

    for(int iobjective = 0; iobjective < nobjectives; iobjective++){
      Multi_Objective_Terms_[iobjective]->gradient(g, z, tol);
      multi_derivative->update(Multi_Objective_Weights_[iobjective], *gp, 1);
    }

    //copy values into g
    Tpetra::deep_copy(*gp, *multi_derivative);
  }
  
  
  void hessVec( ROL::Vector<real_t> &hv, const ROL::Vector<real_t> &v, const ROL::Vector<real_t> &z, real_t &tol ) {
    // Unwrap hv
    ROL::Ptr<MV> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const MV> vp = getVector(v);
    ROL::Ptr<const MV> zp = getVector(z);
    
    if(!derivative_allocated){
      //derivative storage
      multi_derivative = Teuchos::rcp(new MV(hvp->getMap(), 1));
      derivative_allocated = true;
    }
    multi_derivative->putScalar(0);

    for(int iobjective = 0; iobjective < nobjectives; iobjective++){
      Multi_Objective_Terms_[iobjective]->hessVec(hv, v, z, tol);
      multi_derivative->update(Multi_Objective_Weights_[iobjective], *hvp, 1);
    }

    //copy values into g
    Tpetra::deep_copy(*hvp, *multi_derivative);
  }

};

#endif // end header guard
