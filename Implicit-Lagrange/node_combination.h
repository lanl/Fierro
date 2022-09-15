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
 
#ifndef NODE_COMBINATION_H
#define NODE_COMBINATION_H  

#include "matar.h"
#include "utilities.h"

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Kokkos_View.hpp>

using namespace utils;
class Node_Combination;
bool operator< (const Node_Combination &object1, const Node_Combination &object2);

class Node_Combination {

private:
  
  int num_dim_;

public:
  //Trilinos type definitions
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>::size_type SizeType;
  using traits = Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>;
  
  using array_layout    = typename traits::array_layout;
  using execution_space = typename traits::execution_space;
  using device_type     = typename traits::device_type;
  using memory_traits   = typename traits::memory_traits;
    
  CArray<GO> node_set, sort_set;
  GO patch_id, element_id; 
  LO local_patch_id;

  //Default Constructor
  Node_Combination(){}

  //Constructor with initialization
  Node_Combination(CArray<GO> &nodes_init) {
    node_set = nodes_init;
    sort_set = CArray<GO>(node_set.size());
  }

  //Destructor
  ~Node_Combination( ) {}

  //overload = operator
  Node_Combination& operator= (const Node_Combination &not_this){
    node_set = not_this.node_set;
    patch_id = not_this.patch_id;
    element_id = not_this.element_id;
    local_patch_id = not_this.local_patch_id;
    return *this;
  }

  //overload = operator
  bool operator== (Node_Combination &not_this){
    int this_size = this->node_set.size();
    //check if this node combination is identical
    //first check size of the combination
    if(this_size!=not_this.node_set.size())
      return false;

    CArray<GO> sort_set1 = this->sort_set;
    CArray<GO> sort_set2 = not_this.sort_set;
    for(int i = 0; i < this_size; i++){
      sort_set1(i) = this->node_set(i);
      sort_set2(i) = not_this.node_set(i);
    }

    //check if the nodes in the set are the same; sort them to simplify
    std::sort(sort_set1.pointer(),sort_set1.pointer()+sort_set1.size());
    std::sort(sort_set2.pointer(),sort_set2.pointer()+sort_set2.size());\

    //loop through the sorted nodes to check for equivalence
    for(int i = 0; i < this_size; i++)
      if(sort_set1(i)!=sort_set2(i)) return false;

    return true;
    
  }

};

#endif // end STATE_H
