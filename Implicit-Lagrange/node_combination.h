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
  //Trilinos type definitions
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>::size_type SizeType;
  using traits = Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>;
  
  using array_layout    = typename traits::array_layout;
  using execution_space = typename traits::execution_space;
  using device_type     = typename traits::device_type;
  using memory_traits   = typename traits::memory_traits;
  int num_dim_;

public:
    
  CArrayKokkos<GO, array_layout, device_type, memory_traits> node_set;
  GO patch_id, element_id; 
  LO local_patch_id;

  //Default Constructor
  Node_Combination(){}

  //Constructor with initialization
  Node_Combination(CArrayKokkos<GO, array_layout, device_type, memory_traits> &nodes_init) {
    node_set = nodes_init;
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

    //check if the nodes in the set are the same; sort them to simplify
    std::sort(this->node_set.pointer(),this->node_set.pointer()+this->node_set.size());
    std::sort(not_this.node_set.pointer(),not_this.node_set.pointer()+not_this.node_set.size());\

    //loop through the sorted nodes to check for equivalence
    for(int i = 0; i < this_size; i++)
      if(this->node_set(i)!=not_this.node_set(i)) return false;

    return true;
    
  }

};

#endif // end STATE_H
