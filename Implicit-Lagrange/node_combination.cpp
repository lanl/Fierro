
#include "node_combination.h"
#include "utilities.h"
#include "matar.h"

using namespace utils;

//overload < operator
  bool operator< (const Node_Combination &object1, const Node_Combination &object2){
    int size1 = object1.node_set.size();
    //check if this node combination is identical
    //first check size of the combination; if smaller evaluate to true
    //the set using this is then ordered first according to size of the combinations
    if(size1<object2.node_set.size())
      return true;
    
    //This part sorts for segments of the set where combinations have the same size
    //define < using the sort of both combinations. If the first nonequal element of the lhs combination, w.r.t to 
    //the corresponding element of the rhs, is less than the respective element of the rhs < evaluates to true
    std::sort(object1.node_set.get_kokkos_view().data(),object1.node_set.get_kokkos_view().data()+object1.node_set.size());
    std::sort(object2.node_set.get_kokkos_view().data(),object2.node_set.get_kokkos_view().data()+object2.node_set.size());

    //loop through the sorted nodes to check for <
    for(int i = 0; i < size1; i++){
      if(object1.node_set(i)<object2.node_set(i)) return true;
      else if(object1.node_set(i)==object2.node_set(i)) continue;
      else break;
    }

    return false;
    
  }