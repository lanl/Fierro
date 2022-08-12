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
    else if(size1>object2.node_set.size())
      return false;
    
    CArray<Node_Combination::GO> sort_set1(size1);
    CArray<Node_Combination::GO> sort_set2(size1);
    for(int i = 0; i < size1; i++){
      sort_set1(i) = object1.node_set(i);
      sort_set2(i) = object2.node_set(i);
    }
    //This part sorts for segments of the set where combinations have the same size
    //define < using the sort of both combinations. If the first nonequal element of the lhs combination, w.r.t to 
    //the corresponding element of the rhs, is less than the respective element of the rhs < evaluates to true
    std::sort(sort_set1.pointer(),sort_set1.pointer()+sort_set1.size());
    std::sort(sort_set2.pointer(),sort_set2.pointer()+sort_set2.size());

    //loop through the sorted nodes to check for <
    for(int i = 0; i < size1; i++){
      if(sort_set1(i)<sort_set2(i)) return true;
      else if(sort_set1(i)==sort_set2(i)) continue;
      else break;
    }

    return false;
    
  }