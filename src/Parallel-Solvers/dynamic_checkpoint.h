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

#ifndef DYNAMIC_CHECKPOINT_H
#define DYNAMIC_CHECKPOINT_H

#include "matar.h"
#include "utilities.h"
#include <Teuchos_RCP.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Kokkos_Core.hpp>
#include "Tpetra_Details_DefaultTypes.hpp"

using namespace utils;
using namespace mtr;

class Dynamic_Checkpoint;
bool operator<(const Dynamic_Checkpoint& object1, const Dynamic_Checkpoint& object2);

class Dynamic_Checkpoint
{
public:
    // Trilinos type definitions
    typedef Tpetra::Map<>::local_ordinal_type LO;
    typedef Tpetra::Map<>::global_ordinal_type GO;
    typedef Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>::size_type SizeType;
    using traits = Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>;

    typedef Tpetra::MultiVector<real_t, LO, GO> MV;
    typedef Tpetra::MultiVector<GO, LO, GO> MCONN;

    // Default Constructor
    Dynamic_Checkpoint(){}

    // Argument Constructor
    Dynamic_Checkpoint(int vector_count, int timestep, real_t time, real_t dt, int pass_level=0) {
        num_state_vectors = vector_count;
        saved_timestep = timestep;
        saved_time = time;
        saved_dt = dt;
        state_vectors = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(num_state_vectors));
        level = pass_level;
    }

    // Copy Constructor
    Dynamic_Checkpoint(const Dynamic_Checkpoint &copied_checkpoint){
        saved_timestep = copied_checkpoint.saved_timestep;
        num_state_vectors = copied_checkpoint.num_state_vectors;
        state_vectors = copied_checkpoint.state_vectors;
        level = copied_checkpoint.level;
    }

    // Destructor
    ~Dynamic_Checkpoint() {}

    //assignment operator
    Dynamic_Checkpoint& operator= (const Dynamic_Checkpoint &assigned_checkpoint){
        saved_timestep = assigned_checkpoint.saved_timestep;
        num_state_vectors = assigned_checkpoint.num_state_vectors;
        state_vectors = assigned_checkpoint.state_vectors;
        level = assigned_checkpoint.level;
        return *this;
    }

    // overload = operator
    bool operator==(Dynamic_Checkpoint& not_this)
    {
        if(this->saved_timestep!=not_this.saved_timestep){
            return false;
        }
        return true;
    }

    //function to change timestep
    void change_timestep(int new_timestep){
        saved_timestep = new_timestep;
    }

    void change_time(real_t new_time){
        saved_time = new_time;
    }

    void change_level(int new_level){
        level = new_level;
    }

    //function to change stored vectors
    void change_vectors(Teuchos::RCP<std::vector<Teuchos::RCP<MV>>> new_state_vectors){
        state_vectors = new_state_vectors;
    }

    //copy vectors from another checkpoint object
    void copy_vectors(const Dynamic_Checkpoint &copied_checkpoint){
        state_vectors = copied_checkpoint.state_vectors;
    }

    //function to change one of the stored vectors
    void change_vector(int vector_index, Teuchos::RCP<MV> new_vector) const{
        (*state_vectors)[vector_index] = new_vector;
    }

    //function to assign new values to stored vectors
    void assign_vectors(Teuchos::RCP<std::vector<Teuchos::RCP<MV>>> new_state_vectors) const{
        //TODO: add error control for differing sizes of input and member vector
        for(int vector_index = 0; vector_index < (*new_state_vectors).size(); vector_index++){
            (*state_vectors)[vector_index]->assign(*((*new_state_vectors)[vector_index]));
        }
    }

    //function to assign new values to one of the stored vectors
    void assign_vector(int vector_index, Teuchos::RCP<MV> new_vector) const{
        (*state_vectors)[vector_index]->assign(*new_vector);
    }

    //get pointer to contained vector
    Teuchos::RCP<MV> get_vector_pointer(int vector_index) const{
        return (*state_vectors)[vector_index];
    }

    public:
    //checkpoint state data
    int saved_timestep;
    real_t saved_dt;
    real_t saved_time;
    int num_state_vectors;
    int level;

    private:
    Teuchos::RCP<std::vector<Teuchos::RCP<MV>>> state_vectors;

};

#endif // end STATE_H
