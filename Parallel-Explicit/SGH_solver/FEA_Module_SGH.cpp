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

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>  // fmin, fmax, abs note: fminl is long
#include <sys/stat.h>
#include <mpi.h>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Parallel_Reduce.hpp>
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "Tpetra_Details_FixedHashTable.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Import_Util2.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include <set>

#include "elements.h"
#include "swage.h"
#include "matar.h"
#include "utilities.h"
#include "node_combination.h"
#include "Simulation_Parameters_Elasticity.h"
#include "Simulation_Parameters_Topology_Optimization.h"
#include "FEA_Module_SGH.h"
#include "Explicit_Solver_SGH.h"

#define MAX_ELEM_NODES 8
#define STRAIN_EPSILON 0.000000001
#define DENSITY_EPSILON 0.0001
#define BC_EPSILON 1.0e-6

using namespace utils;


FEA_Module_SGH::FEA_Module_SGH(Explicit_Solver_SGH *Solver_Pointer){
  //create parameter object
  simparam = new Simulation_Parameters_Elasticity();
  // ---- Read input file, define state and boundary conditions ---- //
  simparam->input();
  
  Solver_Pointer_ = Solver_Pointer;
  
  //TO parameters
  simparam_TO = dynamic_cast<Simulation_Parameters_Topology_Optimization*>(Solver_Pointer_->simparam);

  //create ref element object
  ref_elem = new elements::ref_element();
  //create mesh objects
  //init_mesh = new swage::mesh_t(simparam);
  //mesh = new swage::mesh_t(simparam);

  //boundary condition data
  max_boundary_sets = 0;
  
  //setup output
  init_output();
}

FEA_Module_SGH::~FEA_Module_SGH(){
   delete simparam;
}

/* ----------------------------------------------------------------------
   Read ANSYS dat format mesh file
------------------------------------------------------------------------- */
void FEA_Module_SGH::read_conditions_ansys_dat(std::ifstream *in, std::streampos before_condition_header){

  char ch;
  int num_dim = simparam->num_dim;
  int buffer_lines = 1000;
  int max_word = 30;
  int p_order = simparam->p_order;
  real_t unit_scaling = simparam->unit_scaling;
  int local_node_index, current_column_index;
  size_t strain_count;
  std::string skip_line, read_line, substring, token;
  std::stringstream line_parse, line_parse2;
  CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;
  CArrayKokkos<long long int, array_layout, HostSpace, memory_traits> read_buffer_indices;
  int buffer_loop, buffer_iteration, buffer_iterations, scan_loop, nodes_per_element, words_per_line;
  size_t read_index_start, node_rid, elem_gid;
  LO local_dof_id;
  GO node_gid;
  real_t dof_value;
  host_vec_array node_densities;
 
} // end read_conditions_ansys_dat

/* ----------------------------------------------------------------------------
   Initialize sets of element boundary surfaces and arrays for input conditions
------------------------------------------------------------------------------- */

void FEA_Module_SGH::init_boundaries(){
  max_boundary_sets = simparam->NB;
  int num_dim = simparam->num_dim;
  
  // set the number of boundary sets
  if(myrank == 0)
    std::cout << "building boundary sets " << std::endl;
  
  //initialize to 1 since there must be at least 1 boundary set anyway; read in may occure later
  if(max_boundary_sets==0) max_boundary_sets = 1;
  //std::cout << "NUM BOUNDARY CONDITIONS ON RANK " << myrank << " FOR INIT " << num_boundary_conditions <<std::endl;
  init_boundary_sets(max_boundary_sets);

  //allocate nodal data
  Node_DOF_Boundary_Condition_Type = CArrayKokkos<int, array_layout, device_type, memory_traits>(nall_nodes*num_dim, "Node_DOF_Boundary_Condition_Type");

  //initialize
  for(int init=0; init < nall_nodes*num_dim; init++)
    Node_DOF_Boundary_Condition_Type(init) = NONE;

  Number_DOF_BCS = 0;
}

/* ----------------------------------------------------------------------
   initialize storage for element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void FEA_Module_SGH::init_boundary_sets (int num_sets){

  if(num_sets == 0){
    std::cout << " Warning: number of boundary conditions = 0";
    return;
  }
  //initialize maximum
  max_boundary_sets = num_sets;
  //std::cout << " DEBUG PRINT "<<num_sets << " " << nboundary_patches << std::endl;
  Boundary_Condition_Type_List = CArrayKokkos<int, array_layout, HostSpace, memory_traits>(num_sets, "Boundary_Condition_Type_List");
  NBoundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, "NBoundary_Condition_Patches");
  //std::cout << "NBOUNDARY PATCHES ON RANK " << myrank << " FOR INIT IS " << nboundary_patches <<std::endl;
  Boundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, nboundary_patches, "Boundary_Condition_Patches");

  //initialize data
  for(int iset = 0; iset < num_sets; iset++) NBoundary_Condition_Patches(iset) = 0;

   //initialize
  for(int ibdy=0; ibdy < num_sets; ibdy++) Boundary_Condition_Type_List(ibdy) = NONE;
}

/* ----------------------------------------------------------------------------
   Grow boundary conditions sets of element boundary surfaces
------------------------------------------------------------------------------- */

void FEA_Module_SGH::grow_boundary_sets(int num_sets){
  int num_dim = simparam->num_dim;

  if(num_sets == 0){
    std::cout << " Warning: number of boundary conditions being set to 0";
    return;
  }

  //std::cout << " DEBUG PRINT "<<num_sets << " " << nboundary_patches << std::endl;
  if(num_sets>max_boundary_sets){
    //temporary storage for previous data
    CArrayKokkos<int, array_layout, HostSpace, memory_traits> Temp_Boundary_Condition_Type_List = Boundary_Condition_Type_List;
    CArrayKokkos<size_t, array_layout, device_type, memory_traits> Temp_NBoundary_Condition_Patches = NBoundary_Condition_Patches;
    CArrayKokkos<size_t, array_layout, device_type, memory_traits> Temp_Boundary_Condition_Patches = Boundary_Condition_Patches;
    
    max_boundary_sets = num_sets + 5; //5 is an arbitrary buffer
    Boundary_Condition_Type_List = CArrayKokkos<int, array_layout, HostSpace, memory_traits>(max_boundary_sets, "Boundary_Condition_Type_List");
    NBoundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(max_boundary_sets, "NBoundary_Condition_Patches");
    //std::cout << "NBOUNDARY PATCHES ON RANK " << myrank << " FOR GROW " << nboundary_patches <<std::endl;
    Boundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(max_boundary_sets, nboundary_patches, "Boundary_Condition_Patches");

    //copy previous data back over
    //std::cout << "NUM BOUNDARY CONDITIONS ON RANK " << myrank << " FOR COPY " << max_boundary_sets <<std::endl;
    for(int iset = 0; iset < num_boundary_conditions; iset++){
      Boundary_Condition_Type_List(iset) = Temp_Boundary_Condition_Type_List(iset);
      NBoundary_Condition_Patches(iset) = Temp_NBoundary_Condition_Patches(iset);
      for(int ipatch = 0; ipatch < nboundary_patches; ipatch++){
        Boundary_Condition_Patches(iset, ipatch) = Temp_Boundary_Condition_Patches(iset, ipatch);
      }
    }
    
    //initialize data
    for(int iset = num_boundary_conditions; iset < max_boundary_sets; iset++) NBoundary_Condition_Patches(iset) = 0;

    //initialize
    for(int ibdy = num_boundary_conditions; ibdy < max_boundary_sets; ibdy++) Boundary_Condition_Type_List(ibdy) = NONE;
  }
}

/* ----------------------------------------------------------------------
   Assign sets of element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void FEA_Module_SGH::generate_bcs(){
  
} // end generate_bcs


/* ----------------------------------------------------------------------
   Loop through applied boundary conditions and tag node ids to remove 
   necessary rows and columns from the assembled linear system
------------------------------------------------------------------------- */

void FEA_Module_SGH::Displacement_Boundary_Conditions(){
 
}

/* ----------------------------------------------------------------------------
   Initialize output data structures
------------------------------------------------------------------------------- */

void FEA_Module_SGH::init_output(){
  //check user parameters for output
  bool output_displacement_flag = simparam->output_displacement_flag;
  displaced_mesh_flag = simparam->displaced_mesh_flag;
  bool output_strain_flag = simparam->output_strain_flag;
  bool output_stress_flag = simparam->output_stress_flag;
  int num_dim = simparam->num_dim;
  int Brows;
  if(num_dim==3) Brows = 6;
  else Brows = 3;

  if(output_displacement_flag){
    //displacement_index is accessed by writers at the solver level for deformed output
    displacement_index = output_displacement_index = noutput;
    noutput += 1;
    module_outputs.resize(noutput);

    vector_style.resize(noutput);
    vector_style[noutput-1] = DOF;

    output_vector_sizes.resize(noutput);
    output_vector_sizes[noutput-1] = num_dim;

    output_dof_names.resize(noutput);
    output_dof_names[noutput-1].resize(num_dim);
    output_dof_names[noutput-1][0] = "vx";
    output_dof_names[noutput-1][1] = "vy";
    output_dof_names[noutput-1][2] = "vz";
  }
  if(output_strain_flag){
    output_strain_index = noutput;
    noutput += 1;
    module_outputs.resize(noutput);

    vector_style.resize(noutput);
    vector_style[noutput-1] = NODAL;

    output_vector_sizes.resize(noutput);
    output_vector_sizes[noutput-1] = Brows;

    output_dof_names.resize(noutput);
    output_dof_names[noutput-1].resize(Brows);
    output_dof_names[noutput-1][0] = "strain_xx";
    output_dof_names[noutput-1][1] = "strain_yy";
    output_dof_names[noutput-1][2] = "strain_zz";
    output_dof_names[noutput-1][3] = "strain_xy";
    output_dof_names[noutput-1][4] = "strain_xz";
    output_dof_names[noutput-1][5] = "strain_yz";
  }
  if(output_stress_flag){
    output_stress_index = noutput;
    noutput += 1;
    module_outputs.resize(noutput);

    vector_style.resize(noutput);
    vector_style[noutput-1] = NODAL;

    output_vector_sizes.resize(noutput);
    output_vector_sizes[noutput-1] = Brows;

    output_dof_names.resize(noutput);
    output_dof_names[noutput-1].resize(Brows);
    output_dof_names[noutput-1][0] = "stress_xx";
    output_dof_names[noutput-1][1] = "stress_yy";
    output_dof_names[noutput-1][2] = "stress_zz";
    output_dof_names[noutput-1][3] = "stress_xy";
    output_dof_names[noutput-1][4] = "stress_xz";
    output_dof_names[noutput-1][5] = "stress_yz";
  }
}

/* -------------------------------------------------------------------------------------------
   Prompts sorting for elastic response output data. For now, nodal strains.
---------------------------------------------------------------------------------------------- */

void FEA_Module_SGH::sort_output(Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > sorted_map){
 
  
}

/* -------------------------------------------------------------------------------------------
   Prompts computation of elastic response output data. For now, nodal strains.
---------------------------------------------------------------------------------------------- */

void FEA_Module_SGH::collect_output(Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > global_reduce_map){
 
}

/* -------------------------------------------------------------------------------------------
   Prompts computation of elastic response output data. For now, nodal strains.
---------------------------------------------------------------------------------------------- */

void FEA_Module_SGH::compute_output(){

}

/* ----------------------------------------------------------------------
   Solve the FEA linear system
------------------------------------------------------------------------- */

int FEA_Module_SGH::solve(){
  //local variable for host view in the dual view
  int num_dim = simparam->num_dim;
  int nodes_per_elem = max_nodes_per_element;
  int local_node_index, current_row, current_column;
  int max_stride = 0;
  size_t access_index, row_access_index, row_counter;
  GO global_index, global_dof_index;
  LO local_dof_index;

 
  
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------------------------
   Communicate ghosts using the current optimization design data
---------------------------------------------------------------------------------------------- */

void FEA_Module_SGH::comm_variables(Teuchos::RCP<const MV> zp){

}


/* -------------------------------------------------------------------------------------------
   enforce constraints on nodes due to BCS
---------------------------------------------------------------------------------------------- */

void FEA_Module_SGH::node_density_constraints(host_vec_array node_densities_lower_bound){

}