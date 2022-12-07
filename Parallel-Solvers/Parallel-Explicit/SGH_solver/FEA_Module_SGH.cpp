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
#include <chrono>
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
#include "Simulation_Parameters_SGH.h"
#include "FEA_Module_SGH.h"
#include "Explicit_Solver_SGH.h"

#define MAX_ELEM_NODES 8
#define STRAIN_EPSILON 0.000000001
#define DENSITY_EPSILON 0.0001
#define BC_EPSILON 1.0e-6

using namespace utils;


FEA_Module_SGH::FEA_Module_SGH(Solver *Solver_Pointer, mesh_t& mesh) :FEA_Module(Solver_Pointer), mesh(mesh){
  //create parameter object
  //recast solver pointer for non-base class access
  Explicit_Solver_Pointer_ = dynamic_cast<Explicit_Solver_SGH*>(Solver_Pointer);
  simparam = dynamic_cast<Simulation_Parameters_SGH*>(Explicit_Solver_Pointer_->simparam);
  // ---- Read input file, define state and boundary conditions ---- //
  //simparam->input();
  
  //TO parameters
  //simparam_TO = dynamic_cast<Simulation_Parameters_Topology_Optimization*>(Explicit_Solver_Pointer_->simparam);

  //create ref element object
  //ref_elem = new elements::ref_element();
  //create mesh objects
  //init_mesh = new swage::mesh_t(simparam);
  //mesh = new swage::mesh_t(simparam);

  //boundary condition data
  max_boundary_sets = 0;
  
  //setup output
  noutput = 0;
  init_output();

}

FEA_Module_SGH::~FEA_Module_SGH(){
   //delete simparam;
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

// -----------------------------------------------------------------------------
// Interfaces read in data with the SGH solver data; currently a hack to streamline
//------------------------------------------------------------------------------
void FEA_Module_SGH::sgh_interface_setup(mesh_t &mesh,
                       node_t &node,
                       elem_t &elem,
                       corner_t &corner){

    const size_t rk_level = 0;
    const size_t num_dims = simparam->num_dim;
    const size_t rk_num_bins = simparam->rk_num_bins;

    size_t num_nodes_in_elem = 1;
    for (int dim=0; dim<num_dims; dim++){
        num_nodes_in_elem *= 2;
    }

    // --- Read in the nodes in the mesh ---

    size_t num_nodes = Explicit_Solver_Pointer_->nall_nodes;
    int myrank = Explicit_Solver_Pointer_->myrank;
    int nranks = Explicit_Solver_Pointer_->nranks;
    //printf("Num nodes assigned to MPI rank %lu is %lu\n" , myrank, num_nodes);

    // intialize node variables
    mesh.initialize_nodes(num_nodes);
    mesh.initialize_local_nodes(Explicit_Solver_Pointer_->nlocal_nodes);
    node.initialize(rk_num_bins, num_nodes, num_dims);
    //std::cout << "Bin counts " << rk_num_bins << " Node counts " << num_nodes << " Num dim " << num_dims << std::endl;

    //view scope
    {
      Explicit_Solver_SGH::host_vec_array interface_node_coords = Explicit_Solver_Pointer_->all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
      //save node data to node.coords
      //std::cout << "NODE DATA ON RANK " << myrank << std::endl;
      if(num_dims==2){
        for(int inode = 0; inode < num_nodes; inode++){
          //std::cout << "Node index " << inode+1 << " ";
          node.coords.host(0,inode,0) = interface_node_coords(inode,0);
          //std::cout << host_node_coords_state(0,inode,0)+1<< " ";
          node.coords.host(0,inode,1) = interface_node_coords(inode,1);
          //std::cout << host_node_coords_state(0,inode,1)+1<< " ";
        }
      }
      else if(num_dims==3){
        for(int inode = 0; inode < num_nodes; inode++){
          //std::cout << "Node index " << inode+1 << " ";
          node.coords.host(0,inode,0) = interface_node_coords(inode,0);
          //std::cout << host_node_coords_state(0,inode,0)+1<< " ";
          node.coords.host(0,inode,1) = interface_node_coords(inode,1);
          //std::cout << host_node_coords_state(0,inode,1)+1<< " ";
        
          node.coords.host(0,inode,2) = interface_node_coords(inode,2);
          //std::cout << host_node_coords_state(0,inode,2)+1<< std::endl;
        }
      }
    } //end view scope
    // --- read in the elements in the mesh ---
    size_t num_elem = 0;
    
    num_elem = Explicit_Solver_Pointer_->rnum_elem;
    //printf("Num elems assigned to MPI rank %lu is %lu\n" , myrank, num_elem);

    // intialize elem variables
    mesh.initialize_elems(num_elem, num_dims);
    elem.initialize(rk_num_bins, num_nodes, 3); // always 3D here, even for 2D

    //save data to mesh.nodes_in_elem.host
    //CArrayKokkos<size_t, DefaultLayout, HostSpace> host_mesh_nodes_in_elem(num_elem, num_nodes_in_elem);
    //view scope
    {
      Explicit_Solver_SGH::host_elem_conn_array interface_nodes_in_elem = Explicit_Solver_Pointer_->nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
      //save node data to node.coords
      //std::cout << "ELEMENT CONNECTIVITY ON RANK " << myrank << std::endl;
      for(int ielem = 0; ielem < num_elem; ielem++){
        //std::cout << "Element index " << ielem+1 << " ";
        for(int inode = 0; inode < num_nodes_in_elem; inode++){
            mesh.nodes_in_elem.host(ielem,inode) = Explicit_Solver_Pointer_->all_node_map->getLocalElement(interface_nodes_in_elem(ielem,inode));
            //debug print
            //std::cout << mesh.nodes_in_elem.get_kokkos_dual_view().h_view(ielem*num_nodes_in_elem + inode)+1<< " ";
        }
        //std::cout << std::endl;
      }
    }
    // update device side
    mesh.nodes_in_elem.update_device();

    //debug print
    
    //CArrayKokkos<size_t> device_mesh_nodes_in_elem(num_elem, num_nodes_in_elem);
    //device_mesh_nodes_in_elem.get_kokkos_view() = mesh.nodes_in_elem.get_kokkos_dual_view().d_view;
    //host_mesh_nodes_in_elem.get_kokkos_view() = mesh.nodes_in_elem.get_kokkos_dual_view().view_host();
    /*
    if(myrank==1){
    std::cout << "ELEMENT CONNECTIVITY ON RANK 1 in LOCAL INDICES" << myrank << std::endl;
    for(int ielem = 0; ielem < num_elem; ielem++){
        std::cout << "Element index " << ielem+1 << " ";
        for(int inode = 0; inode < num_nodes_in_elem; inode++){
            //debug print
            //device_mesh_nodes_in_elem(ielem,inode) = Explicit_Solver_Pointer_->all_node_map->getLocalElement(interface_nodes_in_elem(ielem,inode));
            std::cout << mesh.nodes_in_elem(ielem, inode)+1<< " ";
        }
        std::cout << std::endl;
    }
    }
    */
    /*
    std::cout.flush();
    if(myrank==1){
    std::cout << "ELEMENT CONNECTIVITY ON RANK 1 in GLOBAL INDICES" << myrank << std::endl;
    std::cout << "local node index of global index 275 on rank 1 " << Explicit_Solver_Pointer_->all_node_map->getLocalElement(275) << std::endl;
    for(int ielem = 0; ielem < num_elem; ielem++){
        std::cout << ielem << " ";
        for(int inode = 0; inode < num_nodes_in_elem; inode++){
            //debug print
            //device_mesh_nodes_in_elem(ielem,inode) = Explicit_Solver_Pointer_->all_node_map->getLocalElement(interface_nodes_in_elem(ielem,inode));
            std::cout << Explicit_Solver_Pointer_->all_node_map->getGlobalElement(mesh.nodes_in_elem(ielem, inode))<< " ";
        }
        std::cout << std::endl;
    }
    }
    std::cout.flush();
    */
    /*
    size_t nall_nodes = Explicit_Solver_Pointer_->nall_nodes;
    node.all_coords = DCArrayKokkos <double> (rk_num_bins, nall_nodes, num_dims);
    node.all_vel    = DCArrayKokkos <double> (rk_num_bins, nall_nodes, num_dims);
    node.all_mass   = DCArrayKokkos <double> (nall_nodes);

    //save all data (nlocal +nghost)
    CArrayKokkos<double, DefaultLayout, HostSpace> host_all_node_coords_state(rk_num_bins, nall_nodes, num_dims);
    Explicit_Solver_SGH::host_vec_array interface_all_node_coords = Explicit_Solver_Pointer_->all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    host_all_node_coords_state.get_kokkos_view() = node.all_coords.get_kokkos_dual_view().view_host();
    //host_node_coords_state = CArrayKokkos<double, DefaultLayout, HostSpace>(rk_num_bins, num_nodes, num_dims);
    //host_all_node_coords_state.get_kokkos_view() = Kokkos::View<double*,DefaultLayout, HostSpace>("debug", rk_num_bins*nall_nodes*num_dims);
    //save node data to node.coords
    
    //std::cout << "ALL NODE DATA ON RANK " << myrank << std::endl;
    for(int inode = 0; inode < nall_nodes; inode++){
        //std::cout << "Node index " << inode+1 << " ";
        node.all_coords.host(0,inode,0) = interface_all_node_coords(inode,0);
        //std::cout << host_all_node_coords_state(0,inode,0)+1<< " ";
        node.all_coords.host(0,inode,1) = interface_all_node_coords(inode,1);
        //std::cout << host_all_node_coords_state(0,inode,1)+1<< " ";
        node.all_coords.host(0,inode,2) = interface_all_node_coords(inode,2);
        //std::cout << host_all_node_coords_state(0,inode,2)+1<< std::endl;
    }
    */

    // save the node coords to the current RK value
    for (size_t node_gid=0; node_gid<num_nodes; node_gid++){
        
        for(int rk=1; rk<rk_num_bins; rk++){
            for (int dim = 0; dim < num_dims; dim++){
                node.coords.host(rk, node_gid, dim) = node.coords.host(0, node_gid, dim);
            } // end for dim
        } // end for rk
        
    } // end parallel for
    
    /*
    // save the node coords to the current RK value
    for (size_t node_gid=0; node_gid<nall_nodes; node_gid++){
        
        for(int rk=1; rk<rk_num_bins; rk++){
            for (int dim = 0; dim < num_dims; dim++){
                node.all_coords.host(rk, node_gid, dim) = node.all_coords.host(0, node_gid, dim);
            } // end for dim
        } // end for rk
        
    } // end parallel for
    */
    
    node.coords.update_device();
    //node.all_coords.update_device();

    
    // intialize corner variables
    int num_corners = num_elem*mesh.num_nodes_in_elem;
    mesh.initialize_corners(num_corners);
    corner.initialize(num_corners, num_dims);

    //set Tpetra vector pointers
    node_velocities_distributed = Explicit_Solver_Pointer_->node_velocities_distributed;
    all_node_velocities_distributed = Explicit_Solver_Pointer_->all_node_velocities_distributed;
    all_cached_node_velocities_distributed = Teuchos::rcp(new MV(all_node_map, num_dims));
    
    /*
    for(int inode = 0; inode < num_nodes; inode++){
        std::cout << "Node index " << inode+1 << " ";
        for(int rk=0; rk<rk_num_bins; rk++){
          std::cout << "rk index " << rk+1 << " ";
          std::cout << node.coords(rk,inode,0)+1<< " ";
          std::cout << node.coords(rk,inode,1)+1<< " ";
          std::cout << node.coords(rk,inode,2)+1<< std::endl;
        }
    }
    */
    // Close mesh input file
    //fclose(in);

    return;
    
}

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
  bool output_velocity_flag = simparam->output_velocity_flag;
  displaced_mesh_flag = simparam->displaced_mesh_flag;
  bool output_strain_flag = simparam->output_strain_flag;
  bool output_stress_flag = simparam->output_stress_flag;
  int num_dim = simparam->num_dim;
  int Brows;
  if(num_dim==3) Brows = 6;
  else Brows = 3;

  if(output_velocity_flag){
    //displacement_index is accessed by writers at the solver level for deformed output
    output_velocity_index = noutput;
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
    if(num_dim==3)
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
    if(num_dim==2){
      output_dof_names[noutput-1][0] = "strain_xx";
      output_dof_names[noutput-1][1] = "strain_yy";
      output_dof_names[noutput-1][2] = "strain_xy";
    }
    if(num_dim==3){
      output_dof_names[noutput-1][0] = "strain_xx";
      output_dof_names[noutput-1][1] = "strain_yy";
      output_dof_names[noutput-1][2] = "strain_zz";
      output_dof_names[noutput-1][3] = "strain_xy";
      output_dof_names[noutput-1][4] = "strain_xz";
      output_dof_names[noutput-1][5] = "strain_yz";
    }
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
    if(num_dim==2){
      output_dof_names[noutput-1][0] = "stress_xx";
      output_dof_names[noutput-1][1] = "stress_yy";
      output_dof_names[noutput-1][3] = "stress_xy";
    }
    if(num_dim==3){
      output_dof_names[noutput-1][0] = "stress_xx";
      output_dof_names[noutput-1][1] = "stress_yy";
      output_dof_names[noutput-1][2] = "stress_zz";
      output_dof_names[noutput-1][3] = "stress_xy";
      output_dof_names[noutput-1][4] = "stress_xz";
      output_dof_names[noutput-1][5] = "stress_yz";
    }
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

void FEA_Module_SGH::update_forward_solve(Teuchos::RCP<const MV> zp){
  //local variable for host view in the dual view
  int num_dim = simparam->num_dim;
  int nodes_per_elem = max_nodes_per_element;
  int local_node_index, current_row, current_column;
  int max_stride = 0;
  size_t access_index, row_access_index, row_counter;
  GO global_index, global_dof_index;
  LO local_dof_index;
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

/* ----------------------------------------------------------------------------
   Setup SGH solver data
------------------------------------------------------------------------------- */

void FEA_Module_SGH::setup(){
    
    const size_t num_fills = simparam->num_fills;
    const size_t rk_num_bins = simparam->rk_num_bins;
    const size_t num_bcs = simparam->num_bcs;
    const size_t num_materials = simparam->num_materials;
    const size_t num_state_vars = simparam->max_num_state_vars;

    const CArrayKokkos <mat_fill_t> mat_fill = simparam->mat_fill;
    const CArrayKokkos <boundary_t> boundary = simparam->boundary;
    const CArrayKokkos <material_t> material = simparam->material;
    const CArrayKokkos <double> state_vars = simparam->state_vars; // array to hold init model variables
    
    //--- calculate bdy sets ---//
    mesh.num_nodes_in_patch = 2*(mesh.num_dims-1);  // 2 (2D) or 4 (3D)
    mesh.num_patches_in_elem = 2*mesh.num_dims; // 4 (2D) or 6 (3D)
    mesh.init_bdy_sets(num_bcs);
    printf("Num BC's = %lu\n", num_bcs);
    
    // tag boundary patches in the set
    tag_bdys(boundary, mesh, node_coords);

    build_boundry_node_sets(boundary, mesh);
    
    // loop over BCs
    for (size_t this_bdy = 0; this_bdy < num_bcs; this_bdy++){
        
        RUN_CLASS({
            printf("Boundary Condition number %lu \n", this_bdy);
            printf("  Num bdy patches in this set = %lu \n", mesh.bdy_patches_in_set.stride(this_bdy));
            printf("  Num bdy nodes in this set = %lu \n", mesh.bdy_nodes_in_set.stride(this_bdy));
        });
        Kokkos::fence();

    }// end for
    
    
    // ---- Read model values from a file ----
    // check to see if state_vars come from an external file
    DCArrayKokkos <size_t> read_from_file(num_materials, "read_from_file");
    FOR_ALL_CLASS(mat_id, 0, num_materials, {
        
        read_from_file(mat_id) = material(0).read_state_vars;
        
    }); // end parallel for
    Kokkos::fence();
    
    read_from_file.update_host(); // copy to CPU if code is to read from a file
    Kokkos::fence();
    
    // make memory to store state_vars from an external file
    DCArrayKokkos <double> file_state_vars(num_materials,mesh.num_elems,num_state_vars);
    DCArrayKokkos <size_t> mat_num_state_vars(num_materials); // actual number of state_vars
    FOR_ALL_CLASS(mat_id, 0, num_materials, {
        
        mat_num_state_vars(mat_id) = material(mat_id).num_state_vars;
        
    }); // end parallel for
    Kokkos::fence();
    
    // copy actual number of state_vars to host
    mat_num_state_vars.update_host();
    Kokkos::fence();
    
    for (size_t mat_id=0; mat_id<num_materials; mat_id++){
        
        if (read_from_file.host(mat_id) == 1){
            
            size_t num_vars = mat_num_state_vars.host(mat_id);
            
            user_model_init(file_state_vars,
                            num_vars,
                            mat_id,
                            mesh.num_elems);
            
            // copy the values to the device
            file_state_vars.update_device();
            Kokkos::fence();
            
        } // end if
        
    } // end for
    
    
    //--- apply the fill instructions over the Elements---//
    
    // loop over the fill instructures
    for (int f_id = 0; f_id < num_fills; f_id++){
            
        // parallel loop over elements in mesh
        FOR_ALL_CLASS(elem_gid, 0, mesh.num_elems, {

            const size_t rk_level = 1;

            // calculate the coordinates and radius of the element
            double elem_coords[3]; // note:initialization with a list won't work
            elem_coords[0] = 0.0;
            elem_coords[1] = 0.0;
            elem_coords[2] = 0.0;

            // get the coordinates of the element center
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                elem_coords[0] += node_coords(rk_level, mesh.nodes_in_elem(elem_gid, node_lid), 0);
                elem_coords[1] += node_coords(rk_level, mesh.nodes_in_elem(elem_gid, node_lid), 1);
                if (mesh.num_dims == 3){
                    elem_coords[2] += node_coords(rk_level, mesh.nodes_in_elem(elem_gid, node_lid), 2);
                } else
                {
                    elem_coords[2] = 0.0;
                }
            } // end loop over nodes in element
            elem_coords[0] = elem_coords[0]/mesh.num_nodes_in_elem;
            elem_coords[1] = elem_coords[1]/mesh.num_nodes_in_elem;
            elem_coords[2] = elem_coords[2]/mesh.num_nodes_in_elem;
                
            
            // spherical radius
            double radius = sqrt( elem_coords[0]*elem_coords[0] +
                                  elem_coords[1]*elem_coords[1] +
                                  elem_coords[2]*elem_coords[2] );
                
            // cylinderical radius
            double radius_cyl = sqrt( elem_coords[0]*elem_coords[0] +
                                      elem_coords[1]*elem_coords[1] );   
            
            // default is not to fill the element
            size_t fill_this = 0;
           
            // check to see if this element should be filled
            switch(mat_fill(f_id).volume)
            {
                case region::global:
                {
                    fill_this = 1;
                    break;
                }
                case region::box:
                {
                    if ( elem_coords[0] >= mat_fill(f_id).x1 && elem_coords[0] <= mat_fill(f_id).x2
                      && elem_coords[1] >= mat_fill(f_id).y1 && elem_coords[1] <= mat_fill(f_id).y2
                      && elem_coords[2] >= mat_fill(f_id).z1 && elem_coords[2] <= mat_fill(f_id).z2 )
                        fill_this = 1;
                    break;
                }
                case region::cylinder:
                {
                    if ( radius_cyl >= mat_fill(f_id).radius1
                      && radius_cyl <= mat_fill(f_id).radius2 ) fill_this = 1;
                    break;
                }
                case region::sphere:
                {
                    if ( radius >= mat_fill(f_id).radius1
                      && radius <= mat_fill(f_id).radius2 ) fill_this = 1;
                    break;
                }
            } // end of switch

                 
            // paint the material state on the element
            if (fill_this == 1){
                    
                // density
                elem_den(elem_gid) = mat_fill(f_id).den;
                
                // mass
                elem_mass(elem_gid) = elem_den(elem_gid)*elem_vol(elem_gid);
                
                // specific internal energy
                elem_sie(rk_level, elem_gid) = mat_fill(f_id).sie;
		
                elem_mat_id(elem_gid) = mat_fill(f_id).mat_id;
                size_t mat_id = elem_mat_id(elem_gid); // short name
                
                
                // get state_vars from the input file or read them in
                if (material(mat_id).read_state_vars == 1){
                    
                    // use the values read from a file to get elem state vars
                    for (size_t var=0; var<material(mat_id).num_state_vars; var++){
                        elem_statev(elem_gid,var) = file_state_vars(mat_id,elem_gid,var);
                    } // end for
                    
                }
                else{
                    // use the values in the input file
                    // set state vars for the region where mat_id resides
                    for (size_t var=0; var<material(mat_id).num_state_vars; var++){
                        elem_statev(elem_gid,var) = state_vars(mat_id,var);
                    } // end for
                    
                } // end logical on type
                
                // --- stress tensor ---
                // always 3D even for 2D-RZ
                for (size_t i=0; i<3; i++){
                    for (size_t j=0; j<3; j++){
                        elem_stress(rk_level,elem_gid,i,j) = 0.0;
                    }        
                }  // end for
                
                
                
                // --- Pressure and stress ---
                material(mat_id).eos_model(elem_pres,
                                           elem_stress,
                                           elem_gid,
                                           elem_mat_id(elem_gid),
                                           elem_statev,
                                           elem_sspd,
                                           elem_den(elem_gid),
                                           elem_sie(1,elem_gid));
					    
                
                // loop over the nodes of this element and apply velocity
                for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){

                    // get the mesh node index
                    size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

                
                    // --- Velocity ---
                    switch(mat_fill(f_id).velocity)
                    {
                        case init_conds::cartesian:
                        {
                        
                            node_vel(rk_level, node_gid, 0) = mat_fill(f_id).u;
                            node_vel(rk_level, node_gid, 1) = mat_fill(f_id).v;
                            if (mesh.num_dims == 3) node_vel(rk_level, node_gid, 2) = mat_fill(f_id).w;
                            
                        
                            break;
                        }
                        case init_conds::radial:
                        {
                            // Setting up cylindrical
                            double dir[2]; 
                            dir[0] = 0.0;
                            dir[1] = 0.0;
                            double radius_val = 0.0;
                        
                            for(int dim=0; dim<2; dim++){
                                dir[dim] = node_coords(rk_level, node_gid, dim);
                                radius_val += node_coords(rk_level, node_gid, dim)*node_coords(rk_level, node_gid, dim);
                            } // end for
                            radius_val = sqrt(radius_val);
                        
                            for(int dim=0; dim<2; dim++){
                                if (radius_val > 1.0e-14){
                                    dir[dim] /= (radius_val);
                                }
                                else{
                                    dir[dim] = 0.0;
                                }
                            } // end for
                        
                        
                            node_vel(rk_level, node_gid, 0) = mat_fill(f_id).speed*dir[0];
                            node_vel(rk_level, node_gid, 1) = mat_fill(f_id).speed*dir[1];
                            if (mesh.num_dims == 3) node_vel(rk_level, node_gid, 2) = 0.0;
                            
                            break;
                        }
                        case init_conds::spherical:
                        {
                            
                            // Setting up spherical
                            double dir[3];
                            dir[0] = 0.0;
                            dir[1] = 0.0;
                            dir[2] = 0.0;
                            double radius_val = 0.0;
                        
                            for(int dim=0; dim<3; dim++){
                                dir[dim] = node_coords(rk_level, node_gid, dim);
                                radius_val += node_coords(rk_level, node_gid, dim)*node_coords(rk_level, node_gid, dim);
                            } // end for
                            radius_val = sqrt(radius_val);
                        
                            for(int dim=0; dim<3; dim++){
                                if (radius_val > 1.0e-14){
                                    dir[dim] /= (radius_val);
                                }
                                else{
                                    dir[dim] = 0.0;
                                }
                            } // end for
                        
                            node_vel(rk_level, node_gid, 0) = mat_fill(f_id).speed*dir[0];
                            node_vel(rk_level, node_gid, 1) = mat_fill(f_id).speed*dir[1];
                            if (mesh.num_dims == 3) node_vel(rk_level, node_gid, 2) = mat_fill(f_id).speed*dir[2];

                            break;
                        }
                        case init_conds::radial_linear:
                        {
                        
                            break;
                        }
                        case init_conds::spherical_linear:
                        {
                        
                            break;
                        }
                        case init_conds::tg_vortex:
                        {
                        
                            node_vel(rk_level, node_gid, 0) = sin(PI * node_coords(rk_level,node_gid, 0)) * cos(PI * node_coords(rk_level,node_gid, 1)); 
                            node_vel(rk_level, node_gid, 1) =  -1.0*cos(PI * node_coords(rk_level,node_gid, 0)) * sin(PI * node_coords(rk_level,node_gid, 1)); 
                            if (mesh.num_dims == 3) node_vel(rk_level, node_gid, 2) = 0.0;

                            break;
                        }
                    } // end of switch

                }// end loop over nodes of element
                
                
                if(mat_fill(f_id).velocity == init_conds::tg_vortex)
                {
                    elem_pres(elem_gid) = 0.25*( cos(2.0*PI*elem_coords[0]) + cos(2.0*PI*elem_coords[1]) ) + 1.0;
                
                    // p = rho*ie*(gamma - 1)
                    size_t mat_id = f_id;
                    double gamma = elem_statev(elem_gid,4); // gamma value
                    elem_sie(rk_level, elem_gid) =
                                    elem_pres(elem_gid)/(mat_fill(f_id).den*(gamma - 1.0));
                } // end if

            } // end if fill
          
        }); // end FOR_ALL_CLASS element loop
        Kokkos::fence();
        
  
    } // end for loop over fills
    
   
    
    // apply BC's to velocity
    FEA_Module_SGH::boundary_velocity(mesh, boundary, node_vel);
    
    
    // calculate the corner massess if 2D
    if(mesh.num_dims==2){
        
        FOR_ALL_CLASS(elem_gid, 0, mesh.num_elems, {
            
            // facial area of the corners
            double corner_areas_array[4];
            
            ViewCArrayKokkos <double> corner_areas(&corner_areas_array[0],4);
            ViewCArrayKokkos <size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 4);
            
            get_area_weights2D(corner_areas,
                               elem_gid,
                               node_coords,
                               elem_node_gids);
            
            // loop over the corners of the element and calculate the mass
            for (size_t corner_lid=0; corner_lid<4; corner_lid++){
                
                size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);
                corner_mass(corner_gid) = corner_areas(corner_lid)*elem_den(elem_gid); // node radius is added later
                
            } // end for over corners
        });
    
    } // end of
    
    
    // calculate the nodal mass
    FOR_ALL_CLASS(node_gid, 0, mesh.num_nodes, {
        
        node_mass(node_gid) = 0.0;
        
        if(mesh.num_dims==3){
            
            for(size_t elem_lid=0; elem_lid<mesh.num_corners_in_node(node_gid); elem_lid++){
                size_t elem_gid = mesh.elems_in_node(node_gid,elem_lid);
                node_mass(node_gid) += 1.0/8.0*elem_mass(elem_gid);
            } // end for elem_lid
            
        }// end if dims=3
        else {
            
            // 2D-RZ
            for(size_t corner_lid=0; corner_lid<mesh.num_corners_in_node(node_gid); corner_lid++){
                
                size_t corner_gid = mesh.corners_in_node(node_gid, corner_lid);
                node_mass(node_gid) += corner_mass(corner_gid);  // sans the radius so it is areal node mass
                
                corner_mass(corner_gid) *= node_coords(1,node_gid,1); // true corner mass now
            } // end for elem_lid
            
        } // end else
        
    }); // end FOR_ALL_CLASS
    Kokkos::fence();
    
    return;
    
} // end of setup


/* ----------------------------------------------------------------------------
    set planes for tagging sub sets of boundary patches
    bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
    val = plane value, cyl radius, sphere radius
------------------------------------------------------------------------------- */

void FEA_Module_SGH::tag_bdys(const CArrayKokkos <boundary_t> &boundary,
              mesh_t &mesh,
              const DViewCArrayKokkos <double> &node_coords){

    size_t num_dims = mesh.num_dims;
    //int nboundary_patches = Explicit_Solver_Pointer_->nboundary_patches;
    int nboundary_patches = Explicit_Solver_Pointer_->nboundary_patches;
    
    //if (bdy_set == mesh.num_bdy_sets){
    //    printf(" ERROR: number of boundary sets must be increased by %zu",
    //              bdy_set-mesh.num_bdy_sets+1);
    //    exit(0);
    //} // end if

    //error and debug flag
    //DCArrayKokkos<bool> print_flag(1, "print_flag");
    //print_flag.host(0) = false;
    //print_flag.update_device();
    
    FOR_ALL_CLASS(bdy_set, 0, mesh.num_bdy_sets, {
        
        // tag boundaries
        int bc_tag_id = boundary(bdy_set).surface;
        double val = boundary(bdy_set).value;
        
        // save the boundary patches to this set that are on the plane, spheres, etc.
        for (size_t bdy_patch_lid=0; bdy_patch_lid < nboundary_patches; bdy_patch_lid++){
            
            // save the patch index
            size_t bdy_patch_gid = bdy_patch_lid;
            
            
            // check to see if this patch is on the specified plane
            size_t is_on_bdy = check_bdy(bdy_patch_gid,
                                         bc_tag_id,
                                         val,
                                         mesh,
                                         node_coords); // no=0, yes=1
            
            //debug check
            /*
            for (size_t patch_node_lid=0; patch_node_lid<mesh.num_nodes_in_patch; patch_node_lid++){
              size_t node_gid = mesh.nodes_in_patch(bdy_patch_gid, patch_node_lid);
              //if(bdy_node_gid==549412) print_flag(0) = true;
            }
            */

            if (is_on_bdy == 1){
                
                size_t index = mesh.bdy_patches_in_set.stride(bdy_set);
                
                // increment the number of boundary patches saved
                mesh.bdy_patches_in_set.stride(bdy_set) ++;
                
                
                mesh.bdy_patches_in_set(bdy_set, index) = bdy_patch_gid;
            } // end if
            
            
        } // end for bdy_patch
        
    });  // end FOR_ALL_CLASS bdy_sets
    
    //debug check
    //print_flag.update_host();
    //if(print_flag.host(0)) std::cout << "found boundary node with id 549412" << std::endl;

    return;
} // end tag


/* ----------------------------------------------------------------------------
    routine for checking to see if a vertex is on a boundary
    bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
    val = plane value, radius, radius
------------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
size_t FEA_Module_SGH::check_bdy(const size_t patch_gid,
                 const int this_bc_tag,
                 const double val,
                 const mesh_t &mesh,
                 const DViewCArrayKokkos <double> &node_coords) const {
    
    
    size_t num_dims = mesh.num_dims;
    
    // default bool is not on the boundary
    size_t is_on_bdy = 0;
    
    // the patch coordinates
    double these_patch_coords[3];  // Note: cannot allocated array with num_dims
    
    // loop over the nodes on the patch
    for (size_t patch_node_lid=0; patch_node_lid<mesh.num_nodes_in_patch; patch_node_lid++){
        
        // get the nodal_gid for this node in the patch
        //size_t node_gid = mesh.nodes_in_patch(patch_gid, patch_node_lid);
        size_t node_gid = mesh.Local_Index_Boundary_Patches(patch_gid, patch_node_lid);

        for (size_t dim = 0; dim < num_dims; dim++){
            these_patch_coords[dim] = node_coords(1, node_gid, dim);  // (rk, node_gid, dim)
        } // end for dim
        
        
        // a x-plane
        if (this_bc_tag == 0){
            
            if ( fabs(these_patch_coords[0] - val) <= 1.0e-7 ) is_on_bdy += 1;
            
        }// end if on type
        
        // a y-plane
        else if (this_bc_tag == 1){
            
            if ( fabs(these_patch_coords[1] - val) <= 1.0e-7 ) is_on_bdy += 1;
            
        }// end if on type
        
        // a z-plane
        else if (this_bc_tag == 2){
            
            if ( fabs(these_patch_coords[2] - val) <= 1.0e-7 ) is_on_bdy += 1;
            
        }// end if on type
        
        
        // cylinderical shell where radius = sqrt(x^2 + y^2)
        else if (this_bc_tag == 3){
            
            real_t R = sqrt(these_patch_coords[0]*these_patch_coords[0] +
                            these_patch_coords[1]*these_patch_coords[1]);
            
            if ( fabs(R - val) <= 1.0e-7 ) is_on_bdy += 1;
            
            
        }// end if on type
        
        // spherical shell where radius = sqrt(x^2 + y^2 + z^2)
        else if (this_bc_tag == 4){
            
            real_t R = sqrt(these_patch_coords[0]*these_patch_coords[0] +
                            these_patch_coords[1]*these_patch_coords[1] +
                            these_patch_coords[2]*these_patch_coords[2]);
            
            if ( fabs(R - val) <= 1.0e-7 ) is_on_bdy += 1;
            
        } // end if on type
        
    } // end for nodes in the patch
    
    // if all nodes in the patch are on the surface
    if (is_on_bdy == mesh.num_nodes_in_patch){
        is_on_bdy = 1;
    }
    else {
        is_on_bdy = 0;
    }
    
    
    return is_on_bdy;
    
} // end method to check bdy

/* ----------------------------------------------------------------------------
   Build set of nodes assigned to each boundary condition
------------------------------------------------------------------------------- */

void FEA_Module_SGH::build_boundry_node_sets(const CArrayKokkos <boundary_t> &boundary, mesh_t &mesh){
    
    // build boundary nodes in each boundary set
    int nboundary_patches = Explicit_Solver_Pointer_->nboundary_patches;
    mesh.num_bdy_nodes_in_set = DCArrayKokkos <size_t> (mesh.num_bdy_sets, "num_bdy_nodes_in_set");
    CArrayKokkos <long long int> temp_count_num_bdy_nodes_in_set(mesh.num_bdy_sets, mesh.num_nodes, "temp_count_num_bdy_nodes_in_set");
    
    DynamicRaggedRightArrayKokkos <size_t> temp_nodes_in_set (mesh.num_bdy_sets, nboundary_patches*mesh.num_nodes_in_patch, "temp_nodes_in_set");
    
    // Parallel loop over boundary sets on device
    FOR_ALL_CLASS(bdy_set, 0, mesh.num_bdy_sets, {
	
        // finde the number of patches_in_set
        size_t num_bdy_patches_in_set = mesh.bdy_patches_in_set.stride(bdy_set);

        mesh.num_bdy_nodes_in_set(bdy_set) = 0;
        
        // Loop over boundary patches in boundary set
        for (size_t bdy_patch_gid = 0; bdy_patch_gid<num_bdy_patches_in_set; bdy_patch_gid++){
            
                // get the global id for this boundary patch
                size_t patch_gid = mesh.bdy_patches_in_set(bdy_set, bdy_patch_gid);
                
                // apply boundary condition at nodes on boundary
                for(size_t node_lid = 0; node_lid < mesh.num_nodes_in_patch; node_lid++){
                    
                    size_t node_gid = mesh.Local_Index_Boundary_Patches(patch_gid, node_lid);
                    
                    temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) = -1;
                        
                } // end for node_lid
            
        } // end for bdy_patch_gid
        
        
        // Loop over boundary patches in boundary set
        for (size_t bdy_patch_gid = 0; bdy_patch_gid<num_bdy_patches_in_set; bdy_patch_gid++){
            
                // get the global id for this boundary patch
                size_t patch_gid = mesh.bdy_patches_in_set(bdy_set, bdy_patch_gid);
                
                // apply boundary condition at nodes on boundary
                for(size_t node_lid = 0; node_lid < mesh.num_nodes_in_patch; node_lid++){
                    
                    size_t node_gid = mesh.Local_Index_Boundary_Patches(patch_gid, node_lid);
                    
                    if (temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) == -1){
                        
                        size_t num_saved = mesh.num_bdy_nodes_in_set(bdy_set);
                        
                        mesh.num_bdy_nodes_in_set(bdy_set)++;
                        
                        // replace -1 with node_gid to denote the node was already saved
                        temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) = node_gid;
                        
                        // increment the number of saved nodes, create memory
                        temp_nodes_in_set.stride(bdy_set)++;
                        temp_nodes_in_set(bdy_set, num_saved) = node_gid;
                        
                    } // end if
                    
                } // end for node_lid
            
        } // end for bdy_patch_gid
        
    }); // end FOR_ALL_CLASS bdy_set
    Kokkos::fence();
    
   
    // allocate the RaggedRight bdy_nodes_in_set array
    mesh.bdy_nodes_in_set = RaggedRightArrayKokkos <size_t> (mesh.num_bdy_nodes_in_set, "bdy_nodes_in_set");

    FOR_ALL_CLASS (bdy_set, 0, mesh.num_bdy_sets, {
	
        // Loop over boundary patches in boundary set
        for (size_t bdy_node_lid=0; bdy_node_lid<mesh.num_bdy_nodes_in_set(bdy_set); bdy_node_lid++){

            // save the bdy_node_gid
            mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid) = temp_nodes_in_set(bdy_set, bdy_node_lid);
            
        } // end for
        
    }); // end FOR_ALL_CLASS bdy_set
    
    // update the host side for the number nodes in a bdy_set
    mesh.num_bdy_nodes_in_set.update_host();
    
    return;
} // end method to build boundary nodes

/* ----------------------------------------------------------------------------
   SGH solver loop
------------------------------------------------------------------------------- */

void FEA_Module_SGH::sgh_solve(){

    double time_value = simparam->time_value;
    const double time_final = simparam->time_final;
    const double dt_max = simparam->dt_max;
    const double dt_min = simparam->dt_min;
    const double dt_cfl = simparam->dt_cfl;
    double graphics_time = simparam->graphics_time;
    size_t graphics_cyc_ival = simparam->graphics_cyc_ival;
    double graphics_dt_ival = simparam->graphics_dt_ival;
    const size_t cycle_stop = simparam->cycle_stop;
    const size_t rk_num_stages = simparam->rk_num_stages;
    double dt = simparam->dt;
    const double fuzz = simparam->fuzz;
    const double tiny = simparam->tiny;
    const double small = simparam->small;
    CArray <double> graphics_times = simparam->graphics_times;
    size_t graphics_id = simparam->graphics_id;
    const CArrayKokkos <boundary_t> boundary = simparam->boundary;
    const CArrayKokkos <material_t> material = simparam->material;

    int myrank = Explicit_Solver_Pointer_->myrank;
    if(myrank==0)
      printf("Writing outputs to file at %f \n", time_value);
    /*
    write_outputs(mesh,
                  Explicit_Solver_Pointer_,
                  node_coords,
                  node_vel,
                  node_mass,
                  elem_den,
                  elem_pres,
                  elem_stress,
                  elem_sspd,
                  elem_sie,
                  elem_vol,
                  elem_mass,
                  elem_mat_id,
                  graphics_times,
                  graphics_id,
                  time_value);
      */
    
    
    CArrayKokkos <double> node_extensive_mass(mesh.num_nodes, "node_extensive_mass");
    
    // extensive energy tallies over the mesh elements local to this MPI rank
    double IE_t0 = 0.0;
    double KE_t0 = 0.0;
    double TE_t0 = 0.0;
    
    double IE_sum = 0.0;
    double KE_sum = 0.0;
    
    double IE_loc_sum = 0.0;
    double KE_loc_sum = 0.0;

    // extensive energy tallies over the entire mesh
    double global_IE_t0 = 0.0;
    double global_KE_t0 = 0.0;
    double global_TE_t0 = 0.0;
    int nlocal_elem_non_overlapping = Explicit_Solver_Pointer_->nlocal_elem_non_overlapping;
    const int num_dims = mesh.num_dims;
    
    // extensive IE
    REDUCE_SUM_CLASS(elem_gid, 0, nlocal_elem_non_overlapping, IE_loc_sum, {
        IE_loc_sum += elem_mass(elem_gid)*elem_sie(1,elem_gid);
        
    }, IE_sum);
    IE_t0 = IE_sum;

    MPI_Allreduce(&IE_t0,&global_IE_t0,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    
    // extensive KE
    REDUCE_SUM_CLASS(node_gid, 0, mesh.num_local_nodes, KE_loc_sum, {
        
        double ke = 0;
        for (size_t dim=0; dim<mesh.num_dims; dim++){
            ke += node_vel(1,node_gid,dim)*node_vel(1,node_gid,dim); // 1/2 at end
        } // end for
        
        if(mesh.num_dims==2){
            KE_loc_sum += node_mass(node_gid)*node_coords(1,node_gid,1)*ke;
        }
        else{
            KE_loc_sum += node_mass(node_gid)*ke;
        }
        
    }, KE_sum);
    Kokkos::fence();
    KE_t0 = 0.5*KE_sum;
    
    MPI_Allreduce(&KE_t0,&global_KE_t0,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    // extensive TE
    global_TE_t0 = global_IE_t0 + global_KE_t0;
    TE_t0 = global_TE_t0;
    KE_t0 = global_KE_t0;
    IE_t0 = global_IE_t0;
    
    
    // save the nodal mass
    FOR_ALL_CLASS(node_gid, 0, mesh.num_nodes, {
        
        double radius = 1.0;
        if(mesh.num_dims == 2){
            radius = node_coords(1,node_gid,1);
        }
        node_extensive_mass(node_gid) = node_mass(node_gid)*radius;
        
    }); // end parallel for
    
    
    // a flag to exit the calculation
    size_t stop_calc=0;
    
    auto time_1 = std::chrono::high_resolution_clock::now();
    
	// loop over the max number of time integration cycles
	for (size_t cycle = 0; cycle < cycle_stop; cycle++) {

	    // stop calculation if flag
	    if (stop_calc == 1) break;
        

	    // get the step
        if(mesh.num_dims==2){
            get_timestep2D(mesh,
                           node_coords,
                           node_vel,
                           elem_sspd,
                           elem_vol,
                           time_value,
                           graphics_time,
                           time_final,
                           dt_max,
                           dt_min,
                           dt_cfl,
                           dt,
                           fuzz);
        }
        else {
            get_timestep(mesh,
                         node_coords,
                         node_vel,
                         elem_sspd,
                         elem_vol,
                         time_value,
                         graphics_time,
                         time_final,
                         dt_max,
                         dt_min,
                         dt_cfl,
                         dt,
                         fuzz);
        } // end if 2D
        
        double global_dt;
        MPI_Allreduce(&dt,&global_dt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        dt = global_dt;

        if (cycle==0){
            if(myrank==0)
              printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        }
        // print time step every 10 cycles
        else if (cycle%20==0){
            if(myrank==0)
              printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        } // end if
        
        
        // ---------------------------------------------------------------------
        //  integrate the solution forward to t(n+1) via Runge Kutta (RK) method
        // ---------------------------------------------------------------------
	     
        // save the values at t_n
        rk_init(node_coords,
                node_vel,
                elem_sie,
                elem_stress,
                mesh.num_dims,
                mesh.num_elems,
                mesh.num_nodes);
	    
        
        
        // integrate solution forward in time
        for (size_t rk_stage = 0; rk_stage < rk_num_stages; rk_stage++){

            
            // ---- RK coefficient ----
            double rk_alpha = 1.0/((double)rk_num_stages - (double)rk_stage);
            
            // ---- Calculate velocity diveregence for the element ----
            if(mesh.num_dims==2){
                get_divergence2D(elem_div,
                                 mesh,
                                 node_coords,
                                 node_vel,
                                 elem_vol);
            }
            else {
                get_divergence(elem_div,
                               mesh,
                               node_coords,
                               node_vel,
                               elem_vol);
            } // end if 2D
            
            // ---- calculate the forces on the vertices and evolve stress (hypo model) ----
            if(mesh.num_dims==2){
                get_force_sgh2D(material,
                                mesh,
                                node_coords,
                                node_vel,
                                elem_den,
                                elem_sie,
                                elem_pres,
                                elem_stress,
                                elem_sspd,
                                elem_vol,
                                elem_div,
                                elem_mat_id,
                                corner_force,
                                fuzz,
                                small,
                                elem_statev,
                                dt,
                                rk_alpha);
            }
            else {
                get_force_sgh(material,
                              mesh,
                              node_coords,
                              node_vel,
                              elem_den,
                              elem_sie,
                              elem_pres,
                              elem_stress,
                              elem_sspd,
                              elem_vol,
                              elem_div,
                              elem_mat_id,
                              corner_force,
                              fuzz,
                              small,
                              elem_statev,
                              dt,
                              rk_alpha);
            }

            /*
            debug block
            if(myrank==1){
             std::cout << rk_alpha << " " << dt << std::endl;
             for(int i = 0; i < mesh.num_nodes; i++){
               double node_force[3];
               for (size_t dim = 0; dim < num_dims; dim++){
                 node_force[dim] = 0.0;
               } // end for dim
        
               // loop over all corners around the node and calculate the nodal force
               for (size_t corner_lid=0; corner_lid<mesh.num_corners_in_node(i); corner_lid++){
        
                 // Get corner gid
                 size_t corner_gid = mesh.corners_in_node(i, corner_lid);
                 std::cout << Explicit_Solver_Pointer_->all_node_map->getGlobalElement(i) << " " << corner_gid << " " << corner_force(corner_gid, 0) << " " << corner_force(corner_gid, 1) << " " << corner_force(corner_gid, 2) << std::endl;
                 // loop over dimension
                 for (size_t dim = 0; dim < num_dims; dim++){
                   node_force[dim] += corner_force(corner_gid, dim);
                 } // end for dim
            
               } // end for corner_lid
               //std::cout << Explicit_Solver_Pointer_->all_node_map->getGlobalElement(i) << " " << node_force[0] << " " << node_force[1] << " " << node_force[2] << std::endl;
               //std::cout << Explicit_Solver_Pointer_->all_node_map->getGlobalElement(i) << " " << node_mass(i) << std::endl;
             }
            }
            /*
            //debug print vector values on a rank
            /*
            if(myrank==0)
             for(int i = 0; i < mesh.num_nodes; i++){
               std::cout << Explicit_Solver_Pointer_->all_node_map->getGlobalElement(i) << " " << node_vel(1,i,0) << " " << node_vel(1,i,1) << " " << node_vel(1,i,2) << std::endl;
             }
            */

            // ---- Update nodal velocities ---- //
            update_velocity_sgh(rk_alpha,
                                dt,
                                mesh,
                                node_vel,
                                node_mass,
                                corner_force);
            
            // ---- apply force boundary conditions to the boundary patches----
            boundary_velocity(mesh, boundary, node_vel);

            //current interface has differing velocity arrays; this equates them until we unify memory
            //view scope
            {
              Explicit_Solver_SGH::vec_array node_velocities_interface = Explicit_Solver_Pointer_->node_velocities_distributed->getLocalView<Explicit_Solver_SGH::device_type> (Tpetra::Access::ReadWrite);
              FOR_ALL_CLASS(node_gid, 0, mesh.num_local_nodes, {
                for (int idim = 0; idim < num_dims; idim++){
                  node_velocities_interface(node_gid,idim) = node_vel(1,node_gid,idim);
                }
              }); // end parallel for
            } //end view scope
            Kokkos::fence();
            //communicate ghost velocities
            Explicit_Solver_Pointer_->comm_velocities();

            //this is forcing a copy to the device
            //view scope
            {
              Explicit_Solver_SGH::vec_array all_node_velocities_interface = Explicit_Solver_Pointer_->all_node_velocities_distributed->getLocalView<Explicit_Solver_SGH::device_type> (Tpetra::Access::ReadWrite);

              FOR_ALL_CLASS(node_gid, mesh.num_local_nodes, mesh.num_nodes, {
                for (int idim = 0; idim < num_dims; idim++){
                  node_vel(1,node_gid,idim) = all_node_velocities_interface(node_gid,idim);
                }
        
              }); // end parallel for
            } //end view scope
            Kokkos::fence();
            
            //debug print vector values on a rank
            /*
            if(myrank==0)
             for(int i = 0; i < mesh.num_nodes; i++){
               std::cout << Explicit_Solver_Pointer_->all_node_map->getGlobalElement(i) << " " << node_vel(1,i,0) << " " << node_vel(1,i,1) << " " << node_vel(1,i,2) << std::endl;
             }
            */ 
            // ---- Update specific internal energy in the elements ----
            update_energy_sgh(rk_alpha,
                              dt,
                              mesh,
                              node_vel,
                              node_coords,
                              elem_sie,
                              elem_mass,
                              corner_force);
            
            
            // ---- Update nodal positions ----
            update_position_sgh(rk_alpha,
                                dt,
                                mesh.num_dims,
                                mesh.num_nodes,
                                node_coords,
                                node_vel);
            
            
            // ---- Calculate cell volume for next time step ----
            get_vol();
            
            
            
            // ---- Calculate elem state (den, pres, sound speed, stress) for next time step ----
            if(mesh.num_dims==2){
                update_state2D(material,
                               mesh,
                               node_coords,
                               node_vel,
                               elem_den,
                               elem_pres,
                               elem_stress,
                               elem_sspd,
                               elem_sie,
                               elem_vol,
                               elem_mass,
                               elem_mat_id,
                               elem_statev,
                               dt,
                               rk_alpha);
            }
            else{
                update_state(material,
                             mesh,
                             node_coords,
                             node_vel,
                             elem_den,
                             elem_pres,
                             elem_stress,
                             elem_sspd,
                             elem_sie,
                             elem_vol,
                             elem_mass,
                             elem_mat_id,
                             elem_statev,
                             dt,
                             rk_alpha);
            }
            // ----
            // Notes on strength:
            //    1) hyper-elastic strength models are called in update_state
            //    2) hypo-elastic strength models are called in get_force
            //    3) strength models must be added by the user in user_mat.cpp

            
            // calculate the new corner masses if 2D
            if(mesh.num_dims==2){
                
                // calculate the nodal areal mass
                FOR_ALL_CLASS(node_gid, 0, mesh.num_nodes, {
                    
                    node_mass(node_gid) = 0.0;
                    
                    if (node_coords(1,node_gid,1) > tiny){
                        node_mass(node_gid) = node_extensive_mass(node_gid)/node_coords(1,node_gid,1);
                    }
                    //if(cycle==0&&node_gid==1&&myrank==0)
                      //std::cout << "index " << node_gid << " on rank " << myrank << " node vel " << node_vel(1,node_gid,0) << "  " << node_mass(node_gid) << std::endl << std::flush;

                }); // end parallel for over node_gid
                Kokkos::fence();

                //current interface has differing density arrays; this equates them until we unify memory
                //view scope
                {
                  Explicit_Solver_SGH::vec_array node_densities_interface = Explicit_Solver_Pointer_->design_node_densities_distributed->getLocalView<Explicit_Solver_SGH::device_type> (Tpetra::Access::ReadWrite);
                  FOR_ALL_CLASS(node_gid, 0, mesh.num_local_nodes, {
                    node_densities_interface(node_gid,0) = node_mass(node_gid);
                  }); // end parallel for
                } //end view scope
                Kokkos::fence();
                //communicate ghost densities
                Explicit_Solver_Pointer_->comm_densities();

                //this is forcing a copy to the device
                //view scope
                {
                  Explicit_Solver_SGH::vec_array all_node_densities_interface = Explicit_Solver_Pointer_->all_node_densities_distributed->getLocalView<Explicit_Solver_SGH::device_type> (Tpetra::Access::ReadWrite);

                  FOR_ALL_CLASS(node_gid, mesh.num_local_nodes, mesh.num_nodes, {
                    node_mass(node_gid) = all_node_densities_interface(node_gid,0);
                  }); // end parallel for
                } //end view scope
                Kokkos::fence();
                
                
                // -----------------------------------------------
                // Calcualte the areal mass for nodes on the axis
                // -----------------------------------------------
                // The node order of the 2D element is
                //
                //   J
                //   |
                // 3---2
                // |   |  -- I
                // 0---1
                /*
                FOR_ALL_CLASS(elem_gid, 0, mesh.num_elems, {
                    
                    // loop over the corners of the element and calculate the mass
                    for (size_t node_lid=0; node_lid<4; node_lid++){
                        
                        size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                        size_t node_minus_gid;
                        size_t node_plus_gid;
                        
                        
                        if (node_coords(1,node_gid,1) < tiny){
                            // node is on the axis
                            
                            // minus node
                            if (node_lid==0){
                                node_minus_gid = mesh.nodes_in_elem(elem_gid, 3);
                            } else {
                                node_minus_gid = mesh.nodes_in_elem(elem_gid, node_lid-1);
                            }
                            
                            // plus node
                            if (node_lid==3){
                                node_plus_gid = mesh.nodes_in_elem(elem_gid, 0);
                            } else {
                                node_plus_gid = mesh.nodes_in_elem(elem_gid, node_lid+1);
                            }
                            
                            node_mass(node_gid) = fmax(node_mass(node_plus_gid), node_mass(node_minus_gid))/2.0;
                            
                        } // end if
                        
                    } // end for over corners
                    
                }); // end parallel for over elem_gid
                Kokkos::fence();
                 */
                
                FOR_ALL_CLASS(node_bdy_gid, 0, mesh.num_bdy_nodes, {
                //FOR_ALL_CLASS(node_gid, 0, mesh.num_local_nodes, {    
                    size_t node_gid = mesh.bdy_nodes(node_bdy_gid);
                    
                    if (node_coords(1,node_gid,1) < tiny){
                        // node is on the axis
                        
                        for(size_t node_lid=0; node_lid<mesh.num_nodes_in_node(node_gid); node_lid++){
                            
                            size_t node_neighbor_gid = mesh.nodes_in_node(node_gid, node_lid);
                            
                            // if the node is off the axis, use it's areal mass on the boundary
                            if (node_coords(1,node_neighbor_gid,1) > tiny){
                                node_mass(node_gid) = fmax(node_mass(node_gid), node_mass(node_neighbor_gid)/2.0);
                            }

                        } // end for over neighboring nodes
                        
                    } // end if
                    
                }); // end parallel for over elem_gid
                
            } // end of if 2D-RZ

        } // end of RK loop

	    // increment the time
	    time_value+=dt;
    	
        
        size_t write = 0;
        if ((cycle+1)%graphics_cyc_ival == 0 && cycle>0){
            write = 1;
        }
        else if (cycle == cycle_stop) {
            write = 1;
        }
        else if (time_value >= time_final){
            write = 1;
        }
        else if (time_value >= graphics_time){
            write = 1;
        }
            
        // write outputs
        if (write == 1){
            //interface nodal coordinate data
            //view scope
            {
              Explicit_Solver_SGH::vec_array node_coords_interface = Explicit_Solver_Pointer_->node_coords_distributed->getLocalView<Explicit_Solver_SGH::device_type> (Tpetra::Access::ReadWrite);
              FOR_ALL_CLASS(node_gid, 0, mesh.num_local_nodes, {
                for (int idim = 0; idim < num_dims; idim++){
                  node_coords_interface(node_gid,idim) = node_coords(1,node_gid,idim);
                }
              }); // end parallel for
            } //end view scope

            if(myrank==0)
              printf("Writing outputs to file at %f \n", graphics_time);
              /*
            write_outputs(mesh,
                          Explicit_Solver_Pointer_,
                          node_coords,
                          node_vel,
                          node_mass,
                          elem_den,
                          elem_pres,
                          elem_stress,
                          elem_sspd,
                          elem_sie,
                          elem_vol,
                          elem_mass,
                          elem_mat_id,
                          graphics_times,
                          graphics_id,
                          time_value);
            */
            graphics_time = time_value + graphics_dt_ival;
        } // end if
        
        
        // end of calculation
        if (time_value>=time_final) break;
        
    } // end for cycle loop
    
    
    auto time_2 = std::chrono::high_resolution_clock::now();
    auto calc_time = std::chrono::duration_cast
                           <std::chrono::nanoseconds>(time_2 - time_1).count();
    
    if(myrank==0)
      printf("\nCalculation time in seconds: %f \n", calc_time * 1e-9);
    
    // ---- Calculate energy tallies ----
    double IE_tend = 0.0;
    double KE_tend = 0.0;
    double TE_tend = 0.0;

    double global_IE_tend = 0.0;
    double global_KE_tend = 0.0;
    double global_TE_tend = 0.0;
    
    IE_loc_sum = 0.0;
    KE_loc_sum = 0.0;
    IE_sum = 0.0;
    KE_sum = 0.0;
    
    // extensive IE
    REDUCE_SUM_CLASS(elem_gid, 0, nlocal_elem_non_overlapping, IE_loc_sum, {
        
        IE_loc_sum += elem_mass(elem_gid)*elem_sie(1,elem_gid);
        
    }, IE_sum);
    IE_tend = IE_sum;

    //reduce over MPI ranks
    MPI_Allreduce(&IE_tend,&global_IE_tend,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    // extensive KE
    REDUCE_SUM_CLASS(node_gid, 0, mesh.num_local_nodes, KE_loc_sum, {
        
        double ke = 0;
        for (size_t dim=0; dim<mesh.num_dims; dim++){
            ke += node_vel(1,node_gid,dim)*node_vel(1,node_gid,dim); // 1/2 at end
        } // end for
        
        if(mesh.num_dims==2){
            KE_loc_sum += node_mass(node_gid)*node_coords(1,node_gid,1)*ke;
        }
        else{
            KE_loc_sum += node_mass(node_gid)*ke;
        }
            
        
    }, KE_sum);
    Kokkos::fence();
    KE_tend = 0.5*KE_sum;
    
    //reduce over MPI ranks
    MPI_Allreduce(&KE_tend,&global_KE_tend,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    // extensive TE
    TE_tend = IE_tend + KE_tend;
    KE_tend = global_KE_tend;
    IE_tend = global_IE_tend;
    
    // extensive TE
    TE_tend = IE_tend + KE_tend;

    //reduce over MPI ranks
    
    if(myrank==0)
      printf("Time=0:   KE = %20.15f, IE = %20.15f, TE = %20.15f \n", KE_t0, IE_t0, TE_t0);
    if(myrank==0)
      printf("Time=End: KE = %20.15f, IE = %20.15f, TE = %20.15f \n", KE_tend, IE_tend, TE_tend);
    if(myrank==0)
      printf("total energy conservation error %= %e \n\n", 100*(TE_tend - TE_t0)/TE_t0);
    
    return;
    
} // end of SGH solve
