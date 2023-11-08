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
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "Tpetra_Details_FixedHashTable.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Import_Util2.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include <set>

#include "elements.h"
#include "swage.h"
#include "mesh.h"
#include "ref_elem.h"
#include "matar.h"
#include "utilities.h"
#include "node_combination.h"
#include "Simulation_Parameters_RDH.h"
#include "Simulation_Parameters_Dynamic_Optimization.h"
#include "Simulation_Parameters_Elasticity.h"
#include "FEA_Module_RDH.h"
#include "Explicit_Solver.h"

//optimization
#include "ROL_Algorithm.hpp"
#include "ROL_Solver.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Types.hpp"
#include "ROL_Elementwise_Reduce.hpp"


#define MAX_ELEM_NODES 8
#define STRAIN_EPSILON 0.000000001
#define DENSITY_EPSILON 0.0001
#define BC_EPSILON 1.0e-6
#define BUFFER_GROW 100

using namespace utils;

FEA_Module_RDH::FEA_Module_RDH(Solver *Solver_Pointer, std::shared_ptr<mesh_t> mesh_in, const int my_fea_module_index) :FEA_Module(Solver_Pointer){

  //assign interfacing index
  my_fea_module_index_ = my_fea_module_index;
  
  //recast solver pointer for non-base class access
  Explicit_Solver_Pointer_ = dynamic_cast<Explicit_Solver*>(Solver_Pointer);

  //create parameter object
  simparam = Simulation_Parameters_RDH();
  simparam = Yaml::from_file<Simulation_Parameters_RDH>(Explicit_Solver_Pointer_->filename);
  // ---- Read input file, define state and boundary conditions ---- //
  //simparam->input();
  
  //TO parameters
  simparam_dynamic_opt = Explicit_Solver_Pointer_->simparam_dynamic_opt;

  //create ref element object
  //ref_elem = new elements::ref_element();
  //create mesh objects
  //init_mesh = new swage::mesh_t(simparam);
  //mesh = new swage::mesh_t(simparam);

  // WARNING WARNING WARNING //
  mesh = mesh_in;

  //boundary condition data
  max_boundary_sets = 0;
  Local_Index_Boundary_Patches = Explicit_Solver_Pointer_->Local_Index_Boundary_Patches;

  //set Tpetra vector pointers
  initial_node_velocities_distributed = Explicit_Solver_Pointer_->initial_node_velocities_distributed;
  initial_node_coords_distributed = Explicit_Solver_Pointer_->initial_node_coords_distributed;
  all_initial_node_coords_distributed = Explicit_Solver_Pointer_->all_initial_node_coords_distributed;
  node_coords_distributed = Explicit_Solver_Pointer_->node_coords_distributed;
  node_velocities_distributed = Explicit_Solver_Pointer_->node_velocities_distributed;
  all_node_velocities_distributed = Explicit_Solver_Pointer_->all_node_velocities_distributed;
  if(simparam_dynamic_opt.topology_optimization_on||simparam_dynamic_opt.shape_optimization_on){
    all_cached_node_velocities_distributed = Teuchos::rcp(new MV(all_node_map, simparam.num_dims));
    force_gradient_velocity = Teuchos::rcp(new MV(all_node_map, simparam.num_dims));
    force_gradient_position = Teuchos::rcp(new MV(all_node_map, simparam.num_dims));
    force_gradient_design = Teuchos::rcp(new MV(all_node_map, 1));
    corner_value_storage = Solver_Pointer->corner_value_storage;
    corner_vector_storage = Solver_Pointer->corner_vector_storage;
    relative_element_densities = DCArrayKokkos<double>(rnum_elem, "relative_element_densities");
  }

  if(simparam_dynamic_opt.topology_optimization_on||simparam_dynamic_opt.shape_optimization_on||simparam.num_dims==2){
    node_masses_distributed = Teuchos::rcp(new MV(map, 1));
    ghost_node_masses_distributed = Teuchos::rcp(new MV(ghost_node_map, 1));
    adjoint_vector_distributed = Teuchos::rcp(new MV(map, simparam.num_dims));
    phi_adjoint_vector_distributed = Teuchos::rcp(new MV(map, simparam.num_dims));
    psi_adjoint_vector_distributed = Teuchos::rcp(new MV(all_element_map, 1));
  }
  
  //setup output
  noutput = 0;
  //init_output();

  //optimization flags
  kinetic_energy_objective = false;
  

  //set parameters
  Time_Variables tv = simparam.time_variables;
  time_value = simparam.time_value;
  time_final = tv.time_final;
  dt_max = tv.dt_max;
  dt_min = tv.dt_min;
  dt_cfl = tv.dt_cfl;
  graphics_time = simparam.graphics_options.graphics_time;
  graphics_cyc_ival = simparam.graphics_options.graphics_cyc_ival;
  graphics_dt_ival = simparam.graphics_options.graphics_dt_ival;
  cycle_stop = tv.cycle_stop;
  rk_num_stages = simparam.rk_num_stages;
  dt = tv.dt;
  fuzz = tv.fuzz;
  tiny = tv.tiny;
  small = tv.small;
  graphics_times = simparam.graphics_options.graphics_times;
  graphics_id = simparam.graphics_options.graphics_id;

  if(simparam_dynamic_opt.topology_optimization_on){
    max_time_steps = BUFFER_GROW;
    forward_solve_velocity_data = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps+1));
    time_data.resize(max_time_steps+1);
    forward_solve_coordinate_data = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps+1));
    adjoint_vector_data = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps+1));
    phi_adjoint_vector_data = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps+1));
    psi_adjoint_vector_data = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps+1));
    //assign a multivector of corresponding size to each new timestep in the buffer
    for(int istep = 0; istep < max_time_steps+1; istep++){
      (*forward_solve_velocity_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam.num_dims));
      (*forward_solve_coordinate_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam.num_dims));
      (*adjoint_vector_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam.num_dims));
      (*phi_adjoint_vector_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam.num_dims));
      (*psi_adjoint_vector_data)[istep] = Teuchos::rcp(new MV(all_element_map, 1));
    }
    
  }

  have_loading_conditions = false;

}

FEA_Module_RDH::~FEA_Module_RDH(){
   //delete simparam;
}


// -----------------------------------------------------------------------------
// Interfaces read in data with the RDH solver data; currently a hack to streamline
//------------------------------------------------------------------------------
void FEA_Module_RDH::rdh_interface_setup(node_t &node, elem_t &elem, mesh_t &mesh, corner_t &corner){

    const size_t num_dim = simparam.num_dims;
    const size_t rk_num_bins = simparam.rk_num_bins;

    num_nodes_in_elem = 1;
    for (int dim=0; dim<num_dim; dim++){
        num_nodes_in_elem *= 2;
    }

    // --- Read in the nodes in the mesh ---

    nall_nodes = Explicit_Solver_Pointer_->nall_nodes;
    int myrank = Explicit_Solver_Pointer_->myrank;
    int nranks = Explicit_Solver_Pointer_->nranks;
    //printf("Num nodes assigned to MPI rank %lu is %lu\n" , myrank, nall_nodes);

    // intialize node variables
    mesh.initialize_nodes(nall_nodes);
    // WARNING WARNING WARNING //
    //mesh.initialize_local_nodes(Explicit_Solver_Pointer_->nlocal_nodes);
    node.initialize(rk_num_bins, nall_nodes, num_dim);

    //view scope
    {
      host_vec_array interface_node_coords = Explicit_Solver_Pointer_->all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
      //save node data to node.coords
      if(num_dim==2){
        for(int inode = 0; inode < nall_nodes; inode++){
          //std::cout << "Node index " << inode+1 << " ";
          node.coords.host(0,inode,0) = interface_node_coords(inode,0);
          //std::cout << host_node_coords_state(0,inode,0)+1<< " ";
          node.coords.host(0,inode,1) = interface_node_coords(inode,1);
          //std::cout << host_node_coords_state(0,inode,1)+1<< " ";
        }
      }
      else if(num_dim==3){
        for(int inode = 0; inode < nall_nodes; inode++){
          node.coords.host(0,inode,0) = interface_node_coords(inode,0);
          node.coords.host(0,inode,1) = interface_node_coords(inode,1);
          node.coords.host(0,inode,2) = interface_node_coords(inode,2);
        }
      }
    } //end view scope
    // --- read in the elements in the mesh ---
    
    rnum_elem = Explicit_Solver_Pointer_->rnum_elem;

// intialize elem variables
    mesh.initialize_elems(rnum_elem, num_dim);
    elem.initialize(rk_num_bins, nall_nodes, 3); // always 3D here, even for 2D
    nodes_in_elem = mesh.nodes_in_elem;
    {
      host_elem_conn_array interface_nodes_in_elem = Explicit_Solver_Pointer_->global_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
      for(int ielem = 0; ielem < rnum_elem; ielem++){
        //std::cout << "Element index " << ielem+1 << " ";
        for(int inode = 0; inode < num_nodes_in_elem; inode++){
            nodes_in_elem.host(ielem,inode) = Explicit_Solver_Pointer_->all_node_map->getLocalElement(interface_nodes_in_elem(ielem,inode));
        }
        //std::cout << std::endl;
      }
    }
    // update device side
    nodes_in_elem.update_device();
        
    // save the node coords to the current RK value
    for (size_t node_gid=0; node_gid<nall_nodes; node_gid++){
        
        for(int rk=1; rk<rk_num_bins; rk++){
            for (int dim = 0; dim < num_dim; dim++){
                node.coords.host(rk, node_gid, dim) = node.coords.host(0, node_gid, dim);
            } // end for dim
        } // end for rk
        
    } // end parallel for
    

    node.coords.update_device();
    
    // intialize corner variables
    int num_corners = rnum_elem*num_nodes_in_elem;
    mesh.initialize_corners(num_corners);
    corner.initialize(num_corners, num_dim);
    
    return;
    
}


/* ----------------------------------------------------------------------------
    Deallocate memory used for  material models
------------------------------------------------------------------------------- */
void FEA_Module_RDH::module_cleanup(){
  cleanup_material_models();
}

void FEA_Module_RDH::cleanup_material_models() {

    const DCArrayKokkos <material_t> material = simparam.material;

    // destroy strength model
    destroy_strength_model(elem_strength,
                           material,
                           elem_mat_id,
                           global_vars,
                           elem_user_output_vars,
                           rnum_elem);

    // destroy eos model
    destroy_eos_model(elem_eos,
                      material,
                      elem_mat_id,
                      global_vars,
                      elem_user_output_vars,
                      rnum_elem);
    return;

} // end cleanup_user_strength_model;

/* ----------------------------------------------------------------------------
   Setup RDH solver data
------------------------------------------------------------------------------- */

void FEA_Module_RDH::setup(mesh_t &mesh){

    const size_t rk_level = simparam.rk_num_bins - 1;   
    const size_t num_fills = simparam.region_options.size();
    const size_t rk_num_bins = simparam.rk_num_bins;
    const size_t num_bcs = simparam.boundary_conditions.size();
    const size_t num_materials = simparam.material_options.size();
    const int num_dim = simparam.num_dims;
    const size_t num_lcs = simparam.loading.size();
    if(num_lcs)
      have_loading_conditions = true;


    // ---------------------------------------------------------------------
    //    obtain mesh data
    // --------------------------------------------------------------------- 
    rdh_interface_setup(node_interface, elem_interface, mesh_interface, corner_interface);
    mesh.build_corner_connectivity();
    mesh.build_elem_elem_connectivity();
    mesh.num_bdy_patches = nboundary_patches;
    if(num_dim==2){
      mesh.build_patch_connectivity();
      mesh.build_node_node_connectivity();
    }
        
      // ---------------------------------------------------------------------
      //    allocate memory
      // ---------------------------------------------------------------------

      // shorthand names
    const size_t num_nodes = mesh.num_nodes;
    const size_t num_elems = mesh.num_elems;
    const size_t num_corners = mesh.num_corners;
        
    // create Dual Views of the individual node struct variables
    node_coords = DViewCArrayKokkos<double>(node_interface.coords.get_kokkos_dual_view().view_host().data(),rk_num_bins,num_nodes,num_dim);

    node_vel = DViewCArrayKokkos<double>(node_interface.vel.get_kokkos_dual_view().view_host().data(),rk_num_bins,num_nodes,num_dim);

    node_mass = DViewCArrayKokkos<double>(node_interface.mass.get_kokkos_dual_view().view_host().data(),num_nodes);
        
        
    // create Dual Views of the individual elem struct variables
    elem_den = DViewCArrayKokkos<double>(&elem_interface.den(0),
                                            num_elems);

    elem_pres = DViewCArrayKokkos<double>(&elem_interface.pres(0),
                                              num_elems);

    elem_stress = DViewCArrayKokkos<double>(&elem_interface.stress(0,0,0,0),
                                                rk_num_bins,
                                                num_elems,
                                                3,
                                                3); // always 3D even in 2D-RZ
    
    elem_sspd = DViewCArrayKokkos<double>(&elem_interface.sspd(0),
                                              num_elems);

    elem_sie = DViewCArrayKokkos<double>(&elem_interface.sie(0,0),
                                            rk_num_bins,
                                            num_elems);

    elem_vol = DViewCArrayKokkos<double>(&elem_interface.vol(0),
                                            num_elems);
        
    elem_div = DViewCArrayKokkos<double>(&elem_interface.div(0),
                                            num_elems);
        

    elem_mass = DViewCArrayKokkos<double>(&elem_interface.mass(0),
                                              num_elems);

    elem_mat_id = DViewCArrayKokkos<size_t>(&elem_interface.mat_id(0),
                                                num_elems);
      
    // create Dual Views of the corner struct variables
    corner_force = DViewCArrayKokkos <double>(&corner_interface.force(0,0),
                                                num_corners, 
                                                num_dim);

    corner_mass = DViewCArrayKokkos <double>(&corner_interface.mass(0),
                                                num_corners);
        
    // allocate elem_vel_grad
    elem_vel_grad = DCArrayKokkos <double> (num_elems,3,3);

    // allocate material models
    elem_eos = DCArrayKokkos <eos_t> (num_elems);
    elem_strength = DCArrayKokkos <strength_t> (num_elems); 
      
    // ---------------------------------------------------------------------
    //   calculate geometry
    // ---------------------------------------------------------------------
    node_coords.update_device();
    Kokkos::fence();

    //get_vol();

    //FEA_Module bc variable
    num_boundary_conditions = num_bcs;

    const DCArrayKokkos <mat_fill_t> mat_fill = simparam.mat_fill;
    const DCArrayKokkos <boundary_t> boundary = simparam.boundary;
    const DCArrayKokkos <material_t> material = simparam.material;
    global_vars = simparam.global_vars;
    elem_user_output_vars = DCArrayKokkos <double> (rnum_elem, simparam.output_options.max_num_user_output_vars); 
 
    //--- calculate bdy sets ---//
    mesh.num_nodes_in_patch = 2*(num_dim-1);  // 2 (2D) or 4 (3D)
    mesh.num_patches_in_elem = 2*num_dim; // 4 (2D) or 6 (3D)
    mesh.init_bdy_sets(num_bcs);
    num_bdy_sets = mesh.num_bdy_sets;
    printf("Num BC's = %lu\n", num_bcs);

    // patch ids in bdy set
    bdy_patches_in_set = mesh.bdy_patches_in_set;
    if(num_dim==2)
      bdy_nodes = mesh.bdy_nodes;


    // node ids in bdy_patch set
    bdy_nodes_in_set = mesh.bdy_nodes_in_set;
    num_bdy_nodes_in_set = mesh.num_bdy_nodes_in_set;
    
    //assign mesh views needed by the FEA module

    // elem ids in elem
    elems_in_elem = mesh.elems_in_elem;
    num_elems_in_elem = mesh.num_elems_in_elem;

    //corners
    num_corners_in_node = mesh.num_corners_in_node;
    corners_in_node = mesh.corners_in_node;
    corners_in_elem = mesh.corners_in_elem;

    //elem-node conn & node-node conn
    elems_in_node = mesh.elems_in_node;
    if(num_dim==2){
      nodes_in_node = mesh.nodes_in_node;
      num_nodes_in_node = mesh.num_nodes_in_node;
      //patch conn
    
      patches_in_elem = mesh.patches_in_elem;
      nodes_in_patch = mesh.nodes_in_patch;
      elems_in_patch = mesh.elems_in_patch;
    }

    // loop over BCs
    for (size_t this_bdy = 0; this_bdy < num_bcs; this_bdy++){
        
        RUN_CLASS({
            printf("  Boundary Condition number %lu \n", this_bdy);
            printf("  Num bdy patches in this set = %lu \n", bdy_patches_in_set.stride(this_bdy));
            printf("  Num bdy nodes in this set = %lu \n", bdy_nodes_in_set.stride(this_bdy));
        });
        Kokkos::fence();

    }// end for
    for (int f_id = 0; f_id < num_fills; f_id++){
      FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
        elem_mat_id(elem_gid) = mat_fill(f_id).mat_id;
      });
    }
    elem_mat_id.update_host();
 
    // initialize strength model
    init_strength_model(elem_strength,
                        material,
                        elem_mat_id,
                        global_vars,
                        elem_user_output_vars,
                        rnum_elem);

    // initialize eos model
    init_eos_model(elem_eos,
                   material,
                   elem_mat_id,
                   global_vars,
                   elem_user_output_vars,
                   rnum_elem);
    
//--- apply the fill instructions over each of the Elements---//

    // loop over the fill instructures
    for (int f_id = 0; f_id < num_fills; f_id++){
            
        // parallel loop over elements in mesh
        FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {

            // calculate the coordinates and radius of the element
            double elem_coords[3]; // note:initialization with a list won't work
            elem_coords[0] = 0.0;
            elem_coords[1] = 0.0;
            elem_coords[2] = 0.0;

            // get the coordinates of the element center
            for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
                elem_coords[0] += node_coords(rk_level, nodes_in_elem(elem_gid, node_lid), 0);
                elem_coords[1] += node_coords(rk_level, nodes_in_elem(elem_gid, node_lid), 1);
                if (num_dim == 3){
                    elem_coords[2] += node_coords(rk_level, nodes_in_elem(elem_gid, node_lid), 2);
                } else
                {
                    elem_coords[2] = 0.0;
                }
            } // end loop over nodes in element
            elem_coords[0] = elem_coords[0]/num_nodes_in_elem;
            elem_coords[1] = elem_coords[1]/num_nodes_in_elem;
            elem_coords[2] = elem_coords[2]/num_nodes_in_elem;

            // default is not to fill the element
            bool fill_this = mat_fill(f_id).contains(elem_coords);

            // paint the material state on the element
            if (fill_this){
                    
                // density
                elem_den(elem_gid) = mat_fill(f_id).den;
                
                // mass
                elem_mass(elem_gid) = elem_den(elem_gid)*elem_vol(elem_gid);
                
                // specific internal energy
                elem_sie(rk_level, elem_gid) = mat_fill(f_id).sie;
		
                size_t mat_id = elem_mat_id(elem_gid); // short name
                
                // --- stress tensor ---
                // always 3D even for 2D-RZ
                for (size_t i=0; i<3; i++){
                    for (size_t j=0; j<3; j++){
                        elem_stress(rk_level,elem_gid,i,j) = 0.0;
                    }        
                }  // end for
                
                // short form for clean code
                EOSParent * eos_model = elem_eos(elem_gid).model;

                // --- Pressure ---
                eos_model->calc_pressure(elem_pres,
                                         elem_stress,
                                         elem_gid,
                                         elem_mat_id(elem_gid),
                                         global_vars,
                                         elem_user_output_vars,
                                         elem_sspd,
                                         elem_den(elem_gid),
                                         elem_sie(rk_level,elem_gid));

                // --- Sound speed ---
                eos_model->calc_sound_speed(elem_pres,
                                            elem_stress,
                                            elem_gid,
                                            elem_mat_id(elem_gid),
                                            global_vars,
                                            elem_user_output_vars,
                                            elem_sspd,
                                            elem_den(elem_gid),
                                            elem_sie(rk_level,elem_gid));
                
                // loop over the nodes of this element and apply velocity
                for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){

                    // get the mesh node index
                    size_t node_gid = nodes_in_elem(elem_gid, node_lid);

                
                    // --- Velocity ---
                    switch(mat_fill(f_id).velocity)
                    {
                        case VELOCITY_TYPE::cartesian:
                        {
                        
                            node_vel(rk_level, node_gid, 0) = mat_fill(f_id).u;
                            node_vel(rk_level, node_gid, 1) = mat_fill(f_id).v;
                            if (num_dim == 3) node_vel(rk_level, node_gid, 2) = mat_fill(f_id).w;
                            
                        
                            break;
                        }
                        case VELOCITY_TYPE::radial:
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
                            if (num_dim == 3) node_vel(rk_level, node_gid, 2) = 0.0;
                            
                            break;
                        }
                        case VELOCITY_TYPE::spherical:
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
                            if (num_dim == 3) node_vel(rk_level, node_gid, 2) = mat_fill(f_id).speed*dir[2];

                            break;
                        }
                        case VELOCITY_TYPE::radial_linear:
                        {
                        
                            break;
                        }
                        case VELOCITY_TYPE::spherical_linear:
                        {
                        
                            break;
                        }
                        case VELOCITY_TYPE::tg_vortex:
                        {
                        
                            node_vel(rk_level, node_gid, 0) = sin(PI * node_coords(rk_level,node_gid, 0)) * cos(PI * node_coords(rk_level,node_gid, 1)); 
                            node_vel(rk_level, node_gid, 1) =  -1.0*cos(PI * node_coords(rk_level,node_gid, 0)) * sin(PI * node_coords(rk_level,node_gid, 1)); 
                            if (num_dim == 3) node_vel(rk_level, node_gid, 2) = 0.0;

                            break;
                        }
                    } // end of switch

                }// end loop over nodes of element
                
                
                if(mat_fill(f_id).velocity == VELOCITY_TYPE::tg_vortex)
                {
                    elem_pres(elem_gid) = 0.25*( cos(2.0*PI*elem_coords[0]) + cos(2.0*PI*elem_coords[1]) ) + 1.0;
                
                    // p = rho*ie*(gamma - 1)
                    size_t mat_id = f_id;
                    double gamma = global_vars(mat_id,0); // gamma value
                    elem_sie(rk_level, elem_gid) =
                                    elem_pres(elem_gid)/(mat_fill(f_id).den*(gamma - 1.0));
                } // end if

            } // end if fill
          
        }); // end FOR_ALL_CLASS element loop
        Kokkos::fence();
        
  
    } // end for loop over fills
    Kokkos::fence();
    
    
    // calculate the nodal mass
    // FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
        
    //     //node_mass(node_gid) = 0.0;
        
    //     if(num_dim==3){
            
            
    //     }// end if dims=3
    //     else {
            
            
    //     } // end else
        
    // }); // end FOR_ALL_CLASS
    //Kokkos::fence();
    
    // update host copies of arrays modified in this function
    elem_den.update_host();
    elem_mass.update_host();
    elem_sie.update_host();
    elem_stress.update_host();
    elem_pres.update_host();
    elem_sspd.update_host(); 

    return;
    
} // end of setup



/* ----------------------------------------------------------------------------
   solve function called by solver
------------------------------------------------------------------------------- */

int FEA_Module_RDH::solve(){
  rdh_solve();

  return 0;
}

/* ----------------------------------------------------------------------------
   RDH solver loop
------------------------------------------------------------------------------- */
void FEA_Module_RDH::rdh_solve(){
    Time_Variables tv = simparam.time_variables;
   
    const size_t rk_level = simparam.rk_num_bins - 1; 
    time_value = tv.time_initial;
    time_final = tv.time_final;
    dt_max = tv.dt_max;
    dt_min = tv.dt_min;
    dt_cfl = tv.dt_cfl;
    graphics_time = simparam.output_options.graphics_step;
    graphics_cyc_ival = simparam.graphics_options.graphics_cyc_ival;
    graphics_dt_ival = simparam.output_options.graphics_step;
    cycle_stop = tv.cycle_stop;
    rk_num_stages = simparam.rk_num_stages;
    dt = tv.dt;
    fuzz = tv.fuzz;
    tiny = tv.tiny;
    small = tv.small;
    graphics_times = simparam.graphics_options.graphics_times;
    graphics_id = simparam.graphics_options.graphics_id;
    size_t num_bdy_nodes = mesh->num_bdy_nodes;
    const DCArrayKokkos <boundary_t> boundary = simparam.boundary;
    const DCArrayKokkos <material_t> material = simparam.material;
    int old_max_forward_buffer;
    size_t cycle;
    const int num_dim = simparam.num_dims;

    int myrank = Explicit_Solver_Pointer_->myrank;
    CArrayKokkos <double> node_extensive_mass(nall_nodes, "node_extensive_mass");
    
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

    // ---- Calculate energy tallies ----
    double IE_tend = 0.0;
    double KE_tend = 0.0;
    double TE_tend = 0.0;

    double global_IE_tend = 0.0;
    double global_KE_tend = 0.0;
    double global_TE_tend = 0.0;

    int nlocal_elem_non_overlapping = Explicit_Solver_Pointer_->nlocal_elem_non_overlapping;
    
    // extensive IE
    REDUCE_SUM_CLASS(elem_gid, 0, nlocal_elem_non_overlapping, IE_loc_sum, {
        IE_loc_sum += elem_mass(elem_gid)*elem_sie(rk_level,elem_gid);
        
    }, IE_sum);
    IE_t0 = IE_sum;

    MPI_Allreduce(&IE_t0,&global_IE_t0,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
 
    // extensive KE
    REDUCE_SUM_CLASS(node_gid, 0, nlocal_nodes, KE_loc_sum, {
        
        double ke = 0;
        for (size_t dim=0; dim<num_dim; dim++){
            ke += node_vel(rk_level,node_gid,dim)*node_vel(rk_level,node_gid,dim); // 1/2 at end
        } // end for
        
        if(num_dim==2){
            KE_loc_sum += node_mass(node_gid)*node_coords(rk_level,node_gid,1)*ke;
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
    FOR_ALL_CLASS(node_gid, 0, nall_nodes, {
        
        double radius = 1.0;
        if(num_dim == 2){
            radius = node_coords(rk_level,node_gid,1);
        }
        node_extensive_mass(node_gid) = node_mass(node_gid)*radius;
        
    }); // end parallel for
    
    
    // a flag to exit the calculation
    size_t stop_calc=0;
    
    auto time_1 = std::chrono::high_resolution_clock::now();


    //communicate ghosts
    double comm_time1 = Explicit_Solver_Pointer_->CPU_Time();
            
    //active view scope; triggers host comms from updated data on device
    {
      const_host_vec_array node_velocities_host = Explicit_Solver_Pointer_->node_velocities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
      const_host_vec_array node_coords_host = Explicit_Solver_Pointer_->node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    }
    double comm_time2 = Explicit_Solver_Pointer_->CPU_Time();
    Explicit_Solver_Pointer_->dev2host_time += comm_time2-comm_time1;

    //communicate ghost velocities
    Explicit_Solver_Pointer_->comm_velocities();
    Explicit_Solver_Pointer_->comm_coordinates();
        
            
    double comm_time3 = Explicit_Solver_Pointer_->CPU_Time();

    //view scope
    {
      const_vec_array node_velocities_interface = Explicit_Solver_Pointer_->node_velocities_distributed->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      const_vec_array ghost_node_velocities_interface = Explicit_Solver_Pointer_->ghost_node_velocities_distributed->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      vec_array all_node_velocities_interface = Explicit_Solver_Pointer_->all_node_velocities_distributed->getLocalView<device_type> (Tpetra::Access::ReadWrite);
      const_vec_array node_coords_interface = Explicit_Solver_Pointer_->node_coords_distributed->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      const_vec_array ghost_node_coords_interface = Explicit_Solver_Pointer_->ghost_node_coords_distributed->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      vec_array all_node_coords_interface = Explicit_Solver_Pointer_->all_node_coords_distributed->getLocalView<device_type> (Tpetra::Access::ReadWrite);
      FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
        for (int idim = 0; idim < num_dim; idim++){
          all_node_velocities_interface(node_gid,idim) = node_velocities_interface(node_gid,idim);
          all_node_coords_interface(node_gid,idim) = node_coords_interface(node_gid,idim);
        }
      }); // end parallel for
      Kokkos::fence();

      FOR_ALL_CLASS(node_gid, nlocal_nodes, nlocal_nodes+nghost_nodes, {
        for (int idim = 0; idim < num_dim; idim++){
          all_node_velocities_interface(node_gid,idim) = ghost_node_velocities_interface(node_gid-nlocal_nodes,idim);
          all_node_coords_interface(node_gid,idim) = ghost_node_coords_interface(node_gid-nlocal_nodes,idim);
        }
      }); // end parallel for
      Kokkos::fence();
    } //end view scope
        
    
    auto time_2 = std::chrono::high_resolution_clock::now();
    auto time_difference = time_2 - time_1;
    //double calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count();
    double calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(time_difference).count();
    if(myrank==0)
      printf("\nCalculation time in seconds: %f \n", calc_time*1e-09);

    return;

} // end RDH solve
