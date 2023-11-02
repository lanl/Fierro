
#include "mesh.h"
#include "state.h"
#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>  // fmin, fmax, abs note: fminl is long
#include <mpi.h>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

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
#include "matar.h"
#include "utilities.h"
#include "node_combination.h"
#include "FEA_Module_SGH.h"
#include "Explicit_Solver.h"

//optimization
#include "ROL_Algorithm.hpp"
#include "ROL_Solver.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Types.hpp"
#include "ROL_Elementwise_Reduce.hpp"
#include "ROL_Stream.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_ParameterList.hpp"
#include <ROL_TpetraMultiVector.hpp>
#include "Kinetic_Energy_Minimize.h"

/* ----------------------------------------------------------------------
   Compute new system response due to the design variable update
------------------------------------------------------------------------- */

void FEA_Module_SGH::update_forward_solve(Teuchos::RCP<const MV> zp){
  const size_t rk_level = simparam.dynamic_options.rk_num_bins - 1;
  //local variable for host view in the dual view
  int num_dim = simparam.num_dims;
  int nodes_per_elem = max_nodes_per_element;
  int local_node_index, current_row, current_column;
  int max_stride = 0;
  int current_module_index;
  size_t access_index, row_access_index, row_counter;
  GO global_index, global_dof_index;
  LO local_dof_index;
  const size_t num_fills = simparam.regions.size();
  const size_t rk_num_bins = simparam.dynamic_options.rk_num_bins;
  const size_t num_bcs = fea_params.boundary_conditions.size();
  const size_t num_materials = simparam.materials.size();
  real_t objective_accumulation;

  // --- Read in the nodes in the mesh ---
  int myrank = Explicit_Solver_Pointer_->myrank;
  int nranks = Explicit_Solver_Pointer_->nranks;

  const DCArrayKokkos <mat_fill_t> mat_fill = simparam.mat_fill;
  const DCArrayKokkos <boundary_t> boundary = fea_params.boundary;
  const DCArrayKokkos <material_t> material = simparam.material;
  CArray<double> current_element_nodal_densities = CArray<double>(num_nodes_in_elem);
  
  std::vector<std::vector<int>> FEA_Module_My_TO_Modules = simparam.FEA_Module_My_TO_Modules;
  problem = Explicit_Solver_Pointer_->problem; //Pointer to ROL optimization problem object
  ROL::Ptr<ROL::Objective<real_t>> obj_pointer;

  //compute element averaged density ratios corresponding to nodal density design variables
  {//view scope
    const_host_vec_array all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    //debug print
    //std::cout << "NODE DENSITY TEST " << all_node_densities(0,0) << std::endl;
    for(int elem_id = 0; elem_id < rnum_elem; elem_id++){
      for(int inode = 0; inode < num_nodes_in_elem; inode++){
        current_element_nodal_densities(inode) = all_node_densities(nodes_in_elem(elem_id,inode),0);
      }
      relative_element_densities.host(elem_id) = average_element_density(num_nodes_in_elem, current_element_nodal_densities);
    }//for
  } //view scope
  //debug print
  //std::cout << "ELEMENT RELATIVE DENSITY TEST " << relative_element_densities.host(0) << std::endl;
  relative_element_densities.update_device();

  //set density vector to the current value chosen by the optimizer
  test_node_densities_distributed = zp;

  //reset nodal coordinates to initial values
  node_coords_distributed->assign(*initial_node_coords_distributed);

  //comms for ghosts
  Explicit_Solver_Pointer_->comm_coordinates();

  //view scope
  {
    const_vec_array node_coords_interface = Explicit_Solver_Pointer_->node_coords_distributed->getLocalView<device_type> (Tpetra::Access::ReadOnly);
    const_vec_array ghost_node_coords_interface = Explicit_Solver_Pointer_->ghost_node_coords_distributed->getLocalView<device_type> (Tpetra::Access::ReadOnly);
    vec_array all_node_coords_interface = Explicit_Solver_Pointer_->all_node_coords_distributed->getLocalView<device_type> (Tpetra::Access::ReadWrite);
    FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
      for (int idim = 0; idim < num_dim; idim++){
        all_node_coords_interface(node_gid,idim) = node_coords_interface(node_gid,idim);
      }
    }); // end parallel for
    Kokkos::fence();

    FOR_ALL_CLASS(node_gid, nlocal_nodes, nlocal_nodes+nghost_nodes, {
      for (int idim = 0; idim < num_dim; idim++){
        all_node_coords_interface(node_gid,idim) = ghost_node_coords_interface(node_gid-nlocal_nodes,idim);
      }
    }); // end parallel for
    Kokkos::fence();
  } //end view scope

  //reset velocities to initial conditions
  node_velocities_distributed->assign(*initial_node_velocities_distributed);

  //reset time accumulating objective and constraints
  /*
  for(int imodule = 0 ; imodule < FEA_Module_My_TO_Modules[my_fea_module_index_].size(); imodule++){
    current_module_index = FEA_Module_My_TO_Modules[my_fea_module_index_][imodule];
    //test if module needs reset
    if(){
      
    }
  }
  */
  //simple setup to just request KE for now; above loop to be expanded and used later for scanning modules
  obj_pointer = problem->getObjective();
  KineticEnergyMinimize_TopOpt& kinetic_energy_minimize_function = dynamic_cast<KineticEnergyMinimize_TopOpt&>(*obj_pointer);
  kinetic_energy_minimize_function.objective_accumulation = 0;

  //interface trial density vector

  //interfacing of vectors(should be removed later once made compatible)
  //view scope
  {
    host_vec_array interface_node_coords = Explicit_Solver_Pointer_->all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    for(size_t ibin = 0; ibin < rk_num_bins; ibin++){
      //save node data to node.coords
      //std::cout << "NODE DATA ON RANK " << myrank << std::endl;
      if(num_dim==2){
        for(int inode = 0; inode < nall_nodes; inode++){
          //std::cout << "Node index " << inode+1 << " ";
          node_coords.host(ibin,inode,0) = interface_node_coords(inode,0);
          //std::cout << host_node_coords_state(0,inode,0)+1<< " ";
          node_coords.host(ibin,inode,1) = interface_node_coords(inode,1);
          //std::cout << host_node_coords_state(0,inode,1)+1<< " ";
        }
      }
      else if(num_dim==3){
        for(int inode = 0; inode < nall_nodes; inode++){
          //std::cout << "Node index " << inode+1 << " ";
          node_coords.host(ibin,inode,0) = interface_node_coords(inode,0);
          //std::cout << host_node_coords_state(0,inode,0)+1<< " ";
          node_coords.host(ibin,inode,1) = interface_node_coords(inode,1);
          //std::cout << host_node_coords_state(0,inode,1)+1<< " ";
        
          node_coords.host(ibin,inode,2) = interface_node_coords(inode,2);
          //std::cout << host_node_coords_state(0,inode,2)+1<< std::endl;
        }
      }
    }
  } //end view scope

    // save the node coords to the current RK value
    for (size_t node_gid=0; node_gid < nall_nodes; node_gid++){
        
      for(int rk=1; rk<rk_num_bins; rk++){
        for (int dim = 0; dim < num_dim; dim++){
          node_coords.host(rk, node_gid, dim) = node_coords.host(0, node_gid, dim);
        } // end for dim
      } // end for rk
        
    } // end parallel for
    
    node_coords.update_device();

  //setup that needs repeating
  get_vol();
  //--- apply the fill instructions over the Elements---//
    
    // loop over the fill instructures
    //view scope
    {
    
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
            bool fill_this = mat_fill(f_id).volume.contains(elem_coords);
                 
            // paint the material state on the element
            if (fill_this){
                    
                // density
                elem_den(elem_gid) = mat_fill(f_id).den;

                //compute element average density from initial nodal density variables used as TO design variables
                elem_den(elem_gid) = elem_den(elem_gid)*relative_element_densities(elem_gid);
                
                // mass
                elem_mass(elem_gid) = elem_den(elem_gid)*elem_vol(elem_gid);
                
                // specific internal energy
                elem_sie(rk_level, elem_gid) = mat_fill(f_id).sie;
		
                elem_mat_id(elem_gid) = mat_fill(f_id).material_id;
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
    }//end view scope
    
   
    
    // apply BC's to velocity
    FEA_Module_SGH::boundary_velocity(*mesh, boundary, node_vel);
    
    
    // calculate the corner massess if 2D
    if(num_dim==2){
        
        FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
            
            // facial area of the corners
            double corner_areas_array[4];
            
            ViewCArrayKokkos <double> corner_areas(&corner_areas_array[0],4);
            ViewCArrayKokkos <size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 4);
            
            get_area_weights2D(corner_areas,
                               elem_gid,
                               node_coords,
                               elem_node_gids,
                               rk_level);
            
            // loop over the corners of the element and calculate the mass
            for (size_t corner_lid=0; corner_lid<4; corner_lid++){
                
                size_t corner_gid = corners_in_elem(elem_gid, corner_lid);
                corner_mass(corner_gid) = corner_areas(corner_lid)*elem_den(elem_gid); // node radius is added later
                
            } // end for over corners
        });
    
    } // end of
    
    
    // calculate the nodal mass
    FOR_ALL_CLASS(node_gid, 0, nall_nodes, {
        
        node_mass(node_gid) = 0.0;
        
        if(num_dim==3){
            
            for(size_t elem_lid=0; elem_lid < num_corners_in_node(node_gid); elem_lid++){
                size_t elem_gid = elems_in_node(node_gid,elem_lid);
                node_mass(node_gid) += 1.0/8.0*elem_mass(elem_gid);
            } // end for elem_lid
            
        }// end if dims=3
        else {
            
            // 2D-RZ
            for(size_t corner_lid=0; corner_lid < num_corners_in_node(node_gid); corner_lid++){
                
                size_t corner_gid = corners_in_node(node_gid, corner_lid);
                node_mass(node_gid) += corner_mass(corner_gid);  // sans the radius so it is areal node mass
                
                corner_mass(corner_gid) *= node_coords(rk_level,node_gid,1); // true corner mass now
            } // end for elem_lid
            
        } // end else
        
    }); // end FOR_ALL_CLASS
    Kokkos::fence();


    //current interface has differing mass arrays; this equates them until we unify memory
    //view scope
    {
      vec_array node_mass_interface = node_masses_distributed->getLocalView<device_type> (Tpetra::Access::ReadWrite);
      FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
        node_mass_interface(node_gid,0) = node_mass(node_gid);
      }); // end parallel for
    } //end view scope
    Kokkos::fence();
    //communicate ghost densities
    comm_node_masses();

    //this is forcing a copy to the device
    //view scope
    {
      vec_array ghost_node_mass_interface = ghost_node_masses_distributed->getLocalView<device_type> (Tpetra::Access::ReadWrite);

      FOR_ALL_CLASS(node_gid, nlocal_nodes, nall_nodes, {
        node_mass(node_gid) = ghost_node_mass_interface(node_gid-nlocal_nodes,0);
      }); // end parallel for
    } //end view scope
    Kokkos::fence();

    //update stiffness matrix
    if(simparam.topology_optimization_on||simparam.shape_optimization_on){
      //assemble_matrix();
    }
    
    //execute solve
    sgh_solve();

}

/* -------------------------------------------------------------------------------------------
   Compute average density of an element from nodal densities
---------------------------------------------------------------------------------------------- */

double FEA_Module_SGH::average_element_density(const int nodes_per_elem, const CArray<double> current_element_densities) const
{
  double result = 0;
  for(int i=0; i < nodes_per_elem; i++){
    result += current_element_densities(i)/nodes_per_elem;
  }

  return result;
}


/* ------------------------------------------------------------------------------
  Coupled adjoint problem for the kinetic energy minimization problem
--------------------------------------------------------------------------------- */

void FEA_Module_SGH::compute_topology_optimization_adjoint_full(){
  const size_t rk_level = simparam.dynamic_options.rk_num_bins - 1;
  size_t num_bdy_nodes = mesh->num_bdy_nodes;
  const DCArrayKokkos <boundary_t> boundary = fea_params.boundary;
  const DCArrayKokkos <material_t> material = simparam.material;
  const int num_dim = simparam.num_dims;
  real_t global_dt;
  size_t current_data_index, next_data_index;
  Teuchos::RCP<MV> previous_adjoint_vector_distributed, current_adjoint_vector_distributed, previous_velocity_vector_distributed, current_velocity_vector_distributed;
  Teuchos::RCP<MV> previous_phi_adjoint_vector_distributed, current_phi_adjoint_vector_distributed;
  //initialize first adjoint vector at last_time_step to 0 as the terminal value
  (*adjoint_vector_data)[last_time_step+1]->putScalar(0);
  (*phi_adjoint_vector_data)[last_time_step+1]->putScalar(0);

  //solve terminal value problem, proceeds in time backward. For simplicity, we use the same timestep data from the forward solve.
  //A linear interpolant is assumed between velocity data points; velocity midpoint is used to update the adjoint.
  if(myrank==0)
    std::cout << "Computing adjoint vector " << time_data.size() << std::endl;

  for (int cycle = last_time_step; cycle >= 0; cycle--) {
    //compute timestep from time data
    global_dt = time_data[cycle+1] - time_data[cycle];
    //print
    if(simparam.dynamic_options.output_time_sequence_level==TIME_OUTPUT_LEVEL::extreme){
      if (cycle==last_time_step){
        if(myrank==0)
          printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_data[cycle], global_dt);
      }
          // print time step every 10 cycles
      else if (cycle%1==0){
        if(myrank==0)
          printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_data[cycle], global_dt);
      } // end if
    }
    //else if (cycle==1){
      //if(myrank==0)
        //printf("cycle = %lu, time = %f, time step = %f \n", cycle-1, time_data[cycle-1], global_dt);
    //} // end if

    //compute adjoint vector for this data point; use velocity midpoint
    //view scope
    {
      
      //set velocity, internal energy, and position for this timestep
      const_vec_array previous_velocity_vector = (*forward_solve_velocity_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      const_vec_array current_velocity_vector = (*forward_solve_velocity_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      
      const_vec_array previous_coordinate_vector = (*forward_solve_coordinate_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      const_vec_array current_coordinate_vector = (*forward_solve_coordinate_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      
      const_vec_array previous_element_internal_energy = (*forward_solve_internal_energy_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      const_vec_array current_element_internal_energy = (*forward_solve_internal_energy_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);

      //interface of arrays for current implementation of force calculation

      FOR_ALL_CLASS(node_gid, 0, nlocal_nodes+nghost_nodes, {
        for (int idim = 0; idim < num_dim; idim++){
          node_vel(rk_level,node_gid,idim) = previous_velocity_vector(node_gid,idim);
          node_coords(rk_level,node_gid,idim) = previous_coordinate_vector(node_gid,idim);
        }
      });
      Kokkos::fence();

      FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
        elem_sie(rk_level,elem_gid) = previous_element_internal_energy(elem_gid,0);
      });
      Kokkos::fence();

      //set state according to phase data at this timestep

      get_vol();

      // ---- Calculate velocity diveregence for the element ----
      if(num_dim==2){
          get_divergence2D(elem_div,
                          *mesh,
                          node_coords,
                          node_vel,
                          elem_vol);
      }
      else {
          get_divergence(elem_div,
                        *mesh,
                        node_coords,
                        node_vel,
                        elem_vol);
      } // end if 2D

      // ---- Calculate elem state (den, pres, sound speed, stress) for next time step ----
      if(num_dim==2){
          update_state2D(material,
                          *mesh,
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
                          1.0,
                          cycle);
      }
      else{
          update_state(material,
                        *mesh,
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
                        1.0,
                        cycle);
      }

      //compute gradient matrices
      get_force_egradient_sgh(material,
                              *mesh,
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
                              1.0,
                              cycle);

      get_power_egradient_sgh(1.0,
                              *mesh,
                              node_vel,
                              node_coords,
                              elem_sie,
                              elem_mass,
                              corner_force);

      get_force_vgradient_sgh(material,
                              *mesh,
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
                              1.0,
                              cycle);
      
      get_power_vgradient_sgh(1.0,
                              *mesh,
                              node_vel,
                              node_coords,
                              elem_sie,
                              elem_mass,
                              corner_force);

      get_force_ugradient_sgh(material,
                              *mesh,
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
                              1.0,
                              cycle);

      get_power_ugradient_sgh(1.0,
                              *mesh,
                              node_vel,
                              node_coords,
                              elem_sie,
                              elem_mass,
                              corner_force);
      

      //force_gradient_velocity->describe(*fos,Teuchos::VERB_EXTREME);
      const_vec_array previous_force_gradient_position = force_gradient_position->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      //const_vec_array current_force_gradient_position = force_gradient_position->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      const_vec_array previous_force_gradient_velocity = force_gradient_velocity->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      //const_vec_array current_force_gradient_velocity = force_gradient_velocity->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      //compute gradient of force with respect to velocity
  
      const_vec_array previous_adjoint_vector = (*adjoint_vector_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      const_vec_array phi_previous_adjoint_vector =  (*phi_adjoint_vector_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      const_vec_array psi_previous_adjoint_vector =  (*psi_adjoint_vector_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      vec_array midpoint_adjoint_vector = adjoint_vector_distributed->getLocalView<device_type> (Tpetra::Access::ReadWrite);
      vec_array phi_midpoint_adjoint_vector =  phi_adjoint_vector_distributed->getLocalView<device_type> (Tpetra::Access::ReadWrite);
      vec_array psi_midpoint_adjoint_vector =  psi_adjoint_vector_distributed->getLocalView<device_type> (Tpetra::Access::ReadWrite);
      
      //half step update for RK2 scheme
      FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
        real_t rate_of_change;
        real_t matrix_contribution;
        size_t dof_id;
        size_t elem_id;
        for (int idim = 0; idim < num_dim; idim++){

          //EQUATION 1
          matrix_contribution = 0;
          //compute resulting row of force velocity gradient matrix transpose right multiplied by adjoint vector
          for(int idof = 0; idof < Gradient_Matrix_Strides(node_gid*num_dim+idim); idof++){
            dof_id = DOF_Graph_Matrix(node_gid*num_dim+idim,idof);
            matrix_contribution += previous_adjoint_vector(dof_id/num_dim,dof_id%num_dim)*Force_Gradient_Velocities(node_gid*num_dim+idim,idof);
          }

          //compute resulting row of transpose of power gradient w.r.t velocity matrix right multiplied by psi adjoint vector
          for(int ielem = 0; ielem < DOF_to_Elem_Matrix_Strides(node_gid*num_dim+idim); ielem++){
            elem_id = elems_in_node(node_gid,ielem);
            matrix_contribution += psi_previous_adjoint_vector(elem_id,0)*Power_Gradient_Velocities(node_gid*num_dim+idim,ielem);
          }
          
          rate_of_change = previous_velocity_vector(node_gid,idim)- 
                            matrix_contribution/node_mass(node_gid)-
                            phi_previous_adjoint_vector(node_gid,idim)/node_mass(node_gid);
          midpoint_adjoint_vector(node_gid,idim) = -rate_of_change*global_dt/2 + previous_adjoint_vector(node_gid,idim);

          //EQUATION 2
          matrix_contribution = 0;
          //compute resulting row of force displacement gradient matrix transpose right multiplied by adjoint vector
          for(int idof = 0; idof < Gradient_Matrix_Strides(node_gid*num_dim+idim); idof++){
            dof_id = DOF_Graph_Matrix(node_gid*num_dim+idim,idof);
            matrix_contribution += previous_adjoint_vector(dof_id/num_dim,dof_id%num_dim)*Force_Gradient_Positions(node_gid*num_dim+idim,idof);
          }

          //compute resulting row of transpose of power gradient w.r.t displacement matrix right multiplied by psi adjoint vector
          for(int ielem = 0; ielem < DOF_to_Elem_Matrix_Strides(node_gid*num_dim+idim); ielem++){
            elem_id = elems_in_node(node_gid,ielem);
            matrix_contribution += psi_previous_adjoint_vector(elem_id,0)*Power_Gradient_Positions(node_gid*num_dim+idim,ielem);
          }
          
          rate_of_change = -matrix_contribution;
          //rate_of_change = -0.0000001*previous_adjoint_vector(node_gid,idim);
          phi_midpoint_adjoint_vector(node_gid,idim) = -rate_of_change*global_dt/2 + phi_previous_adjoint_vector(node_gid,idim);

        } 
      }); // end parallel for
      Kokkos::fence();

      //half step update for RK2 scheme
      FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
        real_t rate_of_change;
        real_t matrix_contribution;
        size_t dof_id;
        size_t elem_id;
        //EQUATION 3
        matrix_contribution = 0;
        //compute resulting row of force displacement gradient matrix transpose right multiplied by adjoint vector
        for(int idof = 0; idof < num_nodes_in_elem*num_dim; idof++){
          dof_id = nodes_in_elem(elem_gid,idof/num_dim)*num_dim + idof%num_dim;
          matrix_contribution += previous_adjoint_vector(dof_id/num_dim,dof_id%num_dim)*Force_Gradient_Energies(elem_gid,idof);
        }
        rate_of_change = -(matrix_contribution + psi_previous_adjoint_vector(elem_gid,0)*Power_Gradient_Energies(elem_gid))/elem_mass(elem_gid);
        //rate_of_change = -0.0000001*previous_adjoint_vector(node_gid,idim);
        psi_midpoint_adjoint_vector(elem_gid,0) = -rate_of_change*global_dt/2 + psi_previous_adjoint_vector(elem_gid,0);
      }); // end parallel for
      Kokkos::fence();

      boundary_adjoint(*mesh, boundary,midpoint_adjoint_vector, phi_midpoint_adjoint_vector, psi_midpoint_adjoint_vector);
      comm_adjoint_vectors(cycle);
      //swap names to get ghost nodes for the midpoint vectors
      vec_array current_adjoint_vector = adjoint_vector_distributed->getLocalView<device_type> (Tpetra::Access::ReadWrite);
      vec_array phi_current_adjoint_vector = phi_adjoint_vector_distributed->getLocalView<device_type> (Tpetra::Access::ReadWrite);
      vec_array psi_current_adjoint_vector = psi_adjoint_vector_distributed->getLocalView<device_type> (Tpetra::Access::ReadWrite);
      midpoint_adjoint_vector = (*adjoint_vector_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadWrite);
      phi_midpoint_adjoint_vector =  (*phi_adjoint_vector_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadWrite);
      psi_midpoint_adjoint_vector =  (*psi_adjoint_vector_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadWrite);

      
      //full step update with midpoint gradient for RK2 scheme
      FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
        real_t rate_of_change;
        real_t matrix_contribution;
        size_t dof_id;
        size_t elem_id;
        for (int idim = 0; idim < num_dim; idim++){
          //EQUATION 1
          matrix_contribution = 0;
          //compute resulting row of force velocity gradient matrix transpose right multiplied by adjoint vector
          
          for(int idof = 0; idof < Gradient_Matrix_Strides(node_gid*num_dim+idim); idof++){
            dof_id = DOF_Graph_Matrix(node_gid*num_dim+idim,idof);
            matrix_contribution += midpoint_adjoint_vector(dof_id/num_dim,dof_id%num_dim)*Force_Gradient_Velocities(node_gid*num_dim+idim,idof);
          }

          //compute resulting row of transpose of power gradient w.r.t velocity matrix right multiplied by psi adjoint vector
          for(int ielem = 0; ielem < DOF_to_Elem_Matrix_Strides(node_gid*num_dim+idim); ielem++){
            elem_id = elems_in_node(node_gid,ielem);
            matrix_contribution += psi_midpoint_adjoint_vector(elem_id,0)*Power_Gradient_Velocities(node_gid*num_dim+idim,ielem);
          }
          
          rate_of_change =  (previous_velocity_vector(node_gid,idim) + current_velocity_vector(node_gid,idim))/2- 
                            matrix_contribution/node_mass(node_gid)-
                            phi_midpoint_adjoint_vector(node_gid,idim)/node_mass(node_gid);
          current_adjoint_vector(node_gid,idim) = -rate_of_change*global_dt + previous_adjoint_vector(node_gid,idim);

          //EQUATION 2
          matrix_contribution = 0;
          //compute resulting row of force displacement gradient matrix transpose right multiplied by adjoint vector
          for(int idof = 0; idof < Gradient_Matrix_Strides(node_gid*num_dim+idim); idof++){
            dof_id = DOF_Graph_Matrix(node_gid*num_dim+idim,idof);
            matrix_contribution += midpoint_adjoint_vector(dof_id/num_dim,dof_id%num_dim)*Force_Gradient_Positions(node_gid*num_dim+idim,idof);
          }

          //compute resulting row of transpose of power gradient w.r.t displacement matrix right multiplied by psi adjoint vector
          for(int ielem = 0; ielem < DOF_to_Elem_Matrix_Strides(node_gid*num_dim+idim); ielem++){
            elem_id = elems_in_node(node_gid,ielem);
            matrix_contribution += psi_midpoint_adjoint_vector(elem_id,0)*Power_Gradient_Positions(node_gid*num_dim+idim,ielem);
          }

          rate_of_change = -matrix_contribution;
          //rate_of_change = -0.0000001*midpoint_adjoint_vector(node_gid,idim);
          phi_current_adjoint_vector(node_gid,idim) = -rate_of_change*global_dt + phi_previous_adjoint_vector(node_gid,idim);
        }
      }); // end parallel for
      Kokkos::fence();
      //half step update for RK2 scheme
      FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
        real_t rate_of_change;
        real_t matrix_contribution;
        size_t dof_id;
        size_t elem_id;
        //EQUATION 3
        matrix_contribution = 0;
        //compute resulting row of force displacement gradient matrix transpose right multiplied by adjoint vector
        for(int idof = 0; idof < num_nodes_in_elem*num_dim; idof++){
          dof_id = nodes_in_elem(elem_gid,idof/num_dim)*num_dim + idof%num_dim;
          matrix_contribution += midpoint_adjoint_vector(dof_id/num_dim,dof_id%num_dim)*Force_Gradient_Energies(elem_gid,idof);
        }
        rate_of_change = -(matrix_contribution + psi_midpoint_adjoint_vector(elem_gid,0)*Power_Gradient_Energies(elem_gid))/elem_mass(elem_gid);
        //rate_of_change = -0.0000001*previous_adjoint_vector(node_gid,idim);
        psi_current_adjoint_vector(elem_gid,0) = -rate_of_change*global_dt/2 + psi_midpoint_adjoint_vector(elem_gid,0);
      }); // end parallel for
      Kokkos::fence();

      
      boundary_adjoint(*mesh, boundary,current_adjoint_vector, phi_current_adjoint_vector, psi_midpoint_adjoint_vector);

    } //end view scope
    

    comm_adjoint_vectors(cycle);
    //phi_adjoint_vector_distributed->describe(*fos,Teuchos::VERB_EXTREME);
    
      //debug
      if(cycle==0)
      std::cout << "REACHED THIS STEP OF ADJOINT" << std::endl;
  }
}

/* ----------------------------------------------------------------------------
   Gradient for the (unsimplified) kinetic energy minimization problem
------------------------------------------------------------------------------- */

void FEA_Module_SGH::compute_topology_optimization_gradient_full(Teuchos::RCP<const MV> design_densities_distributed, Teuchos::RCP<MV> design_gradients_distributed){

  size_t num_bdy_nodes = mesh->num_bdy_nodes;
  const DCArrayKokkos <boundary_t> boundary = fea_params.boundary;
  const DCArrayKokkos <material_t> material = simparam.material;
  const int num_dim = simparam.num_dims;
  int num_corners = rnum_elem*num_nodes_in_elem;
  real_t global_dt;
  bool element_constant_density = true;
  size_t current_data_index, next_data_index;
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_element_velocities = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(num_nodes_in_elem,num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_element_adjoint = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(num_nodes_in_elem,num_dim);

  if(myrank==0)
    std::cout << "Computing accumulated kinetic energy gradient" << std::endl;

  compute_topology_optimization_adjoint_full();

  

  { //view scope
    vec_array design_gradients = design_gradients_distributed->getLocalView<device_type> (Tpetra::Access::ReadWrite);
    const_vec_array design_densities = design_densities_distributed->getLocalView<device_type> (Tpetra::Access::ReadOnly);
    //compute design gradients
    FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
      design_gradients(node_id,0) *=0.5/(double)num_nodes_in_elem/(double)num_nodes_in_elem;
      //design_gradients(node_id,0) =0.00001;
    }); // end parallel for
    Kokkos::fence();

    //gradient contribution from time derivative of adjoint \dot{lambda}(dM/drho)v product.
    if(simparam.dynamic_options.output_time_sequence_level==TIME_OUTPUT_LEVEL::extreme){
        if(myrank==0){
          std::cout << "gradient term involving adjoint derivative" << std::endl;
        }
    }

    for (unsigned long cycle = 0; cycle < last_time_step+1; cycle++) {
      //compute timestep from time data
      global_dt = time_data[cycle+1] - time_data[cycle];
      //print
      if(simparam.dynamic_options.output_time_sequence_level==TIME_OUTPUT_LEVEL::extreme){

        if (cycle==0){
          if(myrank==0)
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_data[cycle], global_dt);
        }
            // print time step every 10 cycles
        else if (cycle%20==0){
          if(myrank==0)
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_data[cycle], global_dt);
        } // end if
      }

      //compute adjoint vector for this data point; use velocity midpoint
        //view scope
        {
          const_vec_array current_velocity_vector = (*forward_solve_velocity_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
          const_vec_array current_adjoint_vector = (*adjoint_vector_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
          const_vec_array next_velocity_vector = (*forward_solve_velocity_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
          const_vec_array next_adjoint_vector = (*adjoint_vector_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
          
          FOR_ALL_CLASS(elem_id, 0, rnum_elem, {
            real_t lambda_dot_current;
            real_t lambda_dot_next;
            size_t node_id;
            size_t corner_id;
            real_t inner_product;
            //std::cout << elem_mass(elem_id) <<std::endl;
            //current_nodal_velocities
            for (int inode = 0; inode < num_nodes_in_elem; inode++){
              node_id = nodes_in_elem(elem_id, inode);
              //midpoint rule for integration being used; add velocities and divide by 2
              current_element_velocities(inode,0) = (current_velocity_vector(node_id,0) + next_velocity_vector(node_id,0))/2;
              current_element_velocities(inode,1) = (current_velocity_vector(node_id,1) + next_velocity_vector(node_id,1))/2;
              if(num_dim==3)
              current_element_velocities(inode,2) = (current_velocity_vector(node_id,2) + next_velocity_vector(node_id,2))/2;
            }

            inner_product = 0;
            for(int ifill=0; ifill < num_nodes_in_elem; ifill++){
              node_id = nodes_in_elem(elem_id, ifill);
              for(int idim=0; idim < num_dim; idim++){
                lambda_dot_current = lambda_dot_next = (next_adjoint_vector(node_id,idim)-current_adjoint_vector(node_id,idim))/global_dt;
                //lambda_dot_current = current_velocity_vector(node_id,idim) + damping_constant*current_adjoint_vector(node_id,idim)/node_mass(node_id) - current_phi_adjoint_vector(node_id,idim)/node_mass(node_id);
                //lambda_dot_next = next_velocity_vector(node_id,idim) + damping_constant*next_adjoint_vector(node_id,idim)/node_mass(node_id) - next_phi_adjoint_vector(node_id,idim)/node_mass(node_id);
                inner_product += elem_mass(elem_id)*(lambda_dot_current+lambda_dot_current)*current_element_velocities(ifill,idim)/2;
              }
            }

            for (int inode = 0; inode < num_nodes_in_elem; inode++){
              //compute gradient of local element contribution to v^t*M*v product
              corner_id = elem_id*num_nodes_in_elem + inode;
              corner_value_storage(corner_id) = inner_product*global_dt/relative_element_densities(elem_id);
            }
            
          }); // end parallel for
          Kokkos::fence();
          
          //accumulate node values from corner storage
          //multiply
          FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
            size_t corner_id;
            for(int icorner=0; icorner < num_corners_in_node(node_id); icorner++){
              corner_id = corners_in_node(node_id,icorner);
              design_gradients(node_id,0) += -corner_value_storage(corner_id)/(double)num_nodes_in_elem/(double)num_nodes_in_elem;
            }
          }); // end parallel for
          Kokkos::fence();
          
        } //end view scope
    }

    //gradient contribution from time derivative of psi_adjoint \dot{psi}(dM_E/drho)e product.
    if(simparam.dynamic_options.output_time_sequence_level==TIME_OUTPUT_LEVEL::extreme){
        if(myrank==0){
          std::cout << "gradient term involving adjoint derivative" << std::endl;
        }
    }

    for (unsigned long cycle = 0; cycle < last_time_step+1; cycle++) {
      //compute timestep from time data
      global_dt = time_data[cycle+1] - time_data[cycle];
      //print
      if(simparam.dynamic_options.output_time_sequence_level==TIME_OUTPUT_LEVEL::extreme){

        if (cycle==0){
          if(myrank==0)
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_data[cycle], global_dt);
        }
            // print time step every 10 cycles
        else if (cycle%20==0){
          if(myrank==0)
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_data[cycle], global_dt);
        } // end if
      }

      //compute adjoint vector for this data point; use velocity midpoint
        //view scope
        {
          const_vec_array current_element_internal_energy = (*forward_solve_internal_energy_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
          const_vec_array current_psi_adjoint_vector = (*psi_adjoint_vector_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
          const_vec_array next_element_internal_energy = (*forward_solve_internal_energy_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
          const_vec_array next_psi_adjoint_vector = (*psi_adjoint_vector_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
          
          FOR_ALL_CLASS(elem_id, 0, rnum_elem, {
            real_t psi_dot_current;
            real_t psi_dot_next;
            size_t node_id;
            size_t corner_id;
            real_t inner_product;
            //std::cout << elem_mass(elem_id) <<std::endl;

            psi_dot_current = (next_psi_adjoint_vector(elem_id,0) - current_psi_adjoint_vector(elem_id,0))/global_dt;
            inner_product = elem_mass(elem_id)*(psi_dot_current+psi_dot_current)*current_element_internal_energy(elem_id,0)/2;

            for (int inode = 0; inode < num_nodes_in_elem; inode++){
              //compute gradient of local element contribution to v^t*M*v product
              corner_id = elem_id*num_nodes_in_elem + inode;
              corner_value_storage(corner_id) = inner_product*global_dt/relative_element_densities(elem_id);
            }
            
          }); // end parallel for
          Kokkos::fence();
          
          //accumulate node values from corner storage
          //multiply
          FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
            size_t corner_id;
            for(int icorner=0; icorner < num_corners_in_node(node_id); icorner++){
              corner_id = corners_in_node(node_id,icorner);
              design_gradients(node_id,0) += -corner_value_storage(corner_id)/(double)num_nodes_in_elem;
            }
          }); // end parallel for
          Kokkos::fence();
          
        } //end view scope
    }

    //compute initial condition contribution from velocities
    //view scope
    {
      const_vec_array current_velocity_vector = (*forward_solve_velocity_data)[0]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      const_vec_array current_adjoint_vector = (*adjoint_vector_data)[0]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      
      FOR_ALL_CLASS(elem_id, 0, rnum_elem, {
        real_t lambda_dot;
        size_t node_id;
        size_t corner_id;
        real_t inner_product;
        //std::cout << elem_mass(elem_id) <<std::endl;
        //current_nodal_velocities
        for (int inode = 0; inode < num_nodes_in_elem; inode++){
          node_id = nodes_in_elem(elem_id, inode);
          //midpoint rule for integration being used; add velocities and divide by 2
          current_element_velocities(inode,0) = current_velocity_vector(node_id,0);
          current_element_velocities(inode,1) = current_velocity_vector(node_id,1);
          if(num_dim==3)
          current_element_velocities(inode,2) = current_velocity_vector(node_id,2);
        }

        inner_product = 0;
        for(int ifill=0; ifill < num_nodes_in_elem; ifill++){
          node_id = nodes_in_elem(elem_id, ifill);
          for(int idim=0; idim < num_dim; idim++){
            inner_product += elem_mass(elem_id)*current_adjoint_vector(node_id,idim)*current_element_velocities(ifill,idim);
          }
        }

        for (int inode = 0; inode < num_nodes_in_elem; inode++){
          //compute gradient of local element contribution to v^t*M*v product
          corner_id = elem_id*num_nodes_in_elem + inode;
          corner_value_storage(corner_id) = inner_product/relative_element_densities(elem_id);
        }
        
      }); // end parallel for
      Kokkos::fence();
      
      //accumulate node values from corner storage
      //multiply
      FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
        size_t corner_id;
        for(int icorner=0; icorner < num_corners_in_node(node_id); icorner++){
          corner_id = corners_in_node(node_id,icorner);
          design_gradients(node_id,0) += -corner_value_storage(corner_id)/(double)num_nodes_in_elem/(double)num_nodes_in_elem;
        }
      }); // end parallel for
      Kokkos::fence();
      
    } //end view scope

    //compute initial condition contribution from internal energies
    //view scope
    {
      const_vec_array current_element_internal_energy = (*forward_solve_internal_energy_data)[0]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      const_vec_array current_psi_adjoint_vector = (*psi_adjoint_vector_data)[0]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      
      FOR_ALL_CLASS(elem_id, 0, rnum_elem, {
        real_t lambda_dot;
        size_t node_id;
        size_t corner_id;
        real_t inner_product;
        //std::cout << elem_mass(elem_id) <<std::endl;
        inner_product = elem_mass(elem_id)*current_psi_adjoint_vector(elem_id,0)*current_element_internal_energy(elem_id,0);

        for (int inode = 0; inode < num_nodes_in_elem; inode++){
          //compute gradient of local element contribution to v^t*M*v product
          corner_id = elem_id*num_nodes_in_elem + inode;
          corner_value_storage(corner_id) = inner_product/relative_element_densities(elem_id);
        }
        
      }); // end parallel for
      Kokkos::fence();
      
      //accumulate node values from corner storage
      //multiply
      FOR_ALL_CLASS(node_id, 0, nlocal_nodes, {
        size_t corner_id;
        for(int icorner=0; icorner < num_corners_in_node(node_id); icorner++){
          corner_id = corners_in_node(node_id,icorner);
          design_gradients(node_id,0) += -corner_value_storage(corner_id)/(double)num_nodes_in_elem;
        }
      }); // end parallel for
      Kokkos::fence();
      
    } //end view scope

  }//end view scope design gradients
  
  //view scope
  {
    //host_vec_array host_design_gradients = design_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    //const_host_vec_array host_design_variables = design_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    vec_array design_gradients = design_gradients_distributed->getLocalView<device_type> (Tpetra::Access::ReadWrite);
    const_vec_array design_variables = design_densities_distributed->getLocalView<device_type> (Tpetra::Access::ReadOnly);
    force_design_gradient_term(design_variables, design_gradients);
    power_design_gradient_term(design_variables, design_gradients);
  }//end view scope

}

/* ----------------------------------------------------------------------
   Initialize global vectors and array maps needed for matrix assembly
------------------------------------------------------------------------- */
void FEA_Module_SGH::init_assembly(){
  int num_dim = simparam.num_dims;
  //const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  Gradient_Matrix_Strides = DCArrayKokkos<size_t, array_layout, device_type, memory_traits> (nlocal_nodes*num_dim, "Gradient_Matrix_Strides");
  DOF_to_Elem_Matrix_Strides = DCArrayKokkos<size_t, array_layout, device_type, memory_traits> (nlocal_nodes*num_dim, "Gradient_Matrix_Strides");
  Elem_to_Elem_Matrix_Strides = DCArrayKokkos<size_t, array_layout, device_type, memory_traits> (rnum_elem, "Gradient_Matrix_Strides");
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Graph_Fill(nall_nodes, "nall_nodes");
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> current_row_nodes_scanned;
  CArrayKokkos <size_t> count_saved_corners_in_node(nall_nodes, "count_saved_corners_in_node");
  int local_node_index, current_column_index;
  size_t max_stride = 0;
  size_t nodes_per_element;
  nodal_density_flag = simparam.nodal_density_flag;
  penalty_power = simparam.optimization_options.simp_penalty_power;
  
  //allocate stride arrays
  CArrayKokkos <size_t, array_layout, device_type, memory_traits> Graph_Matrix_Strides_initial(nlocal_nodes, "Graph_Matrix_Strides_initial");
  DCArrayKokkos <size_t, array_layout, device_type, memory_traits> Dual_Graph_Matrix_Strides_initial(nlocal_nodes, "Host_Graph_Matrix_Strides_initial");
  Graph_Matrix_Strides = DCArrayKokkos<size_t, array_layout, device_type, memory_traits>(nlocal_nodes, "Graph_Matrix_Strides");

  //allocate storage for the sparse gradient matrix map for node to node connectivity
  Global_Gradient_Matrix_Assembly_Map = DCArrayKokkos<size_t, array_layout, device_type, memory_traits>(rnum_elem,
                                         max_nodes_per_element,max_nodes_per_element, "Global_Gradient_Matrix_Assembly_Map");
  
  //allocate storage for the sparse gradient matrix map for node to element connectivity
  Element_Gradient_Matrix_Assembly_Map = DCArrayKokkos<size_t, array_layout, device_type, memory_traits>(rnum_elem,
                                         max_nodes_per_element, "Element_Gradient_Matrix_Assembly_Map");

  //allocate array used to determine global node repeats in the sparse graph later
  DCArrayKokkos <int, array_layout, device_type, memory_traits> node_indices_used(nall_nodes, "node_indices_used");

  /*allocate array that stores which column the node index occured on for the current row
    when removing repeats*/
  DCArrayKokkos <size_t, array_layout, device_type, memory_traits> column_index(nall_nodes, "column_index");
  
  //initialize nlocal arrays
  FOR_ALL_CLASS(inode, 0, nlocal_nodes, {
    Graph_Matrix_Strides_initial(inode) = 0;
    Graph_Matrix_Strides(inode) = 0;
    Graph_Fill(inode) = 0;
  }); // end parallel for
  Kokkos::fence();

  //initialize nall arrays
  //initialize nlocal arrays
  FOR_ALL_CLASS(inode, 0, nall_nodes, {
    node_indices_used(inode) = 0;
    column_index(inode) = 0;
    count_saved_corners_in_node(inode) = 0;
  }); // end parallel for
  Kokkos::fence();
  
  //count upper bound of strides for Sparse Pattern Graph by allowing repeats due to connectivity
    if(num_dim == 2)
    for (int ielem = 0; ielem < rnum_elem; ielem++){
      element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
      nodes_per_element = elem2D->num_nodes();
      for (int lnode = 0; lnode < nodes_per_element; lnode++){
        local_node_index = nodes_in_elem(ielem, lnode);
        if(local_node_index < nlocal_nodes){
          Dual_Graph_Matrix_Strides_initial.host(local_node_index) += nodes_per_element;
        }
      }
    }

    if(num_dim == 3)
    for (int ielem = 0; ielem < rnum_elem; ielem++){
      element_select->choose_3Delem_type(Element_Types(ielem), elem);
      nodes_per_element = elem->num_nodes();
      for (int lnode = 0; lnode < nodes_per_element; lnode++){
        local_node_index = nodes_in_elem(ielem, lnode);
        if(local_node_index < nlocal_nodes){
          Dual_Graph_Matrix_Strides_initial.host(local_node_index) += nodes_per_element;
        }
      }
    }
  
  Dual_Graph_Matrix_Strides_initial.update_device();

  //equate strides for later
  FOR_ALL_CLASS(inode, 0, nlocal_nodes, {
    Graph_Matrix_Strides(inode) = Graph_Matrix_Strides_initial(inode) = Dual_Graph_Matrix_Strides_initial(inode);
  }); // end parallel for
  
  //for (int inode = 0; inode < nlocal_nodes; inode++)
    //std::cout << Graph_Matrix_Strides_initial(inode) << std::endl;

  //compute maximum stride
  size_t update = 0;
  REDUCE_MAX_CLASS(inode, 0, nlocal_nodes, update, {
      if(update < Graph_Matrix_Strides_initial(inode))
        update = Graph_Matrix_Strides_initial(inode);
  }, max_stride);
  
  //std::cout << "THE MAX STRIDE" << max_stride << std::endl;
  //allocate array used in the repeat removal process
  current_row_nodes_scanned = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(max_stride, "current_row_nodes_scanned");

  //allocate sparse graph with node repeats
  RaggedRightArrayKokkos<size_t, array_layout, device_type, memory_traits> Repeat_Graph_Matrix(Graph_Matrix_Strides_initial);
  RaggedRightArrayofVectorsKokkos<size_t, array_layout, device_type, memory_traits> Element_local_indices(Graph_Matrix_Strides_initial,num_dim);
  
  //Fill the initial Graph with repeats
  if(num_dim == 2){
    for (int ielem = 0; ielem < rnum_elem; ielem++){
      element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
      nodes_per_element = elem2D->num_nodes();
      for (int lnode = 0; lnode < nodes_per_element; lnode++){
        local_node_index = nodes_in_elem(ielem, lnode);
        if(local_node_index < nlocal_nodes){
          for (int jnode = 0; jnode < nodes_per_element; jnode++){
            current_column_index = Graph_Fill(local_node_index)+jnode;
            Repeat_Graph_Matrix(local_node_index, current_column_index) = nodes_in_elem(ielem,jnode);

            //fill inverse map
            Element_local_indices(local_node_index,current_column_index,0) = ielem;
            Element_local_indices(local_node_index,current_column_index,1) = lnode;
            Element_local_indices(local_node_index,current_column_index,2) = jnode;

            //fill forward map
            Global_Gradient_Matrix_Assembly_Map(ielem,lnode,jnode) = current_column_index;
          }
          Graph_Fill(local_node_index) += nodes_per_element;
        }
      }
    }
  }
  
  if(num_dim == 3){
    for (int ielem = 0; ielem < rnum_elem; ielem++){
      element_select->choose_3Delem_type(Element_Types(ielem), elem);
      nodes_per_element = elem->num_nodes();
      for (int lnode = 0; lnode < nodes_per_element; lnode++){
        local_node_index = nodes_in_elem(ielem, lnode);
        if(local_node_index < nlocal_nodes){
          for (int jnode = 0; jnode < nodes_per_element; jnode++){
            current_column_index = Graph_Fill(local_node_index)+jnode;
            Repeat_Graph_Matrix(local_node_index, current_column_index) = nodes_in_elem(ielem,jnode);

            //fill inverse map
            Element_local_indices(local_node_index,current_column_index,0) = ielem;
            Element_local_indices(local_node_index,current_column_index,1) = lnode;
            Element_local_indices(local_node_index,current_column_index,2) = jnode;

            //fill forward map
            Global_Gradient_Matrix_Assembly_Map(ielem,lnode,jnode) = current_column_index;
          }
          Graph_Fill(local_node_index) += nodes_per_element;
        }
      }
    }
  }
    
  //debug statement
  //std::cout << "started run" << std::endl;
  //std::cout << "Graph Matrix Strides Repeat on task " << myrank << std::endl;
  //for (int inode = 0; inode < nlocal_nodes; inode++)
    //std::cout << Graph_Matrix_Strides(inode) << std::endl;
  RUN_CLASS({
  //remove repeats from the inital graph setup
    int current_node;
    //remove repeats from the inital graph setup
    int current_element_index;
    int element_row_index;
    int element_column_index;
    int current_stride;
    int current_row_n_nodes_scanned;
    for (int inode = 0; inode < nlocal_nodes; inode++){
      current_row_n_nodes_scanned = 0;
      for (int istride = 0; istride < Graph_Matrix_Strides(inode); istride++){
        //convert global index in graph to its local index for the flagging array
        current_node = Repeat_Graph_Matrix(inode,istride);
        //debug
        //if(current_node==-1)
        //std::cout << "Graph Matrix node access on task " << myrank << std::endl;
        //std::cout << Repeat_Graph_Matrix(inode,istride) << std::endl;
        if(node_indices_used(current_node)){
          //set global assembly map index to the location in the graph matrix where this global node was first found
          current_element_index = Element_local_indices(inode,istride,0);
          element_row_index = Element_local_indices(inode,istride,1);
          element_column_index = Element_local_indices(inode,istride,2);
          Global_Gradient_Matrix_Assembly_Map(current_element_index,element_row_index, element_column_index) 
              = column_index(current_node);   

          
          //swap current node with the end of the current row and shorten the stride of the row
          //first swap information about the inverse and forward maps

          current_stride = Graph_Matrix_Strides(inode);
          if(istride!=current_stride-1){
          Element_local_indices(inode,istride,0) = Element_local_indices(inode,current_stride-1,0);
          Element_local_indices(inode,istride,1) = Element_local_indices(inode,current_stride-1,1);
          Element_local_indices(inode,istride,2) = Element_local_indices(inode,current_stride-1,2);
          current_element_index = Element_local_indices(inode,istride,0);
          element_row_index = Element_local_indices(inode,istride,1);
          element_column_index = Element_local_indices(inode,istride,2);

          Global_Gradient_Matrix_Assembly_Map(current_element_index,element_row_index, element_column_index) 
              = istride;

          //now that the element map information has been copied, copy the global node index and delete the last index

          Repeat_Graph_Matrix(inode,istride) = Repeat_Graph_Matrix(inode,current_stride-1);
          }
          istride--;
          Graph_Matrix_Strides(inode)--;
        }
        else{
          /*this node hasn't shown up in the row before; add it to the list of nodes
            that have been scanned uniquely. Use this list to reset the flag array
            afterwards without having to loop over all the nodes in the system*/
          node_indices_used(current_node) = 1;
          column_index(current_node) = istride;
          current_row_nodes_scanned(current_row_n_nodes_scanned) = current_node;
          current_row_n_nodes_scanned++;
        }
      }
      //reset nodes used list for the next row of the sparse list
      for(int node_reset = 0; node_reset < current_row_n_nodes_scanned; node_reset++)
        node_indices_used(current_row_nodes_scanned(node_reset)) = 0;

    }
  });
  Kokkos::fence();  
    
  Graph_Matrix_Strides.update_host();
  //copy reduced content to non_repeat storage
  Graph_Matrix = RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits>(Graph_Matrix_Strides);

  FOR_ALL_CLASS(inode, 0, nlocal_nodes, {
    for(int istride = 0; istride < Graph_Matrix_Strides(inode); istride++){
      Graph_Matrix(inode,istride) = Repeat_Graph_Matrix(inode,istride);
    }
  }); // end parallel for

  //deallocate repeat matrix
  
  /*At this stage the sparse graph should have unique global indices on each row.
    The constructed Assembly map (to the global sparse matrix)
    is used to loop over each element's local stiffness matrix in the assembly process.*/
  
  //expand strides for stiffness matrix by multipling by dim
  FOR_ALL_CLASS(idof, 0, num_dim*nlocal_nodes, {
    Gradient_Matrix_Strides(idof) = num_dim*Graph_Matrix_Strides(idof/num_dim);
  }); // end parallel for

  Gradient_Matrix_Strides.update_host();
  
  DOF_Graph_Matrix = RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits>(Gradient_Matrix_Strides);
  Force_Gradient_Positions = RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout>(Gradient_Matrix_Strides);
  Force_Gradient_Velocities = RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout>(Gradient_Matrix_Strides);
  //needs different graph of node to elem rather than node to node
  //DOF_to_Elem_Matrix_Strides.get_kokkos_dual_view().d_view = elems_in_node.mystrides_;
  Kokkos::View<size_t *,array_layout, device_type, memory_traits> node_to_elem_strides = elems_in_node.mystrides_;
  DOF_to_Elem_Matrix_Strides = DCArrayKokkos<size_t, array_layout, device_type, memory_traits>(nlocal_nodes*num_dim);
  FOR_ALL_CLASS(inode, 0, nlocal_nodes, {
    for(int idim = 0; idim < num_dim; idim++){
      DOF_to_Elem_Matrix_Strides(inode*num_dim + idim) = node_to_elem_strides(inode);
    }
  }); // end parallel for
  DOF_to_Elem_Matrix_Strides.update_host();
  Force_Gradient_Energies = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(rnum_elem, num_nodes_in_elem*num_dim);
  Power_Gradient_Energies = CArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits>(rnum_elem);
  Power_Gradient_Positions = RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout>(DOF_to_Elem_Matrix_Strides);
  Power_Gradient_Velocities = RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout>(DOF_to_Elem_Matrix_Strides);

  //set stiffness Matrix Graph
  //debug print
    //std::cout << "DOF GRAPH MATRIX ENTRIES ON TASK " << myrank << std::endl;
  FOR_ALL_CLASS(idof, 0, num_dim*nlocal_nodes, {
    for (int istride = 0; istride < Gradient_Matrix_Strides(idof); istride++){
      DOF_Graph_Matrix(idof,istride) = Graph_Matrix(idof/num_dim,istride/num_dim)*num_dim + istride%num_dim;
    }
  }); // end parallel for

  
  //build inverse map for element gradient assembly
  for (size_t elem_gid = 0; elem_gid < rnum_elem; elem_gid++){
      FOR_ALL_CLASS(node_lid, 0, num_nodes_in_elem, {
          
          // get the global_id of the node
          size_t node_gid = nodes_in_elem(elem_gid, node_lid);
          
          // the column index is the num corners saved
          size_t j = count_saved_corners_in_node(node_gid);
          
          Element_Gradient_Matrix_Assembly_Map(elem_gid, node_lid) = j;

          // increment the number of corners saved to this node_gid
          count_saved_corners_in_node(node_gid)++;

      });  // end FOR_ALL over nodes in element
  } // end for elem_gid
  Kokkos::fence();  

  /*
  //construct distributed gradient matrix from local kokkos data
  //build column map for the global gradient matrix
  Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > colmap;
  const Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > dommap = local_dof_map;
  
  Tpetra::Details::makeColMap<LO,GO,node_type>(colmap,dommap,DOF_Graph_Matrix.get_kokkos_view(), nullptr);

  size_t nnz = DOF_Graph_Matrix.size();

  //debug print
  //std::cout << "DOF GRAPH SIZE ON RANK " << myrank << " IS " << nnz << std::endl;
  
  //local indices in the graph using the constructed column map
  CArrayKokkos<LO, array_layout, device_type, memory_traits> gradient_local_indices(nnz, "gradient_local_indices");
  
  //row offsets with compatible template arguments
  Kokkos::View<size_t *,array_layout, device_type, memory_traits> row_offsets = DOF_Graph_Matrix.start_index_;
  row_pointers row_offsets_pass("row_offsets", nlocal_nodes*num_dim+1);
  for(int ipass = 0; ipass < nlocal_nodes*num_dim + 1; ipass++){
    row_offsets_pass(ipass) = row_offsets(ipass);
  }

  size_t entrycount = 0;
  for(int irow = 0; irow < nlocal_nodes*num_dim; irow++){
    for(int istride = 0; istride < Gradient_Matrix_Strides(irow); istride++){
      gradient_local_indices(entrycount) = colmap->getLocalElement(DOF_Graph_Matrix(irow,istride));
      entrycount++;
    }
  }
  
  
  //sort values and indices
  Tpetra::Import_Util::sortCrsEntries<row_pointers, indices_array, values_array>(row_offsets_pass, gradient_local_indices.get_kokkos_view(), Force_Gradient_Positions.get_kokkos_view());
  Tpetra::Import_Util::sortCrsEntries<row_pointers, indices_array, values_array>(row_offsets_pass, gradient_local_indices.get_kokkos_view(), Force_Gradient_Velocities.get_kokkos_view());
  
  //Teuchos::RCP<Teuchos::ParameterList> crs_matrix_params = Teuchos::rcp(new Teuchos::ParameterList("crsmatrix"));
  //crs_matrix_params->set("sorted", false);
  distributed_force_gradient_positions = Teuchos::rcp(new MAT(local_dof_map, colmap, row_offsets_pass, gradient_local_indices.get_kokkos_view(), Force_Gradient_Positions.get_kokkos_view()));
  distributed_force_gradient_positions->fillComplete();
  distributed_force_gradient_velocities = Teuchos::rcp(new MAT(local_dof_map, colmap, row_offsets_pass, gradient_local_indices.get_kokkos_view(), Force_Gradient_Velocities.get_kokkos_view()));
  distributed_force_gradient_velocities->fillComplete();
  */
  //distributed_force_gradient_positions->describe(*fos,Teuchos::VERB_EXTREME);
  //distributed_force_gradient_velocities->describe(*fos,Teuchos::VERB_EXTREME);
}

/* ----------------------------------------------------------------------
   Enforce boundary conditions on the adjoint vectors
------------------------------------------------------------------------- */

void FEA_Module_SGH::boundary_adjoint(const mesh_t &mesh,
                       const DCArrayKokkos <boundary_t> &boundary,
                       vec_array &node_adjoint,
                       vec_array &node_phi_adjoint,
                       vec_array &node_psi_adjoint){

    //error and debug flag
    //DCArrayKokkos<bool> print_flag(1, "print_flag");
    //print_flag.host(0) = false;
    //print_flag.update_device();
   
    const size_t rk_level = simparam.dynamic_options.rk_num_bins - 1; 
    int num_dims = simparam.num_dims;
    // Loop over boundary sets
    for (size_t bdy_set=0; bdy_set<num_bdy_sets; bdy_set++){
        
        // Loop over boundary nodes in a boundary set
        FOR_ALL_CLASS(bdy_node_lid, 0, num_bdy_nodes_in_set.host(bdy_set), {
                
            // reflected (boundary array is on the device)
            if (boundary(bdy_set).type == BOUNDARY_CONDITION_TYPE::reflected){
            
                // directions with hydro_bc:
                // x_plane  = 0,
                // y_plane  = 1,
                // z_plane  = 2,
                size_t direction = boundary(bdy_set).surface.planar_surface_index();
                
                size_t bdy_node_gid = bdy_nodes_in_set(bdy_set, bdy_node_lid);
                    
                // Set velocity to zero in that directdion
                if(bdy_node_gid<nlocal_nodes)
                  node_adjoint(bdy_node_gid, direction) = 0.0;
                //node_phi_adjoint(bdy_node_gid, direction) = 0.0;
                        
            }
            else if (boundary(bdy_set).type == BOUNDARY_CONDITION_TYPE::fixed_position){
                
                size_t bdy_node_gid = bdy_nodes_in_set(bdy_set, bdy_node_lid);
                
                //debug clause
                //if(bdy_node_gid==549412) print_flag(0) = true;

                for(size_t dim=0; dim < num_dims; dim++){
                    // Set velocity to zero
                  if(bdy_node_gid<nlocal_nodes)
                    node_adjoint(bdy_node_gid, dim) = 0.0;
                  //node_phi_adjoint(bdy_node_gid, dim) = 0.0;
                }
                
            }
            
            
                
        }); // end for bdy_node_lid
	    
    } // end for bdy_set
    
    //debug check
    //print_flag.update_host();
    //if(print_flag.host(0)) std::cout << "found boundary node with id 549412" << std::endl;
    return;
} // end boundary_velocity function

/* ----------------------------------------------------------------------
   Communicate updated nodal adjoint vectors to ghost nodes
------------------------------------------------------------------------- */

void FEA_Module_SGH::comm_adjoint_vectors(int cycle){
  
  //debug print of design vector
      //std::ostream &out = std::cout;
      //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
      //if(myrank==0)
      //*fos << "Density data :" << std::endl;
      //node_densities_distributed->describe(*fos,Teuchos::VERB_EXTREME);
      //*fos << std::endl;
      //std::fflush(stdout);

  //communicate design densities
  //create import object using local node indices map and all indices map
  //Tpetra::Import<LO, GO> importer(map, all_node_map);
  
  //comms to get ghosts
  (*adjoint_vector_data)[cycle]->doImport(*adjoint_vector_distributed, *importer, Tpetra::INSERT);
  (*phi_adjoint_vector_data)[cycle]->doImport(*phi_adjoint_vector_distributed, *importer, Tpetra::INSERT);
  //all_node_map->describe(*fos,Teuchos::VERB_EXTREME);
  //all_node_velocities_distributed->describe(*fos,Teuchos::VERB_EXTREME);
  
  //update_count++;
  //if(update_count==1){
      //MPI_Barrier(world);
      //MPI_Abort(world,4);
  //}
}
