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
#include "mesh_writer.hpp"
#include "mesh_reader.hpp"   // PointIndexFromIJK

#include "region.hpp"
#include "string_utils.hpp"

#include <map>
#include <memory>
#include <cstring>
#include <sys/stat.h>
#include <iostream>
#include <regex>
#include <fstream>
#include <sstream>
#include <mpi.h>

namespace mesh_io_mpi_detail {

/// @brief Returns (0,1) if MPI is not initialized (e.g. some unit paths).
static void query_world_rank_size(int& rank, int& world_size) noexcept
{
    rank        = 0;
    world_size  = 1;
    int init = 0;
    if (MPI_Initialized(&init) == MPI_SUCCESS && init) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    }
}

} // namespace mesh_io_mpi_detail


/////////////////////////////////////////////////////////////////////////////
///
/// \fn writes mesh with the format given in the input.yaml file
///
/// \param Simulation mesh
/// \param Element related state
/// \param Node related state
/// \param Corner related state
/// \param Simulation input parameters
///
/////////////////////////////////////////////////////////////////////////////
void MeshWriter::write_mesh(swage::Mesh_t& mesh,
    State_t& State,
    SimulationParameters_t& SimulationParamaters,
    double dt,
    double time_value,
    CArray<double> graphics_times,
    std::vector<node_state> node_states,
    std::vector<gauss_pt_state> gauss_pt_states,
    std::vector<material_pt_state> material_pt_states,
    const size_t solver_id)
{


    // node_state is an enum for possible fields (e.g., coords, velocity, etc.), see state.h
    // gauss_pt_state is an enum for possible fields (e.g., vol, divergence, etc.)
    // material_pt_state is an enum for possible fields (e.g., den, pres, etc.)


    // *******************
    //  Update host 
    // *******************

    const size_t num_mats = State.MaterialPoints.num_material_points.size();

    // material point values
        
    //  Update host data for mat_pt state
    for (auto field : material_pt_states){
        switch(field){
            // scalar vars to write out
            case material_pt_state::density:
                State.MaterialPoints.den.update_host();
                break;
            case material_pt_state::pressure:
                State.MaterialPoints.pres.update_host();
                break;
            case material_pt_state::specific_internal_energy:
                State.MaterialPoints.sie.update_host();
                break;
            case material_pt_state::sound_speed:
                State.MaterialPoints.sspd.update_host();
                break;
            case material_pt_state::mass:
                State.MaterialPoints.mass.update_host();
                break;
            case material_pt_state::volume_fraction:
                State.MaterialPoints.mat_volfrac.update_host();
                State.MaterialPoints.geo_volfrac.update_host();
                break;
            case material_pt_state::eroded_flag:
                State.MaterialPoints.eroded.update_host();
                break;
            // tensor vars to write out
            case material_pt_state::stress:
                State.MaterialPoints.stress.update_host();
                break;
            case material_pt_state::strain:
                State.MaterialPoints.strain.update_host();
                break;
            
            // additional vars for thermal-mechanical solver
            case material_pt_state::thermal_conductivity:
                State.MaterialPoints.conductivity.update_host();
                break;
            
            case material_pt_state::specific_heat:
                State.MaterialPoints.specific_heat.update_host();
                break;

            case material_pt_state::heat_flux:
                State.MaterialPoints.q_flux.update_host();
                break;

            // add other variables here
            
            // not used
            case material_pt_state::elastic_modulii:
                break;
            case material_pt_state::shear_modulii:
                break;
            case material_pt_state::poisson_ratios:
                break;
            
            default:
                std::cout<<"Desired material point state not understood in outputs"<<std::endl;
        } // end switch
    } // end for over mat_pt_states



    // update gauss point values
    for (auto field : gauss_pt_states){
        switch(field){
            // scalar vars to write out
            case gauss_pt_state::volume:
                State.GaussPoints.vol.update_host();
                break;
            case gauss_pt_state::divergence_velocity:
                State.GaussPoints.div.update_host();
                break;
            case gauss_pt_state::level_set:
                State.GaussPoints.level_set.update_host();
                break;      
            // tensor vars to write out
            case gauss_pt_state::gradient_velocity:
                State.GaussPoints.vel_grad.update_host();
                break;
            case gauss_pt_state::shock_detector:
                State.GaussPoints.shock_detector.update_host();
                break;
            default:
                std::cout<<"Desired Gauss point state not understood in vtk outputs"<<std::endl;

        } // end switch
    } // end loop

    // nodal values
    for (auto field : node_states){
        switch(field){
            case node_state::mass:
                State.node.mass.update_host();
                break;
            case node_state::temp:
                State.node.temp.update_host();
                break;
            case node_state::coords:
                State.node.coords.update_host();
                break;
            case node_state::velocity:
                State.node.vel.update_host();
                break;
            case node_state::gradient_level_set:
                State.node.gradient_level_set.update_host();
                break;
            case node_state::displacement:
                State.node.displacement.update_host();
                break;
            case node_state::coords_t0:
                // blank because never changes and update host called in the driver
                break;

            case node_state::force:
                break;

            // heat transer vars
            case node_state::heat_transfer:
                break;

        } // end switch
    } // end for over 
    Kokkos::fence();


    // ******************************************
    //  Build Material and Element state outputs
    // ******************************************

    size_t num_mat_pt_scalar_vars = 0;
    size_t num_mat_pt_tensor_vars = 0;
    size_t num_mat_pt_vector_vars = 0;
        
    // count the number of material point state vars to write out
    for (auto field : SimulationParamaters.OutputOptions.output_mat_pt_state){
        switch(field){
            // scalar vars to write out
            case material_pt_state::density:
                num_mat_pt_scalar_vars ++;
                break;
            case material_pt_state::pressure:
                num_mat_pt_scalar_vars ++;
                break;
            case material_pt_state::specific_internal_energy:
                num_mat_pt_scalar_vars ++;
                break;
            case material_pt_state::sound_speed:
                num_mat_pt_scalar_vars ++;
                break;
            case material_pt_state::mass:
                num_mat_pt_scalar_vars ++;
                break;
            case material_pt_state::volume_fraction:
                num_mat_pt_scalar_vars ++; // mat volfrac
                num_mat_pt_scalar_vars ++; // geometric volfrac
                break;
            case material_pt_state::eroded_flag:
                num_mat_pt_scalar_vars ++;
                break;
            // tensor vars to write out
            case material_pt_state::stress:
                num_mat_pt_tensor_vars ++;
                break;
            case material_pt_state::strain:
                num_mat_pt_tensor_vars ++;
                break;
            
            // additional vars for thermal-mechanical solver
            case material_pt_state::thermal_conductivity:
                num_mat_pt_scalar_vars ++;
                break;
            
            case material_pt_state::specific_heat:
                num_mat_pt_scalar_vars ++;
                break;

            case material_pt_state::heat_flux:
                num_mat_pt_vector_vars ++;
                break;

            // add other variables here

            // not used
            case material_pt_state::elastic_modulii:
                break;
            case material_pt_state::shear_modulii:
                break;
            case material_pt_state::poisson_ratios:
                break;
            
            default:
                std::cout<<"Desired material point state not understood in outputs"<<std::endl;
        } // end switch
    } // end for over mat_pt_states



    size_t num_elem_scalar_vars = 0;
    size_t num_elem_vector_vars = 0;
    size_t num_elem_tensor_vars = 0;

    // count the number of element average fields to write out
    for (auto field : SimulationParamaters.OutputOptions.output_elem_state){
        switch(field){
            // scalar vars to write out
            case material_pt_state::density:
                num_elem_scalar_vars ++;
                break;
            case material_pt_state::pressure:
                num_elem_scalar_vars ++;
                break;
            case material_pt_state::specific_internal_energy:
                num_elem_scalar_vars ++;
                break;
            case material_pt_state::sound_speed:
                num_elem_scalar_vars ++;
                break;
            case material_pt_state::mass:
                num_elem_scalar_vars ++;
                break;
            // tensor vars to write out
            case material_pt_state::stress:
                num_elem_tensor_vars ++;
                break;

            // additional vars for thermal-mechanical solver
            case material_pt_state::thermal_conductivity:
                num_elem_scalar_vars ++;
                break;
            
            case material_pt_state::specific_heat:
                num_elem_scalar_vars ++;
                break;

            // add other variables here

            // not used
            case material_pt_state::volume_fraction:
                break;
            case material_pt_state::eroded_flag:
                break;
            case material_pt_state::elastic_modulii:
                break;
            case material_pt_state::shear_modulii:
                break;
            case material_pt_state::poisson_ratios:
                break;
            case material_pt_state::heat_flux:
                break;
            default:
                std::cout<<"Desired material point state not understood in outputs"<<std::endl;
        } // end switch
    } // end for over mat_pt_states


    size_t num_gauss_pt_scalar_vars = 0;
    size_t num_gauss_pt_tensor_vars = 0;

    // gauss point values to ouptput
    for (auto field : SimulationParamaters.OutputOptions.output_gauss_pt_state){
        switch(field){
            // scalar vars to write out
            case gauss_pt_state::volume:
                num_gauss_pt_scalar_vars ++;
                break;
            case gauss_pt_state::level_set:
                num_gauss_pt_scalar_vars ++;
                break;
            case gauss_pt_state::divergence_velocity:
                num_gauss_pt_scalar_vars ++;
                break;
            case gauss_pt_state::shock_detector:
                num_gauss_pt_scalar_vars ++;
                break;
            // tensor vars to write out
            case gauss_pt_state::gradient_velocity:
                num_gauss_pt_tensor_vars ++;
                break;
            
            default:
                std::cout<<"Desired Gauss point state not understood in vtk outputs"<<std::endl;

        } // end switch
    } // end loop

    // add the Gauss point state to the element state
    num_elem_scalar_vars += num_gauss_pt_scalar_vars;
    num_elem_tensor_vars += num_gauss_pt_tensor_vars;


    // Scalar, vector, and tensor value names associated with a elem
    std::vector<std::string> elem_scalar_var_names(num_elem_scalar_vars);
    std::vector<std::string> elem_tensor_var_names(num_elem_tensor_vars);
    std::vector<std::string> elem_vector_var_names(num_elem_vector_vars);

    // Scalar, vector, and tensor values associated with a material in part elems
    std::vector<std::string> mat_elem_scalar_var_names(num_mat_pt_scalar_vars);
    std::vector<std::string> mat_elem_tensor_var_names(num_mat_pt_tensor_vars);
    std::vector<std::string> mat_elem_vector_var_names(num_mat_pt_vector_vars);


    // the ids to access a variable in the mat_scalar_var_name or tensor list
    int mat_den_id = -1;
    int mat_pres_id = -1;
    int mat_sie_id = -1;
    int mat_sspd_id = -1;
    int mat_mass_id = -1;
    int mat_mat_volfrac_id = -1;  
    int mat_geo_volfrac_id = -1;  // geometric volume fraction of part
    int mat_eroded_id = -1;
    int mat_stress_id = -1;
    int mat_strain_id = -1;
    int mat_heat_flux_id = -1;

    int mat_conductivity_id = -1;
    int mat_specific_heat_id = -1;

    // the index for the scalar, vector, and tensor fields
    size_t var = 0;
    size_t vector_var = 0;
    size_t tensor_var = 0;

    // material point state to output
    for (auto field : SimulationParamaters.OutputOptions.output_mat_pt_state){
        switch(field){
            // scalar vars
            case material_pt_state::density:
                mat_elem_scalar_var_names[var] = "mat_den";
                mat_den_id = var;
                var++;
                break;
            case material_pt_state::pressure:
                mat_elem_scalar_var_names[var] = "mat_pres";
                mat_pres_id = var;
                var++;
                break;
            case material_pt_state::specific_internal_energy:
                mat_elem_scalar_var_names[var] = "mat_sie";
                mat_sie_id = var;
                var++;
                break;
            case material_pt_state::sound_speed:
                mat_elem_scalar_var_names[var] = "mat_sspd";
                mat_sspd_id = var;
                var++;
                break;
            case material_pt_state::mass:
                mat_elem_scalar_var_names[var] = "mat_mass";
                mat_mass_id = var;
                var++;
                break;
            case material_pt_state::volume_fraction:
                mat_elem_scalar_var_names[var] = "mat_volfrac";  // mat_mat_volfrac
                mat_mat_volfrac_id = var; 
                var++;

                mat_elem_scalar_var_names[var] = "mat_geo_volfrac";
                mat_geo_volfrac_id = var; 
                var++;
                break;
            case material_pt_state::eroded_flag:
                mat_elem_scalar_var_names[var] = "mat_eroded";
                mat_eroded_id = var;
                var++;
                break;
            // tensor vars
            case material_pt_state::stress:
                mat_elem_tensor_var_names[tensor_var] = "mat_stress";
                mat_stress_id = tensor_var;
                tensor_var++;
                break;
            case material_pt_state::strain:
                mat_elem_tensor_var_names[tensor_var] = "mat_strain";
                mat_strain_id = tensor_var;
                tensor_var++;
                break;


            // additional vars for thermal-mechanical solver
            case material_pt_state::thermal_conductivity:
                mat_elem_scalar_var_names[var] = "mat_thermal_K";
                mat_conductivity_id = var;
                var++;
                break;
            
            case material_pt_state::specific_heat:
                mat_elem_scalar_var_names[var] = "mat_Cp";
                mat_specific_heat_id = var;
                var++;
                break;

            case material_pt_state::heat_flux:
                mat_elem_vector_var_names[var] = "mat_heat_flux";
                mat_heat_flux_id = vector_var;
                vector_var++;
                break;


            // add other variables here

            // not used
            case material_pt_state::elastic_modulii:
                break;
            case material_pt_state::shear_modulii:
                break;
            case material_pt_state::poisson_ratios:
                break;
            
        } // end switch
    } // end for over mat_pt_states


    // element average fields to output

    // the ids to access a variable in the elem_scalar_var_name or tensor list
    int den_id = -1;
    int pres_id = -1;
    int sie_id = -1;
    int sspd_id = -1;
    int mass_id = -1; 
    int stress_id = -1;
    int strain_id = -1;

    int conductivity_id = -1;
    int specific_heat_id = -1;

    // reset the counters
    var = 0;
    vector_var = 0;
    tensor_var = 0;

    // element state to output
    for (auto field : SimulationParamaters.OutputOptions.output_elem_state){
        switch(field){
            // scalar vars
            case material_pt_state::density:
                elem_scalar_var_names[var] = "den";
                den_id = var;
                var++;
                break;
            case material_pt_state::pressure:
                elem_scalar_var_names[var] = "pres";
                pres_id = var;
                var++;
                break;
            case material_pt_state::specific_internal_energy:
                elem_scalar_var_names[var] = "sie";
                sie_id = var;
                var++;
                break;
            case material_pt_state::sound_speed:
                elem_scalar_var_names[var] = "sspd";
                sspd_id = var;
                var++;
                break;
            case material_pt_state::mass:
                elem_scalar_var_names[var] = "mass";
                mass_id = var;
                var++;
                break;
            // tensor vars
            case material_pt_state::stress:
                elem_tensor_var_names[tensor_var] = "stress";
                stress_id = tensor_var;
                tensor_var++;
                break;
            case material_pt_state::strain:
                elem_tensor_var_names[tensor_var] = "strain";
                strain_id = tensor_var;
                tensor_var++;
                break;

            // heat transfer variables
            case material_pt_state::thermal_conductivity:
                elem_scalar_var_names[var] = "thermal_K";
                conductivity_id = var;
                var++;
                break;
            
            case material_pt_state::specific_heat:
                elem_scalar_var_names[var] = "Cp";
                specific_heat_id = var;
                var++;
                break;

            // add other variables here

            // not used
            case material_pt_state::volume_fraction:
                break;
            case material_pt_state::eroded_flag:
                break;
            case material_pt_state::elastic_modulii:
                break;
            case material_pt_state::shear_modulii:
                break;
            case material_pt_state::poisson_ratios:
                break;
            case material_pt_state::heat_flux:
                break;
        } // end switch
    } // end for over mat_pt_states

    // append Gauss point vars to the element arrays
    int vol_id = -1;
    int div_id = -1;
    int level_set_id = -1;
    int vel_grad_id = -1;
    int shock_detector_id = -1;

    for (auto field : SimulationParamaters.OutputOptions.output_gauss_pt_state){
        switch(field){
            // scalars
            case gauss_pt_state::volume:
                elem_scalar_var_names[var] = "vol";
                vol_id = var;
                var++;
                break;
            case gauss_pt_state::divergence_velocity:
                elem_scalar_var_names[var] = "div";
                div_id = var;
                var++;
                break;

            case gauss_pt_state::level_set:
                elem_scalar_var_names[var] = "level_set";
                level_set_id = var;
                var++;
                break;
            
            case gauss_pt_state::shock_detector:
                elem_scalar_var_names[var] = "shock_detector";
                shock_detector_id = var;
                var++;
                break;

            // tensors
            case gauss_pt_state::gradient_velocity:
                elem_tensor_var_names[tensor_var] = "vel_grad";
                vel_grad_id = tensor_var;
                tensor_var++;
                break;
        } // end switch
    } // end loop over gauss_pt_states


    // *******************
    //  nodal values
    // *******************

    size_t num_node_scalar_vars = 0;
    size_t num_node_vector_vars = 0;

    for (auto field : SimulationParamaters.OutputOptions.output_node_state){
        switch(field){
            // --- scalars
            case node_state::mass:
                num_node_scalar_vars ++;
                break;
            case node_state::temp:
                num_node_scalar_vars ++;
                break;
            // -- vectors
            case node_state::coords:
                num_node_vector_vars ++;
                break;
            case node_state::velocity:
                num_node_vector_vars ++; // for velocity
                num_node_vector_vars ++; // for acceleration
                break;
            case node_state::gradient_level_set:
                num_node_vector_vars ++;
                break;
            case node_state::displacement:
                num_node_vector_vars ++;
                break;
            case node_state::coords_t0:
                // not necessary to output as coords at t=0 will be coords_t0
                break;                
            case node_state::force:
                break;
            
            // heat transer vars
            case node_state::heat_transfer:
                break;
        } // end switch
    } // end for over 
    Kokkos::fence();


    // Scalar and vector values associated with a node
    std::vector<std::string> node_scalar_var_names(num_node_scalar_vars);
    std::vector<std::string> node_vector_var_names(num_node_vector_vars);

    int node_mass_id = -1;
    int node_vel_id = -1;
    int node_accel_id = -1;
    int node_coord_id = -1;
    int node_temp_id = -1;
    int node_grad_level_set_id = -1;
    int node_disp_id = -1;

    // reset counters for node fields
    var = 0;
    vector_var = 0;
    tensor_var = 0;

    for (auto field : SimulationParamaters.OutputOptions.output_node_state){
        switch(field){
            // scalars
            case node_state::mass:
                node_scalar_var_names[var] = "node_mass";
                node_mass_id = var;
                var++;
                break;
            case node_state::temp:
                node_scalar_var_names[var] = "node_temp";
                node_temp_id = var;
                var++;
                break;

            // vector fields

            case node_state::coords:
                node_vector_var_names[vector_var] = "node_coords";
                node_coord_id = vector_var;
                vector_var++;
                break;

            case node_state::velocity:
                node_vector_var_names[vector_var] = "node_vel";
                node_vel_id = vector_var;
                vector_var++;

                node_vector_var_names[vector_var] = "node_accel";
                node_accel_id = vector_var;
                vector_var++;
                break;

            case node_state::gradient_level_set:
                node_vector_var_names[vector_var] = "node_grad_lvlset";
                node_grad_level_set_id = vector_var;
                vector_var++;
                break;
            case node_state::displacement:
                node_vector_var_names[vector_var] = "node_disp";
                node_disp_id = vector_var;
                vector_var ++;
                break;

            // -- not used vars
            case node_state::coords_t0:
                break;

            case node_state::force:
                break;

            // heat transer vars
            case node_state::heat_transfer:
                break;

            // tensors

        } // end switch
    } // end for over 


    // **************************************
    //  build and save element average fields
    // **************************************

    // short hand
    const size_t num_nodes = mesh.num_nodes;
    const size_t num_elems = mesh.num_elems;
    const size_t num_dims  = mesh.num_dims;
    const size_t num_nodes_in_elem = mesh.num_nodes_in_elem;
    const int Pn_order = mesh.Pn;

    // save the elem state to an array for exporting to graphics files
    DCArrayKokkos<double> elem_scalar_fields(num_elem_scalar_vars, num_elems, "elem_scalars");
    DCArrayKokkos<double> elem_tensor_fields(num_elem_tensor_vars, num_elems, 3, 3, "elem_tensors");
    DCArrayKokkos<double> elem_vector_fields(num_elem_vector_vars, num_elems, 3, "elem_vectors");
    elem_scalar_fields.set_values(0.0);
    elem_tensor_fields.set_values(0.0);


    // -----------------------------------------------------------------------
    // save the output fields to a single element average array for all state
    // -----------------------------------------------------------------------
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {

        // material point and guass point state are concatenated together
        concatenate_elem_fields(State.MaterialPoints,
                                State.GaussPoints,
                                elem_scalar_fields,
                                elem_tensor_fields,
                                State.MaterialToMeshMaps.elem_in_mat_elem,
                                SimulationParamaters.OutputOptions.output_elem_state,
                                SimulationParamaters.OutputOptions.output_gauss_pt_state,
                                State.MaterialToMeshMaps.num_mat_elems.host(mat_id),
                                mat_id,
                                num_elems,
                                den_id,
                                pres_id,
                                sie_id,
                                sspd_id,
                                mass_id,
                                stress_id,
                                strain_id,
                                vol_id,
                                div_id,
                                shock_detector_id,
                                level_set_id,
                                vel_grad_id,
                                conductivity_id,
                                specific_heat_id);
    } // end for mats

    // make specific fields for the element average
    if (sie_id>=0){
        FOR_ALL(elem_gid, 0, num_elems, {
            // get sie by dividing by the mass
            elem_scalar_fields(sie_id, elem_gid) /= (elem_scalar_fields(mass_id, elem_gid)+1.e-20); 
        });
    } // end if

    Kokkos::fence();
    elem_scalar_fields.update_host();
    elem_tensor_fields.update_host();
    

    // ************************
    //  Build the nodal fields 
    // ************************

    // save the nodal fields to an array for exporting to graphics files
    DCArrayKokkos<double> node_scalar_fields(num_node_scalar_vars, num_nodes, "node_scalars");
    DCArrayKokkos<double> node_vector_fields(num_node_vector_vars, num_nodes, 3, "node_tenors");
  
    concatenate_nodal_fields(State.node,
                             node_scalar_fields,
                             node_vector_fields,
                             SimulationParamaters.OutputOptions.output_node_state,
                             dt,
                             num_nodes,
                             num_dims,
                             node_mass_id,
                             node_vel_id,
                             node_accel_id,
                             node_coord_id,
                             node_grad_level_set_id,
                             node_temp_id,
                             node_disp_id);
                             

    Kokkos::fence();
    node_scalar_fields.update_host();
    node_vector_fields.update_host();


    // ********************************
    //  Write the nodal and elem fields 
    // ********************************

    if (SimulationParamaters.OutputOptions.format == output_options::viz ||
        SimulationParamaters.OutputOptions.format == output_options::viz_and_state) {

        int mpi_rank = 0;
        int mpi_size = 1;
        mesh_io_mpi_detail::query_world_rank_size(mpi_rank, mpi_size);

        MPI_Barrier(MPI_COMM_WORLD);

        if (mpi_rank == 0 && solver_id == 0 && graphics_id == 0) {
            struct stat st_rm;
            if (stat("vtk", &st_rm) == 0) {
                (void)system("rm -f vtk/Fierro*");
            }
            if (stat("vtk/data", &st_rm) == 0) {
                (void)system("rm -f vtk/data/Fierro*");
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        struct stat st;
        if (stat("vtk", &st) != 0) {
            if (system("mkdir vtk") == 1) {
                std::cout << "Unable to make vtk directory" << std::endl;
            }
        }
        if (stat("vtk/data", &st) != 0) {
            if (system("mkdir vtk/data") == 1) {
                std::cout << "Unable to make vtk/data directory" << std::endl;
            }
        }

        // Per-rank VTU: write owned elements only (rows 0..num_owned_elems-1 of nodes_in_elem), but
        // Points must list all local nodes (0..num_nodes-1). Connectivity is full local indexing;
        // boundary elements reference ghost nodes — truncating coords to num_owned_nodes (the prior
        // bug) misaligned indices and produced inverted/wrong cells in ParaView.
        const size_t n_owned_elems = mesh.num_owned_elems; // num_elems; 


        const std::string elem_fields_name = "fields";

        ViewCArray<double> node_coords_host(&State.node.coords.host(0, 0), num_nodes, num_dims);
        ViewCArray<size_t> nodes_in_elem_host(&mesh.nodes_in_elem.host(0, 0), n_owned_elems, num_nodes_in_elem);

        // VTK diagnostics (optional CellData / PointData in write_vtu): host-side scratch arrays,
        // then const pointers so write_vtu can emit one value per owned cell or per owned node.
        //   diag_mpi_rank_elem — CellData "mpi_rank" (which MPI rank owns each exported element).
        //   diag_global_elem   — CellData "global_elem_id" from mesh.local_to_global_elem_mapping (or local e).
        //   diag_global_node   — PointData "global_node_id" from mesh.local_to_global_node_mapping (or local n).
        //   p_*                — nullptr until filled; passed to write_vtu (nullptr disables that array).
        std::vector<double> diag_mpi_rank_elem;
        std::vector<double> diag_global_elem;
        std::vector<double> diag_global_node;
        const double*       p_rank_elem = nullptr;
        const double*       p_glob_elem = nullptr;
        const double*       p_glob_node = nullptr;

        if (n_owned_elems > 0 && num_nodes > 0) {
            diag_mpi_rank_elem.assign(n_owned_elems, static_cast<double>(mpi_rank));
            p_rank_elem = diag_mpi_rank_elem.data();

            diag_global_elem.resize(n_owned_elems);
            if (mpi_size > 1 || mesh.num_elems > mesh.num_owned_elems) {
                mesh.local_to_global_elem_mapping.update_host();
                for (size_t e = 0; e < n_owned_elems; e++) {
                    diag_global_elem[e] =
                        static_cast<double>(mesh.local_to_global_elem_mapping.host(e));
                }
            }
            else {
                for (size_t e = 0; e < n_owned_elems; e++) {
                    diag_global_elem[e] = static_cast<double>(e);
                }
            }
            p_glob_elem = diag_global_elem.data();

            diag_global_node.resize(num_nodes);
            if (mpi_size > 1 || mesh.num_nodes > mesh.num_owned_nodes) {
                mesh.local_to_global_node_mapping.update_host();
                for (size_t n = 0; n < num_nodes; n++) {
                    diag_global_node[n] =
                        static_cast<double>(mesh.local_to_global_node_mapping.host(n));
                }
            }
            else {
                for (size_t n = 0; n < num_nodes; n++) {
                    diag_global_node[n] = static_cast<double>(n);
                }
            }
            p_glob_node = diag_global_node.data();

            write_vtu(node_coords_host,
                      nodes_in_elem_host,
                      elem_scalar_fields,
                      elem_tensor_fields,
                      elem_vector_fields,
                      node_scalar_fields,
                      node_vector_fields,
                      elem_scalar_var_names,
                      elem_tensor_var_names,
                      elem_vector_var_names,
                      node_scalar_var_names,
                      node_vector_var_names,
                      elem_fields_name,
                      graphics_id,
                      num_nodes,
                      n_owned_elems,
                      num_nodes_in_elem,
                      Pn_order,
                      num_dims,
                      solver_id,
                      mpi_rank,
                      mpi_size,
                      p_rank_elem,
                      p_glob_elem,
                      p_glob_node);
        }

        // ********************************
        //  Build and write the mat fields 
        // ********************************

        if(num_mat_pt_scalar_vars > 0 || num_mat_pt_tensor_vars >0 || num_mat_pt_vector_vars >0){

            for (int mat_id = 0; mat_id < num_mats; mat_id++) {

                const size_t num_mat_elems = State.MaterialToMeshMaps.num_mat_elems.host(mat_id);

                // only save material data if the mat lives on the mesh, ie. has state allocated
                if (num_mat_elems>0){

                    // set the nodal vars to zero size, we don't write these fields again
                    node_scalar_var_names.clear();
                    node_vector_var_names.clear();

                    // the arrays storing all the material field data
                    DCArrayKokkos<double> mat_elem_scalar_fields(num_mat_pt_scalar_vars, num_mat_elems, "mat_pt_scalars");
                    DCArrayKokkos<double> mat_elem_tensor_fields(num_mat_pt_tensor_vars, num_mat_elems, 3, 3, "mat_pt_tensors");
                    DCArrayKokkos<double> mat_elem_vector_fields(num_mat_pt_vector_vars, num_mat_elems, 3, "mat_pt_vectors");


                    // concatenate material fields into a single array
                    concatenate_mat_fields(State.MaterialPoints,
                                           mat_elem_scalar_fields,
                                           mat_elem_tensor_fields,
                                           mat_elem_vector_fields,
                                           State.MaterialToMeshMaps.elem_in_mat_elem,
                                           SimulationParamaters.OutputOptions.output_mat_pt_state,
                                           num_mat_elems,
                                           mat_id,
                                           mat_den_id,
                                           mat_pres_id,
                                           mat_sie_id,
                                           mat_sspd_id,
                                           mat_mass_id,
                                           mat_mat_volfrac_id,
                                           mat_geo_volfrac_id,  
                                           mat_eroded_id,
                                           mat_stress_id,
                                           mat_strain_id,
                                           mat_conductivity_id,
                                           mat_specific_heat_id,
                                           mat_heat_flux_id);
                    Kokkos::fence();
                    mat_elem_scalar_fields.update_host();
                    mat_elem_tensor_fields.update_host();
                    mat_elem_vector_fields.update_host();


                    std::string str_mat_val = std::to_string(mat_id);                       
                    std::string mat_fields_name = "mat";
                    mat_fields_name += str_mat_val;  // add the mat number

                    // save the nodes belonging to this part (i.e., the material)
                    DCArrayKokkos <double> mat_node_coords(num_nodes,num_dims, "mat_node_coords");
                    DCArrayKokkos <size_t> mat_nodes_in_mat_elem(num_mat_elems, num_nodes_in_elem, "mat_nodes_in_mat_elem");

                    // the number of actual nodes belonging to the part (i.e., the material)
                    size_t num_mat_nodes = 0;

                    // build a unique mesh (element and nodes) for the material (i.e., the part)
                    build_material_elem_node_lists(mesh,
                                                   State.node.coords,
                                                   mat_node_coords,
                                                   mat_nodes_in_mat_elem,
                                                   State.MaterialToMeshMaps.elem_in_mat_elem,
                                                   mat_id,
                                                   num_mat_nodes,
                                                   num_mat_elems,
                                                   num_nodes_in_elem,
                                                   num_dims);

                    ViewCArray <double> mat_node_coords_host(&mat_node_coords.host(0,0), num_mat_nodes, num_dims);
                    ViewCArray <size_t> mat_nodes_in_elem_host(&mat_nodes_in_mat_elem.host(0,0), num_mat_elems, num_nodes_in_elem);
                    
                    // write out a vtu file this
                    write_vtu(mat_node_coords_host,
                              mat_nodes_in_elem_host,
                              mat_elem_scalar_fields,
                              mat_elem_tensor_fields,
                              mat_elem_vector_fields,
                              node_scalar_fields,
                              node_vector_fields,
                              mat_elem_scalar_var_names,
                              mat_elem_tensor_var_names,
                              mat_elem_vector_var_names,
                              node_scalar_var_names,
                              node_vector_var_names,
                              mat_fields_name,
                              graphics_id,
                              num_mat_nodes,
                              num_mat_elems,
                              num_nodes_in_elem,
                              Pn_order,
                              num_dims,
                              solver_id,
                              mpi_rank,
                              mpi_size,
                              nullptr,
                              nullptr,
                              nullptr);

                } // end for mat_id

            } // end if material is on the mesh

        } // end if mat variables are to be written


        // *************************************************
        //  write Paraview files to open the graphics files
        // *************************************************

        // save the graphics time
        graphics_times(graphics_id) = time_value;

        std::vector<unsigned long long> local_mat_counts(num_mats, 0ULL);
        for (size_t mi = 0; mi < num_mats; mi++) {
            local_mat_counts[mi] =
                static_cast<unsigned long long>(State.MaterialToMeshMaps.num_mat_elems.host(mi));
        }
        std::vector<unsigned long long> gathered_mat_elems;
        if (mpi_rank == 0) {
            gathered_mat_elems.assign(static_cast<size_t>(mpi_size) * num_mats, 0ULL);
        }
        MPI_Gather(local_mat_counts.data(),
                   static_cast<int>(num_mats),
                   MPI_UNSIGNED_LONG_LONG,
                   mpi_rank == 0 ? gathered_mat_elems.data() : nullptr,
                   static_cast<int>(num_mats),
                   MPI_UNSIGNED_LONG_LONG,
                   0,
                   MPI_COMM_WORLD);

        const unsigned long long      local_owned_elems_ull = static_cast<unsigned long long>(n_owned_elems);
        std::vector<unsigned long long> gathered_owned_elems;
        if (mpi_rank == 0) {
            gathered_owned_elems.assign(static_cast<size_t>(mpi_size), 0ULL);
        }
        MPI_Gather(&local_owned_elems_ull,
                   1,
                   MPI_UNSIGNED_LONG_LONG,
                   mpi_rank == 0 ? gathered_owned_elems.data() : nullptr,
                   1,
                   MPI_UNSIGNED_LONG_LONG,
                   0,
                   MPI_COMM_WORLD);

        // check to see if an mesh state was written 
        bool write_mesh_state = false;
        if( num_elem_scalar_vars > 0 ||
            num_elem_tensor_vars > 0 ||
            num_node_scalar_vars > 0 ||
            num_node_vector_vars > 0)
        {
            write_mesh_state = true;
        }

        // check to see if a mat state was written
        bool write_mat_pt_state = false;
        if( num_mat_pt_scalar_vars > 0 ||
            num_mat_pt_tensor_vars > 0 ||
            num_mat_pt_vector_vars > 0)
        {
             write_mat_pt_state = true;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (mpi_rank == 0) {
            const std::string mat_fields_name = "mat";
            write_vtm(graphics_times,
                      elem_fields_name,
                      mat_fields_name,
                      time_value,
                      graphics_id,
                      num_mats,
                      write_mesh_state,
                      write_mat_pt_state,
                      solver_id,
                      mpi_size,
                      gathered_owned_elems.empty() ? nullptr : gathered_owned_elems.data(),
                      gathered_mat_elems.empty() ? nullptr : gathered_mat_elems.data());

            write_pvd(graphics_times,
                      time_value,
                      graphics_id,
                      solver_id,
                      mpi_rank);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // increment graphics id counter
        graphics_id++; // this is private variable in the class

    } // end if viz paraview output is to be written


    // STATE
    if (SimulationParamaters.OutputOptions.format == output_options::state ||
        SimulationParamaters.OutputOptions.format == output_options::viz_and_state) {

        /* write_material_point_state(mesh,
                                   State,
                                   SimulationParamaters,
                                   time_value,
                                   graphics_times,
                                   node_states,
                                   gauss_pt_states,
                                   material_pt_states); */

        write_text_state(mesh,
                         State,
                         mat_elem_scalar_var_names,
                         mat_elem_tensor_var_names,
                         mat_elem_vector_var_names,
                         node_scalar_var_names,
                         node_vector_var_names,
                         mesh.num_nodes_in_elem,
                         mesh.num_dims,
                         time_value,
                         mat_den_id,
                         mat_pres_id,
                         mat_sie_id,
                         mat_sspd_id,
                         mat_mass_id,
                         mat_mat_volfrac_id,
                         mat_geo_volfrac_id,
                         mat_eroded_id,
                         mat_stress_id,
                         mat_strain_id,
                         mat_conductivity_id,
                         mat_specific_heat_id,
                         mat_heat_flux_id,
                         node_mass_id,
                         node_vel_id,
                         node_coord_id,
                         node_temp_id,
                         node_grad_level_set_id,
                         node_disp_id);

    } // end if state is to be written


    // will drop ensight outputs in the near future
    if (SimulationParamaters.OutputOptions.format == output_options::ensight){
       write_ensight(mesh,
                     State,
                     SimulationParamaters,
                     dt,
                     time_value,
                     graphics_times,
                     node_states,
                     gauss_pt_states,
                     material_pt_states);
    }

    return;

} // end write_mesh


/////////////////////////////////////////////////////////////////////////////
///
/// \fn write_ensight
///
/// \brief Writes an ensight output file
///
/// \param Simulation mesh
/// \param State data
/// \param Simulation parameters
/// \param current time value
/// \param Vector of all graphics output times
///
/////////////////////////////////////////////////////////////////////////////
void MeshWriter::write_ensight(swage::Mesh_t& mesh,
    State_t& State,
    SimulationParameters_t& SimulationParamaters,
    double dt,
    double time_value,
    CArray<double> graphics_times,
    std::vector<node_state> node_states,
    std::vector<gauss_pt_state> gauss_pt_states,
    std::vector<material_pt_state> material_pt_states)
{
    size_t num_mats = State.MaterialPoints.num_material_points.size();

    // ---- Update host data ----

    // material point values
    State.MaterialPoints.den.update_host();
    State.MaterialPoints.pres.update_host();
    State.MaterialPoints.stress.update_host();
    State.MaterialPoints.strain.update_host();
    State.MaterialPoints.sspd.update_host();
    State.MaterialPoints.sie.update_host();
    State.MaterialPoints.mass.update_host();
    State.MaterialPoints.eroded.update_host();


    // gauss point values
    State.GaussPoints.vol.update_host();

    // nodal values
    State.node.coords.update_host();
    State.node.vel.update_host();
    State.node.mass.update_host();

    Kokkos::fence();

    // --------------------------

    const int num_scalar_vars = 10;
    const int num_vec_vars    = 3;

    std::string name_tmp;
    name_tmp = "Outputs_SGH";

    char* name = new char [name_tmp.length() + 1];
    std::strcpy(name, name_tmp.c_str());

    const char scalar_var_names[num_scalar_vars][15] = {
        "den", "pres", "sie", "vol", "mass", "sspd", "speed", "mat_id", "elem_switch", "eroded"
    };

    const char vec_var_names[num_vec_vars][15] = {
        "pos", "vel", "accel"
    };

    // short hand
    const size_t num_nodes = mesh.num_nodes;
    const size_t num_elems = mesh.num_elems;
    const size_t num_dims  = mesh.num_dims;

    // save the cell state to an array for exporting to graphics files
    auto elem_fields = CArray<double>(num_elems, num_scalar_vars);
    int  elem_switch = 1;


    DCArrayKokkos<double> speed(num_elems, "speed");
    FOR_ALL(elem_gid, 0, num_elems, {
        double elem_vel[3]; // note:initialization with a list won't work
        elem_vel[0] = 0.0;
        elem_vel[1] = 0.0;
        elem_vel[2] = 0.0;
        // get the coordinates of the element center
        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
            elem_vel[0] += State.node.vel(mesh.nodes_in_elem(elem_gid, node_lid), 0);
            elem_vel[1] += State.node.vel(mesh.nodes_in_elem(elem_gid, node_lid), 1);
            if (mesh.num_dims == 3) {
                elem_vel[2] += State.node.vel(mesh.nodes_in_elem(elem_gid, node_lid), 2);
            }
            else{
                elem_vel[2] = 0.0;
            }
        } // end loop over nodes in element
        elem_vel[0] = elem_vel[0] / mesh.num_nodes_in_elem;
        elem_vel[1] = elem_vel[1] / mesh.num_nodes_in_elem;
        elem_vel[2] = elem_vel[2] / mesh.num_nodes_in_elem;

        double speed_sqrd = 0.0;
        for (int dim = 0; dim < num_dims; dim++) {
            speed_sqrd += elem_vel[dim] * elem_vel[dim];
        }
        speed(elem_gid) = sqrt(speed_sqrd);
    }); // end parallel for
    speed.update_host();
    Kokkos::fence();

    // save the output scale fields to a single 2D array

    // export material centeric data to the elements
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {
        size_t num_mat_elems = State.MaterialToMeshMaps.num_mat_elems.host(mat_id);

        for (size_t mat_elem_sid = 0; mat_elem_sid < num_mat_elems; mat_elem_sid++) {
            // 1 material per element

            // get elem gid
            size_t elem_gid = State.MaterialToMeshMaps.elem_in_mat_elem.host(mat_id, mat_elem_sid);

            // save outputs
            elem_fields(elem_gid, 0) = State.MaterialPoints.den.host(mat_id, mat_elem_sid);
            elem_fields(elem_gid, 1) = State.MaterialPoints.pres.host(mat_id, mat_elem_sid);
            elem_fields(elem_gid, 2) = State.MaterialPoints.sie.host(mat_id, mat_elem_sid);
            // 3 is guass point vol
            elem_fields(elem_gid, 4) = State.MaterialPoints.mass.host(mat_id, mat_elem_sid);
            elem_fields(elem_gid, 5) = State.MaterialPoints.sspd.host(mat_id, mat_elem_sid);
            // 6 is elem speed
            elem_fields(elem_gid, 7) = (double)mat_id;
            // 8 is the e_switch
            elem_fields(elem_gid, 9) = (double)State.MaterialPoints.eroded.host(mat_id, mat_elem_sid);
        } // end for mat elems storage
    } // end parallel loop over materials

    // export element centric data
    double e_switch = 1;
    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
        elem_fields(elem_gid, 3) = State.GaussPoints.vol.host(elem_gid);
        elem_fields(elem_gid, 6) = speed.host(elem_gid);
        elem_fields(elem_gid, 8) = e_switch;
        elem_switch *= -1;
    } // end for elem_gid

    // save the vertex vector fields to an array for exporting to graphics files
    CArray<double> vec_fields(num_nodes, num_vec_vars, 3);

    for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
        // position, var 0
        vec_fields(node_gid, 0, 0) = State.node.coords.host(node_gid, 0);
        vec_fields(node_gid, 0, 1) = State.node.coords.host(node_gid, 1);
        if (num_dims == 2) {
            vec_fields(node_gid, 0, 2) = 0.0;
        }
        else{
            vec_fields(node_gid, 0, 2) = State.node.coords.host(node_gid, 2);
        }

        // velocity, var 1
        vec_fields(node_gid, 1, 0) = State.node.vel.host(node_gid, 0);
        vec_fields(node_gid, 1, 1) = State.node.vel.host(node_gid, 1);
        if (num_dims == 2) {
            vec_fields(node_gid, 1, 2) = 0.0;
        }
        else{
            vec_fields(node_gid, 1, 2) = State.node.vel.host(node_gid, 2);
        }

        // accelleration, var 2
        vec_fields(node_gid, 2, 0) = (State.node.vel.host(node_gid, 0) - State.node.vel_n0.host(node_gid, 0))/dt;
        vec_fields(node_gid, 2, 1) = (State.node.vel.host(node_gid, 1) - State.node.vel_n0.host(node_gid, 1))/dt;
        if (num_dims == 2) {
            vec_fields(node_gid, 2, 2) = 0.0;
        }
        else{
            vec_fields(node_gid, 2, 2) = (State.node.vel.host(node_gid, 2) - State.node.vel_n0.host(node_gid, 2))/dt;
        }


    } // end for loop over vertices


    //  ---------------------------------------------------------------------------
    //  Setup of file and directoring for exporting
    //  ---------------------------------------------------------------------------
    FILE* out[20];   // the output files that are written to
    char  filename[128];
    int   max_len = sizeof filename;
    int   str_output_len;

    struct stat st;

    if (stat("ensight", &st) != 0) {
        system("mkdir ensight");
    }

    if (stat("ensight/data", &st) != 0) {
        system("mkdir ensight/data");
    }

    //  ---------------------------------------------------------------------------
    //  Write the Geometry file
    //  ---------------------------------------------------------------------------
    // sprintf(filename, "ensight/data/%s.%05d.geo", name, graphics_id);
    str_output_len = snprintf(filename, max_len, "ensight/data/%s.%05d.geo", name, graphics_id);
    // filename has the full string
    if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

    out[0] = fopen(filename, "w");

    fprintf(out[0], "A graphics dump by Fierro \n");

    fprintf(out[0], "%s", "EnSight Gold geometry\n");
    fprintf(out[0], "%s", "node id assign\n");
    fprintf(out[0], "%s", "element id assign\n");

    fprintf(out[0], "part\n");
    fprintf(out[0], "%10d\n", 1);
    fprintf(out[0], "Mesh\n");

    // --- vertices ---
    fprintf(out[0], "coordinates\n");
    fprintf(out[0], "%10lu\n", num_nodes);

    // write all components of the point coordinates
    for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
        fprintf(out[0], "%12.5e\n", State.node.coords.host(node_gid, 0));
    }

    for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
        fprintf(out[0], "%12.5e\n", State.node.coords.host(node_gid, 1));
    }

    for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
        if (num_dims == 3) {
            fprintf(out[0], "%12.5e\n", State.node.coords.host(node_gid, 2));
        }
        else{
            fprintf(out[0], "%12.5e\n", 0.0);
        }
    }

    // --- elements ---
    if (num_dims == 3) {
        fprintf(out[0], "hexa8\n");
    }
    else{
        fprintf(out[0], "quad4\n");
    }
    fprintf(out[0], "%10lu\n", num_elems);


    int convert_ijk_to_ensight[8];
    if(mesh.num_dims==3){
        convert_ijk_to_ensight[0] = 0;
        convert_ijk_to_ensight[1] = 1;
        convert_ijk_to_ensight[2] = 3;
        convert_ijk_to_ensight[3] = 2;
        convert_ijk_to_ensight[4] = 4;
        convert_ijk_to_ensight[5] = 5;
        convert_ijk_to_ensight[6] = 7;
        convert_ijk_to_ensight[7] = 6;
    }
    else{
    
        convert_ijk_to_ensight[0] = 0;
        convert_ijk_to_ensight[1] = 1;
        convert_ijk_to_ensight[2] = 2;
        convert_ijk_to_ensight[3] = 3;
        convert_ijk_to_ensight[4] = 4;
        convert_ijk_to_ensight[5] = 5;
        convert_ijk_to_ensight[6] = 6;
        convert_ijk_to_ensight[7] = 7;
    } // end if


    // write all global point numbers for this cell
    for (int elem_gid = 0; elem_gid < num_elems; elem_gid++) {
        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
            fprintf(out[0], "%10lu\t", mesh.nodes_in_elem.host(elem_gid, convert_ijk_to_ensight[node_lid]) + 1); // note: node_gid starts at 1
        }
        fprintf(out[0], "\n");
    }

    fclose(out[0]);

    // ---------------------------------------------------------------------------
    // Write the Scalar variable files
    // ---------------------------------------------------------------------------

    // ensight_vars = (den, pres,...)
    for (int var = 0; var < num_scalar_vars; var++) {
        // write a scalar value
        // sprintf(filename, "ensight/data/%s.%05d.%s", name, graphics_id, scalar_var_names[var]);
        str_output_len = snprintf(filename, max_len, "ensight/data/%s.%05d.%s", name, graphics_id, scalar_var_names[var]);
        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

        out[0] = fopen(filename, "w");

        fprintf(out[0], "Per_elem scalar values\n");
        fprintf(out[0], "part\n");
        fprintf(out[0], "%10d\n", 1);
        if (num_dims == 3) {
            fprintf(out[0], "hexa8\n");
        }
        else{
            fprintf(out[0], "quad4\n");
        }

        for (int elem_id = 0; elem_id < num_elems; elem_id++) {
            fprintf(out[0], "%12.5e\n", elem_fields(elem_id, var));
        }

        fclose(out[0]);
    } // end for var

    //  ---------------------------------------------------------------------------
    //  Write the Vector variable files
    //  ---------------------------------------------------------------------------

    // ensight vector vars = (position, velocity, force)
    for (int var = 0; var < num_vec_vars; var++) {
        // sprintf(filename, "ensight/data/%s.%05d.%s", name, graphics_id, vec_var_names[var]);
        str_output_len = snprintf(filename, max_len, "ensight/data/%s.%05d.%s", name, graphics_id, vec_var_names[var]);
        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

        out[0] = fopen(filename, "w");
        // fprintf(out[0],"Per_node vector values\n");
        // fprintf(out[0],"part\n");
        // fprintf(out[0],"%10d \n",1);
        // fprintf(out[0],"hexa8\n"); // WARNING, maybe bug here?

        fprintf(out[0], "Per_node vector values\n");
        fprintf(out[0], "part\n");
        fprintf(out[0], "%10d\n", 1);
        fprintf(out[0], "block\n");

        for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
            fprintf(out[0], "%12.5e\n", vec_fields(node_gid, var, 0));
        }

        for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
            fprintf(out[0], "%12.5e\n", vec_fields(node_gid, var, 1));
        }

        for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
            fprintf(out[0], "%12.5e\n", vec_fields(node_gid, var, 2));
        }

        fclose(out[0]);
    } // end for var

    // ---------------------------------------------------------------------------
    // Write the case file
    // ---------------------------------------------------------------------------

    // sprintf(filename, "ensight/%s.case", name);
    str_output_len = snprintf(filename, max_len, "ensight/%s.case", name);
    if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

    out[0] = fopen(filename, "w");

    fprintf(out[0], "FORMAT\n");
    fprintf(out[0], "type: ensight gold\n");
    fprintf(out[0], "GEOMETRY\n");

    // sprintf(filename, "model: data/%s.*****.geo\n", name);
    str_output_len = snprintf(filename, max_len, "model: data/%s.*****.geo\n", name);
    if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

    fprintf(out[0], "%s", filename);
    fprintf(out[0], "VARIABLE\n");

    for (int var = 0; var < num_scalar_vars; var++) {
        // sprintf(filename, "scalar per element: %s data/%s.*****.%s\n", scalar_var_names[var], name, scalar_var_names[var]);
        str_output_len = snprintf(filename, max_len, "scalar per element: %s data/%s.*****.%s\n", scalar_var_names[var], name, scalar_var_names[var]);
        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

        fprintf(out[0], "%s", filename);
    }

    for (int var = 0; var < num_vec_vars; var++) {
        // sprintf(filename, "vector per node: %s data/%s.*****.%s\n", vec_var_names[var], name, vec_var_names[var]);
        str_output_len = snprintf(filename, max_len, "vector per node: %s data/%s.*****.%s\n", vec_var_names[var], name, vec_var_names[var]);
        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }
        fprintf(out[0], "%s", filename);
    }

    fprintf(out[0], "TIME\n");
    fprintf(out[0], "time set: 1\n");
    fprintf(out[0], "number of steps: %4d\n", graphics_id + 1);
    fprintf(out[0], "filename start number: 0\n");
    fprintf(out[0], "filename increment: 1\n");
    fprintf(out[0], "time values: \n");

    graphics_times(graphics_id) = time_value;

    for (int i = 0; i <= graphics_id; i++) {
        fprintf(out[0], "%12.5e\n", graphics_times(i));
    }
    fclose(out[0]);

    // ---------------------------------------------------------------------------
    // Done writing the graphics dump
    // ---------------------------------------------------------------------------

    // increment graphics id counter
    graphics_id++;

    delete[] name;


    return;
}


/////////////////////////////////////////////////////////////////////////////
///
/// \fn write_vtk_old
///
/// \brief Writes a vtk output file
///
/// \param Simulation mesh
/// \param State data
/// \param Simulation parameters
/// \param current time value
/// \param Vector of all graphics output times
///
/////////////////////////////////////////////////////////////////////////////
void MeshWriter::write_vtk_old(swage::Mesh_t& mesh,
    State_t& State,
    SimulationParameters_t& SimulationParamaters,
    double dt,
    double time_value,
    CArray<double> graphics_times,
    std::vector<node_state> node_states,
    std::vector<gauss_pt_state> gauss_pt_states,
    std::vector<material_pt_state> material_pt_states)
{

    size_t num_mats = State.MaterialPoints.num_material_points.size();

    // ---- Update host data ----

    // material point values
    State.MaterialPoints.den.update_host();
    State.MaterialPoints.pres.update_host();
    State.MaterialPoints.stress.update_host();
    State.MaterialPoints.strain.update_host();
    State.MaterialPoints.sspd.update_host();
    State.MaterialPoints.sie.update_host();
    State.MaterialPoints.mass.update_host();
    State.MaterialPoints.conductivity.update_host();
    State.MaterialPoints.temp_grad.update_host();
    State.MaterialPoints.eroded.update_host();


    // gauss point values
    State.GaussPoints.vol.update_host();

    // nodal values
    State.node.coords.update_host();
    State.node.vel.update_host();
    State.node.mass.update_host();
    State.node.temp.update_host();

    Kokkos::fence();


    const int num_cell_scalar_vars = 13;
    const int num_cell_vec_vars    = 0;
    const int num_cell_tensor_vars = 0;

    const int num_point_scalar_vars = 1;
    const int num_point_vec_vars = 2;


    // Scalar values associated with a cell
    const char cell_scalar_var_names[num_cell_scalar_vars][15] = {
        "den", "pres", "sie", "vol", "mass", "sspd", "speed", "mat_id", "elem_switch","eroded", "temp_grad_x", "temp_grad_y", "temp_grad_z"
    };
    
    const char cell_vec_var_names[num_cell_vec_vars][15] = {
        
    };

    const char point_scalar_var_names[num_point_scalar_vars][15] = {
        "temp"
    };

    const char point_vec_var_names[num_point_vec_vars][15] = {
        "pos", "vel" 
    };

    // short hand
    const size_t num_nodes = mesh.num_nodes;
    const size_t num_elems = mesh.num_elems;
    const size_t num_dims  = mesh.num_dims;

    // save the cell state to an array for exporting to graphics files
    auto elem_fields = CArray<double>(num_elems, num_cell_scalar_vars);
    int  elem_switch = 1;

    DCArrayKokkos<double> speed(num_elems, "speed");
    FOR_ALL(elem_gid, 0, num_elems, {
        double elem_vel[3]; // note:initialization with a list won't work
        elem_vel[0] = 0.0;
        elem_vel[1] = 0.0;
        elem_vel[2] = 0.0;
        // get the coordinates of the element center
        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
            elem_vel[0] += State.node.vel(mesh.nodes_in_elem(elem_gid, node_lid), 0);
            elem_vel[1] += State.node.vel(mesh.nodes_in_elem(elem_gid, node_lid), 1);
            if (mesh.num_dims == 3) {
                elem_vel[2] += State.node.vel(mesh.nodes_in_elem(elem_gid, node_lid), 2);
            }
            else{
                elem_vel[2] = 0.0;
            }
        } // end loop over nodes in element
        elem_vel[0] = elem_vel[0] / mesh.num_nodes_in_elem;
        elem_vel[1] = elem_vel[1] / mesh.num_nodes_in_elem;
        elem_vel[2] = elem_vel[2] / mesh.num_nodes_in_elem;

        double speed_sqrd = 0.0;
        for (int dim = 0; dim < num_dims; dim++) {
            speed_sqrd += elem_vel[dim] * elem_vel[dim];
        }
        speed(elem_gid) = sqrt(speed_sqrd);
    }); // end parallel for
    speed.update_host();
    Kokkos::fence();

    // save the output scale fields to a single 2D array


    // export material centeric data to the elements
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {
        size_t num_mat_elems = State.MaterialToMeshMaps.num_mat_elems.host(mat_id);

        for (size_t mat_elem_sid = 0; mat_elem_sid < num_mat_elems; mat_elem_sid++) {
            // 1 material per element

            // get elem gid
            size_t elem_gid = State.MaterialToMeshMaps.elem_in_mat_elem.host(mat_id, mat_elem_sid);

            // save outputs
            elem_fields(elem_gid, 0) = State.MaterialPoints.den.host(mat_id,mat_elem_sid);
            elem_fields(elem_gid, 1) = State.MaterialPoints.pres.host(mat_id, mat_elem_sid);
            elem_fields(elem_gid, 2) = State.MaterialPoints.sie.host(mat_id, mat_elem_sid);
            // 3 is guass point vol
            elem_fields(elem_gid, 4) = State.MaterialPoints.mass.host(mat_id, mat_elem_sid);
            elem_fields(elem_gid, 5) = State.MaterialPoints.sspd.host(mat_id, mat_elem_sid);
            // 6 is elem speed
            elem_fields(elem_gid, 7) = (double)mat_id;
            // 8 is the e_switch
            elem_fields(elem_gid, 9) = (double)State.MaterialPoints.eroded.host(mat_id, mat_elem_sid);
            elem_fields(elem_gid, 10) = (double)State.MaterialPoints.temp_grad.host(mat_id, elem_gid,0);
            elem_fields(elem_gid, 11) = (double)State.MaterialPoints.temp_grad.host(mat_id, elem_gid,1);
            elem_fields(elem_gid, 12) = (double)State.MaterialPoints.temp_grad.host(mat_id, elem_gid,2);
        } // end for mat elems storage
    } // end parallel loop over materials

    // export element centric data
    double e_switch = 1;
    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
        elem_fields(elem_gid, 3) = State.GaussPoints.vol.host(elem_gid);
        elem_fields(elem_gid, 6) = speed.host(elem_gid);
        elem_fields(elem_gid, 8) = State.GaussPoints.div.host(elem_gid);
        elem_switch *= -1;
    } // end for elem_gid

    // save the vertex vector fields to an array for exporting to graphics files
    CArray<double> vec_fields(num_nodes, num_point_vec_vars, 3);
    CArray<double> point_scalar_fields(num_nodes, num_point_scalar_vars);

    for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
        // position, var 0
        vec_fields(node_gid, 0, 0) = State.node.coords.host(node_gid, 0);
        vec_fields(node_gid, 0, 1) = State.node.coords.host(node_gid, 1);
        if (num_dims == 2) {
            vec_fields(node_gid, 0, 2) = 0.0;
        }
        else{
            vec_fields(node_gid, 0, 2) = State.node.coords.host(node_gid, 2);
        }

        // position, var 1
        vec_fields(node_gid, 1, 0) = State.node.vel.host(node_gid, 0);
        vec_fields(node_gid, 1, 1) = State.node.vel.host(node_gid, 1);
        if (num_dims == 2) {
            vec_fields(node_gid, 1, 2) = 0.0;
        }
        else{
            vec_fields(node_gid, 1, 2) = State.node.vel.host(node_gid, 2);
        }

        point_scalar_fields(node_gid, 0) = State.node.temp.host(node_gid);
    } // end for loop over vertices


    FILE* out[20];   // the output files that are written to
    char  filename[100]; // char string
    int   max_len = sizeof filename;
    int   str_output_len;

    struct stat st;

    if (stat("vtk", &st) != 0) {
        system("mkdir vtk");
    }

    // snprintf(filename, max_len, "ensight/data/%s.%05d.%s", name, graphics_id, vec_var_names[var]);

    //sprintf(filename, "vtk/Fierro.%05d.vtk", graphics_id);  // mesh file
    str_output_len = snprintf(filename, max_len, "vtk/Fierro.%05d.vtk", graphics_id);
    if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }
     // mesh file
    
    out[0] = fopen(filename, "w");

    fprintf(out[0], "# vtk DataFile Version 2.0\n");  // part 2
    fprintf(out[0], "Mesh for Fierro\n");             // part 2
    fprintf(out[0], "ASCII \n");                      // part 3
    fprintf(out[0], "DATASET UNSTRUCTURED_GRID\n\n"); // part 4

    fprintf(out[0], "POINTS %zu float\n", mesh.num_nodes);

    // write all components of the point coordinates
    for (size_t node_gid = 0; node_gid < mesh.num_nodes; node_gid++) {
        fprintf(out[0],
                "%f %f %f\n",
                State.node.coords.host(node_gid, 0),
                State.node.coords.host(node_gid, 1),
                State.node.coords.host(node_gid, 2));
    } // end for

    /*
    ---------------------------------------------------------------------------
    Write the elems
    ---------------------------------------------------------------------------
    */

    fprintf(out[0], "\n");
    fprintf(out[0], "CELLS %lu %lu\n", mesh.num_elems, mesh.num_elems + mesh.num_elems * mesh.num_nodes_in_elem);  // size=all printed values

    int Pn_order   = mesh.Pn;
    int order[3]   = { Pn_order, Pn_order, Pn_order };

    // const int num_1D_points = Pn_order+1;

    // write all global point numbers for this elem
    for (size_t elem_gid = 0; elem_gid < mesh.num_elems; elem_gid++) {
        fprintf(out[0], "%lu ", mesh.num_nodes_in_elem); // num points in this elem

        for (int k = 0; k <= Pn_order; k++) {
            for (int j = 0; j <= Pn_order; j++) {
                for (int i = 0; i <= Pn_order; i++) {
                    size_t node_lid = PointIndexFromIJK(i, j, k, order);
                    fprintf(out[0], "%lu ", mesh.nodes_in_elem.host(elem_gid, node_lid));
                }
            }
        }

        fprintf(out[0], "\n");
    } // end for

    // Write the element types
    fprintf(out[0], "\n");
    fprintf(out[0], "CELL_TYPES %zu \n", mesh.num_elems);
    // VTK_LAGRANGE_HEXAHEDRON: 72,
    // VTK_HIGHER_ORDER_HEXAHEDRON: 67
    // VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33
    // element types: https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
    // element types: https://kitware.github.io/vtk-js/api/Common_DataModel_CellTypes.html
    // vtk format: https://www.kitware.com//modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/
    for (size_t elem_gid = 0; elem_gid < mesh.num_elems; elem_gid++) {
        fprintf(out[0], "%d \n", 72);
    }

    /*
    ---------------------------------------------------------------------------
    Write the nodal vector variables to file
    ---------------------------------------------------------------------------
    */

    fprintf(out[0], "\n");
    fprintf(out[0], "POINT_DATA %zu \n", mesh.num_nodes);

    // vtk vector vars = (position, velocity)
    for (int var = 0; var < num_point_vec_vars; var++) {
        fprintf(out[0], "VECTORS %s double \n", point_vec_var_names[var]);
        for (size_t node_gid = 0; node_gid < mesh.num_nodes; node_gid++) {
            fprintf(out[0], "%.15e %.15e %.15e\n",
                    vec_fields(node_gid, var, 0),
                    vec_fields(node_gid, var, 1),
                    vec_fields(node_gid, var, 2));
        } // end for nodes
    } // end for vec_vars


    // vtk scalar vars = (temp)
    for (int var = 0; var < num_point_scalar_vars; var++) {
        fprintf(out[0], "SCALARS %s double 1\n", point_scalar_var_names[var]);
        fprintf(out[0], "LOOKUP_TABLE default\n");
        for (size_t node_gid = 0; node_gid < mesh.num_nodes; node_gid++) {
            fprintf(out[0], "%.15e\n",
                    point_scalar_fields(node_gid, 0));
        } // end for nodes
    } // end for vec_vars

    /*
    ---------------------------------------------------------------------------
    Write the scalar elem variable to file
    ---------------------------------------------------------------------------
    */
    fprintf(out[0], "\n");
    fprintf(out[0], "CELL_DATA %zu \n", mesh.num_elems);

    for (int var = 0; var < num_cell_scalar_vars; var++) {
        fprintf(out[0], "SCALARS %s double 1\n", cell_scalar_var_names[var]); // the 1 is number of scalar components [1:4]
        fprintf(out[0], "LOOKUP_TABLE default\n");
        for (size_t elem_gid = 0; elem_gid < mesh.num_elems; elem_gid++) {
            fprintf(out[0], "%.15e\n", elem_fields(elem_gid, var));
        } // end for elem
    } // end for cell scalar_vars

    fclose(out[0]);

    graphics_times(graphics_id) = time_value;

    // Write time series metadata
    //sprintf(filename, "vtk/Fierro.vtk.series", graphics_id);  // mesh file
    str_output_len = snprintf(filename, max_len, "vtk/Fierro.vtk.series"); 
    if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }
    // mesh file

    out[0] = fopen(filename, "w");

    fprintf(out[0], "{\n");
    fprintf(out[0], "  \"file-series-version\" : \"1.0\",\n");
    fprintf(out[0], "  \"files\" : [\n");

    for (int i = 0; i <= graphics_id; i++) {
        fprintf(out[0], "    { \"name\" : \"Fierro.%05d.vtk\", \"time\" : %12.5e },\n", i, graphics_times(i) );
    }

    // fprintf(out[0], "%12.5e\n", graphics_times(i));
    fprintf(out[0], "  ]\n"); // part 4
    fprintf(out[0], "}"); // part 4

    fclose(out[0]);

    // increment graphics id counter
    graphics_id++;


} // end write vtk old


/////////////////////////////////////////////////////////////////////////////
///
/// \fn concatenate_elem_fields
///
/// \brief A function to calculate the average of elem fields and concatentate into 1 array
///
///
/// \param MaterialPoints a struct containing the material point state arrays
/// \param elem_scalar_fields the scalar fields
/// \param elem_tensor_fields the tensor fields
/// \param elem_in_mat_elem a listing of the element ids the material resides in
/// \param output_elem_state a std::vector of enums specifying the elem avg outputs
/// \param num_mat_elems the number of elements the material resides in
/// \param mat_id the index for the material
///
/////////////////////////////////////////////////////////////////////////////
void MeshWriter::concatenate_elem_fields(const MaterialPoint_t& MaterialPoints,
                             const GaussPoint_t& GaussPoints,
                             DCArrayKokkos<double>& elem_scalar_fields,
                             DCArrayKokkos<double>& elem_tensor_fields,
                             const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
                             const std::vector<material_pt_state>& output_elem_state,
                             const std::vector<gauss_pt_state>& output_gauss_pt_states,
                             const size_t num_mat_elems,
                             const size_t mat_id,
                             const size_t num_elems,
                             const int den_id,
                             const int pres_id,
                             const int sie_id,
                             const int sspd_id,
                             const int mass_id,
                             const int stress_id,
                             const int strain_id,
                             const int vol_id,
                             const int div_id,
                             const int shock_detector_id,
                             const int level_set_id,
                             const int vel_grad_id,
                             const int conductivity_id,
                             const int specific_heat_id)
{

    // --- loop over the material point states

    for (auto field : output_elem_state){
        switch(field){
            // scalar vars
            case material_pt_state::density:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // get elem gid
                    size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

                    // field
                    elem_scalar_fields(den_id, elem_gid) += MaterialPoints.den(mat_id, mat_elem_sid)*
                                                            MaterialPoints.mat_volfrac(mat_id, mat_elem_sid)*
                                                            MaterialPoints.geo_volfrac(mat_id, mat_elem_sid);
                });
                break;
            case material_pt_state::pressure:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // get elem gid
                    size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

                    // field
                    elem_scalar_fields(pres_id, elem_gid) += MaterialPoints.pres(mat_id, mat_elem_sid)*
                                                            MaterialPoints.mat_volfrac(mat_id, mat_elem_sid)*
                                                            MaterialPoints.geo_volfrac(mat_id, mat_elem_sid);
                });
                break;
            case material_pt_state::specific_internal_energy:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // get elem gid
                    size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

                    // field
                    // extensive ie here, but after this function, it will become specific ie
                    elem_scalar_fields(sie_id, elem_gid) += MaterialPoints.mass(mat_id, mat_elem_sid)*
                                                            MaterialPoints.sie(mat_id, mat_elem_sid);
                });
                break;
            case material_pt_state::sound_speed:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // get elem gid
                    size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

                    // field
                    elem_scalar_fields(sspd_id, elem_gid) += MaterialPoints.sspd(mat_id, mat_elem_sid)*
                                                            MaterialPoints.mat_volfrac(mat_id, mat_elem_sid)*
                                                            MaterialPoints.geo_volfrac(mat_id, mat_elem_sid);
                });
                break;
            case material_pt_state::mass:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // get elem gid
                    size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

                    // field
                    elem_scalar_fields(mass_id, elem_gid) += MaterialPoints.mass(mat_id, mat_elem_sid);
                });
                break;
            // ---------------    
            // tensor vars
            // ---------------
            case material_pt_state::stress:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // get elem gid
                    size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

                    // field
                    // average tensor fields, it is always 3D
                    // note: paraview is row-major, CArray convention
                    for (size_t i=0; i<3; i++){
                        for(size_t j=0; j<3; j++){

                            // stress tensor 
                            elem_tensor_fields(stress_id, elem_gid, i, j) +=
                                            MaterialPoints.stress(mat_id, mat_elem_sid,i,j) *
                                            MaterialPoints.mat_volfrac(mat_id, mat_elem_sid)*
                                            MaterialPoints.geo_volfrac(mat_id, mat_elem_sid);
                        } // end for
                    } // end for
                });
                break;

            case material_pt_state::strain:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // get elem gid
                    size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

                    // field
                    // average tensor fields, it is always 3D
                    // note: paraview is row-major, CArray convention
                    for (size_t i=0; i<3; i++){
                        for(size_t j=0; j<3; j++){

                            // strain tensor 
                            elem_tensor_fields(strain_id, elem_gid, i, j) +=
                                            MaterialPoints.strain(mat_id, mat_elem_sid,i,j) *
                                            MaterialPoints.mat_volfrac(mat_id, mat_elem_sid)*
                                            MaterialPoints.geo_volfrac(mat_id, mat_elem_sid);
                        } // end for
                    } // end for
                });
                break;

            // thermal solver vars
            case material_pt_state::thermal_conductivity:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // get elem gid
                    size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

                    // field
                    elem_scalar_fields(conductivity_id, elem_gid) += MaterialPoints.conductivity(mat_id, mat_elem_sid)*
                                                                         MaterialPoints.mat_volfrac(mat_id, mat_elem_sid)*
                                                                         MaterialPoints.geo_volfrac(mat_id, mat_elem_sid);
                });
                break;

            case material_pt_state::specific_heat:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // get elem gid
                    size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

                    // field
                    elem_scalar_fields(specific_heat_id, elem_gid) += MaterialPoints.specific_heat(mat_id, mat_elem_sid)*
                                                                          MaterialPoints.mat_volfrac(mat_id, mat_elem_sid)*
                                                                          MaterialPoints.geo_volfrac(mat_id, mat_elem_sid);
                });
                break;


            // add other variables here

            // not used variables
            case material_pt_state::volume_fraction:
                break;
            case material_pt_state::eroded_flag:
                break;
            case material_pt_state::elastic_modulii:
                break;
            case material_pt_state::shear_modulii:
                break;
            case material_pt_state::poisson_ratios:
                break;
            case material_pt_state::heat_flux:
                break;
        } // end switch
    }// end for over mat point state

    
    // --- add loop over gauss points ---

    // export element centric data
    for (auto field : output_gauss_pt_states){
        switch(field){
            // scalars
            case gauss_pt_state::volume:

                FOR_ALL(elem_gid, 0, num_elems, {
                    elem_scalar_fields(vol_id, elem_gid) = GaussPoints.vol(elem_gid);
                });

                break;
            case gauss_pt_state::divergence_velocity:

                FOR_ALL(elem_gid, 0, num_elems, {
                    elem_scalar_fields(div_id, elem_gid) = GaussPoints.div(elem_gid);
                });

                break;

            case gauss_pt_state::level_set:

                FOR_ALL(elem_gid, 0, num_elems, {
                    elem_scalar_fields(level_set_id, elem_gid) = GaussPoints.level_set(elem_gid);
                });

                break;

            case gauss_pt_state::shock_detector:
                FOR_ALL(elem_gid, 0, num_elems, {
                    elem_scalar_fields(shock_detector_id, elem_gid) = GaussPoints.shock_detector(elem_gid);
                });
                break;

            // tensors
            case gauss_pt_state::gradient_velocity:
                // note: paraview is row-major, CArray convention
                FOR_ALL(elem_gid, 0, num_elems, {
                    for (size_t i=0; i<3; i++){
                        for(size_t j=0; j<3; j++){
                            elem_tensor_fields(vel_grad_id, elem_gid, i, j) = 
                                                GaussPoints.vel_grad(elem_gid, i, j);
                        }
                    } // end for
                });

                break;

            // add other gauss variables here

        } // end switch
    } // end loop over gauss_pt_states


    // --- add end gauss point loop --

} // end of function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn concatenate_mat_fields
///
/// \brief A function to concatentate material fields into 1 array
///
///
/// \param MaterialPoints a struct containing the material point state arrays
/// \param elem_scalar_fields the scalar fields
/// \param elem_tensor_fields the tensor fields
/// \param elem_in_mat_elem a listing of the element ids the material resides in
/// \param output_material_pt_states a std::vector of enums specifying the model
/// \param num_mat_elems the number of elements the material resides in
/// \param mat_id the index for the material
///
/////////////////////////////////////////////////////////////////////////////
void MeshWriter::concatenate_mat_fields(const MaterialPoint_t& MaterialPoints,
                            DCArrayKokkos<double>& mat_elem_scalar_fields,
                            DCArrayKokkos<double>& mat_elem_tensor_fields,
                            DCArrayKokkos<double>& mat_elem_vector_fields,
                            const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
                            const std::vector<material_pt_state>& output_material_pt_states,
                            const size_t num_mat_elems,
                            const size_t mat_id,
                            const int mat_den_id,
                            const int mat_pres_id,
                            const int mat_sie_id,
                            const int mat_sspd_id,
                            const int mat_mass_id,
                            const int mat_mat_volfrac_id,  
                            const int mat_geo_volfrac_id,  
                            const int mat_eroded_id,
                            const int mat_stress_id,
                            const int mat_strain_id,
                            const int mat_conductivity_id,
                            const int mat_specific_heat_id,
                            const int mat_heat_flux_id)
{
  
    // --- loop over the material point states

    for (auto field : output_material_pt_states){
        switch(field){
            // scalar vars
            case material_pt_state::density:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // field
                    mat_elem_scalar_fields(mat_den_id, mat_elem_sid) = MaterialPoints.den(mat_id, mat_elem_sid);
                });
                break;
            case material_pt_state::pressure:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // field
                    mat_elem_scalar_fields(mat_pres_id, mat_elem_sid) = MaterialPoints.pres(mat_id, mat_elem_sid);
                });
                break;
            case material_pt_state::specific_internal_energy:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // field
                    // extensive ie here, but after this function, it will become specific ie
                    mat_elem_scalar_fields(mat_sie_id, mat_elem_sid) = MaterialPoints.sie(mat_id, mat_elem_sid);
                });
                break;
            case material_pt_state::sound_speed:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // field
                    mat_elem_scalar_fields(mat_sspd_id, mat_elem_sid) = MaterialPoints.sspd(mat_id, mat_elem_sid);
                });
                break;
            case material_pt_state::mass:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // field
                    mat_elem_scalar_fields(mat_mass_id, mat_elem_sid) = MaterialPoints.mass(mat_id, mat_elem_sid);
                });
                break;
            case material_pt_state::volume_fraction:
                // material volume fraction
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // field
                    // this is the volume fraction of a material within a part
                    mat_elem_scalar_fields(mat_mat_volfrac_id, mat_elem_sid) = MaterialPoints.mat_volfrac(mat_id, mat_elem_sid);
                });

                // geometric volume fraction
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // field
                    // this is the geometric volume fraction (interface reconstruction)
                    mat_elem_scalar_fields(mat_geo_volfrac_id, mat_elem_sid) = MaterialPoints.geo_volfrac(mat_id, mat_elem_sid);
                });
                break;
            case material_pt_state::eroded_flag:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // field
                    mat_elem_scalar_fields(mat_eroded_id, mat_elem_sid) = (double)MaterialPoints.eroded(mat_id, mat_elem_sid);
                });
                break;
            // ---------------    
            // tensor vars
            // ---------------
            case material_pt_state::stress:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // field
                    // average tensor fields, it is always 3D
                    // note: paraview is row-major, CArray convention
                    for (size_t i=0; i<3; i++){
                        for(size_t j=0; j<3; j++){

                            // stress tensor 
                            mat_elem_tensor_fields(mat_stress_id, mat_elem_sid, i, j) =
                                            MaterialPoints.stress(mat_id, mat_elem_sid,i,j);
                        } // end for
                    } // end for
                });
                break;
            case material_pt_state::strain:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // field
                    // average tensor fields, it is always 3D
                    // note: paraview is row-major, CArray convention
                    for (size_t i=0; i<3; i++){
                        for(size_t j=0; j<3; j++){

                            // strain tensor 
                            mat_elem_tensor_fields(mat_strain_id, mat_elem_sid, i, j) =
                                            MaterialPoints.strain(mat_id, mat_elem_sid,i,j);
                        } // end for
                    } // end for
                });
                break;

            // thermal solver vars
            case material_pt_state::thermal_conductivity:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // get elem gid
                    size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

                    // field
                    mat_elem_scalar_fields(mat_conductivity_id, elem_gid) += MaterialPoints.conductivity(mat_id, mat_elem_sid);
                });
                break;

            case material_pt_state::specific_heat:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // get elem gid
                    size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

                    // field
                    mat_elem_scalar_fields(mat_specific_heat_id, elem_gid) += MaterialPoints.specific_heat(mat_id, mat_elem_sid);
                });
                break;

            case material_pt_state::heat_flux:
                FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

                    // field
                    for (size_t i=0; i<3; i++) {
                        mat_elem_vector_fields(mat_heat_flux_id, mat_elem_sid) += MaterialPoints.q_flux(mat_id, mat_elem_sid,i);
                    }
                });
                break;

            // add other variables here

            // not used variables
            case material_pt_state::elastic_modulii:
                break;
            case material_pt_state::shear_modulii:
                break;
            case material_pt_state::poisson_ratios:
                break;
            
        } // end switch
    }// end for over mat point state


} // end of function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn concatenate_nodal_fields
///
/// \brief A function to calculate the average of elem fields
///
///
/// \param Node a struct containing the material point state arrays
/// \param elem_scalar_fields the scalar fields
/// \param elem_tensor_fields the tensor fields
/// \param elem_in_mat_elem a listing of the element ids the material resides in
/// \param output_node_states a std::vector of enums specifying the model
/// \param num_mat_elems the number of elements the material resides in
/// \param mat_id the index for the material
///
/////////////////////////////////////////////////////////////////////////////
void MeshWriter::concatenate_nodal_fields(const node_t& Node,
                                          DCArrayKokkos<double>& node_scalar_fields,
                                          const DCArrayKokkos<double>& node_vector_fields,
                                          std::vector<node_state>& output_node_states,
                                          double dt,
                                          const size_t num_nodes,
                                          const size_t num_dims,
                                          const int node_mass_id,
                                          const int node_vel_id,
                                          const int node_accel_id,
                                          const int node_coord_id,
                                          const int node_grad_level_set_id,
                                          const int node_temp_id,
                                          const int node_disp_id)
{
    for (auto field : output_node_states){
        switch(field){
            // scalars
            case node_state::mass:

                FOR_ALL(node_gid, 0, num_nodes, {
                    node_scalar_fields(node_mass_id, node_gid) = Node.mass(node_gid);
                });

                break;
            case node_state::temp:
                FOR_ALL(node_gid, 0, num_nodes, {
                    node_scalar_fields(node_temp_id, node_gid) = Node.temp(node_gid);
                });

                break;

            // vector fields

            case node_state::coords:

                FOR_ALL(node_gid, 0, num_nodes, {

                    node_vector_fields(node_coord_id, node_gid, 0) = Node.coords(node_gid, 0);
                    node_vector_fields(node_coord_id, node_gid, 1) = Node.coords(node_gid, 1);
                    if (num_dims == 2) {
                        node_vector_fields(node_coord_id, node_gid, 2) = 0.0;
                    }
                    else{
                        node_vector_fields(node_coord_id, node_gid, 2) = Node.coords(node_gid, 2);
                    } // end if

                }); // end parallel for

                break;
            case node_state::velocity:

                FOR_ALL(node_gid, 0, num_nodes, {

                    // velocity, var is node_vel_id 
                    node_vector_fields(node_vel_id, node_gid, 0) = Node.vel(node_gid, 0);
                    node_vector_fields(node_vel_id, node_gid, 1) = Node.vel(node_gid, 1);
                    if (num_dims == 2) {
                        node_vector_fields(node_vel_id, node_gid, 2) = 0.0;
                    }
                    else{
                        node_vector_fields(node_vel_id, node_gid, 2) = Node.vel(node_gid, 2);
                    } // end if

                    // accellerate, var is node_accel_id            
                    node_vector_fields(node_accel_id, node_gid, 0) = (Node.vel(node_gid, 0) - Node.vel_n0(node_gid, 0))/dt;
                    node_vector_fields(node_accel_id, node_gid, 1) = (Node.vel(node_gid, 1) - Node.vel_n0(node_gid, 1))/dt;
                    if (num_dims == 2) {
                        node_vector_fields(node_accel_id, node_gid, 2) = 0.0;
                    }
                    else{
                        node_vector_fields(node_accel_id, node_gid, 2) = (Node.vel(node_gid, 2) - Node.vel_n0(node_gid, 2))/dt;
                    } // end if

                }); // end parallel for

                break;

            case node_state::displacement:
                FOR_ALL(node_gid, 0, num_nodes, {
                    node_vector_fields(node_disp_id,node_gid,0) = Node.displacement(node_gid,0);
                    node_vector_fields(node_disp_id,node_gid,1) = Node.displacement(node_gid,1);
                    node_vector_fields(node_disp_id,node_gid,2) = Node.displacement(node_gid,2);
                });
                break;

            case node_state::coords_t0:
                break;
                
                
            case node_state::gradient_level_set:

                FOR_ALL(node_gid, 0, num_nodes, {

                    // velocity, var is node_vel_id 
                    node_vector_fields(node_grad_level_set_id, node_gid, 0) = Node.gradient_level_set(node_gid, 0);
                    node_vector_fields(node_grad_level_set_id, node_gid, 1) = Node.gradient_level_set(node_gid, 1);
                    if (num_dims == 2) {
                        node_vector_fields(node_grad_level_set_id, node_gid, 2) = 0.0;
                    }
                    else{
                        node_vector_fields(node_grad_level_set_id, node_gid, 2) = Node.gradient_level_set(node_gid, 2);
                    } // end if

                }); // end parallel for

                break;                
            
            // -- not used vars
            case node_state::force:
                break;

            // heat transer vars
            case node_state::heat_transfer:
                break;
            // tensors
        } // end switch
    } // end for over

    

} // end function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn write_vtu
///
/// \brief Writes a vtu ASCII output file
///
/// \param Simulation mesh
/// \param State data
/// \param Simulation parameters
/// \param current time value
/// \param Vector of all graphics output times
///
/////////////////////////////////////////////////////////////////////////////
void MeshWriter::write_vtu(
    const ViewCArray<double>& node_coords_host,
    const ViewCArray<size_t>& nodes_in_elem_host,
    const DCArrayKokkos<double>& elem_scalar_fields,
    const DCArrayKokkos<double>& elem_tensor_fields,
    const DCArrayKokkos<double>& elem_vector_fields,
    const DCArrayKokkos<double>& node_scalar_fields,
    const DCArrayKokkos<double>& node_vector_fields,
    const std::vector<std::string>& elem_scalar_var_names,
    const std::vector<std::string>& elem_tensor_var_names,
    const std::vector<std::string>& elem_vector_var_names,
    const std::vector<std::string>& node_scalar_var_names,
    const std::vector<std::string>& node_vector_var_names,
    const std::string partname,
    const int graphics_id,
    const size_t num_nodes,
    const size_t num_elems,
    const size_t num_nodes_in_elem,
    const int Pn_order,
    const size_t num_dims,
    const size_t solver_id,
    int mpi_rank,
    int mpi_size,
    const double* diag_mpi_rank_per_elem,
    const double* diag_global_elem_id,
    const double* diag_global_node_id)
{
    FILE* out[20];
    char  filename[512];
    int   max_len = static_cast<int>(sizeof filename);
    int   str_output_len;

    const size_t num_elem_scalar_vars = elem_scalar_var_names.size();
    const size_t num_elem_tensor_vars = elem_tensor_var_names.size();
    const size_t num_elem_vector_vars = elem_vector_var_names.size();

    const size_t num_node_scalar_vars = node_scalar_var_names.size();
    const size_t num_node_vector_vars = node_vector_var_names.size();


    if (mpi_size > 1) {
        str_output_len = snprintf(filename,
                                    static_cast<size_t>(max_len),
                                    "vtk/data/Fierro.solver%zu.%s.%05d_r%04d.vtu",
                                    solver_id,
                                    partname.c_str(),
                                    graphics_id,
                                    mpi_rank);
    }
    else {
        str_output_len = snprintf(filename,
                                    static_cast<size_t>(max_len),
                                    "vtk/data/Fierro.solver%zu.%s.%05d.vtu",
                                    solver_id,
                                    partname.c_str(),
                                    graphics_id);
    }

    if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

    out[0] = fopen(filename, "w");
    if (!out[0]) {
        std::cerr << "[MeshWriter] Failed to open VTU file: " << filename << std::endl;
        return;
    }

    fprintf(out[0], "<?xml version=\"1.0\"?>\n");  
    fprintf(out[0], "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"); 
    fprintf(out[0], "  <UnstructuredGrid> \n");
    fprintf(out[0], "    <Piece NumberOfPoints=\"%zu\" NumberOfCells=\"%zu\">\n", num_nodes, num_elems); 

    /*
    ---------------------------------------------------------------------------
    Write the mesh points
    ---------------------------------------------------------------------------
    */
    fprintf(out[0], "\n");
    fprintf(out[0], "      <!-- Define the mesh nodes -->\n");
    fprintf(out[0], "      <Points>\n");
    fprintf(out[0], "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");

    // write all components of the point coordinates
    for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
        double coord_z = 0.0;
        if(num_dims==3){
            coord_z = node_coords_host(node_gid, 2);
        } 
        fprintf(out[0],
                "          %.15e %.15e %.15e\n",
                node_coords_host(node_gid, 0),
                node_coords_host(node_gid, 1),
                coord_z);
    } // end for
    fprintf(out[0], "        </DataArray>\n");
    fprintf(out[0], "      </Points>\n");

    /*
    ---------------------------------------------------------------------------
    Write the elems
    ---------------------------------------------------------------------------
    */
    fprintf(out[0], "\n");
    fprintf(out[0], "      <!-- Define the elements -->\n");
    fprintf(out[0], "      <Cells>\n");
    fprintf(out[0], "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");  

    // WARNING: look into high-order Pn 2D elements with paraview
    int Pn_order_z = 0;
    if (num_dims == 3){
        Pn_order_z = Pn_order;
    }
    int order[3] = {Pn_order, Pn_order, Pn_order_z};

    // const int num_1D_points = Pn_order+1;

    // write all global point numbers for this elem
    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
        fprintf(out[0], "          ");  // adding indentation before printing nodes in element
        if (num_dims==3 && Pn_order>1){
            CArray<int> convert_fierro_to_vtk((Pn_order+1)*(Pn_order+1)*(Pn_order+1));
            int this_node_lid = 0;
            for (int k = 0; k <= Pn_order; k++) {
                for (int j = 0; j <= Pn_order; j++) {
                    for (int i = 0; i <= Pn_order; i++) {
                        // convert this_point index to the FE index convention
                        size_t vtk_index = PointIndexFromIJK(i, j, k, order);

                        // store the points in this elem according the the finite
                        // element numbering convention
                        convert_fierro_to_vtk(vtk_index) = this_node_lid;
                        // increment the point counting index
                        this_node_lid ++;

                    } // end for icount
                } // end for jcount
            }  // end for kcount
            // Write connectivity: all node IDs for all elements, space-separated
            for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {

                int this_node_lid = 0;
                for (int k = 0; k <= Pn_order; k++) {
                    for (int j = 0; j <= Pn_order; j++) {
                        for (int i = 0; i <= Pn_order; i++) {
                            // size_t node_lid = PointIndexFromIJK(i, j, k, order);
                            int node_lid = convert_fierro_to_vtk(this_node_lid);
                            fprintf(out[0], "%lu ", nodes_in_elem_host(elem_gid, node_lid));
                            this_node_lid ++;
                        }
                    }
                }
            } // end for
        }
        else if (num_dims == 3 && Pn_order == 1){
           // 3D linear hexahedral elements
            for (int node_lid = 0; node_lid < 8; node_lid++) {
                fprintf(out[0], "%lu ", nodes_in_elem_host(elem_gid, node_lid));
            } // end for
        }
        else if (num_dims == 2){
            // 2D linear is the only supported option
            for (int node_lid = 0; node_lid < 4; node_lid++) {
                fprintf(out[0], "%lu ", nodes_in_elem_host(elem_gid, node_lid));
            } // end for
        }
        else {
            std::cout << "ERROR: outputs failed, dimensions and element types are not compatible \n";
        } // end if
        fprintf(out[0], "\n");
    } // end for
    fprintf(out[0], "        </DataArray>\n");

    // Write the element offsets
    fprintf(out[0], "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");  
    size_t count=0;
    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
        count += num_nodes_in_elem;
        fprintf(out[0], "          %lu\n", count); // num points in this elem + all others before it
    } // end for
    fprintf(out[0], "        </DataArray>\n");


    // Write the element types
    fprintf(out[0], "        <DataArray type=\"Int8\" Name=\"types\" format=\"ascii\">\n"); 
    // ----
    // linear element types
    //   VTK_PIXEL = 8,   linear 2D quad with i,j,k indexing (future format for 2D solver)
    //   VTK_Quad = 9,    linear 2D quad with ensight index ordering (current 2D rz convention)
    //   VTK_VOXEL = 11,  linear 3D hex with i,j,k indexing (current format)
    // arbitrary order types
    //   VTK_LAGRANGE_QUADRILATERAL = 70, use this type when a 2D high-order scheme exists
    //   VTK_LAGRANGE_HEXAHEDRON: 72, this is the current 3D high-order 
    //   VTK_HIGHER_ORDER_HEXAHEDRON: 67
    //   VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33
    // element types: https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
    // element types: https://kitware.github.io/vtk-js/api/Common_DataModel_CellTypes.html
    // vtk format: https://www.kitware.com//modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/
    for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
        if (num_dims==3 && Pn_order>1){
            fprintf(out[0], "          %d \n", 72);
        }
        else if (num_dims == 3 && Pn_order == 1){
            // 3D linear hex
            fprintf(out[0], "          %d \n", 11);
        }
        else {
            // 2D ensight mesh ordering
            fprintf(out[0], "          %d \n", 9);
        }
    }
    fprintf(out[0], "        </DataArray>\n");
    fprintf(out[0], "      </Cells>\n");


    /*
    ---------------------------------------------------------------------------
    Write the nodal variables to file
    ---------------------------------------------------------------------------
    */
    // vtk vector vars = (position, velocity)
    fprintf(out[0], "\n");
    fprintf(out[0], "      <!-- Define the node vector data -->\n");
    if (num_node_vector_vars > 0 || num_node_scalar_vars > 0 || diag_global_node_id != nullptr) {

        fprintf(out[0], "      <PointData>\n");

        for (int a_var = 0; a_var < num_node_vector_vars; a_var++) {
            fprintf(out[0], "        <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\">\n", node_vector_var_names[a_var].c_str());
           
            for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
                fprintf(out[0], "          %.15e %.15e %.15e\n",
                        node_vector_fields.host(a_var, node_gid, 0),
                        node_vector_fields.host(a_var, node_gid, 1),
                        node_vector_fields.host(a_var, node_gid, 2));
            }
            fprintf(out[0], "        </DataArray>\n");
        }

        for (int a_var = 0; a_var < num_node_scalar_vars; a_var++) {
            fprintf(out[0], "        <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", node_scalar_var_names[a_var].c_str());
            for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
                fprintf(out[0], "          %.15e\n", node_scalar_fields.host(a_var, node_gid));
            }
            fprintf(out[0], "        </DataArray>\n");
        }

        if (diag_global_node_id != nullptr) {
            fprintf(out[0], "        <DataArray type=\"Float64\" Name=\"global_node_id\" format=\"ascii\">\n");
            for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
                fprintf(out[0], "          %.15e\n", diag_global_node_id[node_gid]);
            }
            fprintf(out[0], "        </DataArray>\n");
        }

        fprintf(out[0], "      </PointData>\n");

    } // end if

    /*
    ---------------------------------------------------------------------------
    Write the elem variables to file
    ---------------------------------------------------------------------------
    */
    fprintf(out[0], "\n");
    fprintf(out[0], "      <!-- Define the cell data -->\n");
    if (num_elem_scalar_vars > 0 || num_elem_tensor_vars > 0 || num_elem_vector_vars > 0 ||
        diag_mpi_rank_per_elem != nullptr || diag_global_elem_id != nullptr) {

        fprintf(out[0], "      <CellData>\n");

        for (int a_var = 0; a_var < num_elem_scalar_vars; a_var++) {

            fprintf(out[0], "        <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", elem_scalar_var_names[a_var].c_str());

            for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
                fprintf(out[0], "          %.15e\n", elem_scalar_fields.host(a_var, elem_gid));
            }
            fprintf(out[0], "        </DataArray>\n");
        }

        for (int a_var = 0; a_var < num_elem_tensor_vars; a_var++) {
            fprintf(out[0], "        <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"9\" format=\"ascii\">\n", elem_tensor_var_names[a_var].c_str());
            
            for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
                for (size_t i = 0; i < 3; i++) {
                    for (size_t j = 0; j < 3; j++) {
                        fprintf(out[0], "          %.15e ", elem_tensor_fields.host(a_var, elem_gid, i, j));
                    }
                }
            }
            fprintf(out[0], "\n");
            fprintf(out[0], "        </DataArray>\n");
        }

        for (int a_var = 0; a_var < num_elem_vector_vars; a_var++) {
            fprintf(out[0], "        <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\">\n", elem_vector_var_names[a_var].c_str());

            for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
                fprintf(out[0], "          %.15e %.15e %.15e\n",
                        elem_vector_fields.host(a_var, elem_gid, 0),
                        elem_vector_fields.host(a_var, elem_gid, 1),
                        elem_vector_fields.host(a_var, elem_gid, 2));
            }
            fprintf(out[0], "        </DataArray>\n");
        }

        if (diag_mpi_rank_per_elem != nullptr) {
            fprintf(out[0], "        <DataArray type=\"Float64\" Name=\"mpi_rank\" format=\"ascii\">\n");
            for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
                fprintf(out[0], "          %.15e\n", diag_mpi_rank_per_elem[elem_gid]);
            }
            fprintf(out[0], "        </DataArray>\n");
        }
        if (diag_global_elem_id != nullptr) {
            fprintf(out[0], "        <DataArray type=\"Float64\" Name=\"global_elem_id\" format=\"ascii\">\n");
            for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
                fprintf(out[0], "          %.15e\n", diag_global_elem_id[elem_gid]);
            }
            fprintf(out[0], "        </DataArray>\n");
        }

        fprintf(out[0], "      </CellData>\n");
    } // end if

    // end of the vtu file
    fprintf(out[0], "    </Piece>\n");
    fprintf(out[0], "  </UnstructuredGrid>\n");
    fprintf(out[0], "</VTKFile>\n");
    
    //-----------------
    // close the vtu file for element fields
    //-----------------
    fclose(out[0]);

} // end write vtu


/////////////////////////////////////////////////////////////////////////////
///
/// \fn writes mesh with the format given in the input.yaml file
///
/// \param Simulation mesh
/// \param Element related state
/// \param Node related state
/// \param Corner related state
/// \param Simulation input parameters
///
/////////////////////////////////////////////////////////////////////////////
void MeshWriter::write_mesh_Pn(swage::Mesh_t& mesh,
    State_t& State,
    SimulationParameters_t& SimulationParamaters,
    double dt,
    double time_value,
    CArray<double> graphics_times,
    std::vector<node_state> node_states,
    std::vector<gauss_pt_state> gauss_pt_states,
    std::vector<material_pt_state> material_pt_states,
    const size_t solver_id,
    elements::ReferenceElement_t& ref_elem)
{


    // node_state is an enum for possible fields (e.g., coords, velocity, etc.), see state.h
    // gauss_pt_state is an enum for possible fields (e.g., vol, divergence, etc.)
    // material_pt_state is an enum for possible fields (e.g., den, pres, etc.)


    // *******************
    //  Update host 
    // *******************

    const size_t num_mats = State.MaterialPoints.num_material_points.size();

    // material point values

    //  Update host data for mat_pt state
    for (auto field : material_pt_states){
        switch(field){
            // scalar vars to write out
            case material_pt_state::density:
                State.MaterialPoints.den.update_host();
                break;
            case material_pt_state::pressure:
                State.MaterialPoints.pres.update_host();
                break;
            case material_pt_state::specific_internal_energy:
                State.MaterialPoints.sie.update_host();
                break;
            case material_pt_state::sound_speed:
                State.MaterialPoints.sspd.update_host();
                break;
            case material_pt_state::mass:
                State.MaterialPoints.mass.update_host();
                break;
            case material_pt_state::volume_fraction:
                State.MaterialPoints.mat_volfrac.update_host();
                State.MaterialPoints.geo_volfrac.update_host();
                break;
            case material_pt_state::eroded_flag:
                State.MaterialPoints.eroded.update_host();
                break;
            // tensor vars to write out
            case material_pt_state::stress:
                State.MaterialPoints.stress.update_host();
                break;

            case material_pt_state::strain:
                State.MaterialPoints.strain.update_host();
                break;

            // additional vars for thermal-mechanical solver
            case material_pt_state::thermal_conductivity:
                State.MaterialPoints.conductivity.update_host();
                break;

            case material_pt_state::specific_heat:
                State.MaterialPoints.specific_heat.update_host();
                break;

            case material_pt_state::heat_flux:
                State.MaterialPoints.q_flux.update_host();
                break;

            // add other variables here

            // not used
            case material_pt_state::elastic_modulii:
                break;
            case material_pt_state::shear_modulii:
                break;
            case material_pt_state::poisson_ratios:
                break;
            
            default:
                std::cout<<"Desired material point state not understood in outputs"<<std::endl;
        } // end switch
    } // end for over mat_pt_states



    // update gauss point values
    for (auto field : gauss_pt_states){
        switch(field){
            // scalar vars to write out
            case gauss_pt_state::volume:
                State.GaussPoints.vol.update_host();
                break;
            case gauss_pt_state::divergence_velocity:
                State.GaussPoints.div.update_host();
                break;
            case gauss_pt_state::level_set:
                State.GaussPoints.level_set.update_host();
                break;      

            // tensor vars to write out
            case gauss_pt_state::gradient_velocity:
                State.GaussPoints.vel_grad.update_host();
                break;
            default:
                std::cout<<"Desired Gauss point state not understood in vtk outputs"<<std::endl;

        } // end switch
    } // end loop

    // nodal values
    for (auto field : node_states){
        switch(field){
            case node_state::mass:
                State.node.mass.update_host();
                break;
            case node_state::temp:
                State.node.temp.update_host();
                break;
            case node_state::coords:
                State.node.coords.update_host();
                break;
            case node_state::velocity:
                State.node.vel.update_host();
                break;
            case node_state::gradient_level_set:
                State.node.gradient_level_set.update_host();
                break;  
            case node_state::displacement:
                State.node.displacement.update_host();
                break;
            case node_state::coords_t0:
                // blank because never changes and update host called in the driver
                break;

            case node_state::force:
                break;

            // heat transer vars
            case node_state::heat_transfer:
                break;

        } // end switch
    } // end for over 
    Kokkos::fence();


    // ******************************************
    //  Build Material and Element state outputs
    // ******************************************

    size_t num_mat_pt_scalar_vars = 0;
    size_t num_mat_pt_tensor_vars = 0;
    size_t num_mat_pt_vector_vars = 0;

    // count the number of material point state vars to write out
    for (auto field : SimulationParamaters.OutputOptions.output_mat_pt_state){
        switch(field){
            // scalar vars to write out
            case material_pt_state::density:
                num_mat_pt_scalar_vars ++;
                break;
            case material_pt_state::pressure:
                num_mat_pt_scalar_vars ++;
                break;
            case material_pt_state::specific_internal_energy:
                num_mat_pt_scalar_vars ++;
                break;
            case material_pt_state::sound_speed:
                num_mat_pt_scalar_vars ++;
                break;
            case material_pt_state::mass:
                num_mat_pt_scalar_vars ++;
                break;
            case material_pt_state::volume_fraction:
                num_mat_pt_scalar_vars ++; // mat volfrac
                num_mat_pt_scalar_vars ++; // geometric volfrac
                break;
            case material_pt_state::eroded_flag:
                num_mat_pt_scalar_vars ++;
                break;
            // tensor vars to write out
            case material_pt_state::stress:
                num_mat_pt_tensor_vars ++;
                break;

            case material_pt_state::strain:
                num_mat_pt_tensor_vars ++;
                break;

            // additional vars for thermal-mechanical solver
            case material_pt_state::thermal_conductivity:
                num_mat_pt_scalar_vars ++;
                break;

            case material_pt_state::specific_heat:
                num_mat_pt_scalar_vars ++;
                break;

            case material_pt_state::heat_flux:
                num_mat_pt_vector_vars ++;
                break;

            // add other variables here

            // not used
            case material_pt_state::elastic_modulii:
                break;
            case material_pt_state::shear_modulii:
                break;
            case material_pt_state::poisson_ratios:
                break;
            
            default:
                std::cout<<"Desired material point state not understood in outputs"<<std::endl;
        } // end switch
    } // end for over mat_pt_states



    size_t num_elem_scalar_vars = 0;
    size_t num_elem_vector_vars = 0;
    size_t num_elem_tensor_vars = 0;

    // count the number of element average fields to write out
    for (auto field : SimulationParamaters.OutputOptions.output_elem_state){
        switch(field){
            // scalar vars to write out
            case material_pt_state::density:
                num_elem_scalar_vars ++;
                break;
            case material_pt_state::pressure:
                num_elem_scalar_vars ++;
                break;
            case material_pt_state::specific_internal_energy:
                num_elem_scalar_vars ++;
                break;
            case material_pt_state::sound_speed:
                num_elem_scalar_vars ++;
                break;
            case material_pt_state::mass:
                num_elem_scalar_vars ++;
                break;
            // tensor vars to write out
            case material_pt_state::stress:
                num_elem_tensor_vars ++;
                break;

            case material_pt_state::strain:
                num_elem_tensor_vars ++;
                break;

            // additional vars for thermal-mechanical solver
            case material_pt_state::thermal_conductivity:
                num_elem_scalar_vars ++;
                break;

            case material_pt_state::specific_heat:
                num_elem_scalar_vars ++;
                break;

            // add other variables here

            // not used
            case material_pt_state::volume_fraction:
                break;
            case material_pt_state::eroded_flag:
                break;
            case material_pt_state::elastic_modulii:
                break;
            case material_pt_state::shear_modulii:
                break;
            case material_pt_state::poisson_ratios:
                break;
            case material_pt_state::heat_flux:
                break;
            default:
                std::cout<<"Desired material point state not understood in outputs"<<std::endl;
        } // end switch
    } // end for over mat_pt_states


    size_t num_gauss_pt_scalar_vars = 0;
    size_t num_gauss_pt_tensor_vars = 0;

    // gauss point values to ouptput
    for (auto field : SimulationParamaters.OutputOptions.output_gauss_pt_state){
        switch(field){
            // scalar vars to write out
            case gauss_pt_state::volume:
                num_gauss_pt_scalar_vars ++;
                break;
            case gauss_pt_state::level_set:
                num_gauss_pt_scalar_vars ++;
                break;
            case gauss_pt_state::divergence_velocity:
                num_gauss_pt_scalar_vars ++;
                break;

            // tensor vars to write out
            case gauss_pt_state::gradient_velocity:
                num_gauss_pt_tensor_vars ++;
                break;
            default:
                std::cout<<"Desired Gauss point state not understood in vtk outputs"<<std::endl;

        } // end switch
    } // end loop

    // add the Gauss point state to the element state
    num_elem_scalar_vars += num_gauss_pt_scalar_vars;
    num_elem_tensor_vars += num_gauss_pt_tensor_vars;


    // Scalar, vector, and tensor value names associated with a elem
    std::vector<std::string> elem_scalar_var_names(num_elem_scalar_vars);
    std::vector<std::string> elem_tensor_var_names(num_elem_tensor_vars);
    std::vector<std::string> elem_vector_var_names(num_elem_vector_vars);

    // Scalar, vector, and tensor values associated with a material in part elems
    std::vector<std::string> mat_elem_scalar_var_names(num_mat_pt_scalar_vars);
    std::vector<std::string> mat_elem_tensor_var_names(num_mat_pt_tensor_vars);
    std::vector<std::string> mat_elem_vector_var_names(num_mat_pt_vector_vars);


    // the ids to access a variable in the mat_scalar_var_name or tensor list
    int mat_den_id = -1;
    int mat_pres_id = -1;
    int mat_sie_id = -1;
    int mat_sspd_id = -1;
    int mat_mass_id = -1;
    int mat_mat_volfrac_id = -1;  
    int mat_geo_volfrac_id = -1;  // geometric volume fraction of part
    int mat_eroded_id = -1;
    int mat_stress_id = -1;
    int mat_strain_id = -1;
    int mat_heat_flux_id = -1;

    int mat_conductivity_id = -1;
    int mat_specific_heat_id = -1;

    // the index for the scalar, vector, and tensor fields
    size_t var = 0;
    size_t vector_var = 0;
    size_t tensor_var = 0;

    // material point state to output
    for (auto field : SimulationParamaters.OutputOptions.output_mat_pt_state){
        switch(field){
            // scalar vars
            case material_pt_state::density:
                mat_elem_scalar_var_names[var] = "mat_den";
                mat_den_id = var;
                var++;
                break;
            case material_pt_state::pressure:
                mat_elem_scalar_var_names[var] = "mat_pres";
                mat_pres_id = var;
                var++;
                break;
            case material_pt_state::specific_internal_energy:
                mat_elem_scalar_var_names[var] = "mat_sie";
                mat_sie_id = var;
                var++;
                break;
            case material_pt_state::sound_speed:
                mat_elem_scalar_var_names[var] = "mat_sspd";
                mat_sspd_id = var;
                var++;
                break;
            case material_pt_state::mass:
                mat_elem_scalar_var_names[var] = "mat_mass";
                mat_mass_id = var;
                var++;
                break;
            case material_pt_state::volume_fraction:
                mat_elem_scalar_var_names[var] = "mat_volfrac";
                mat_mat_volfrac_id = var; 
                var++;

                mat_elem_scalar_var_names[var] = "mat_geo_volfrac";
                mat_geo_volfrac_id = var; 
                var++;
                break;
            case material_pt_state::eroded_flag:
                mat_elem_scalar_var_names[var] = "mat_eroded";
                mat_eroded_id = var;
                var++;
                break;
            // tensor vars
            case material_pt_state::stress:
                mat_elem_tensor_var_names[tensor_var] = "mat_stress";
                mat_stress_id = tensor_var;
                tensor_var++;
                break;

            case material_pt_state::strain:
                mat_elem_tensor_var_names[tensor_var] = "mat_strain";
                mat_strain_id = tensor_var;
                tensor_var++;
                break;


            // additional vars for thermal-mechanical solver
            case material_pt_state::thermal_conductivity:
                mat_elem_scalar_var_names[var] = "mat_thermal_K";
                mat_conductivity_id = var;
                var++;
                break;

            case material_pt_state::specific_heat:
                mat_elem_scalar_var_names[var] = "mat_Cp";
                mat_specific_heat_id = var;
                var++;
                break;

            case material_pt_state::heat_flux:
                mat_elem_vector_var_names[var] = "mat_heat_flux";
                mat_heat_flux_id = vector_var;
                vector_var++;
                break;


            // add other variables here

            // not used
            case material_pt_state::elastic_modulii:
                break;
            case material_pt_state::shear_modulii:
                break;
            case material_pt_state::poisson_ratios:
                break;
            
        } // end switch
    } // end for over mat_pt_states


    // element average fields to output

    // the ids to access a variable in the elem_scalar_var_name or tensor list
    int den_id = -1;
    int pres_id = -1;
    int sie_id = -1;
    int sspd_id = -1;
    int mass_id = -1; 
    int stress_id = -1;
    int strain_id = -1;

    int conductivity_id = -1;
    int specific_heat_id = -1;

    // reset the counters
    var = 0;
    vector_var = 0;
    tensor_var = 0;

    // element state to output
    for (auto field : SimulationParamaters.OutputOptions.output_elem_state){
        switch(field){
            // scalar vars
            case material_pt_state::density:
                elem_scalar_var_names[var] = "den";
                den_id = var;
                var++;
                break;
            case material_pt_state::pressure:
                elem_scalar_var_names[var] = "pres";
                pres_id = var;
                var++;
                break;
            case material_pt_state::specific_internal_energy:
                elem_scalar_var_names[var] = "sie";
                sie_id = var;
                var++;
                break;
            case material_pt_state::sound_speed:
                elem_scalar_var_names[var] = "sspd";
                sspd_id = var;
                var++;
                break;
            case material_pt_state::mass:
                elem_scalar_var_names[var] = "mass";
                mass_id = var;
                var++;
                break;
            // tensor vars
            case material_pt_state::stress:
                elem_tensor_var_names[tensor_var] = "stress";
                stress_id = tensor_var;
                tensor_var++;
                break;

            case material_pt_state::strain:
                elem_tensor_var_names[tensor_var] = "strain";
                strain_id = tensor_var;
                tensor_var++;
                break;

            // heat transfer variables
            case material_pt_state::thermal_conductivity:
                elem_scalar_var_names[var] = "thermal_K";
                conductivity_id = var;
                var++;
                break;

            case material_pt_state::specific_heat:
                elem_scalar_var_names[var] = "Cp";
                specific_heat_id = var;
                var++;
                break;

            // add other variables here

            // not used
            case material_pt_state::volume_fraction:
                break;
            case material_pt_state::eroded_flag:
                break;
            case material_pt_state::elastic_modulii:
                break;
            case material_pt_state::shear_modulii:
                break;
            case material_pt_state::poisson_ratios:
                break;
            case material_pt_state::heat_flux:
                break;
        } // end switch
    } // end for over mat_pt_states

    // append Gauss point vars to the element arrays
    int vol_id = -1;
    int div_id = -1;
    int level_set_id = -1;
    int vel_grad_id = -1;
    int shock_detector_id = -1;


    for (auto field : SimulationParamaters.OutputOptions.output_gauss_pt_state){
        switch(field){
            // scalars
            case gauss_pt_state::volume:
                elem_scalar_var_names[var] = "vol";
                vol_id = var;
                var++;
                break;
            case gauss_pt_state::divergence_velocity:
                elem_scalar_var_names[var] = "div";
                div_id = var;
                var++;
                break;

            case gauss_pt_state::shock_detector:
                elem_scalar_var_names[var] = "shock_detector";
                shock_detector_id = var;
                var++;
                break;

            case gauss_pt_state::level_set:
                elem_scalar_var_names[var] = "level_set";
                level_set_id = var;
                var++;
                break;

            // tensors
            case gauss_pt_state::gradient_velocity:
                elem_tensor_var_names[tensor_var] = "vel_grad";
                vel_grad_id = tensor_var;
                tensor_var++;
                break;
        } // end switch
    } // end loop over gauss_pt_states


    // *******************
    //  nodal values
    // *******************

    size_t num_node_scalar_vars = 0;
    size_t num_node_vector_vars = 0;

    for (auto field : SimulationParamaters.OutputOptions.output_node_state){
        switch(field){
            // --- scalars
            case node_state::mass:
                num_node_scalar_vars ++;
                break;
            case node_state::temp:
                num_node_scalar_vars ++;
                break;
            // -- vectors
            case node_state::coords:
                num_node_vector_vars ++;
                break;
            case node_state::velocity:
                num_node_vector_vars ++; // for velocity
                num_node_vector_vars ++; // for acceleration
                break;
            case node_state::gradient_level_set:
                num_node_vector_vars ++;
                break;
            case node_state::displacement:
                num_node_vector_vars ++;
                break;
            case node_state::coords_t0:
                // not necessary to output as coords at t=0 will be coords_t0
                break;                    
            case node_state::force:
                break;

            // heat transer vars
            case node_state::heat_transfer:
                break;
        } // end switch
    } // end for over 
    Kokkos::fence();


    // Scalar and vector values associated with a node
    std::vector<std::string> node_scalar_var_names(num_node_scalar_vars);
    std::vector<std::string> node_vector_var_names(num_node_vector_vars);

    int node_mass_id = -1;
    int node_vel_id = -1;
    int node_accel_id = -1;
    int node_coord_id = -1;
    int node_temp_id = -1;
    int node_grad_level_set_id = -1;
    int node_disp_id = -1;

    // reset counters for node fields
    var = 0;
    vector_var = 0;
    tensor_var = 0;

    for (auto field : SimulationParamaters.OutputOptions.output_node_state){
        switch(field){
            // scalars
            case node_state::mass:
                node_scalar_var_names[var] = "node_mass";
                node_mass_id = var;
                var++;
                break;
            case node_state::temp:
                node_scalar_var_names[var] = "node_temp";
                node_temp_id = var;
                var++;
                break;

            // vector fields

            case node_state::coords:
                node_vector_var_names[vector_var] = "node_coords";
                node_coord_id = vector_var;
                vector_var++;
                break;

            case node_state::velocity:
                node_vector_var_names[vector_var] = "node_vel";
                node_vel_id = vector_var;
                vector_var++;

                node_vector_var_names[vector_var] = "node_accel";
                node_accel_id = vector_var;
                vector_var++;
                break;

            case node_state::gradient_level_set:
                node_vector_var_names[vector_var] = "node_grad_lvlset";
                node_grad_level_set_id = vector_var;
                vector_var++;
                break;

            case node_state::displacement:
                node_vector_var_names[vector_var] = "node_disp";
                node_disp_id = vector_var;
                vector_var ++;
                break;

            // -- not used vars
            case node_state::coords_t0:
                break;

            case node_state::force:
                break;

            // heat transer vars
            case node_state::heat_transfer:
                break;

            // tensors

        } // end switch
    } // end for over 


    // **************************************
    //  build and save element average fields
    // **************************************

    // short hand
    const size_t num_nodes = mesh.num_nodes;
    const size_t num_elems = mesh.num_elems;
    const size_t num_dims  = mesh.num_dims;
    const size_t num_nodes_in_elem = mesh.num_nodes_in_elem;
    const int Pn_order = mesh.Pn;

    // save the elem state to an array for exporting to graphics files
    DCArrayKokkos<double> elem_scalar_fields(num_elem_scalar_vars, num_elems, "elem_scalars");
    DCArrayKokkos<double> elem_tensor_fields(num_elem_tensor_vars, num_elems, 3, 3, "elem_tensors");
    DCArrayKokkos<double> elem_vector_fields(num_elem_vector_vars, num_elems, 3, "elem_vectors");

    elem_scalar_fields.set_values(0.0);
    elem_tensor_fields.set_values(0.0);


    // -----------------------------------------------------------------------
    // save the output fields to a single element average array for all state
    // -----------------------------------------------------------------------
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {

        // material point and guass point state are concatenated together
        concatenate_elem_fields(State.MaterialPoints,
                                State.GaussPoints,
                                elem_scalar_fields,
                                elem_tensor_fields,
                                State.MaterialToMeshMaps.elem_in_mat_elem,
                                SimulationParamaters.OutputOptions.output_elem_state,
                                SimulationParamaters.OutputOptions.output_gauss_pt_state,
                                State.MaterialToMeshMaps.num_mat_elems.host(mat_id),
                                mat_id,
                                num_elems,
                                den_id,
                                pres_id,
                                sie_id,
                                sspd_id,
                                mass_id,
                                stress_id,
                                strain_id,
                                vol_id,
                                div_id,
                                shock_detector_id,
                                level_set_id,
                                vel_grad_id,
                                conductivity_id,
                                specific_heat_id);
    } // end for mats

    // make specific fields for the element average
    if (sie_id>=0){
        FOR_ALL(elem_gid, 0, num_elems, {
            // get sie by dividing by the mass
            elem_scalar_fields(sie_id, elem_gid) /= (elem_scalar_fields(mass_id, elem_gid)+1.e-20); 
        });
    } // end if

    Kokkos::fence();
    elem_scalar_fields.update_host();
    elem_tensor_fields.update_host();

    // ************************
    //  Build the nodal fields 
    // ************************

    // save the nodal fields to an array for exporting to graphics files
    DCArrayKokkos<double> node_scalar_fields(num_node_scalar_vars, num_nodes, "node_scalars");
    DCArrayKokkos<double> node_vector_fields(num_node_vector_vars, num_nodes, 3, "node_vectors");
    Kokkos::fence();

    concatenate_nodal_fields(State.node,
                             node_scalar_fields,
                             node_vector_fields,
                             SimulationParamaters.OutputOptions.output_node_state,
                             dt,
                             num_nodes,
                             num_dims,
                             node_mass_id,
                             node_vel_id,
                             node_accel_id,
                             node_coord_id,
                             node_grad_level_set_id,
                             node_temp_id,
                             node_disp_id);


    Kokkos::fence();
    node_scalar_fields.update_host();
    node_vector_fields.update_host();


    // ********************************
    //  Write the nodal and elem fields 
    // ********************************

    if (SimulationParamaters.OutputOptions.format == output_options::viz ||
        SimulationParamaters.OutputOptions.format == output_options::viz_and_state) {

        int mpi_rank = 0;
        int mpi_size = 1;
        mesh_io_mpi_detail::query_world_rank_size(mpi_rank, mpi_size);

        MPI_Barrier(MPI_COMM_WORLD);

        if (mpi_rank == 0 && solver_id == 0 && graphics_id == 0) {
            struct stat st_rm;
            if (stat("vtk", &st_rm) == 0) {
                (void)system("rm -f vtk/Fierro*");
            }
            if (stat("vtk/data", &st_rm) == 0) {
                (void)system("rm -f vtk/data/Fierro*");
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        struct stat st;
        if (stat("vtk", &st) != 0) {
            if (system("mkdir vtk") == 1) {
                std::cout << "Unable to make vtk directory" << std::endl;
            }
        }
        if (stat("vtk/data", &st) != 0) {
            if (system("mkdir vtk/data") == 1) {
                std::cout << "Unable to make vtk/data directory" << std::endl;
            }
        }

        // Per-rank VTU: write owned elements only (rows 0..num_owned_elems-1 of nodes_in_elem), but
        // Points must list all local nodes (0..num_nodes-1). Connectivity is full local indexing;
        // boundary elements reference ghost nodes — truncating coords to num_owned_nodes (the prior
        // bug) misaligned indices and produced inverted/wrong cells in ParaView.
        const size_t n_owned_elems = mesh.num_owned_elems; // num_elems; 


        const std::string elem_fields_name = "fields";

        ViewCArray<double> node_coords_host(&State.node.coords.host(0, 0), num_nodes, num_dims);
        ViewCArray<size_t> nodes_in_elem_host(&mesh.nodes_in_elem.host(0, 0), n_owned_elems, num_nodes_in_elem);

        // VTK diagnostics (optional CellData / PointData in write_vtu): host-side scratch arrays,
        // then const pointers so write_vtu can emit one value per owned cell or per owned node.
        //   diag_mpi_rank_elem — CellData "mpi_rank" (which MPI rank owns each exported element).
        //   diag_global_elem   — CellData "global_elem_id" from mesh.local_to_global_elem_mapping (or local e).
        //   diag_global_node   — PointData "global_node_id" from mesh.local_to_global_node_mapping (or local n).
        //   p_*                — nullptr until filled; passed to write_vtu (nullptr disables that array).
        std::vector<double> diag_mpi_rank_elem;
        std::vector<double> diag_global_elem;
        std::vector<double> diag_global_node;
        const double*       p_rank_elem = nullptr;
        const double*       p_glob_elem = nullptr;
        const double*       p_glob_node = nullptr;

        if (n_owned_elems > 0 && num_nodes > 0) {
            diag_mpi_rank_elem.assign(n_owned_elems, static_cast<double>(mpi_rank));
            p_rank_elem = diag_mpi_rank_elem.data();

            diag_global_elem.resize(n_owned_elems);
            if (mpi_size > 1 || mesh.num_elems > mesh.num_owned_elems) {
                mesh.local_to_global_elem_mapping.update_host();
                for (size_t e = 0; e < n_owned_elems; e++) {
                    diag_global_elem[e] =
                        static_cast<double>(mesh.local_to_global_elem_mapping.host(e));
                }
            }
            else {
                for (size_t e = 0; e < n_owned_elems; e++) {
                    diag_global_elem[e] = static_cast<double>(e);
                }
            }
            p_glob_elem = diag_global_elem.data();

            diag_global_node.resize(num_nodes);
            if (mpi_size > 1 || mesh.num_nodes > mesh.num_owned_nodes) {
                mesh.local_to_global_node_mapping.update_host();
                for (size_t n = 0; n < num_nodes; n++) {
                    diag_global_node[n] =
                        static_cast<double>(mesh.local_to_global_node_mapping.host(n));
                }
            }
            else {
                for (size_t n = 0; n < num_nodes; n++) {
                    diag_global_node[n] = static_cast<double>(n);
                }
            }
            p_glob_node = diag_global_node.data();

            write_vtu(node_coords_host,
                      nodes_in_elem_host,
                      elem_scalar_fields,
                      elem_tensor_fields,
                      elem_vector_fields,
                      node_scalar_fields,
                      node_vector_fields,
                      elem_scalar_var_names,
                      elem_tensor_var_names,
                      elem_vector_var_names,
                      node_scalar_var_names,
                      node_vector_var_names,
                      elem_fields_name,
                      graphics_id,
                      num_nodes,
                      n_owned_elems,
                      num_nodes_in_elem,
                      Pn_order,
                      num_dims,
                      solver_id,
                      mpi_rank,
                      mpi_size,
                      p_rank_elem,
                      p_glob_elem,
                      p_glob_node);
        }


        // ********************************
        //  Build and write the mat fields 
        // ********************************

        // note: the file path and folder was created in the elem and node outputs
        size_t num_mat_files_written = 0;
        if(num_mat_pt_scalar_vars > 0 || num_mat_pt_tensor_vars >0){

            for (int mat_id = 0; mat_id < num_mats; mat_id++) {

                const size_t num_mat_elems = State.MaterialToMeshMaps.num_mat_elems.host(mat_id);

                // only save material data if the mat lives on the mesh, ie. has state allocated
                if (num_mat_elems>0){

                    // set the nodal vars to zero size, we don't write these fields again
                    node_scalar_var_names.clear();
                    node_vector_var_names.clear();

                    // the arrays storing all the material field data
                    DCArrayKokkos<double> mat_elem_scalar_fields(num_mat_pt_scalar_vars, num_mat_elems, "mat_pt_scalars");
                    DCArrayKokkos<double> mat_elem_tensor_fields(num_mat_pt_tensor_vars, num_mat_elems, 3, 3, "mat_pt_tensors");
                    DCArrayKokkos<double> mat_elem_vector_fields(num_mat_pt_vector_vars, num_mat_elems, 3, "mat_pt_vectors");


                    // concatenate material fields into a single array
                    concatenate_mat_fields(State.MaterialPoints,
                                           mat_elem_scalar_fields,
                                           mat_elem_tensor_fields,
                                           mat_elem_vector_fields,
                                           State.MaterialToMeshMaps.elem_in_mat_elem,
                                           SimulationParamaters.OutputOptions.output_mat_pt_state,
                                           num_mat_elems,
                                           mat_id,
                                           mat_den_id,
                                           mat_pres_id,
                                           mat_sie_id,
                                           mat_sspd_id,
                                           mat_mass_id,
                                           mat_mat_volfrac_id,
                                           mat_geo_volfrac_id,  
                                           mat_eroded_id,
                                           mat_stress_id,
                                           mat_strain_id,
                                           mat_conductivity_id,
                                           mat_specific_heat_id,
                                           mat_heat_flux_id);
                    Kokkos::fence();
                    mat_elem_scalar_fields.update_host();
                    mat_elem_tensor_fields.update_host();
                    mat_elem_vector_fields.update_host();


                    std::string str_mat_val = std::to_string(mat_id);                       
                    std::string mat_fields_name = "mat";
                    mat_fields_name += str_mat_val;  // add the mat number

                    // save the nodes belonging to this part (i.e., the material)
                    DCArrayKokkos <double> mat_node_coords(num_nodes,num_dims, "mat_node_coords");
                    DCArrayKokkos <size_t> mat_nodes_in_mat_elem(num_mat_elems, num_nodes_in_elem, "mat_nodes_in_mat_elem");

                    // the number of actual nodes belonging to the part (i.e., the material)
                    size_t num_mat_nodes = 0;

                    // build a unique mesh (element and nodes) for the material (i.e., the part)
                    build_material_elem_node_lists(mesh,
                                                   State.node.coords,
                                                   mat_node_coords,
                                                   mat_nodes_in_mat_elem,
                                                   State.MaterialToMeshMaps.elem_in_mat_elem,
                                                   mat_id,
                                                   num_mat_nodes,
                                                   num_mat_elems,
                                                   num_nodes_in_elem,
                                                   num_dims);

                    ViewCArray <double> mat_node_coords_host(&mat_node_coords.host(0,0), num_mat_nodes, num_dims);
                    ViewCArray <size_t> mat_nodes_in_elem_host(&mat_nodes_in_mat_elem.host(0,0), num_mat_elems, num_nodes_in_elem);

                    // write out a vtu file this 
                    write_vtu(mat_node_coords_host,
                              mat_nodes_in_elem_host,
                              mat_elem_scalar_fields,
                              mat_elem_tensor_fields,
                              mat_elem_vector_fields,
                              node_scalar_fields,
                              node_vector_fields,
                              mat_elem_scalar_var_names,
                              mat_elem_tensor_var_names,
                              mat_elem_vector_var_names,
                              node_scalar_var_names,
                              node_vector_var_names,
                              mat_fields_name,
                              graphics_id,
                              num_mat_nodes,
                              num_mat_elems,
                              num_nodes_in_elem,
                              Pn_order,
                              num_dims,
                              solver_id,
                              mpi_rank,
                              mpi_size,
                              nullptr,
                              nullptr,
                              nullptr);

                    // --- Dataset B: Gauss-point cloud for ParaView interpolation ---
                    if (num_mat_pt_scalar_vars > 0 || num_mat_pt_tensor_vars > 0) {

                        std::string pn_name = "mat" + str_mat_val + "_Pn";

                        write_vtu_Pn(mesh,
                                     State,
                                     ref_elem,
                                     mat_elem_scalar_var_names,
                                     mat_elem_tensor_var_names,
                                     pn_name,
                                     graphics_id,
                                     mat_id,
                                     num_mat_elems,
                                     num_nodes_in_elem,
                                     num_dims,
                                     solver_id,
                                     mpi_rank,
                                     mpi_size,
                                     mat_den_id,
                                     mat_pres_id,
                                     mat_sie_id,
                                     mat_sspd_id,
                                     mat_mass_id,
                                     mat_mat_volfrac_id,
                                     mat_geo_volfrac_id,
                                     mat_eroded_id,
                                     mat_stress_id,
                                     mat_strain_id,
                                     mat_conductivity_id,
                                     mat_specific_heat_id);

                    }


                    num_mat_files_written++;

                } // end for mat_id

            } // end if material is on the mesh

        } // end if mat variables are to be written


        // *************************************************
        //  write Paraview files to open the graphics files
        // *************************************************

        // save the graphics time
        graphics_times(graphics_id) = time_value;

        std::vector<unsigned long long> local_mat_counts(num_mats, 0ULL);
        for (size_t mi = 0; mi < num_mats; mi++) {
            local_mat_counts[mi] =
                static_cast<unsigned long long>(State.MaterialToMeshMaps.num_mat_elems.host(mi));
        }
        std::vector<unsigned long long> gathered_mat_elems;
        if (mpi_rank == 0) {
            gathered_mat_elems.assign(static_cast<size_t>(mpi_size) * num_mats, 0ULL);
        }
        MPI_Gather(local_mat_counts.data(),
                   static_cast<int>(num_mats),
                   MPI_UNSIGNED_LONG_LONG,
                   mpi_rank == 0 ? gathered_mat_elems.data() : nullptr,
                   static_cast<int>(num_mats),
                   MPI_UNSIGNED_LONG_LONG,
                   0,
                   MPI_COMM_WORLD);

        const unsigned long long      local_owned_elems_ull = static_cast<unsigned long long>(n_owned_elems);
        std::vector<unsigned long long> gathered_owned_elems;
        if (mpi_rank == 0) {
            gathered_owned_elems.assign(static_cast<size_t>(mpi_size), 0ULL);
        }
        MPI_Gather(&local_owned_elems_ull,
                   1,
                   MPI_UNSIGNED_LONG_LONG,
                   mpi_rank == 0 ? gathered_owned_elems.data() : nullptr,
                   1,
                   MPI_UNSIGNED_LONG_LONG,
                   0,
                   MPI_COMM_WORLD);

        // check to see if an mesh state was written 
        bool write_mesh_state = false;
        if( num_elem_scalar_vars > 0 ||
            num_elem_tensor_vars > 0 ||
            num_node_scalar_vars > 0 ||
            num_node_vector_vars > 0)
        {
            write_mesh_state = true;
        }

        // check to see if a mat state was written
        bool write_mat_pt_state = false;
        if( num_mat_pt_scalar_vars > 0 ||
            num_mat_pt_tensor_vars > 0)
        {
             write_mat_pt_state = true;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (mpi_rank == 0) {
            // call the vtm file writer
            std::string mat_fields_name = "mat";
            //write_vtm_Pn(graphics_times,
            //        elem_fields_name,
            //        mat_fields_name,
            //        time_value,
            //        graphics_id,
            //        num_mats,
            //        write_mesh_state,
            //        write_mat_pt_state,
            //        solver_id);
            write_vtm_Pn(graphics_times,
                         elem_fields_name,
                         mat_fields_name,
                         time_value,
                         graphics_id,
                         num_mats,
                         write_mesh_state,
                         write_mat_pt_state,
                         solver_id,
                         mpi_size,
                         gathered_owned_elems.empty() ? nullptr : gathered_owned_elems.data(),
                         gathered_mat_elems.empty() ? nullptr : gathered_mat_elems.data());

            // call the pvd file writer
            write_pvd(graphics_times,
                        time_value,
                        graphics_id,
                        solver_id,
                        mpi_rank);
        }


        // increment graphics id counter
        graphics_id++; // this is private variable in the class

    } // end if viz paraview output is to be written


    // STATE
    if (SimulationParamaters.OutputOptions.format == output_options::state ||
        SimulationParamaters.OutputOptions.format == output_options::viz_and_state) {

        /* write_material_point_state(mesh,
                                   State,
                                   SimulationParamaters,
                                   time_value,
                                   graphics_times,
                                   node_states,
                                   gauss_pt_states,
                                   material_pt_states); */
        write_text_state_Pn(mesh,
                            State,
                            ref_elem,
                            mat_elem_scalar_var_names,
                            mat_elem_tensor_var_names,
                            node_scalar_var_names,
                            node_vector_var_names,
                            mesh.num_nodes_in_elem,
                            mesh.num_dims,
                            time_value,
                            mat_den_id,
                            mat_pres_id,
                            mat_sie_id,
                            mat_sspd_id,
                            mat_mass_id,
                            mat_mat_volfrac_id,
                            mat_geo_volfrac_id,
                            mat_eroded_id,
                            mat_stress_id,
                            mat_strain_id,
                            mat_conductivity_id,
                            mat_specific_heat_id,
                            node_mass_id,
                            node_vel_id,
                            node_coord_id,
                            node_temp_id,
                            node_grad_level_set_id,
                            node_disp_id);

    } // end if state is to be written


    // will drop ensight outputs in the near future
    if (SimulationParamaters.OutputOptions.format == output_options::ensight){
       write_ensight(mesh,
                     State,
                     SimulationParamaters,
                     dt,
                     time_value,
                     graphics_times,
                     node_states,
                     gauss_pt_states,
                     material_pt_states);
    }

    return;

} // end write_mesh_Pn


/////////////////////////////////////////////////////////////////////////////
///
/// \fn write_vtu_Pn
///
/// \brief Writes a vtu ASCII output file for gauss points
///
/// \param Simulation mesh
/// \param State data
/// \param Simulation parameters
/// \param current time value
/// \param Vector of all graphics output times
///
/////////////////////////////////////////////////////////////////////////////

void MeshWriter::write_vtu_Pn(
    const swage::Mesh_t&                   mesh,
    const State_t&                       State,
    elements::ReferenceElement_t&             ref_elem,
    const std::vector<std::string>&      mat_scalar_var_names,   // length == num_mat_pt_scalar_vars
    const std::vector<std::string>&      mat_tensor_var_names,   // length == num_mat_pt_tensor_vars
    const std::string&                   partname,               // e.g. "mat0_Pn"
    const int                            graphics_id,
    const int                            mat_id,
    const size_t                         num_mat_elems,
    const size_t                         num_nodes_in_elem,
    const size_t                         num_dims,
    const size_t                         solver_id,
    int                                  mpi_rank,
    int                                  mpi_size,
    // field slot IDs (-1 means "not requested")
    const int mat_den_id,
    const int mat_pres_id,
    const int mat_sie_id,
    const int mat_sspd_id,
    const int mat_mass_id,
    const int mat_mat_volfrac_id,
    const int mat_geo_volfrac_id,
    const int mat_eroded_id,
    const int mat_stress_id,
    const int mat_strain_id,
    const int mat_conductivity_id,
    const int mat_specific_heat_id
)
{
    // -----------------------------------------------------------------------
    //  Derived sizes
    // -----------------------------------------------------------------------

    // Number of Gauss points per element: (2*Pn)^num_dims
    // Read directly from the reference element so we stay consistent with
    // whatever quadrature order was set up for this run.
    const size_t num_gp_per_elem   = ref_elem.qpt_grad_basis.dims(0);
    const size_t num_total_gp      = num_mat_elems * num_gp_per_elem;

    const size_t num_scalar_vars   = mat_scalar_var_names.size();
    const size_t num_tensor_vars   = mat_tensor_var_names.size();

    // -----------------------------------------------------------------------
    //  Open file
    // -----------------------------------------------------------------------

    FILE* fp;
    char  filename[512];
    int   max_len = static_cast<int>(sizeof filename);
    int   str_output_len;

    // File lives next to the standard per-material .vtu files.
    // The "_Pn" suffix (already part of partname) distinguishes this
    // point-cloud dataset from the element-averaged dataset written by
    // write_vtu. Filename convention matches write_vtu exactly so the
    // .vtm DataSet references (Fierro.solverN.matM_Pn.FFFFF_rRRRR.vtu)
    // resolve correctly under MPI.
    if (mpi_size > 1) {
        str_output_len = snprintf(filename,
                                   static_cast<size_t>(max_len),
                                   "vtk/data/Fierro.solver%zu.%s.%05d_r%04d.vtu",
                                   solver_id,
                                   partname.c_str(),
                                   graphics_id,
                                   mpi_rank);
    }
    else {
        str_output_len = snprintf(filename,
                                   static_cast<size_t>(max_len),
                                   "vtk/data/Fierro.solver%zu.%s.%05d.vtu",
                                   solver_id,
                                   partname.c_str(),
                                   graphics_id);
    }

    if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

    fp = fopen(filename, "w");
    if (!fp) {
        std::cerr << "write_vtu_Pn: could not open " << filename << std::endl;
        return;
    }

    // -----------------------------------------------------------------------
    //  VTK file header
    // -----------------------------------------------------------------------

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" "
                "byte_order=\"LittleEndian\">\n");
    fprintf(fp, "  <UnstructuredGrid>\n");
    fprintf(fp, "    <Piece NumberOfPoints=\"%zu\" NumberOfCells=\"%zu\">\n",
            num_total_gp, num_total_gp);

    // -----------------------------------------------------------------------
    //  Points  –  physical coordinates of every Gauss point
    //
    //  Mapping:  x_phys[dim] = Σ_a  N_a(ξ_gp) * x_node_a[dim]
    // -----------------------------------------------------------------------

    fprintf(fp, "\n");
    fprintf(fp, "      <!-- Gauss point physical coordinates "
                "(isoparametric mapping) -->\n");
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float64\" "
                "NumberOfComponents=\"3\" format=\"ascii\">\n");

    DCArrayKokkos <double> x_phys(num_mat_elems, num_gp_per_elem, 3);
    x_phys.set_values(0);
    FOR_ALL(elem, 0, num_mat_elems, {
        const size_t elem_id =
            State.MaterialToMeshMaps.elem_in_mat_elem(mat_id, elem);

        for (size_t gp = 0; gp < num_gp_per_elem; gp++) {

            for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                const size_t node_gid = mesh.nodes_in_elem(elem_id, node_lid);
                const double N        = ref_elem.qpt_basis(gp, node_lid);

                for (size_t dim = 0; dim < num_dims; dim++) {
                    x_phys(elem, gp, dim) += N * State.node.coords(node_gid, dim);
                }
            }

        } // end gp
    });
    x_phys.update_host();
    for (size_t elem = 0; elem < num_mat_elems; elem++) {

        for (size_t gp = 0; gp < num_gp_per_elem; gp++) {

            // VTK always expects 3 components; pad Z to 0 for 2-D runs.
            fprintf(fp, "          %.15e %.15e %.15e\n",
                    x_phys.host(elem, gp, 0),
                    x_phys.host(elem, gp, 1),
                    (num_dims == 3) ? x_phys.host(elem, gp, 2) : 0.0);

        } // end gp
    } // end elem

    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");

    // -----------------------------------------------------------------------
    //  Cells  –  one VTK_VERTEX (type = 1) per Gauss point
    //
    //  connectivity[i] = i   (each vertex cell references exactly one point)
    //  offsets[i]      = i+1
    //  types[i]        = 1   (VTK_VERTEX)
    // -----------------------------------------------------------------------

    fprintf(fp, "\n");
    fprintf(fp, "      <!-- One VTK_VERTEX cell per Gauss point -->\n");
    fprintf(fp, "      <Cells>\n");

    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" "
                "format=\"ascii\">\n");
    for (size_t pt = 0; pt < num_total_gp; pt++) {
        fprintf(fp, "          %zu\n", pt);
    }
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" "
                "format=\"ascii\">\n");
    for (size_t pt = 0; pt < num_total_gp; pt++) {
        fprintf(fp, "          %zu\n", pt + 1);
    }
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "        <DataArray type=\"Int8\" Name=\"types\" "
                "format=\"ascii\">\n");
    for (size_t pt = 0; pt < num_total_gp; pt++) {
        fprintf(fp, "          1\n");  // VTK_VERTEX
    }
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "      </Cells>\n");

    // -----------------------------------------------------------------------
    //  PointData  –  material-point fields at every Gauss point
    //
    //  The flat material-point index is obtained the same way the device
    //  kernels do it:  pt_id = State.points_in_mat_elem.host(elem, gp)
    // -----------------------------------------------------------------------

    if (num_scalar_vars > 0 || num_tensor_vars > 0) {

        fprintf(fp, "\n");
        fprintf(fp, "      <!-- Material-point field data at each "
                    "Gauss point -->\n");
        fprintf(fp, "      <PointData>\n");

        // -------------------------------------------------------------------
        //  Scalar fields
        // -------------------------------------------------------------------

        // density
        if (mat_den_id >= 0) {
            fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
                        "format=\"ascii\">\n",
                    mat_scalar_var_names[mat_den_id].c_str());
            for (size_t elem = 0; elem < num_mat_elems; elem++) {
                for (size_t gp = 0; gp < num_gp_per_elem; gp++) {
                    const size_t pt_id =
                        State.points_in_mat_elem.host(elem, gp);
                    fprintf(fp, "          %.15e\n",
                            State.MaterialPoints.den.host(mat_id, pt_id));
                }
            }
            fprintf(fp, "        </DataArray>\n");
        }

        // pressure
        if (mat_pres_id >= 0) {
            fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
                        "format=\"ascii\">\n",
                    mat_scalar_var_names[mat_pres_id].c_str());
            for (size_t elem = 0; elem < num_mat_elems; elem++) {
                for (size_t gp = 0; gp < num_gp_per_elem; gp++) {
                    const size_t pt_id =
                        State.points_in_mat_elem.host(elem, gp);
                    fprintf(fp, "          %.15e\n",
                            State.MaterialPoints.pres.host(mat_id, pt_id));
                }
            }
            fprintf(fp, "        </DataArray>\n");
        }

        // specific internal energy
        if (mat_sie_id >= 0) {
            fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
                        "format=\"ascii\">\n",
                    mat_scalar_var_names[mat_sie_id].c_str());
            for (size_t elem = 0; elem < num_mat_elems; elem++) {
                for (size_t gp = 0; gp < num_gp_per_elem; gp++) {
                    const size_t pt_id =
                        State.points_in_mat_elem.host(elem, gp);
                    fprintf(fp, "          %.15e\n",
                            State.MaterialPoints.sie.host(mat_id, pt_id));
                }
            }
            fprintf(fp, "        </DataArray>\n");
        }

        // sound speed
        if (mat_sspd_id >= 0) {
            fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
                        "format=\"ascii\">\n",
                    mat_scalar_var_names[mat_sspd_id].c_str());
            for (size_t elem = 0; elem < num_mat_elems; elem++) {
                for (size_t gp = 0; gp < num_gp_per_elem; gp++) {
                    const size_t pt_id =
                        State.points_in_mat_elem.host(elem, gp);
                    fprintf(fp, "          %.15e\n",
                            State.MaterialPoints.sspd.host(mat_id, pt_id));
                }
            }
            fprintf(fp, "        </DataArray>\n");
        }

        // mass
        if (mat_mass_id >= 0) {
            fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
                        "format=\"ascii\">\n",
                    mat_scalar_var_names[mat_mass_id].c_str());
            for (size_t elem = 0; elem < num_mat_elems; elem++) {
                for (size_t gp = 0; gp < num_gp_per_elem; gp++) {
                    const size_t pt_id =
                        State.points_in_mat_elem.host(elem, gp);
                    fprintf(fp, "          %.15e\n",
                            State.MaterialPoints.mass.host(mat_id, pt_id));
                }
            }
            fprintf(fp, "        </DataArray>\n");
        }

        // material volume fraction
        if (mat_mat_volfrac_id >= 0) {
            fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
                        "format=\"ascii\">\n",
                    mat_scalar_var_names[mat_mat_volfrac_id].c_str());
            for (size_t elem = 0; elem < num_mat_elems; elem++) {
                for (size_t gp = 0; gp < num_gp_per_elem; gp++) {
                    const size_t pt_id =
                        State.points_in_mat_elem.host(elem, gp);
                    fprintf(fp, "          %.15e\n",
                            State.MaterialPoints.mat_volfrac.host(mat_id, pt_id));
                }
            }
            fprintf(fp, "        </DataArray>\n");
        }

        // geometric volume fraction
        if (mat_geo_volfrac_id >= 0) {
            fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
                        "format=\"ascii\">\n",
                    mat_scalar_var_names[mat_geo_volfrac_id].c_str());
            for (size_t elem = 0; elem < num_mat_elems; elem++) {
                for (size_t gp = 0; gp < num_gp_per_elem; gp++) {
                    const size_t pt_id =
                        State.points_in_mat_elem.host(elem, gp);
                    fprintf(fp, "          %.15e\n",
                            State.MaterialPoints.geo_volfrac.host(mat_id, pt_id));
                }
            }
            fprintf(fp, "        </DataArray>\n");
        }

        // eroded flag (cast to double so VTK Float64 array stays consistent)
        if (mat_eroded_id >= 0) {
            fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
                        "format=\"ascii\">\n",
                    mat_scalar_var_names[mat_eroded_id].c_str());
            for (size_t elem = 0; elem < num_mat_elems; elem++) {
                for (size_t gp = 0; gp < num_gp_per_elem; gp++) {
                    const size_t pt_id =
                        State.points_in_mat_elem.host(elem, gp);
                    fprintf(fp, "          %.15e\n",
                            static_cast<double>(
                                State.MaterialPoints.eroded.host(mat_id, pt_id)));
                }
            }
            fprintf(fp, "        </DataArray>\n");
        }

        // thermal conductivity
        if (mat_conductivity_id >= 0) {
            fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
                        "format=\"ascii\">\n",
                    mat_scalar_var_names[mat_conductivity_id].c_str());
            for (size_t elem = 0; elem < num_mat_elems; elem++) {
                for (size_t gp = 0; gp < num_gp_per_elem; gp++) {
                    const size_t pt_id =
                        State.points_in_mat_elem.host(elem, gp);
                    fprintf(fp, "          %.15e\n",
                            State.MaterialPoints.conductivity.host(mat_id, pt_id));
                }
            }
            fprintf(fp, "        </DataArray>\n");
        }

        // specific heat
        if (mat_specific_heat_id >= 0) {
            fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
                        "format=\"ascii\">\n",
                    mat_scalar_var_names[mat_specific_heat_id].c_str());
            for (size_t elem = 0; elem < num_mat_elems; elem++) {
                for (size_t gp = 0; gp < num_gp_per_elem; gp++) {
                    const size_t pt_id =
                        State.points_in_mat_elem.host(elem, gp);
                    fprintf(fp, "          %.15e\n",
                            State.MaterialPoints.specific_heat.host(mat_id, pt_id));
                }
            }
            fprintf(fp, "        </DataArray>\n");
        }

        // -------------------------------------------------------------------
        //  Tensor fields  (stored as 9-component Float64 arrays, row-major)
        //
        //  ParaView convention:  Txx Txy Txz  Tyx Tyy Tyz  Tzx Tzy Tzz
        //  which matches CArray row-major storage (i outer, j inner).
        // -------------------------------------------------------------------

        // stress
        if (mat_stress_id >= 0) {
            fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
                        "NumberOfComponents=\"9\" format=\"ascii\">\n",
                    mat_tensor_var_names[mat_stress_id].c_str());
            for (size_t elem = 0; elem < num_mat_elems; elem++) {
                for (size_t gp = 0; gp < num_gp_per_elem; gp++) {
                    const size_t pt_id =
                        State.points_in_mat_elem.host(elem, gp);
                    fprintf(fp, "         ");
                    for (size_t i = 0; i < 3; i++) {
                        for (size_t j = 0; j < 3; j++) {
                            fprintf(fp, " %.15e",
                                    State.MaterialPoints.stress.host(
                                        mat_id, pt_id, i, j));
                        }
                    }
                    fprintf(fp, "\n");
                }
            }
            fprintf(fp, "        </DataArray>\n");
        }

        // strain
        if (mat_strain_id >= 0) {
            fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
                        "NumberOfComponents=\"9\" format=\"ascii\">\n",
                    mat_tensor_var_names[mat_strain_id].c_str());
            for (size_t elem = 0; elem < num_mat_elems; elem++) {
                for (size_t gp = 0; gp < num_gp_per_elem; gp++) {
                    const size_t pt_id =
                        State.points_in_mat_elem.host(elem, gp);
                    fprintf(fp, "         ");
                    for (size_t i = 0; i < 3; i++) {
                        for (size_t j = 0; j < 3; j++) {
                            fprintf(fp, " %.15e",
                                    State.MaterialPoints.strain.host(
                                        mat_id, pt_id, i, j));
                        }
                    }
                    fprintf(fp, "\n");
                }
            }
            fprintf(fp, "        </DataArray>\n");
        }

        fprintf(fp, "      </PointData>\n");

    } // end if any fields to write

    // -----------------------------------------------------------------------
    //  Footer
    // -----------------------------------------------------------------------

    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");

    fclose(fp);

} // end write_vtu_Pn


void MeshWriter::write_vtm_Pn(CArray<double>& graphics_times,
                  const std::string& elem_part_name,
                  const std::string& mat_part_name,
                  double time_value,
                  int graphics_id,
                  size_t num_mats_global,
                  bool write_mesh_state,
                  bool write_mat_pt_state,
                  const size_t solver_id,
                  int mpi_size,
                  const unsigned long long* owned_elems_by_rank,
                  const unsigned long long* mat_elem_counts_by_rank)
{
    for (int file_id = 0; file_id <= graphics_id; file_id++) {

        FILE* out[1];
        char  filename[512];
        int   max_len = static_cast<int>(sizeof filename);
        int   str_output_len;

        str_output_len =
            snprintf(filename, static_cast<size_t>(max_len), "vtk/data/Fierro.solver%zu.%05d.vtm", solver_id, file_id);

        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

        out[0] = fopen(filename, "w");
        if (!out[0]) {
            std::cerr << "[MeshWriter] Failed to open VTM file: " << filename << std::endl;
            continue;
        }

        fprintf(out[0], "<?xml version=\"1.0\"?>\n");
        fprintf(out[0], "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
        fprintf(out[0], "  <vtkMultiBlockDataSet>\n");

        size_t block_id = 0;

        // ---------------------------------------------------------------
        //  Block: element-averaged mesh (same schema as write_vtm)
        // ---------------------------------------------------------------
        if (write_mesh_state) {
            int mesh_pieces = 0;
            if (mpi_size > 1 && owned_elems_by_rank != nullptr) {
                for (int r = 0; r < mpi_size; r++) {
                    if (owned_elems_by_rank[static_cast<size_t>(r)] > 0ULL) {
                        mesh_pieces++;
                    }
                }
            }
            else if (owned_elems_by_rank != nullptr && owned_elems_by_rank[0] > 0ULL) {
                mesh_pieces = 1;
            }
            else if (owned_elems_by_rank == nullptr) {
                mesh_pieces = 1;
            }

            if (mesh_pieces > 0) {
                fprintf(out[0], "    <Block index=\"%zu\" name=\"Mesh\">\n", block_id);
                block_id++;

                int ds_index = 0;
                if (mpi_size > 1 && owned_elems_by_rank != nullptr) {
                    for (int r = 0; r < mpi_size; r++) {
                        if (owned_elems_by_rank[static_cast<size_t>(r)] == 0ULL) {
                            continue;
                        }
                        fprintf(out[0],
                                "      <DataSet index=\"%d\" name=\"Field_r%04d\" file=\"Fierro.solver%zu.%s.%05d_r%04d.vtu\" />\n",
                                ds_index,
                                r,
                                solver_id,
                                elem_part_name.c_str(),
                                file_id,
                                r);
                        ds_index++;
                    }
                }
                else {
                    fprintf(out[0],
                            "      <DataSet index=\"0\" name=\"Field\" file=\"Fierro.solver%zu.%s.%05d.vtu\" />\n",
                            solver_id,
                            elem_part_name.c_str(),
                            file_id);
                }

                fprintf(out[0], "    </Block>\n");
            }
        }

        // ---------------------------------------------------------------
        //  Block: per-material data
        //
        //  Each material that has data gets its own "MatN" sub-block
        //  (contiguous mat_block_idx, same reasoning as write_vtm to avoid
        //  duplicate sibling indices). Inside MatN, each rank that owns
        //  data for that material contributes TWO DataSet children:
        //    - "Elements"   -> matN.vtu           (element-averaged field)
        //    - "GaussPoints"-> matN_Pn.vtu         (Gauss-point cloud)
        //  Both are direct children of MatN with unique indices, so
        //  ParaView won't collapse them the way the old <Piece> layout did.
        // ---------------------------------------------------------------
        if (write_mat_pt_state && mat_elem_counts_by_rank != nullptr && num_mats_global > 0) {
            int mat_pieces = 0;
            for (size_t mat_id = 0; mat_id < num_mats_global; mat_id++) {
                if (mpi_size > 1) {
                    for (int r = 0; r < mpi_size; r++) {
                        if (mat_elem_counts_by_rank[static_cast<size_t>(r) * num_mats_global + mat_id] > 0ULL) {
                            mat_pieces++;
                        }
                    }
                }
                else if (mat_elem_counts_by_rank[mat_id] > 0ULL) {
                    mat_pieces++;
                }
            }

            if (mat_pieces > 0) {
                fprintf(out[0], "    <Block index=\"%zu\" name=\"Mat\">\n", block_id);
                size_t mat_block_idx = 0;

                for (size_t mat_id = 0; mat_id < num_mats_global; mat_id++) {
                    int pieces_this_mat = 0;
                    if (mpi_size > 1) {
                        for (int r = 0; r < mpi_size; r++) {
                            if (mat_elem_counts_by_rank[static_cast<size_t>(r) * num_mats_global + mat_id] > 0ULL) {
                                pieces_this_mat++;
                            }
                        }
                    }
                    else if (mat_elem_counts_by_rank[mat_id] > 0ULL) {
                        pieces_this_mat = 1;
                    }
                    if (pieces_this_mat == 0) {
                        continue;
                    }

                    fprintf(out[0], "      <Block index=\"%zu\" name=\"Mat%zu\">\n", mat_block_idx, mat_id);

                    int ds_index = 0;
                    if (mpi_size > 1) {
                        for (int r = 0; r < mpi_size; r++) {
                            const unsigned long long nm =
                                mat_elem_counts_by_rank[static_cast<size_t>(r) * num_mats_global + mat_id];
                            if (nm == 0ULL) {
                                continue;
                            }

                            // Elements (element-averaged material field)
                            fprintf(out[0],
                                    "        <DataSet index=\"%d\" name=\"Mat%zu_Elements_r%04d\" "
                                    "file=\"Fierro.solver%zu.%s%zu.%05d_r%04d.vtu\" />\n",
                                    ds_index,
                                    mat_id,
                                    r,
                                    solver_id,
                                    mat_part_name.c_str(),
                                    mat_id,
                                    file_id,
                                    r);
                            ds_index++;

                            // GaussPoints (Pn point cloud)
                            fprintf(out[0],
                                    "        <DataSet index=\"%d\" name=\"Mat%zu_GaussPoints_r%04d\" "
                                    "file=\"Fierro.solver%zu.%s%zu_Pn.%05d_r%04d.vtu\" />\n",
                                    ds_index,
                                    mat_id,
                                    r,
                                    solver_id,
                                    mat_part_name.c_str(),
                                    mat_id,
                                    file_id,
                                    r);
                            ds_index++;
                        }
                    }
                    else {
                        // Elements
                        fprintf(out[0],
                                "        <DataSet index=\"%d\" name=\"Mat%zu_Elements\" "
                                "file=\"Fierro.solver%zu.%s%zu.%05d.vtu\" />\n",
                                ds_index,
                                mat_id,
                                solver_id,
                                mat_part_name.c_str(),
                                mat_id,
                                file_id);
                        ds_index++;

                        // GaussPoints
                        fprintf(out[0],
                                "        <DataSet index=\"%d\" name=\"Mat%zu_GaussPoints\" "
                                "file=\"Fierro.solver%zu.%s%zu_Pn.%05d.vtu\" />\n",
                                ds_index,
                                mat_id,
                                solver_id,
                                mat_part_name.c_str(),
                                mat_id,
                                file_id);
                        ds_index++;
                    }

                    fprintf(out[0], "      </Block>\n");
                    mat_block_idx++;
                }

                fprintf(out[0], "    </Block>\n");
            }
        }

        fprintf(out[0], "  </vtkMultiBlockDataSet>\n");
        fprintf(out[0], "</VTKFile>");

        fclose(out[0]);

    } // end for file_id

} // end write_vtm_Pn


/////////////////////////////////////////////////////////////////////////////
///
/// \fn write_pvd
///
/// \brief Writes a pvd ASCII output file for the element and nodal fields
///
/// \param Vector of all graphics output times
/// \param element average field names
/// \param current time value
/// \param graphics index
///
/////////////////////////////////////////////////////////////////////////////
void MeshWriter::write_pvd(CArray<double>& graphics_times,
               double time_value,
               int graphics_id,
               const size_t solver_id,
               int mpi_rank)
{
    if (mpi_rank != 0) {
        return;
    }

    FILE* out[20];
    char  filename[512];
    int   max_len = static_cast<int>(sizeof filename);
    int   str_output_len;

    str_output_len = snprintf(filename, static_cast<size_t>(max_len), "vtk/Fierro.solver%zu.pvd", solver_id);

    if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

    out[0] = fopen(filename, "w");

    fprintf(out[0], "<?xml version=\"1.0\"?>\n");
    fprintf(out[0], "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\">\n");
    fprintf(out[0], "  <Collection>\n");

    for (int i = 0; i <= graphics_id; i++) {
        fprintf(out[0], "    <DataSet timestep=\"%.5e\" file=\"data/Fierro.solver%zu.%05d.vtm\" time=\"%.5e\" />\n", 
                                                 graphics_times(i), solver_id, i, graphics_times(i) );
        //fprintf(out[0], "    <DataSet timestep=\"%d\" file=\"data/Fierro.solver%zu.%05d.vtm\" time= \"%12.5e\" />\n", 
        //                                         i, solver_id, i, graphics_times(i) );
    }

    fprintf(out[0], "  </Collection>\n");
    fprintf(out[0], "</VTKFile>"); 

    fclose(out[0]);

} // end pvd


/////////////////////////////////////////////////////////////////////////////
///
/// \fn write_vtm
///
/// \brief Writes a vtm ASCII output file for all fields -- mesh and material
///
/// \param Vector of all graphics output times
/// \param element average field names
/// \param current time value
/// \param graphics index
///
/////////////////////////////////////////////////////////////////////////////
void MeshWriter::write_vtm(CArray<double>& graphics_times,
               const std::string& elem_part_name,
               const std::string& mat_part_name,
               double time_value,
               int graphics_id,
               size_t num_mats_global,
               bool write_mesh_state,
               bool write_mat_pt_state,
               const size_t solver_id,
               int mpi_size,
               const unsigned long long* owned_elems_by_rank,
               const unsigned long long* mat_elem_counts_by_rank)
{
    for (int file_id = 0; file_id <= graphics_id; file_id++) {

        FILE* out[20];
        char  filename[512];
        int   max_len = static_cast<int>(sizeof filename);
        int   str_output_len;

        str_output_len =
            snprintf(filename, static_cast<size_t>(max_len), "vtk/data/Fierro.solver%zu.%05d.vtm", solver_id, file_id);

        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

        out[0] = fopen(filename, "w");
        if (!out[0]) {
            std::cerr << "[MeshWriter] Failed to open VTM file: " << filename << std::endl;
            continue;
        }

        fprintf(out[0], "<?xml version=\"1.0\"?>\n");
        fprintf(out[0], "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
        fprintf(out[0], "  <vtkMultiBlockDataSet>\n");

        size_t block_id = 0;
        if (write_mesh_state) {
            int mesh_pieces = 0;
            if (mpi_size > 1 && owned_elems_by_rank != nullptr) {
                for (int r = 0; r < mpi_size; r++) {
                    if (owned_elems_by_rank[static_cast<size_t>(r)] > 0ULL) {
                        mesh_pieces++;
                    }
                }
            }
            else if (owned_elems_by_rank != nullptr && owned_elems_by_rank[0] > 0ULL) {
                mesh_pieces = 1;
            }
            else if (owned_elems_by_rank == nullptr) {
                mesh_pieces = 1;
            }

            if (mesh_pieces > 0) {
                // vtkMultiBlockDataSet schema: <DataSet> children must be direct
                // children of <Block>, not wrapped in <Piece>. Sibling index attrs
                // must be unique within the parent.
                fprintf(out[0], "    <Block index=\"%zu\" name=\"Mesh\">\n", block_id);
                {
                    block_id++;
                    int ds_index = 0;
                    if (mpi_size > 1 && owned_elems_by_rank != nullptr) {
                        for (int r = 0; r < mpi_size; r++) {
                            if (owned_elems_by_rank[static_cast<size_t>(r)] == 0ULL) {
                                continue;
                            }
                            fprintf(out[0],
                                    "      <DataSet index=\"%d\" name=\"Field_r%04d\" file=\"Fierro.solver%zu.%s.%05d_r%04d.vtu\" />\n",
                                    ds_index,
                                    r,
                                    solver_id,
                                    elem_part_name.c_str(),
                                    file_id,
                                    r);
                            ds_index++;
                        }
                    }
                    else {
                        fprintf(out[0],
                                "      <DataSet index=\"0\" name=\"Field\" file=\"Fierro.solver%zu.%s.%05d.vtu\" />\n",
                                solver_id,
                                elem_part_name.c_str(),
                                file_id);
                    }
                }
                fprintf(out[0], "    </Block>\n");
            }
        }

        if (write_mat_pt_state && mat_elem_counts_by_rank != nullptr && num_mats_global > 0) {
            int mat_pieces = 0;
            for (size_t mat_id = 0; mat_id < num_mats_global; mat_id++) {
                if (mpi_size > 1) {
                    for (int r = 0; r < mpi_size; r++) {
                        if (mat_elem_counts_by_rank[static_cast<size_t>(r) * num_mats_global + mat_id] > 0ULL) {
                            mat_pieces++;
                        }
                    }
                }
                else if (mat_elem_counts_by_rank[mat_id] > 0ULL) {
                    mat_pieces++;
                }
            }
            if (mat_pieces > 0) {
                // Nested layout: <Block name="Mat"> contains one <Block name="MatN"> per
                // material that has data, each with per-rank <DataSet> children. This is
                // required because the prior flat layout reset piece_index to 0 per
                // material, producing duplicate sibling indices (e.g. Mat0_r0000 and
                // Mat1_r0001 both at index 0) which ParaView silently collapses, dropping
                // the first piece. mat_block_idx is a contiguous counter so sibling
                // indices in the outer Mat block remain unique and dense.
                fprintf(out[0], "    <Block index=\"%zu\" name=\"Mat\">\n", block_id);
                size_t mat_block_idx = 0;
                for (size_t mat_id = 0; mat_id < num_mats_global; mat_id++) {
                    int pieces_this_mat = 0;
                    if (mpi_size > 1) {
                        for (int r = 0; r < mpi_size; r++) {
                            if (mat_elem_counts_by_rank[static_cast<size_t>(r) * num_mats_global + mat_id] > 0ULL) {
                                pieces_this_mat++;
                            }
                        }
                    }
                    else if (mat_elem_counts_by_rank[mat_id] > 0ULL) {
                        pieces_this_mat = 1;
                    }
                    if (pieces_this_mat == 0) {
                        continue;
                    }

                    fprintf(out[0], "      <Block index=\"%zu\" name=\"Mat%zu\">\n", mat_block_idx, mat_id);
                    int ds_index = 0;
                    if (mpi_size > 1) {
                        for (int r = 0; r < mpi_size; r++) {
                            const unsigned long long nm =
                                mat_elem_counts_by_rank[static_cast<size_t>(r) * num_mats_global + mat_id];
                            if (nm == 0ULL) {
                                continue;
                            }
                            fprintf(out[0],
                                    "        <DataSet index=\"%d\" name=\"Mat%zu_r%04d\" "
                                    "file=\"Fierro.solver%zu.%s%zu.%05d_r%04d.vtu\" />\n",
                                    ds_index,
                                    mat_id,
                                    r,
                                    solver_id,
                                    mat_part_name.c_str(),
                                    mat_id,
                                    file_id,
                                    r);
                            ds_index++;
                        }
                    }
                    else {
                        fprintf(out[0],
                                "        <DataSet index=\"0\" name=\"Mat%zu\" "
                                "file=\"Fierro.solver%zu.%s%zu.%05d.vtu\" />\n",
                                mat_id,
                                solver_id,
                                mat_part_name.c_str(),
                                mat_id,
                                file_id);
                    }
                    fprintf(out[0], "      </Block>\n");
                    mat_block_idx++;
                }
                fprintf(out[0], "    </Block>\n");
            }
        }

        fprintf(out[0], "  </vtkMultiBlockDataSet>\n");
        fprintf(out[0], "</VTKFile>");

        fclose(out[0]);

    } // end for file_id

} // end vtm


/////////////////////////////////////////////////////////////////////////////
///
/// \fn build_material_elem_node_lists
///
/// \brief Creates elems and nodes for a unique mesh of a material (i.e, a part)
///
/// \param Simulation mesh
/// \param State node data
/// \param Material node coordinates
/// \param Material nodes in the material element
/// \param Material to mesh map for elements
/// \param number of material nodes
/// \param number of material elements
/// \param number of nodes in the element
/// \param number of dimensions 
///
/////////////////////////////////////////////////////////////////////////////
void MeshWriter::build_material_elem_node_lists(
    const swage::Mesh_t& mesh,
    const MPICArrayKokkos<double>& state_node_coords,
    DCArrayKokkos<double>& mat_node_coords,
    DCArrayKokkos <size_t>& mat_nodes_in_mat_elem,
    const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
    const size_t mat_id,
    size_t& num_mat_nodes,
    const size_t num_mat_elems,
    const size_t num_nodes_in_elem,
    const size_t num_dims)
{

    // helper arrays
    DCArrayKokkos <size_t> dummy_counter(mesh.num_nodes, "dummy_counter");
    DCArrayKokkos <size_t> access_mat_node_gids(mesh.num_nodes, "access_mat_node_gids");
    dummy_counter.set_values(0);

    // tag and count the number of nodes in this part
    FOR_ALL (mat_elem_sid, 0, num_mat_elems, {
        // get elem gid
        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);  // WARNING not GPU compatible
        
        // parallel loop over the nodes in the element
        for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++) {
            size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

            Kokkos::atomic_add(&dummy_counter(node_gid), 1); // values in part will be >0

        } // end for nodes in element
        
    }); // end parallel for
    Kokkos::fence();
    dummy_counter.update_host();

    // loop opperation is not thread safe, must be run serially
    size_t mat_node_gid = 0;
    for(size_t node_gid = 0; node_gid<mesh.num_nodes; node_gid ++) {
        
        // save the nodes on the part (i.e., that belong to the material)
        if (dummy_counter.host(node_gid)>0){
            mat_node_coords.host(mat_node_gid, 0) = state_node_coords.host(node_gid, 0);
            mat_node_coords.host(mat_node_gid, 1) = state_node_coords.host(node_gid, 1);
            if (num_dims == 3){ 
                mat_node_coords.host(mat_node_gid, 2) = state_node_coords.host(node_gid, 2);
            } // end if on dims

            access_mat_node_gids.host(node_gid) = mat_node_gid; // the part node id

            mat_node_gid ++;

            dummy_counter.host(node_gid) = 0; // set counter to zero, it was accounted for
        } // end if this node is on the part

    } // end loop over all mesh nodes
    mat_node_coords.update_device();
    access_mat_node_gids.update_device();
    dummy_counter.update_device();
    Kokkos::fence();

    // save the number of nodes defining the material region, i.e., the part
    num_mat_nodes = mat_node_gid;
    
    // save the new node id's
    FOR_ALL (mat_elem_sid, 0, num_mat_elems, {
        // get elem gid
        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);
        
        // parallel loop over the nodes in the element
        for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++) {
            size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

            // save the mat_node to the mat elem list
            mat_nodes_in_mat_elem(mat_elem_sid, node_lid) = access_mat_node_gids(node_gid);

        } // end for nodes in element
        
    }); // end parallel for
    Kokkos::fence();
    mat_nodes_in_mat_elem.update_host();

} // end build part (i.e., material elem and point lists) function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn write_material_point_state
///
/// \brief Writes a state output file at each material point
///
/// \param Simulation mesh
/// \param State data
/// \param Simulation parameters
/// \param current time value
/// \param Vector of all graphics output times
///
/////////////////////////////////////////////////////////////////////////////
void MeshWriter::write_material_point_state(swage::Mesh_t& mesh,
    State_t& State,
    SimulationParameters_t& SimulationParamaters,
    double time_value,
    CArray<double> graphics_times,
    std::vector<node_state> node_states,
    std::vector<gauss_pt_state> gauss_pt_states,
    std::vector<material_pt_state> material_pt_states)
{
    // WARNING WARNING WARNING:
    // This currently assumes the gauss and material point IDs are the same as the element ID
    // This will need to be updated for high order methods

    int mpi_rank = 0;
    int mpi_size = 1;
    mesh_io_mpi_detail::query_world_rank_size(mpi_rank, mpi_size);

    // Update host data
    // ---- Update host data ----
    size_t num_mats = State.MaterialPoints.num_material_points.size();

    State.MaterialPoints.den.update_host();
    State.MaterialPoints.pres.update_host();
    State.MaterialPoints.stress.update_host();
    State.MaterialPoints.strain.update_host();
    State.MaterialPoints.sspd.update_host();
    State.MaterialPoints.sie.update_host();
    State.MaterialPoints.mass.update_host();

    State.GaussPoints.vol.update_host();
    State.GaussPoints.shock_detector.update_host();

    State.node.coords.update_host();
    State.node.vel.update_host();
    State.node.mass.update_host();

    Kokkos::fence();

    struct stat st;
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        if (stat("state", &st) != 0) {
            system("mkdir state");
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    size_t num_dims = mesh.num_dims;

    std::vector<double> local_mat_pt_data;
    std::vector<double> local_node_data;

    // write out values for the elem
    for (size_t mat_id = 0; mat_id < num_mats; mat_id++) {

        size_t num_mat_elems = State.MaterialToMeshMaps.num_mat_elems.host(mat_id);

        for (size_t mat_elem_sid = 0; mat_elem_sid < num_mat_elems; mat_elem_sid++)
        {

            const size_t elem_gid = State.MaterialToMeshMaps.elem_in_mat_elem.host(mat_id, mat_elem_sid);

            // Skip if this element is not owned by the current rank
            if (elem_gid >= mesh.num_owned_elems) {
                continue;
            }

            double elem_coords[3];
            elem_coords[0] = 0.0;
            elem_coords[1] = 0.0;
            elem_coords[2] = 0.0;

            // get the coordinates of the element center
            for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {

                elem_coords[0] += State.node.coords.host(mesh.nodes_in_elem.host(elem_gid, node_lid), 0);
                elem_coords[1] += State.node.coords.host(mesh.nodes_in_elem.host(elem_gid, node_lid), 1);
                if (num_dims == 3) {
                    elem_coords[2] += State.node.coords.host(mesh.nodes_in_elem.host(elem_gid, node_lid), 2);
                }
                else{
                    elem_coords[2] = 0.0;
                }
            } // end loop over nodes in element

            elem_coords[0] = elem_coords[0] / ((double)mesh.num_nodes_in_elem);
            elem_coords[1] = elem_coords[1] / ((double)mesh.num_nodes_in_elem);
            elem_coords[2] = elem_coords[2] / ((double)mesh.num_nodes_in_elem);

            double rad2 = sqrt(elem_coords[0] * elem_coords[0] +
                               elem_coords[1] * elem_coords[1]);

            double rad3 = sqrt(elem_coords[0] * elem_coords[0] +
                               elem_coords[1] * elem_coords[1] +
                               elem_coords[2] * elem_coords[2]);

            local_mat_pt_data.push_back(elem_coords[0]);
            local_mat_pt_data.push_back(elem_coords[1]);
            local_mat_pt_data.push_back(elem_coords[2]);
            local_mat_pt_data.push_back(rad2);
            local_mat_pt_data.push_back(rad3);
            local_mat_pt_data.push_back(State.MaterialPoints.den.host(mat_id, mat_elem_sid));
            local_mat_pt_data.push_back(State.MaterialPoints.pres.host(mat_id, mat_elem_sid));
            local_mat_pt_data.push_back(State.MaterialPoints.sie.host(mat_id, mat_elem_sid));
            local_mat_pt_data.push_back(State.MaterialPoints.sspd.host(mat_id, mat_elem_sid));
            local_mat_pt_data.push_back(State.GaussPoints.vol.host(elem_gid));
            local_mat_pt_data.push_back(State.MaterialPoints.mass.host(mat_id, mat_elem_sid));

        } // end for elements

    } // end for materials

    // get the coordinates of the node
    for (size_t node_gid = 0; node_gid < mesh.num_owned_nodes; node_gid++) {

        double node_coords[3];

        node_coords[0] = State.node.coords.host(node_gid, 0);
        node_coords[1] = State.node.coords.host(node_gid, 1);
        if (num_dims == 3) {
            node_coords[2] = State.node.coords.host(node_gid, 2);
        }
        else{
            node_coords[2] = 0.0;
        }

        double rad2 = sqrt(node_coords[0] * node_coords[0] +
                           node_coords[1] * node_coords[1]);
        double rad3 = sqrt(node_coords[0] * node_coords[0] +
                           node_coords[1] * node_coords[1] +
                           node_coords[2] * node_coords[2]);

        double node_vel[3];

       node_vel[0] = State.node.vel.host(node_gid, 0);
       node_vel[1] = State.node.vel.host(node_gid, 1);
        if (num_dims == 3) {
           node_vel[2] = State.node.vel.host(node_gid, 2);
        }
        else{
           node_vel[2] = 0.0;
        }

        double speed = sqrt(node_vel[0] * node_vel[0] +
                            node_vel[1] * node_vel[1] +
                            node_vel[2] * node_vel[2]);



        // looking at perfect radial motion
        double unit_r_vec[2];
        unit_r_vec[0] = node_coords[0]/rad2;
        unit_r_vec[1] = node_coords[1]/rad2;

        //the radial motion
        double v_dot_r = node_vel[0] * unit_r_vec[0] +
                         node_vel[1] * unit_r_vec[1];
        

        double err_v_dot_r[3]; 
        err_v_dot_r[0] = node_vel[0]-unit_r_vec[0]*v_dot_r;
        err_v_dot_r[1] = node_vel[1]-unit_r_vec[1]*v_dot_r;

        double mag_err_v_dot_r = sqrt(err_v_dot_r[0]*err_v_dot_r[0] + err_v_dot_r[1]*err_v_dot_r[1]);

        local_node_data.push_back(node_coords[0]);
        local_node_data.push_back(node_coords[1]);
        local_node_data.push_back(node_coords[2]);
        local_node_data.push_back(rad2);
        local_node_data.push_back(rad3);
        local_node_data.push_back(node_vel[0]);
        local_node_data.push_back(node_vel[1]);
        local_node_data.push_back(node_vel[2]);
        local_node_data.push_back(speed);
        local_node_data.push_back(mag_err_v_dot_r);


    } // end loop over nodes in element

    int local_mat_pt_rows = static_cast<int>(local_mat_pt_data.size() / 11);
    int local_node_rows = static_cast<int>(local_node_data.size() / 10);

    std::vector<int> mat_pt_row_counts;
    std::vector<int> mat_pt_row_displs;
    std::vector<int> node_row_counts;
    std::vector<int> node_row_displs;
    if (mpi_rank == 0) {
        mat_pt_row_counts.assign(mpi_size, 0);
        mat_pt_row_displs.assign(mpi_size, 0);
        node_row_counts.assign(mpi_size, 0);
        node_row_displs.assign(mpi_size, 0);
    }

    MPI_Gather(&local_mat_pt_rows, 1, MPI_INT,
               mpi_rank == 0 ? mat_pt_row_counts.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_node_rows, 1, MPI_INT,
               mpi_rank == 0 ? node_row_counts.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int total_mat_pt_rows = 0;
    int total_node_rows = 0;
    std::vector<int> mat_pt_counts;
    std::vector<int> mat_pt_displs;
    std::vector<int> node_counts;
    std::vector<int> node_displs;
    std::vector<double> gathered_mat_pt_data;
    std::vector<double> gathered_node_data;
    if (mpi_rank == 0) {
        mat_pt_counts.assign(mpi_size, 0);
        mat_pt_displs.assign(mpi_size, 0);
        node_counts.assign(mpi_size, 0);
        node_displs.assign(mpi_size, 0);

        int mat_disp = 0;
        int node_disp = 0;
        for (int r = 0; r < mpi_size; r++) {
            mat_pt_row_displs[r] = mat_disp;
            mat_disp += mat_pt_row_counts[r];
            total_mat_pt_rows += mat_pt_row_counts[r];

            node_row_displs[r] = node_disp;
            node_disp += node_row_counts[r];
            total_node_rows += node_row_counts[r];
        }

        for (int r = 0; r < mpi_size; r++) {
            mat_pt_counts[r] = mat_pt_row_counts[r] * 11;
            mat_pt_displs[r] = mat_pt_row_displs[r] * 11;
            node_counts[r] = node_row_counts[r] * 10;
            node_displs[r] = node_row_displs[r] * 10;
        }

        gathered_mat_pt_data.assign(total_mat_pt_rows * 11, 0.0);
        gathered_node_data.assign(total_node_rows * 10, 0.0);
    }

    MPI_Gatherv(local_mat_pt_data.empty() ? nullptr : local_mat_pt_data.data(),
                static_cast<int>(local_mat_pt_data.size()), MPI_DOUBLE,
                mpi_rank == 0 ? gathered_mat_pt_data.data() : nullptr,
                mpi_rank == 0 ? mat_pt_counts.data() : nullptr,
                mpi_rank == 0 ? mat_pt_displs.data() : nullptr,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(local_node_data.empty() ? nullptr : local_node_data.data(),
                static_cast<int>(local_node_data.size()), MPI_DOUBLE,
                mpi_rank == 0 ? gathered_node_data.data() : nullptr,
                mpi_rank == 0 ? node_counts.data() : nullptr,
                mpi_rank == 0 ? node_displs.data() : nullptr,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {

    //  ---------------------------------------------------------------------------
    //  Setup of file and directory for exporting
    //  ---------------------------------------------------------------------------

    // output file
    FILE* out_elem_state;  // element average state
    char  filename[128];

    int max_len = sizeof filename;

    //snprintf(filename, max_len, "state/mat_pt_state_rank_%d_t_%6.4e.txt", mpi_rank, time_value);
    snprintf(filename, max_len, "state/mat_pt_state_t_%6.4e.txt", time_value);

    // output files
    out_elem_state = fopen(filename, "w");

    // write state dump
    fprintf(out_elem_state, "# state dump file\n");
    fprintf(out_elem_state, "# x  y  z  radius_2D  radius_3D  den  pres  sie  sspd  vol  mass \n");

    for (int row = 0; row < total_mat_pt_rows; row++) {
        size_t i = static_cast<size_t>(row) * 11;
        fprintf(out_elem_state, "%4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t \n",
                 gathered_mat_pt_data[i],
                 gathered_mat_pt_data[i + 1],
                 gathered_mat_pt_data[i + 2],
                 gathered_mat_pt_data[i + 3],
                 gathered_mat_pt_data[i + 4],
                 gathered_mat_pt_data[i + 5],
                 gathered_mat_pt_data[i + 6],
                 gathered_mat_pt_data[i + 7],
                 gathered_mat_pt_data[i + 8],
                 gathered_mat_pt_data[i + 9],
                 gathered_mat_pt_data[i + 10]);
    }
    fclose(out_elem_state);



    // printing nodal state
        
    FILE* out_point_state;  // element average state

    snprintf(filename, max_len, "state/node_state_t_%6.4e.txt", time_value);

    // output files
    out_point_state = fopen(filename, "w");

    // write state dump
    fprintf(out_point_state, "# state node dump file\n");
    fprintf(out_point_state, "# x  y  z  radius_2D  radius_3D  vel_x  vel_y  vel_z  speed  ||err_v_dot_r|| \n");

    for (int row = 0; row < total_node_rows; row++) {
        size_t i = static_cast<size_t>(row) * 10;
        fprintf(out_point_state, "%4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t  %4.12e\t %4.12e\t \n",
                     gathered_node_data[i],
                     gathered_node_data[i + 1],
                     gathered_node_data[i + 2],
                     gathered_node_data[i + 3],
                     gathered_node_data[i + 4],
                     gathered_node_data[i + 5],
                     gathered_node_data[i + 6],
                     gathered_node_data[i + 7],
                     gathered_node_data[i + 8],
                     gathered_node_data[i + 9]);
    }


    fclose(out_point_state);

    } // mpi_rank == 0

    return;
} // end of state output


void MeshWriter::write_text_state(
    const swage::Mesh_t&                   mesh,
    const State_t&                       State,
    const std::vector<std::string>&      mat_scalar_var_names,
    const std::vector<std::string>&      mat_tensor_var_names,
    const std::vector<std::string>&      mat_vector_var_names,
    const std::vector<std::string>&      node_scalar_var_names,
    const std::vector<std::string>&      node_vector_var_names,
    const size_t                         num_nodes_in_elem,
    const size_t                         num_dims,
    const double                         time_value,
    // material-point field slot IDs (-1 means "not requested")
    const int mat_den_id,
    const int mat_pres_id,
    const int mat_sie_id,
    const int mat_sspd_id,
    const int mat_mass_id,
    const int mat_mat_volfrac_id,
    const int mat_geo_volfrac_id,
    const int mat_eroded_id,
    const int mat_stress_id,
    const int mat_strain_id,
    const int mat_conductivity_id,
    const int mat_specific_heat_id,
    const int mat_heat_flux_id,
    // node field slot IDs (-1 means "not requested")
    const int node_mass_id,
    const int node_vel_id,
    const int node_coord_id,
    const int node_temp_id,
    const int node_grad_level_set_id,
    const int node_disp_id
)
{
    // -----------------------------------------------------------------------
    //  Ensure the output directory exists before opening any files.
    // -----------------------------------------------------------------------
    struct stat st;
    MPI_Barrier(MPI_COMM_WORLD);
    if (stat("state", &st) != 0) {
        system("mkdir state");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    int rank;
    int world_size;
    mesh_io_mpi_detail::query_world_rank_size(rank,world_size);

    const size_t num_mats = State.MaterialPoints.num_material_points.size();

    // sizing owned mat elems if not sized already
    if (num_owned_mat_elems.size() <= 0) {
        num_owned_mat_elems = CArray <long long int> (num_mats);
        for (int mat = 0; mat < num_mats; mat++) {
            num_owned_mat_elems(mat) = -1;
        }
    }

    // -----------------------------------------------------------------------
    //  FILE 1 – Material-point state
    //
    //  One file per material.  One row per element (single integration
    //  point assumed); the element coordinate is the centroid, computed
    //  as the simple average of its node positions.
    //
    //  Filename:  state/mat_pt_state_t_<time>_mat_id_<id>.txt
    // -----------------------------------------------------------------------
    for (int mat_id = 0; mat_id < (int)State.MaterialToMeshMaps.num_mat_elems.dims(0); mat_id++) {

        size_t num_mat_elems = State.MaterialToMeshMaps.num_mat_elems.host(mat_id);

        // populating the counts for num_owned_mat_elems if not already populated
        if (num_owned_mat_elems(mat_id) == -1) {
            int loc_tally = 0;
            int tally = 0;
            FOR_REDUCE_SUM(mat_elem, 0, num_mat_elems, loc_tally, {
                const size_t elem_id = State.MaterialToMeshMaps.elem_in_mat_elem(mat_id, mat_elem);
                if (elem_id < mesh.num_owned_elems) {
                    loc_tally += 1;
                }
            }, tally);
            num_owned_mat_elems(mat_id) = tally;
        }

        // ------------------------------------------------------------------
        //  Element centroid coordinates
        // ------------------------------------------------------------------
        DCArrayKokkos<double> x_phys(num_owned_mat_elems(mat_id), 3);
        x_phys.set_values(0);

        FOR_ALL(elem, 0, num_owned_mat_elems(mat_id), {
            const size_t elem_id =
                State.MaterialToMeshMaps.elem_in_mat_elem(mat_id, elem);

            for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                const size_t node_gid = mesh.nodes_in_elem(elem_id, node_lid);
                for (size_t dim = 0; dim < num_dims; dim++) {
                    x_phys(elem, dim) += State.node.coords(node_gid, dim);
                }
            }
            for (size_t dim = 0; dim < num_dims; dim++) {
                x_phys(elem, dim) /= (double)num_nodes_in_elem;
            }
        });
        x_phys.update_host();

        FILE* out_elem_state;
        char  filename[128];
        int   max_len = sizeof filename;
        if (world_size == 1) {
            snprintf(filename, max_len, "state/mat_pt_state_t_%6.4e_mat_id_%d.txt", time_value, mat_id);
        } else {
            snprintf(filename, max_len, "state/mat_pt_state_t_%6.4e_mat_id_%d_rank_%d.txt", time_value, mat_id, rank);
        }

        out_elem_state = fopen(filename, "w");
        if (!out_elem_state) {
            std::cerr << "write_text_state: could not open "
                      << filename << std::endl;
            return;
        }

        // ---- Metadata comment -------------------------------------------
        fprintf(out_elem_state,
                "# Material-point state mat_id=%d  time=%.10e\n",
                mat_id, time_value);
        fprintf(out_elem_state,
                "# num_mat_elems=%zu  num_gp_per_elem=1\n",
                (size_t)num_owned_mat_elems(mat_id));

        // ---- Column header line -----------------------------------------
        //
        // No 'gp' column: with a single integration point per element
        // it carries no information.
        // -----------------------------------------------------------------
        fprintf(out_elem_state,
                "# %-12s %-8s  %-22s %-22s %-22s %-22s %-22s",
                "elem_gid", "mat_elem_id", "x", "y", "z", "radius_2D", "radius_3D");

        if (mat_den_id           >= 0 && State.MaterialPoints.den.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_den_id].c_str());
        if (mat_pres_id          >= 0 && State.MaterialPoints.pres.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_pres_id].c_str());
        if (mat_sie_id           >= 0 && State.MaterialPoints.sie.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_sie_id].c_str());
        if (mat_sspd_id          >= 0 && State.MaterialPoints.sspd.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_sspd_id].c_str());
        if (mat_mass_id          >= 0 && State.MaterialPoints.mass.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_mass_id].c_str());
        if (mat_mat_volfrac_id   >= 0 && State.MaterialPoints.mat_volfrac.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_mat_volfrac_id].c_str());
        if (mat_geo_volfrac_id   >= 0 && State.MaterialPoints.geo_volfrac.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_geo_volfrac_id].c_str());
        if (mat_eroded_id        >= 0 && State.MaterialPoints.eroded.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_eroded_id].c_str());
        if (mat_conductivity_id  >= 0 && State.MaterialPoints.conductivity.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_conductivity_id].c_str());
        if (mat_specific_heat_id >= 0 && State.MaterialPoints.specific_heat.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_specific_heat_id].c_str());
        if (mat_heat_flux_id >= 0 && State.MaterialPoints.q_flux.size() > 0) {
            fprintf(out_elem_state, "  %-22s", 
                    (mat_vector_var_names[mat_heat_flux_id] + "_x").c_str());

            fprintf(out_elem_state, "  %-22s", 
                    (mat_vector_var_names[mat_heat_flux_id] + "_y").c_str());

            fprintf(out_elem_state, "  %-22s", 
                    (mat_vector_var_names[mat_heat_flux_id] + "_z").c_str());
            }

        // Tensor headers: <name>_ij  (i, j in {x, y, z}), row-major
        const char* comp[3] = {"x", "y", "z"};
        if (mat_stress_id >= 0 && State.MaterialPoints.stress.size() > 0) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    char col[40];
                    snprintf(col, sizeof col, "%s_%s%s",
                             mat_tensor_var_names[mat_stress_id].c_str(),
                             comp[i], comp[j]);
                    fprintf(out_elem_state, "  %-22s", col);
                }
            }
        }
        if (mat_strain_id >= 0 && State.MaterialPoints.strain.size() > 0) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    char col[40];
                    snprintf(col, sizeof col, "%s_%s%s",
                             mat_tensor_var_names[mat_strain_id].c_str(),
                             comp[i], comp[j]);
                    fprintf(out_elem_state, "  %-22s", col);
                }
            }
        }
        fprintf(out_elem_state, "\n");

        // ---- Data rows --------------------------------------------------
        //
        // One row per element.  pt_id uses gp=0 since there is exactly
        // one integration point per element in this formulation.
        // -----------------------------------------------------------------
        for (size_t elem = 0; elem < (size_t)num_owned_mat_elems(mat_id); elem++) {
            const size_t elem_rid =
                State.MaterialToMeshMaps.elem_in_mat_elem.host(mat_id, elem);
            size_t elem_gid;
            if (world_size > 1) {
                elem_gid = mesh.local_to_global_elem_mapping.host(elem_rid);
            }
            else {
                elem_gid = elem_rid;
            }
            const size_t pt_id = State.points_in_mat_elem.host(elem, 0);

            fprintf(out_elem_state,
                    "  %-12zu %-8zu  %22.14e %22.14e %22.14e %22.14e %22.14e",
                    elem_gid, elem,
                    x_phys.host(elem, 0),
                    x_phys.host(elem, 1),
                    (num_dims == 3) ? x_phys.host(elem, 2) : 0.0,
                    sqrt(x_phys.host(elem, 0)*x_phys.host(elem, 0)+x_phys.host(elem, 1)*x_phys.host(elem, 1)),
                    (num_dims == 3) ? sqrt(x_phys.host(elem, 0)*x_phys.host(elem, 0)+x_phys.host(elem, 1)*x_phys.host(elem, 1)+x_phys.host(elem, 2)*x_phys.host(elem, 2)) : 0.0);

            // Scalar fields
            if (mat_den_id           >= 0 && State.MaterialPoints.den.size() > 0)
                fprintf(out_elem_state, "  %22.14e",
                        State.MaterialPoints.den.host(mat_id, pt_id));
            if (mat_pres_id          >= 0 && State.MaterialPoints.pres.size() > 0)
                fprintf(out_elem_state, "  %22.14e",
                        State.MaterialPoints.pres.host(mat_id, pt_id));
            if (mat_sie_id           >= 0 && State.MaterialPoints.sie.size() > 0)
                fprintf(out_elem_state, "  %22.14e",
                        State.MaterialPoints.sie.host(mat_id, pt_id));
            if (mat_sspd_id          >= 0 && State.MaterialPoints.sspd.size() > 0)
                fprintf(out_elem_state, "  %22.14e",
                        State.MaterialPoints.sspd.host(mat_id, pt_id));
            if (mat_mass_id          >= 0 && State.MaterialPoints.mass.size() > 0)
                fprintf(out_elem_state, "  %22.14e",
                        State.MaterialPoints.mass.host(mat_id, pt_id));
            if (mat_mat_volfrac_id   >= 0 && State.MaterialPoints.mat_volfrac.size() > 0)
                fprintf(out_elem_state, "  %22.14e",
                        State.MaterialPoints.mat_volfrac.host(mat_id, pt_id));
            if (mat_geo_volfrac_id   >= 0 && State.MaterialPoints.geo_volfrac.size() > 0)
                fprintf(out_elem_state, "  %22.14e",
                        State.MaterialPoints.geo_volfrac.host(mat_id, pt_id));
            if (mat_eroded_id        >= 0 && State.MaterialPoints.eroded.size() > 0)
                fprintf(out_elem_state, "  %22.14e",
                        static_cast<double>(
                            State.MaterialPoints.eroded.host(mat_id, pt_id)));
            if (mat_conductivity_id  >= 0 && State.MaterialPoints.conductivity.size() > 0)
                fprintf(out_elem_state, "  %22.14e",
                        State.MaterialPoints.conductivity.host(mat_id, pt_id));
            if (mat_specific_heat_id >= 0 && State.MaterialPoints.specific_heat.size() > 0)
                fprintf(out_elem_state, "  %22.14e",
                        State.MaterialPoints.specific_heat.host(mat_id, pt_id));
            if (mat_heat_flux_id >= 0 && State.MaterialPoints.q_flux.size() > 0)
                for (size_t i = 0; i < 3; i++)
                    fprintf(out_elem_state, "  %22.14e",
                            State.MaterialPoints.q_flux.host(
                                mat_id, pt_id, i));

            // Tensor fields: all 9 components, row-major (Txx Txy Txz ...)
            if (mat_stress_id >= 0 && State.MaterialPoints.stress.size() > 0) {
                for (size_t i = 0; i < 3; i++)
                    for (size_t j = 0; j < 3; j++)
                        fprintf(out_elem_state, "  %22.14e",
                                State.MaterialPoints.stress.host(
                                    mat_id, pt_id, i, j));
            }
            if (mat_strain_id >= 0 && State.MaterialPoints.strain.size() > 0) {
                for (size_t i = 0; i < 3; i++)
                    for (size_t j = 0; j < 3; j++)
                        fprintf(out_elem_state, "  %22.14e",
                                State.MaterialPoints.strain.host(
                                    mat_id, pt_id, i, j));
            }

            fprintf(out_elem_state, "\n");

        } // end elem

        fclose(out_elem_state);

    } // end mat_id loop

    // -----------------------------------------------------------------------
    //  FILE 2 – Node state  (identical to write_text_state_Pn)
    //
    //  Filename:  state/node_state_t_<time>.txt
    // -----------------------------------------------------------------------
    {
        FILE* out_point_state;
        char  filename[128];
        int   max_len = sizeof filename;
        if (world_size == 1) {
            snprintf(filename, max_len, "state/node_state_t_%6.4e.txt", time_value);
        } else {
            snprintf(filename, max_len, "state/node_state_t_%6.4e_rank_%d.txt", time_value, rank);
        }

        out_point_state = fopen(filename, "w");
        if (!out_point_state) {
            std::cerr << "write_text_state: could not open "
                      << filename << std::endl;
            return;
        }

        const size_t num_owned_nodes = mesh.num_owned_nodes;

        // ---- Metadata comment -------------------------------------------
        fprintf(out_point_state,
                "# Node state  time=%.10e\n", time_value);
        fprintf(out_point_state,
                "# num_owned_nodes=%zu  num_dims=%zu\n",
                num_owned_nodes, num_dims);

        // ---- Column header line -----------------------------------------
        fprintf(out_point_state, "# %-12s", "node_id");

        auto hdr_vec = [&out_point_state](const std::string& name) {
            char col[48];
            snprintf(col, sizeof col, "%s_x", name.c_str());
            fprintf(out_point_state, "  %-22s", col);
            snprintf(col, sizeof col, "%s_y", name.c_str());
            fprintf(out_point_state, "  %-22s", col);
            snprintf(col, sizeof col, "%s_z", name.c_str());
            fprintf(out_point_state, "  %-22s", col);
            snprintf(col, sizeof col, "|%s|", name.c_str());
            fprintf(out_point_state, "  %-22s", col);
        };

        if (node_coord_id          >= 0)
            hdr_vec(node_vector_var_names[node_coord_id]);
        if (node_vel_id            >= 0 && State.node.vel.size() > 0)
            hdr_vec(node_vector_var_names[node_vel_id]);
        if (node_disp_id           >= 0 && State.node.displacement.size() > 0)
            hdr_vec(node_vector_var_names[node_disp_id]);
        if (node_grad_level_set_id >= 0 && State.node.gradient_level_set.size() > 0)
            hdr_vec(node_vector_var_names[node_grad_level_set_id]);
        if (node_mass_id           >= 0 && State.node.mass.size() > 0)
            fprintf(out_point_state, "  %-22s",
                    node_scalar_var_names[node_mass_id].c_str());
        if (node_temp_id           >= 0 && State.node.temp.size() > 0)
            fprintf(out_point_state, "  %-22s",
                    node_scalar_var_names[node_temp_id].c_str());

        fprintf(out_point_state, "\n");

        // ---- Data rows --------------------------------------------------
        auto dat_vec = [&out_point_state](double vx, double vy, double vz_raw, size_t num_dims) {
            const double vz  = (num_dims == 3) ? vz_raw : 0.0;
            const double mag = sqrt(vx*vx + vy*vy + vz*vz);
            fprintf(out_point_state,
                    "  %22.14e  %22.14e  %22.14e  %22.14e",
                    vx, vy, vz, mag);
        };

        for (size_t node_gid = 0; node_gid < num_owned_nodes; node_gid++) {

            fprintf(out_point_state, "  %-12zu", node_gid);

            if (node_coord_id >= 0)
                dat_vec(State.node.coords.host(node_gid, 0),
                        State.node.coords.host(node_gid, 1),
                        State.node.coords.host(node_gid, 2),
                        mesh.num_dims);

            if (node_vel_id >= 0 && State.node.vel.size() > 0)
                dat_vec(State.node.vel.host(node_gid, 0),
                        State.node.vel.host(node_gid, 1),
                        State.node.vel.host(node_gid, 2),
                        mesh.num_dims);

            if (node_disp_id >= 0 && State.node.displacement.size() > 0)
                dat_vec(State.node.displacement.host(node_gid, 0),
                        State.node.displacement.host(node_gid, 1),
                        State.node.displacement.host(node_gid, 2),
                        mesh.num_dims);

            if (node_grad_level_set_id >= 0 && State.node.gradient_level_set.size() > 0)
                dat_vec(State.node.gradient_level_set.host(node_gid, 0),
                        State.node.gradient_level_set.host(node_gid, 1),
                        State.node.gradient_level_set.host(node_gid, 2),
                        mesh.num_dims);

            if (node_mass_id >= 0 && State.node.mass.size() > 0)
                fprintf(out_point_state, "  %22.14e",
                        State.node.mass.host(node_gid));

            if (node_temp_id >= 0 && State.node.temp.size() > 0)
                fprintf(out_point_state, "  %22.14e",
                        State.node.temp.host(node_gid));

            fprintf(out_point_state, "\n");

        } // end node loop

        fclose(out_point_state);

    } // end FILE 2 scope

    // -----------------------------------------------------------------------
    //  Concatenate per-rank material-point files into one file per mat_id
    // -----------------------------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    CArray <int> num_universal_mat_elems(num_mats);
    for (int mat_id = 0; mat_id < num_mats; mat_id++){
        MPI_Allreduce(&num_owned_mat_elems(mat_id), &num_universal_mat_elems(mat_id), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

    if (world_size > 1 && rank == 0) {
        for (int mat_id = 0; mat_id < num_mats; mat_id++) {

            // open rank 0's file for this mat_id to count columns and capture header
            char probe_filename[128];
            snprintf(probe_filename, sizeof probe_filename,
                    "state/mat_pt_state_t_%6.4e_mat_id_%d_rank_%d.txt",
                    time_value, mat_id, 0);

            std::ifstream probe_in(probe_filename);
            if (!probe_in.is_open()) {
                std::cerr << "concatenate mat_pt state: could not open "
                        << probe_filename << " to count columns" << std::endl;
                continue;
            }

            std::vector<std::string> header_lines;
            int num_columns = 0;
            std::string probe_line;
            while (std::getline(probe_in, probe_line)) {
                if (probe_line.empty()) continue;

                if (probe_line[0] == '#') {
                    header_lines.push_back(probe_line);
                    continue;
                }

                std::istringstream iss(probe_line);
                std::string tok;
                while (iss >> tok) num_columns++;
                break; // only need the first data row for the column count
            }
            probe_in.close();

            if (num_columns == 0) {
                std::cerr << "concatenate mat_pt state: no data rows found in "
                        << probe_filename << std::endl;
                continue;
            }

            // Store rows sequentially as they are read 
            std::vector<std::vector<double>> all_rows;
            // Pre-allocate capacity to prevent reallocation overhead
            all_rows.reserve(num_universal_mat_elems(mat_id));

            for (int rank_count = 0; rank_count < world_size; rank_count++) {

                char filename[128];
                snprintf(filename, sizeof filename,
                        "state/mat_pt_state_t_%6.4e_mat_id_%d_rank_%d.txt",
                        time_value, mat_id, rank_count);

                std::ifstream in(filename);
                if (!in.is_open()) {
                    std::cerr << "concatenate mat_pt state: could not open "
                            << filename << std::endl;
                    continue;
                }

                std::string line;
                while (std::getline(in, line)) {
                    if (line.empty() || line[0] == '#') continue;

                    std::istringstream iss(line);
                    std::vector<std::string> tok;
                    std::string t;
                    while (iss >> t) tok.push_back(t);
                    if ((int)tok.size() != num_columns) continue; // skip malformed rows
                    
                    // Read row sequentially 
                    std::vector<double> current_row(num_columns);
                    for (int k = 0; k < num_columns; k++) {
                        current_row[k] = std::stod(tok[k]);
                    }

                    // Append the row to our list
                    all_rows.push_back(std::move(current_row));
                }
                in.close();

                // delete the per-rank file now that it's been read in
                std::remove(filename);
            }

            // ---- Sort the accumulated rows by elem_gid (column 0) ----
            std::sort(all_rows.begin(), all_rows.end(), 
                [](const std::vector<double>& a, const std::vector<double>& b) {
                    return a[0] < b[0]; 
                }
            );

            // ---- Make the second column sequential from zero ----
            for (size_t i = 0; i < all_rows.size(); ++i) {
                all_rows[i][1] = (i);
            }

            // ---- write out the concatenated file ----
            char out_filename[128];
            snprintf(out_filename, sizeof out_filename,
                    "state/mat_pt_state_t_%6.4e_mat_id_%d.txt",
                    time_value, mat_id);

            FILE* out = fopen(out_filename, "w");
            if (!out) {
                std::cerr << "concatenate mat_pt state: could not open "
                        << out_filename << " for writing" << std::endl;
                continue;
            }

            for (const std::string& hline : header_lines) {
                fprintf(out, "%s\n", hline.c_str());
            }

            // Iterate over the sorted sequential rows
            for (size_t r = 0; r < all_rows.size(); r++) {
                fprintf(out, "  %-12zu %-8zu ",
                        (size_t)all_rows[r][0], (size_t)all_rows[r][1]);
                for (int k = 2; k < num_columns; k++) {
                    fprintf(out, " %22.14e", all_rows[r][k]);
                }
                fprintf(out, "\n");
            }

            fclose(out);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

} // end write_text_state


/////////////////////////////////////////////////////////////////////////////
///
/// \fn write_text_state_Pn
///
/// \brief  ASCII columnar text-dump companion to write_vtu_Pn.
///         Produces two files per call:
///
///           state/mat_pt_state_t_<time>.txt  – one row per Gauss point,
///                                              all material-point fields
///           state/node_state_t_<time>.txt    – one row per mesh node,
///                                              coords + velocity (+ mass)
///
///         Field selection follows the same slot-ID convention as
///         write_vtu_Pn (-1 means "not requested").  Column headers in
///         both files are prefixed with '#' so Gnuplot, NumPy, etc. can
///         load the files directly without skipping.
///
///         Assumes that host mirrors for State.node.coords, State.node.vel,
///         State.node.mass, and all State.MaterialPoints.* arrays are
///         up-to-date at the time of the call (i.e. update_host() was
///         already called by the caller for each device array).
///
/// \param mesh                  Simulation mesh
/// \param State                 State data
/// \param ref_elem              Reference finite element (quadrature layout)
/// \param mat_scalar_var_names  Registered scalar field names (for headers)
/// \param mat_tensor_var_names  Registered tensor field names (for headers)
/// \param node_scalar_var_names    Registered node scalar field names
/// \param node_vector_var_names    Registered node vector field names
/// \param mat_id                Material index
/// \param num_mat_elems         Number of elements containing this material
/// \param num_nodes_in_elem     Nodes per element
/// \param num_dims              Spatial dimension count (2 or 3)
/// \param time_value            Current simulation time (embedded in filename)
/// \param mat_den_id … mat_specific_heat_id   Slot IDs (-1 = not requested)
///
/////////////////////////////////////////////////////////////////////////////
void MeshWriter::write_text_state_Pn(
    const swage::Mesh_t&                   mesh,
    const State_t&                       State,
    elements::ReferenceElement_t&             ref_elem,
    const std::vector<std::string>&      mat_scalar_var_names,
    const std::vector<std::string>&      mat_tensor_var_names,
    const std::vector<std::string>&      node_scalar_var_names,
    const std::vector<std::string>&      node_vector_var_names,
    const size_t                         num_nodes_in_elem,
    const size_t                         num_dims,
    const double                         time_value,
    // field slot IDs (-1 means "not requested")
    const int mat_den_id,
    const int mat_pres_id,
    const int mat_sie_id,
    const int mat_sspd_id,
    const int mat_mass_id,
    const int mat_mat_volfrac_id,
    const int mat_geo_volfrac_id,
    const int mat_eroded_id,
    const int mat_stress_id,
    const int mat_strain_id,
    const int mat_conductivity_id,
    const int mat_specific_heat_id,
    // node field slot IDs (-1 means "not requested")
    const int node_mass_id,
    const int node_vel_id,
    const int node_coord_id,
    const int node_temp_id,
    const int node_grad_level_set_id,
    const int node_disp_id
)
{
    // -----------------------------------------------------------------------
    //  Derived sizes
    // -----------------------------------------------------------------------

    const size_t num_gp_per_elem = ref_elem.qpt_grad_basis.dims(0);
    const size_t num_mats = State.MaterialPoints.num_material_points.size();

    // -----------------------------------------------------------------------
    //  Gauss-point physical coordinates  (isoparametric mapping)
    //
    //  Same kernel as write_vtu_Pn so the two outputs are numerically
    //  identical for the coordinate columns.
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //  Ensure the output directory exists before opening any files.
    //  create_directories is a no-op (and does not error) if the path
    //  already exists.  Requires #include <filesystem> and C++17.
    // -----------------------------------------------------------------------
    struct stat st;
    MPI_Barrier(MPI_COMM_WORLD);
    //if (mpi_rank == 0) {
        if (stat("state", &st) != 0) {
            system("mkdir state");
        }
    //}
    MPI_Barrier(MPI_COMM_WORLD);

    int rank;
    int world_size;
    mesh_io_mpi_detail::query_world_rank_size(rank,world_size);

    // sizing owned mat elems if not sized already
    if (num_owned_mat_elems.size() <= 0) {
        num_owned_mat_elems = CArray <long long int> (num_mats);
        for (int mat = 0; mat < num_mats; mat++) {
            num_owned_mat_elems(mat) = -1;
        }
    }

    // -----------------------------------------------------------------------
    //  FILE 1 – Material-point (Gauss-point) state
    //
    //  Filename:  state/mat_pt_state_t_<time>.txt
    //  Format:    fixed-width columns, one row per Gauss point,
    //             header line prefixed with '#'
    // -----------------------------------------------------------------------
    for (int mat_id = 0; mat_id < (int)State.MaterialToMeshMaps.num_mat_elems.dims(0); mat_id++) {
        // allocating array for storing real space locations of gauss points
        DCArrayKokkos<double> x_phys;

        // populating the counts for num_owned_mat_elems if not already populated
        if (num_owned_mat_elems(mat_id) == -1) {
            int loc_tally = 0;
            int tally = 0;
            FOR_REDUCE_SUM(mat_elem, 0, State.MaterialToMeshMaps.num_mat_elems.host(mat_id), loc_tally, {
                const size_t elem_id = State.MaterialToMeshMaps.elem_in_mat_elem(mat_id, mat_elem);
                if (elem_id < mesh.num_owned_elems) {
                    loc_tally += 1;
                }
            }, tally);
            // sizing array for storing real space locations of gauss points
            num_owned_mat_elems(mat_id) = tally;
        }

        // initializing array for storing real space locations of gauss points
        x_phys = DCArrayKokkos <double> (num_owned_mat_elems(mat_id), num_gp_per_elem, 3);
        x_phys.set_values(0);

        FOR_ALL(elem, 0, num_owned_mat_elems(mat_id), {
            const size_t elem_id =
                State.MaterialToMeshMaps.elem_in_mat_elem(mat_id, elem);
            for (size_t gp = 0; gp < num_gp_per_elem; gp++) {
                for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                    const size_t node_gid = mesh.nodes_in_elem(elem_id, node_lid);
                    const double N        = ref_elem.qpt_basis(gp, node_lid);
                    for (size_t dim = 0; dim < num_dims; dim++) {
                        x_phys(elem, gp, dim) += N * State.node.coords(node_gid, dim);
                    }
                }
            }
        });
        x_phys.update_host();

        FILE* out_elem_state;
        char  filename[128];
        int   max_len = sizeof filename;
        if (world_size == 1) {
            snprintf(filename, max_len, "state/mat_pt_state_t_%6.4e_mat_id_%d.txt", time_value, mat_id);
        } else {
            snprintf(filename, max_len, "state/mat_pt_state_t_%6.4e_mat_id_%d_rank_%d.txt", time_value, mat_id, rank);
        }

        out_elem_state = fopen(filename, "w");
        if (!out_elem_state) {
            std::cerr << "write_text_state_Pn: could not open "
                    << filename << std::endl;
            return;
        }

        // ---- Metadata comment -------------------------------------------
        fprintf(out_elem_state,
                "# Material-point state  mat_id=%d  time=%.10e\n",
                mat_id, time_value);
        fprintf(out_elem_state,
                "# num_mat_elems=%zu  num_gp_per_elem=%zu\n",
                State.MaterialToMeshMaps.num_mat_elems.host(mat_id), num_gp_per_elem);

        // ---- Column header line -----------------------------------------
        //
        // Always-present columns first, then conditionally-present scalars
        // in the same order as the VTK writer, then expanded tensor columns.
        // All headers are left-aligned in a 22-character field to match the
        // data column width below.
        // -----------------------------------------------------------------
        fprintf(out_elem_state,
                "# %-12s %-8s %-8s  %-22s %-22s %-22s %-22s %-22s",
                "elem_gid", "elem", "gp", "x", "y", "z", "radius_2D", "radius_3D");

        if (mat_den_id          >= 0 && State.MaterialPoints.den.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_den_id].c_str());
        if (mat_pres_id         >= 0 && State.MaterialPoints.pres.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_pres_id].c_str());
        if (mat_sie_id          >= 0 && State.MaterialPoints.sie.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_sie_id].c_str());
        if (mat_sspd_id         >= 0 && State.MaterialPoints.sspd.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_sspd_id].c_str());
        if (mat_mass_id         >= 0 && State.MaterialPoints.mass.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_mass_id].c_str());
        if (mat_mat_volfrac_id      >= 0 && State.MaterialPoints.mat_volfrac.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_mat_volfrac_id].c_str());
        if (mat_geo_volfrac_id  >= 0 && State.MaterialPoints.geo_volfrac.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_geo_volfrac_id].c_str());
        if (mat_eroded_id       >= 0 && State.MaterialPoints.eroded.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_eroded_id].c_str());
        if (mat_conductivity_id >= 0 && State.MaterialPoints.conductivity.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_conductivity_id].c_str());
        if (mat_specific_heat_id >= 0 && State.MaterialPoints.specific_heat.size() > 0)
            fprintf(out_elem_state, "  %-22s",
                    mat_scalar_var_names[mat_specific_heat_id].c_str());

        // Tensor headers: <name>_ij  (i, j ∈ {x, y, z}), row-major
        const char* comp[3] = {"x", "y", "z"};
        if (mat_stress_id >= 0 && State.MaterialPoints.stress.size() > 0) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    char col[40];
                    snprintf(col, sizeof col, "%s_%s%s",
                            mat_tensor_var_names[mat_stress_id].c_str(),
                            comp[i], comp[j]);
                    fprintf(out_elem_state, "  %-22s", col);
                }
            }
        }
        if (mat_strain_id >= 0 && State.MaterialPoints.strain.size() > 0) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    char col[40];
                    snprintf(col, sizeof col, "%s_%s%s",
                            mat_tensor_var_names[mat_strain_id].c_str(),
                            comp[i], comp[j]);
                    fprintf(out_elem_state, "  %-22s", col);
                }
            }
        }
        fprintf(out_elem_state, "\n");

        // ---- Data rows --------------------------------------------------
        for (size_t elem = 0; elem < num_owned_mat_elems(mat_id); elem++) {
            const size_t elem_rid =
                State.MaterialToMeshMaps.elem_in_mat_elem.host(mat_id, elem);
            size_t elem_gid;
            if (world_size > 1) {
                elem_gid = mesh.local_to_global_elem_mapping.host(elem_rid);
            }
            else {
                elem_gid = elem_rid;
            }

            for (size_t gp = 0; gp < num_gp_per_elem; gp++) {
                const size_t pt_id = State.points_in_mat_elem.host(elem, gp);

                // Always-present: global elem id, local elem/gp indices, coords
                fprintf(out_elem_state,
                        "  %-12zu %-8zu %-8zu %22.14e %22.14e %22.14e %22.14e %22.14e",
                        elem_gid, elem, gp,
                        x_phys.host(elem, gp, 0),
                        x_phys.host(elem, gp, 1),
                        (num_dims == 3) ? x_phys.host(elem, gp, 2) : 0.0,
                        sqrt(x_phys.host(elem, gp, 0)*x_phys.host(elem, gp, 0)+x_phys.host(elem, gp, 1)*x_phys.host(elem, gp, 1)),
                        (num_dims == 3) ? sqrt(x_phys.host(elem, gp, 0)*x_phys.host(elem, gp, 0)+x_phys.host(elem, gp, 1)*x_phys.host(elem, gp, 1)+x_phys.host(elem, gp, 2)*x_phys.host(elem, gp, 2)) : 0.0);
                
                // Scalar fields (in the same order as the header)
                if (mat_den_id          >= 0 && State.MaterialPoints.den.size() > 0)
                    fprintf(out_elem_state, "  %22.14e",
                            State.MaterialPoints.den.host(mat_id, pt_id));
                if (mat_pres_id         >= 0 && State.MaterialPoints.pres.size() > 0)
                    fprintf(out_elem_state, "  %22.14e",
                            State.MaterialPoints.pres.host(mat_id, pt_id));
                if (mat_sie_id          >= 0 && State.MaterialPoints.sie.size() > 0)
                    fprintf(out_elem_state, "  %22.14e",
                            State.MaterialPoints.sie.host(mat_id, pt_id));
                if (mat_sspd_id         >= 0 && State.MaterialPoints.sspd.size() > 0)
                    fprintf(out_elem_state, "  %22.14e",
                            State.MaterialPoints.sspd.host(mat_id, pt_id));
                if (mat_mass_id         >= 0 && State.MaterialPoints.mass.size() > 0)
                    fprintf(out_elem_state, "  %22.14e",
                            State.MaterialPoints.mass.host(mat_id, pt_id));
                if (mat_mat_volfrac_id      >= 0 && State.MaterialPoints.mat_volfrac.size() > 0)
                    fprintf(out_elem_state, "  %22.14e",
                            State.MaterialPoints.mat_volfrac.host(mat_id, pt_id));
                if (mat_geo_volfrac_id  >= 0 && State.MaterialPoints.geo_volfrac.size() > 0)
                    fprintf(out_elem_state, "  %22.14e",
                            State.MaterialPoints.geo_volfrac.host(mat_id, pt_id));
                if (mat_eroded_id       >= 0 && State.MaterialPoints.eroded.size() > 0)
                    fprintf(out_elem_state, "  %22.14e",
                            static_cast<double>(
                                State.MaterialPoints.eroded.host(mat_id, pt_id)));
                if (mat_conductivity_id >= 0 && State.MaterialPoints.conductivity.size() > 0)
                    fprintf(out_elem_state, "  %22.14e",
                            State.MaterialPoints.conductivity.host(mat_id, pt_id));
                if (mat_specific_heat_id >= 0 && State.MaterialPoints.specific_heat.size() > 0)
                    fprintf(out_elem_state, "  %22.14e",
                            State.MaterialPoints.specific_heat.host(mat_id, pt_id));

                // Tensor fields: all 9 components, row-major (Txx Txy Txz ...)
                if (mat_stress_id >= 0 && State.MaterialPoints.stress.size() > 0) {
                    for (size_t i = 0; i < 3; i++)
                        for (size_t j = 0; j < 3; j++)
                            fprintf(out_elem_state, "  %22.14e",
                                    State.MaterialPoints.stress.host(
                                        mat_id, pt_id, i, j));
                }
                if (mat_strain_id >= 0 && State.MaterialPoints.strain.size() > 0) {
                    for (size_t i = 0; i < 3; i++)
                        for (size_t j = 0; j < 3; j++)
                            fprintf(out_elem_state, "  %22.14e",
                                    State.MaterialPoints.strain.host(
                                        mat_id, pt_id, i, j));
                }

                fprintf(out_elem_state, "\n");

            } // end gp
        } // end elem

        fclose(out_elem_state);

    } // end FILE 1 scope

    // -----------------------------------------------------------------------
    //  FILE 2 – Node state
    //
    //  Filename:  state/node_state_t_<time>.txt
    //  Format:    fixed-width columns, one row per owned mesh node,
    //             header line prefixed with '#'
    //
    //  All fields are driven by their slot ID (-1 = skip), following the
    //  same convention as the mat-pt file above.
    //
    //  Vector fields (coord, vel, accel, disp, grad_level_set) produce
    //  four columns each:  <name>_x  <name>_y  <name>_z  |<name>|
    //  The Z component is 0.0 for 2-D runs.
    //
    //  Scalar fields (mass, temp) produce a single column each.
    //
    //  Uses mesh.num_owned_nodes to skip ghost nodes in MPI runs,
    //  mirroring the behaviour of write_material_point_state.
    // -----------------------------------------------------------------------
    {
        FILE* out_point_state;
        char  filename[128];
        int   max_len = sizeof filename;
        if (world_size == 1) {
            snprintf(filename, max_len, "state/node_state_t_%6.4e.txt", time_value);
        } else {
            snprintf(filename, max_len, "state/node_state_t_%6.4e_rank_%d.txt", time_value, rank);
        }

        out_point_state = fopen(filename, "w");
        if (!out_point_state) {
            std::cerr << "write_text_state_Pn: could not open "
                    << filename << std::endl;
            return;
        }

        const size_t num_owned_nodes = mesh.num_owned_nodes;

        // ---- Metadata comment -------------------------------------------
        fprintf(out_point_state,
                "# Node state  time=%.10e\n", time_value);
        fprintf(out_point_state,
                "# num_owned_nodes=%zu  num_dims=%zu\n",
                num_owned_nodes, num_dims);

        // ---- Column header line -----------------------------------------
        //
        // node_id is always present.  Every requested vector field adds four
        // columns (<name>_x/y/z, |<name>|); every requested scalar adds one.
        // -----------------------------------------------------------------

        fprintf(out_point_state, "# %-12s", "node_id");

        // Lambda: emit the four header columns for one vector field.
        auto hdr_vec = [&out_point_state](const std::string& name) {
            char col[48];
            snprintf(col, sizeof col, "%s_x", name.c_str());
            fprintf(out_point_state, "  %-22s", col);
            snprintf(col, sizeof col, "%s_y", name.c_str());
            fprintf(out_point_state, "  %-22s", col);
            snprintf(col, sizeof col, "%s_z", name.c_str());
            fprintf(out_point_state, "  %-22s", col);
            snprintf(col, sizeof col, "|%s|", name.c_str());
            fprintf(out_point_state, "  %-22s", col);
        };

        if (node_coord_id          >= 0)
            hdr_vec(node_vector_var_names[node_coord_id]);
        if (node_vel_id            >= 0 && State.node.vel.size() > 0)
            hdr_vec(node_vector_var_names[node_vel_id]);
        if (node_disp_id           >= 0 && State.node.displacement.size() > 0)
            hdr_vec(node_vector_var_names[node_disp_id]);
        if (node_grad_level_set_id >= 0 && State.node.gradient_level_set.size() > 0)
            hdr_vec(node_vector_var_names[node_grad_level_set_id]);
        if (node_mass_id           >= 0 && State.node.mass.size() > 0)
            fprintf(out_point_state, "  %-22s",
                    node_scalar_var_names[node_mass_id].c_str());
        if (node_temp_id           >= 0 && State.node.temp.size() > 0)
            fprintf(out_point_state, "  %-22s",
                    node_scalar_var_names[node_temp_id].c_str());

        fprintf(out_point_state, "\n");

        // ---- Data rows --------------------------------------------------
        //
        // Lambda: write 3 components + magnitude for one vector field.
        // Raw Z value is zeroed here for 2-D runs so callers always pass
        // the full [2] index without a conditional at the call site.
        auto dat_vec = [&out_point_state](double vx, double vy, double vz_raw, size_t num_dims) {
            const double vz  = (num_dims == 3) ? vz_raw : 0.0;
            const double mag = sqrt(vx*vx + vy*vy + vz*vz);
            fprintf(out_point_state,
                    "  %22.14e  %22.14e  %22.14e  %22.14e",
                    vx, vy, vz, mag);
        };

        for (size_t node_gid = 0; node_gid < num_owned_nodes; node_gid++) {

            // Always-present row identifier
            fprintf(out_point_state, "  %-12zu", node_gid);

            // ---- Vector fields ------------------------------------------

            if (node_coord_id >= 0)
                dat_vec(State.node.coords.host(node_gid, 0),
                        State.node.coords.host(node_gid, 1),
                        State.node.coords.host(node_gid, 2),
                        mesh.num_dims);

            if (node_vel_id >= 0 && State.node.vel.size() > 0)
                dat_vec(State.node.vel.host(node_gid, 0),
                        State.node.vel.host(node_gid, 1),
                        State.node.vel.host(node_gid, 2),
                        mesh.num_dims);

            if (node_disp_id >= 0 && State.node.displacement.size() > 0)
                dat_vec(State.node.displacement.host(node_gid, 0),
                        State.node.displacement.host(node_gid, 1),
                        State.node.displacement.host(node_gid, 2),
                        mesh.num_dims);

            if (node_grad_level_set_id >= 0 && State.node.gradient_level_set.size() > 0)
                dat_vec(State.node.gradient_level_set.host(node_gid, 0),
                        State.node.gradient_level_set.host(node_gid, 1),
                        State.node.gradient_level_set.host(node_gid, 2),
                        mesh.num_dims);

            // ---- Scalar fields ------------------------------------------

            if (node_mass_id >= 0 && State.node.mass.size() > 0)
                fprintf(out_point_state, "  %22.14e",
                        State.node.mass.host(node_gid));

            if (node_temp_id >= 0 && State.node.temp.size() > 0)
                fprintf(out_point_state, "  %22.14e",
                        State.node.temp.host(node_gid));

            fprintf(out_point_state, "\n");

        } // end node loop

        fclose(out_point_state);

    } // end FILE 2 scope

    MPI_Barrier(MPI_COMM_WORLD);
    CArray <int> num_universal_mat_elems(num_mats);
    for (int mat_id = 0; mat_id < num_mats; mat_id++){
        MPI_Allreduce(&num_owned_mat_elems(mat_id), &num_universal_mat_elems(mat_id), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

    // concatenate material outputs if more than one rank
    if (world_size > 1 && rank == 0) {
        // looping through materials
        for (int mat_id = 0; mat_id < num_mats; mat_id++) {

            // open rank 0's file for this mat_id to count columns and capture header
            char probe_filename[128];
            snprintf(probe_filename, sizeof probe_filename,
                    "state/mat_pt_state_t_%6.4e_mat_id_%d_rank_%d.txt",
                    time_value, mat_id, 0);

            std::ifstream probe_in(probe_filename);
            if (!probe_in.is_open()) {
                std::cerr << "concatenate mat_pt state: could not open "
                        << probe_filename << " to count columns" << std::endl;
                continue;
            }

            std::vector<std::string> header_lines;
            int num_columns = 0;
            std::string probe_line;
            while (std::getline(probe_in, probe_line)) {
                if (probe_line.empty()) continue;

                if (probe_line[0] == '#') {
                    header_lines.push_back(probe_line);
                    continue;
                }

                std::istringstream iss(probe_line);
                std::string tok;
                while (iss >> tok) num_columns++;
                break; // only need the first data row for the column count
            }
            probe_in.close();

            if (num_columns == 0) {
                std::cerr << "concatenate mat_pt state: no data rows found in "
                        << probe_filename << std::endl;
                continue;
            }

            // Store rows sequentially as they are read 
            std::vector<std::vector<double>> all_rows;
            // Pre-allocate capacity to prevent reallocation overhead
            all_rows.reserve(num_universal_mat_elems(mat_id) * num_gp_per_elem);

            for (int rank_count = 0; rank_count < world_size; rank_count++) {

                char filename[128];
                snprintf(filename, sizeof filename,
                        "state/mat_pt_state_t_%6.4e_mat_id_%d_rank_%d.txt",
                        time_value, mat_id, rank_count);

                std::ifstream in(filename);
                if (!in.is_open()) {
                    std::cerr << "concatenate mat_pt state: could not open "
                            << filename << std::endl;
                    continue;
                }

                std::string line;
                while (std::getline(in, line)) {
                    if (line.empty() || line[0] == '#') continue;

                    std::istringstream iss(line);
                    std::vector<std::string> tok;
                    std::string t;
                    while (iss >> t) tok.push_back(t);
                    if ((int)tok.size() != num_columns) continue; // skip malformed rows

                    // Read row sequentially 
                    std::vector<double> current_row(num_columns);
                    for (int k = 0; k < num_columns; k++) {
                        current_row[k] = std::stod(tok[k]);
                    }

                    // Append the row to our list
                    all_rows.push_back(std::move(current_row));
                }
                in.close();

                // delete the per-rank file now that it's been read in
                std::remove(filename);
            }

            // ---- Sort the accumulated rows ----
            // Sort primarily by elem_gid (column 0), secondarily by gp (column 2)
            std::sort(all_rows.begin(), all_rows.end(), 
                [](const std::vector<double>& a, const std::vector<double>& b) {
                    if (a[0] != b[0]) {
                        return a[0] < b[0]; // Sort by elem_gid
                    }
                    return a[2] < b[2];     // Sort by gp if elem_gid is identical
                }
            );

            // ---- Make the second column sequential from zero ----
            for (size_t i = 0; i < all_rows.size(); ++i) {
                all_rows[i][1] = static_cast<double>(all_rows[i][0]);
            }

            // ---- write out the concatenated file ----
            char out_filename[128];
            snprintf(out_filename, sizeof out_filename,
                    "state/mat_pt_state_t_%6.4e_mat_id_%d.txt",
                    time_value, mat_id);

            FILE* out = fopen(out_filename, "w");
            if (!out) {
                std::cerr << "concatenate mat_pt state: could not open "
                        << out_filename << " for writing" << std::endl;
                continue;
            }

            for (const std::string& hline : header_lines) {
                fprintf(out, "%s\n", hline.c_str());
            }

            // Iterate over the sorted sequential rows
            for (size_t r = 0; r < all_rows.size(); r++) {
                fprintf(out, "  %-12zu %-8zu %-8zu",
                        (size_t)all_rows[r][0], (size_t)all_rows[r][1], (size_t)all_rows[r][2]);
                for (int k = 3; k < num_columns; k++) {
                    fprintf(out, "  %22.14e", all_rows[r][k]);
                }
                fprintf(out, "\n");
            }

            fclose(out);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

} // end write_text_state_Pn
MeshWriter::MeshWriter()
{}

MeshWriter::~MeshWriter()
{
}
