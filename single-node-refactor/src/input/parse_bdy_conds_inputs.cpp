/**********************************************************************************************
� 2020. Triad National Security, LLC. All rights reserved.
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
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include <vector>
#include <variant>
#include <algorithm>
#include <map>

#include "string_utils.h"

#include "matar.h"
#include "parse_tools.hpp"
#include "parse_bdy_conds_inputs.hpp"

// simulation parameters contains:
//   mesh_input
//   output_options
//   dynamic_options
//   solver_inputs
//   region_setups
#include "simulation_parameters.h"

// boundary conditions
#include "boundary_conditions.h"

// velocity bc files
#include "constant_velocity_bc.h"
#include "no_velocity_bc.h"
#include "piston_velocity_bc.h"
#include "reflected_velocity_bc.h"
#include "time_varying_velocity_bc.h"
#include "user_defined_velocity_bc.h"
#include "zero_velocity_bc.h"


// temperature bc files
#include "constant_temp_bc.h"

// stress bc files
#include "constant_stress_bc.h"
#include "no_stress_bc.h"
#include "time_varying_stress_bc.h"
#include "user_defined_stress_bc.h"





// ==============================================================================
//   Function Definitions
// ==============================================================================



// =================================================================================
//    Parse Boundary Conditions
// =================================================================================
void parse_bcs(Yaml::Node& root, BoundaryCondition_t& BoundaryConditions, const size_t num_solvers)
{

    Yaml::Node& bc_yaml = root["boundary_conditions"];

    size_t num_bcs = bc_yaml.Size();

    std::cout<<"Number of boundary conditions = " << num_bcs << std::endl;

    BoundaryConditions.num_bcs = num_bcs;

    BoundaryConditions.BoundaryConditionSetup = CArrayKokkos <BoundaryConditionSetup_t>(num_bcs, "bc_setup_vars");

    // device functions
    BoundaryConditions.BoundaryConditionFunctions = CArrayKokkos <BoundaryConditionFunctions_t>(num_bcs, "bc_fcns");

    // enums to select options with boundary conditions
    BoundaryConditions.BoundaryConditionEnums  = DCArrayKokkos<BoundaryConditionEnums_t> (num_bcs,"bc_enums");  

    // --- BC velocity ---
    // stores the velocity bdy node lists per solver, in the future, this needs to be a DualRaggedRight
    BoundaryConditions.vel_bdy_sets_in_solver = DCArrayKokkos<size_t> (num_solvers, num_bcs, "vel_bdy_sets_in_solver");  
    BoundaryConditions.temperature_bdy_sets_in_solver = DCArrayKokkos<size_t> (num_solvers, num_bcs, "temperature_bdy_sets_in_solver");
    // this stores the number of bdy sets for a solver
   

    // this stores the number of vel bdy sets for a solver
    BoundaryConditions.num_vel_bdy_sets_in_solver = DCArrayKokkos<size_t> (num_solvers, "num_vel_bdy_sets_in_solver");   
    BoundaryConditions.num_temperature_bdy_sets_in_solver = DCArrayKokkos<size_t> (num_solvers, "num_temperature_bdy_sets_in_solver");
    
    // set the storage counter to zero
    for(size_t solver_id=0; solver_id<num_solvers; solver_id++){
        BoundaryConditions.num_vel_bdy_sets_in_solver.host(solver_id) = 0;
        BoundaryConditions.num_temperature_bdy_sets_in_solver.host(solver_id) = 0;
    } // end for

    // --- BC stress ---
    // stores the stress bdy node lists per solver, in the future, this needs to be a DualRaggedRight
    BoundaryConditions.stress_bdy_sets_in_solver = DCArrayKokkos<size_t> (num_solvers, num_bcs, "stress_bdy_sets_in_solver");  
    // this stores the number of stess bdy sets for a solver
    BoundaryConditions.num_stress_bdy_sets_in_solver = DCArrayKokkos<size_t> (num_solvers, "num_stress_bdy_sets_in_solver");   
    // set the storage counter to zero
    for(size_t solver_id=0; solver_id<num_solvers; solver_id++){
        BoundaryConditions.num_stress_bdy_sets_in_solver.host(solver_id) = 0;
    } // end for


    // temporary arrays for boundary condition variables
    DCArrayKokkos<double> tempVelocityBCGlobalVars (num_bcs, 100, "temporary_velocity_bc_global_values");

    DCArrayKokkos<double> tempTemperatureBCGlobalVars (num_bcs, 100, "temporary_temperature_bc_global_values");
    DCArrayKokkos<double> tempStressBCGlobalVars (num_bcs, 100, "temporary_stress_bc_global_values");
    // DCArrayKokkos<double> tempHeatFluxBCGlobalVars (num_bcs, 100, "temporary_heat_flux_bc_global_values");
    
    BoundaryConditions.num_velocity_bc_global_vars = CArrayKokkos <size_t>(num_bcs, "BoundaryConditions.num_velocity_bc_global_vars"); 
    BoundaryConditions.num_temperature_bc_global_vars = CArrayKokkos <size_t>(num_bcs, "BoundaryConditions.num_temperature_bc_global_vars");
    BoundaryConditions.num_stress_bc_global_vars = CArrayKokkos <size_t>(num_bcs, "BoundaryConditions.num_stress_bc_global_vars");
    // BoundaryConditions.num_heat_flux_bc_global_vars = CArrayKokkos <size_t>(num_bcs, "BoundaryConditions.num_heat_flux_bc_global_vars"); 

    
    


    // initialize the num of global vars to 0 for all models
    FOR_ALL(bc_id, 0, num_bcs, {
        BoundaryConditions.num_velocity_bc_global_vars(bc_id) = 0;

        BoundaryConditions.num_temperature_bc_global_vars(bc_id) = 0;
        BoundaryConditions.num_stress_bc_global_vars(bc_id) = 0;
        // BoundaryConditions.num_heat_flux_bc_global_vars(bc_id) = 0;
    }); // end parallel for


    // state place holder is here
    BoundaryConditions.bc_state_vars  = DCArrayKokkos<double>(num_bcs, 4, "bc_state_values");  // WARNING a place holder

    // loop over the BC specified
    for (size_t bc_id = 0; bc_id < num_bcs; bc_id++) {
        // read the variables names
        Yaml::Node& inps_yaml = bc_yaml[bc_id]["boundary_condition"];

        // get the boundary condition variables names set by the user
        std::vector<std::string> user_str_bc_inps;

        // extract words from the input file and validate they are correct
        validate_inputs(inps_yaml, user_str_bc_inps, str_bc_inps, bc_required_inps);


        // verify the boundary condition block connects to a solver
        // loop over the words in the boundary input definition and find the solver id
        int solver_id = -1;
        for (auto& a_word : user_str_bc_inps) {

            Yaml::Node& inps_yaml = bc_yaml[bc_id]["boundary_condition"][a_word];

            if (a_word.compare("solver_id") == 0) {
                solver_id = bc_yaml[bc_id]["boundary_condition"][a_word].As<int>();

                if (solver_id<0 || solver_id>=num_solvers){
                    std::cout << "ERROR: invalid solver_id specified in the boundary condition definition. Either negative or >= num_solvers" << std::endl;
            
                    throw std::runtime_error("**** Solver_id is out of bounds ****");
                } // end check on m_id range

            } // end id

            // add other checks here...

        } // end loop over all boundary condition inputs

        if (solver_id<0){
            std::cout << "ERROR: solver_id must be specified in the boundary condition definition " << std::endl;
            
            throw std::runtime_error("**** Solver_id is missing ****");
        } // end check on m_id range


        // loop over the words in the boundary condition input definition
        for (auto& a_word : user_str_bc_inps) {
            
            Yaml::Node& inps_yaml = bc_yaml[bc_id]["boundary_condition"][a_word];

            // get solver for this boundary condition
            if (a_word.compare("solver_id") == 0) {
                // do nothing, I already have solver_id since the check above

            } // solver id
            // get boundary condition type
            else if (a_word.compare("velocity_model") == 0) {
                
                // Note: solver_id was retrieved at the top of the bc_id loop

                // find out how many velocity bdy sets have been saved 
                size_t num_saved = BoundaryConditions.num_vel_bdy_sets_in_solver.host(solver_id);
                BoundaryConditions.vel_bdy_sets_in_solver.host(solver_id, num_saved) = bc_id;
                BoundaryConditions.num_vel_bdy_sets_in_solver.host(solver_id) += 1;  // increment saved counter

                std::string velocity_model = bc_yaml[bc_id]["boundary_condition"][a_word].As<std::string>();

                auto map = bc_velocity_model_map; 

                // set the velocity_model
                if (map.find(velocity_model) != map.end()) {
                    auto bc_velocity_model = map[velocity_model];

                    // bc_velocity_model_map[velocity_model] returns enum value, e.g., boundary_conditions::velocity_constant
                    switch(map[velocity_model]){

                        case boundary_conditions::constantVelocityBC :
                            std::cout << "Setting constant velocity bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCVelocityModel = boundary_conditions::constantVelocityBC ;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).velocity = &ConstantVelocityBC::velocity;
                            });
                            break;

                        case boundary_conditions::timeVaryingVelocityBC:
                            std::cout << "Setting time varying velocity bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCVelocityModel = boundary_conditions::timeVaryingVelocityBC;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).velocity = &TimeVaryingVelocityBC::velocity;
                            });
                            break;
                        
                        case boundary_conditions::reflectedVelocityBC:
                            std::cout << "Setting reflected velocity bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCVelocityModel = boundary_conditions::reflectedVelocityBC;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).velocity = &ReflectedVelocityBC::velocity;
                            });
                            break;

                        case boundary_conditions::zeroVelocityBC:
                            std::cout << "Setting zero velocity bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCVelocityModel = boundary_conditions::zeroVelocityBC;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).velocity = &ZeroVelocityBC::velocity;
                            });
                            break;
                        case boundary_conditions::userDefinedVelocityBC:
                            std::cout << "Setting user defined velocity bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCVelocityModel = boundary_conditions::userDefinedVelocityBC;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).velocity = &UserDefinedVelocityBC::velocity;
                            });
                            break;
                        case boundary_conditions::pistonVelocityBC:
                            std::cout << "Setting piston velocity bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCVelocityModel = boundary_conditions::pistonVelocityBC;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).velocity = &PistonVelocityBC::velocity;
                            });
                            break;                        
                        default:
                            
                            std::cout << "ERROR: invalid velocity boundary condition input: " << velocity_model << std::endl;
                            throw std::runtime_error("**** Velocity BC model Not Understood ****");
                            break;
                        
                    } // end switch

                }
                else{
                    std::cout << "ERROR: invalid boundary condition option input in YAML file: " << velocity_model << std::endl;
                    std::cout << "Valid options are: " << std::endl;

                    for (const auto& pair : map) {
                        std::cout << "\t" << pair.first << std::endl;
                    }

                    throw std::runtime_error("**** Boundary Condition Velocity Model Not Understood ****");
                } // end if
            } // type

            else if (a_word.compare("temperature_model") == 0) {
                std::cout << "Inside temperature_model check" << std::endl;

                // Note: solver_id was retrieved at the top of the bc_id loop

                std::cout<<"Solver id = " << solver_id << std::endl;
                std::cout<<"bc_id = " << bc_id << std::endl;

                // find out how many temperature bdy sets have been saved 
                size_t num_saved = BoundaryConditions.num_temperature_bdy_sets_in_solver.host(solver_id);

                std::cout<<"num_saved = " << num_saved << std::endl;

                BoundaryConditions.temperature_bdy_sets_in_solver.host(solver_id,num_saved) = bc_id;
                BoundaryConditions.num_temperature_bdy_sets_in_solver.host(solver_id) += 1;  // increment saved counter

                std::string temperature_model = bc_yaml[bc_id]["boundary_condition"][a_word].As<std::string>();

                auto map = bc_temperature_model_map; 
                std::cout<<"Before map check" << std::endl;
                // set the temperature_model
                if (map.find(temperature_model) != map.end()) {
                    std::cout<<"Inside map check" << std::endl;
                    auto bc_temperature_model = map[temperature_model];
                    
                    switch(map[temperature_model]){
                        case boundary_conditions::constantTemperatureBC:
                            std::cout << "Setting constant temperature bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCTemperatureModel = boundary_conditions::constantTemperatureBC;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).temperature = &ConstantTemperatureBC::temperature;
                            });
                            break;

                        case boundary_conditions::convectionTemperatureBC:
                            std::cout << "Setting convection bc " << std::endl;
                            
                            RUN({   
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCTemperatureModel = boundary_conditions::convectionTemperatureBC;
                                // BoundaryConditions.BoundaryConditionFunctions(bc_id).temperature = &ConvectionTemperatureBC::temperature;
                            });
                            break;

                        case boundary_conditions::radiationTemperatureBC:
                            std::cout << "Setting radiation bc " << std::endl;
                            
                            RUN({   
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCTemperatureModel = boundary_conditions::radiationTemperatureBC;
                                // BoundaryConditions.BoundaryConditionFunctions(bc_id).temperature = &RadiationTemperatureBC::temperature;
                            });
                            break;

                        default:
                            std::cout << "ERROR: invalid temperature boundary condition input: " << temperature_model << std::endl;
                            throw std::runtime_error("**** Temperature BC model Not Understood ****");
                            break;
                    }
                }

            }   

            // get boundary condition type
            else if (a_word.compare("stress_model") == 0) {
                
                // Note: solver_id was retrieved at the top of the bc_id loop

                // find out how many stress bdy sets have been saved 
                size_t num_saved = BoundaryConditions.num_stress_bdy_sets_in_solver.host(solver_id);
                BoundaryConditions.stress_bdy_sets_in_solver.host(solver_id, num_saved) = bc_id;
                BoundaryConditions.num_stress_bdy_sets_in_solver.host(solver_id) += 1;  // increment saved counter

                std::string stress_model = bc_yaml[bc_id]["boundary_condition"][a_word].As<std::string>();

                auto map = bc_stress_model_map; 

                // set the stress_model
                if (map.find(stress_model) != map.end()) {
                    auto bc_stress_model = map[stress_model];

                    // bc_stress_model_map[stress_model] returns enum value, e.g., boundary_conditions::stress_constant
                    switch(map[stress_model]){

                        case boundary_conditions::constantStressBC :
                            std::cout << "Setting stress bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCStressModel = boundary_conditions::constantStressBC ;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).stress = &ConstantStressBC::stress;
                            });
                            break;

                        case boundary_conditions::timeVaryingStressBC:
                            std::cout << "Setting stress bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCStressModel = boundary_conditions::timeVaryingStressBC;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).stress = &TimeVaryingStressBC::stress;
                            });
                            break;

                        case boundary_conditions::userDefinedStressBC:
                            std::cout << "Setting stress bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCStressModel = boundary_conditions::userDefinedStressBC;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).stress = &UserDefinedStressBC::stress;
                            });
                            break;
                      
                        default:
                            
                            std::cout << "ERROR: invalid stress boundary condition input: " << stress_model << std::endl;
                            throw std::runtime_error("**** stress BC model Not Understood ****");
                            break;
                        
                    } // end switch

                }
                else{
                    std::cout << "ERROR: invalid boundary condition option input in YAML file: " << stress_model << std::endl;
                    std::cout << "Valid options are: " << std::endl;

                    for (const auto& pair : map) {
                        std::cout << "\t" << pair.first << std::endl;
                    }


                    throw std::runtime_error("**** Boundary Condition Stress Model Not Understood ****");
                } // end if
            } // type of stress model

            // get boundary condition location -- host or device
            else if (a_word.compare("location") == 0) {
                std::string location = bc_yaml[bc_id]["boundary_condition"][a_word].As<std::string>();

                auto map = bc_location_map;

                // set the location
                if (map.find(location) != map.end()) {
                    auto bc_location = map[location];
                    RUN({
                        BoundaryConditions.BoundaryConditionEnums(bc_id).Location = bc_location;
                    });

                }
                else{
                    std::cout << "ERROR: invalid boundary condition option input in YAML file: " << location << std::endl;
                    std::cout << "Valid options are: " << std::endl;

                    for (const auto& pair : map) {
                        std::cout << "\t" << pair.first << std::endl;
                    }
                    throw std::runtime_error("**** Boundary Conditions Not Understood ****");
                } // end if
            } // location

            // get boundary condition surface geometry
            else if (a_word.compare("surface") == 0) {

                // -----
                // loop over the sub fields under surface
                // -----
                Yaml::Node& inps_subfields_yaml = bc_yaml[bc_id]["boundary_condition"]["surface"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_bc_surface_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_bc_surface_inputs, str_bc_surface_inps, bc_surface_required_inps);


                // loop over the subfield words
                for(auto& a_subfield_word : user_bc_surface_inputs){ 

                    if (a_subfield_word.compare("type") == 0){
                        std::string surface = bc_yaml[bc_id]["boundary_condition"]["surface"][a_subfield_word].As<std::string>();

                        auto map = bc_surface_map;

                        // set the surface
                        if (map.find(surface) != map.end()) {
                            auto bc_surface = map[surface];
                            RUN({
                                BoundaryConditions.BoundaryConditionSetup(bc_id).surface = bc_surface;
                            });

                        }
                        else{
                            std::cout << "ERROR: invalid boundary condition option input in YAML file: " << surface << std::endl;
                            std::cout << "Valid options are: " << std::endl;

                            for (const auto& pair : map) {
                                std::cout << "\t" << pair.first << std::endl;
                            }
                            throw std::runtime_error("**** Boundary Condition Surface Inputs Not Understood ****");
                        } // end if

                    } // end if type
                    else if (a_subfield_word.compare("plane_position") == 0) {
                        double value = bc_yaml[bc_id]["boundary_condition"]["surface"][a_subfield_word].As<double>();
                        RUN({
                            BoundaryConditions.BoundaryConditionSetup(bc_id).value = value;
                        });
                    } // end if plane position
                    else if (a_subfield_word.compare("radius") == 0) {
                        double value = bc_yaml[bc_id]["boundary_condition"]["surface"][a_subfield_word].As<double>();
                        RUN({
                            BoundaryConditions.BoundaryConditionSetup(bc_id).value = value;
                        });
                    } // end if radius
                    else if (a_subfield_word.compare("tolerance") == 0) {
                        // the tolerance to tag a surface
                        double tolerance = bc_yaml[bc_id]["boundary_condition"]["surface"][a_subfield_word].As<double>();
                        RUN({
                            BoundaryConditions.BoundaryConditionSetup(bc_id).tolerance = tolerance;
                        });
                    } // end if tolerance 
                    else if (a_subfield_word.compare("origin") == 0) {
                        std::string origin = bc_yaml[bc_id]["boundary_condition"]["surface"][a_subfield_word].As<std::string>();

                        // get the origin numbers, values are words
                        std::vector<std::string> numbers = exact_array_values(origin, ",");

                        double x1 = std::stod(numbers[0]);
                        double y1 = std::stod(numbers[1]);
                        double z1;

                        if(numbers.size()==3){ 
                            // 3D
                            z1 = std::stod(numbers[2]);
                        }
                        else {
                            // 2D
                            z1 = 0.0;
                        } //

                        // storing the origin values as (x1,y1,z1)

                        RUN({
                            BoundaryConditions.BoundaryConditionSetup(bc_id).origin[0] = x1;
                            BoundaryConditions.BoundaryConditionSetup(bc_id).origin[1] = y1;
                            BoundaryConditions.BoundaryConditionSetup(bc_id).origin[2] = z1;
                        });

                    } //end origin
                    else {
                        // word is unknown
                        std::cout << "ERROR: invalid input under boundary condition geometery: " << a_subfield_word << std::endl;
                        std::cout << "Valid options are: " << std::endl;
                        for (const auto& element : str_bc_surface_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** Boundary Conditions Not Understood ****");
                    }

                } // end loop over words in the subfield
            } // surface

            
            // Set the global variables for velocity boundary condition models
            else if (a_word.compare("velocity_bc_global_vars") == 0) {
                Yaml::Node & vel_bc_global_vars_yaml = bc_yaml[bc_id]["boundary_condition"][a_word];

                size_t num_global_vars = vel_bc_global_vars_yaml.Size();

                if(num_global_vars > 100){
                    throw std::runtime_error("**** Per boundary condition, the code only supports up to 100 velocity global vars in the input file ****");
                } // end check on num_global_vars

                RUN({ 
                    BoundaryConditions.num_velocity_bc_global_vars(bc_id) = num_global_vars;
                });

                // store the global eos model parameters
                for (int global_var_id = 0; global_var_id < num_global_vars; global_var_id++) {
                    double velocity_bc_var = bc_yaml[bc_id]["boundary_condition"]["velocity_bc_global_vars"][global_var_id].As<double>();
                    
                    RUN({
                        tempVelocityBCGlobalVars(bc_id, global_var_id) = velocity_bc_var;
                    });

                } // end loop over global vars
            } // end else if on velocity_bc_global_vars

            
            // Set the global variables for temperature boundary condition models
            else if (a_word.compare("temperature_bc_global_vars") == 0) {
                std::cout << "Inside temperature_bc_global_vars" << std::endl;
                Yaml::Node & temp_bc_global_vars_yaml = bc_yaml[bc_id]["boundary_condition"][a_word];

                size_t num_global_vars = temp_bc_global_vars_yaml.Size();

                if(num_global_vars > 100){
                    throw std::runtime_error("**** Per boundary condition, the code only supports up to 100 temperature global vars in the input file ****");
                } // end check on num_global_vars

                RUN({ 
                    BoundaryConditions.num_temperature_bc_global_vars(bc_id) = num_global_vars;
                });

                // Store the global temperature boundary condition variables
                for (int global_var_id = 0; global_var_id < num_global_vars; global_var_id++) {
                    double temperature_bc_var = bc_yaml[bc_id]["boundary_condition"]["temperature_bc_global_vars"][global_var_id].As<double>();
                    
                    RUN({
                        tempTemperatureBCGlobalVars(bc_id, global_var_id) = temperature_bc_var;
                    });

                }
            } // end else if on temperature_bc_global_vars

      
            // set the stress global values
            else if (a_word.compare("stress_bc_global_vars") == 0) {

                Yaml::Node & stress_bc_global_vars_yaml = bc_yaml[bc_id]["boundary_condition"][a_word];

                size_t num_global_vars = stress_bc_global_vars_yaml.Size();

                if(num_global_vars>100){
                    throw std::runtime_error("**** Per boundary condition, the code only supports up to 100 velocity global vars in the input file ****");
                } // end check on num_global_vars

                RUN({ 
                    printf("num global stress vars = %zu \n", num_global_vars);
                    BoundaryConditions.num_stress_bc_global_vars(bc_id) = num_global_vars;
                });

                // store the global eos model parameters
                for (int global_var_id = 0; global_var_id < num_global_vars; global_var_id++) {
                    double stress_bc_var = bc_yaml[bc_id]["boundary_condition"]["stress_bc_global_vars"][global_var_id].As<double>();
                    

                    RUN({
                        tempStressBCGlobalVars(bc_id, global_var_id) = stress_bc_var;
                    });

                } // end loop over global vars
            } // end else if on stress_bc_global_vars
            else {
                std::cout << "ERROR: invalid input: " << a_word << std::endl;
                std::cout << "Valid options are: " << std::endl;
                for (const auto& element : str_bc_inps) {
                    std::cout << element << std::endl;
                }
                throw std::runtime_error("**** Boundary Conditions Not Understood ****");
            }
        } // end for words in boundary conditions


        // add checks for velocity vs time boundary condition


    } // end loop over BCs specified



     // allocate ragged right memory to hold the model global variables
    BoundaryConditions.velocity_bc_global_vars = RaggedRightArrayKokkos <double> (BoundaryConditions.num_velocity_bc_global_vars, "BoundaryConditions.velocity_bc_global_vars");

    BoundaryConditions.temperature_bc_global_vars = RaggedRightArrayKokkos <double> (BoundaryConditions.num_temperature_bc_global_vars, "BoundaryConditions.temperature_bc_global_vars");

    BoundaryConditions.stress_bc_global_vars = RaggedRightArrayKokkos <double> (BoundaryConditions.num_stress_bc_global_vars, "BoundaryConditions.stress_bc_global_vars");
   

    // ... allocate other bc global vars here

    // save the global variables
    FOR_ALL(bc_id, 0, num_bcs, {

        for (size_t var_lid = 0; var_lid < BoundaryConditions.num_velocity_bc_global_vars(bc_id); var_lid++){
            BoundaryConditions.velocity_bc_global_vars(bc_id, var_lid) = tempVelocityBCGlobalVars(bc_id, var_lid);
        } // end for eos var_lid

        for (size_t var_lid = 0; var_lid < BoundaryConditions.num_temperature_bc_global_vars(bc_id); var_lid++){
            BoundaryConditions.temperature_bc_global_vars(bc_id, var_lid) = tempTemperatureBCGlobalVars(bc_id, var_lid);
        } // end for var_lid
      
        for (size_t var_lid=0; var_lid<BoundaryConditions.num_stress_bc_global_vars(bc_id); var_lid++){
            BoundaryConditions.stress_bc_global_vars(bc_id, var_lid) = tempStressBCGlobalVars(bc_id, var_lid);
        } // end for eos var_lid


        // ... add other bc global vars here

    }); // end for loop over boundary conditions

    // copy the enum values to the host 
    BoundaryConditions.BoundaryConditionEnums.update_host();



} // end of function to parse bdy conditions