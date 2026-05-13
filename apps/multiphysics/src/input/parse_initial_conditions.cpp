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

#include "string_utils.hpp"

#include "Yaml.hpp"
#include "matar.h"
#include "parse_tools.hpp"
#include "simulation_parameters.hpp"
#include "parse_initial_conditions.hpp"

#include "state.hpp"


// remember: 
//    fill_gauss_state is in state.hpp
//    fill_node_state is in state.hpp

// =================================================================================
//    Parse region ICs
// =================================================================================
void parse_initial_conditions(Yaml::Node& root, 
                              CArrayKokkos<RegionICs_t>& region_ics, 
                              std::vector <fill_gauss_state>& fill_gauss_states,
                              std::vector <fill_node_state>& fill_node_states)

{

    // node state to setup, this array is built based on the fills
    // possible node states:
    //     fill_node_state::velocity,
    //     fill_node_state::temperature
    //  ...
    fill_node_states = std::vector <fill_node_state> (0);

    // mat_pt state to setup, this array is built based on the fills
    // possible gauss states:
    //    fill_gauss_state::density
    //    fill_gauss_state::specific_internal_energy
    //    fill_gauss_state::thermal_conductivity
    //    fill_gauss_state::specific_heat
    // ...
    fill_gauss_states = std::vector <fill_gauss_state> (0);

    // the above fill state vectors are also used with checking solver allocations



    Yaml::Node& region_yaml = root["initial_conditions"];

    size_t num_intial_conditions = region_yaml.Size();

    printf("Number of initial conditions = %zu \n", num_intial_conditions);
   
    region_ics = CArrayKokkos<RegionICs_t>(num_intial_conditions, "sim_param.initial_condition_setup.region_ics");


    // loop over the fill ics specified in yaml file
    for (int ic_id = 0; ic_id < num_intial_conditions; ic_id++) {

        // read the variables names
        Yaml::Node& inps_yaml = root["initial_conditions"][ic_id]["initial_condition"];

        // get the intial condition variables names set by the user
        std::vector<std::string> user_str_ics_inps;

        // extract words from the input file and validate they are correct
        validate_inputs(inps_yaml, user_str_ics_inps, str_ics_inps, ics_required_inps);


        // loop over the words in the material input definition
        for (auto& a_word : user_str_ics_inps) {

            if (a_word.compare("region_id") == 0) {
                int reg_id = root["initial_conditions"][ic_id]["initial_condition"][a_word].As<int>();

                RUN({
                    region_ics(ic_id).region_id = reg_id;
                });
            } // material_id
            else if (a_word.compare("material_id") == 0) {
                int mat_id = root["initial_conditions"][ic_id]["initial_condition"][a_word].As<int>();

                RUN({
                    region_ics(ic_id).material_id = mat_id;
                });
            } // material_id
            else if (a_word.compare("density") == 0) {

                // check to see if density enum was saved
                bool store = true;
                for (auto field : fill_gauss_states){
                    if (field == fill_gauss_state::density){store = false;}
                }
                // store density name if it has not been stored already
                if(store){
                    fill_gauss_states.push_back(fill_gauss_state::density);
                }
                
                // -----
                // loop over the sub fields under den
                // -----
                Yaml::Node& inps_subfields_yaml = root["initial_conditions"][ic_id]["initial_condition"]["density"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_ics_den_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_ics_den_inputs, str_ics_den_inps, ics_den_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_ics_den_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        // density
                        double value = root["initial_conditions"][ic_id]["initial_condition"]["density"]["value"].As<double>();

                        // check for a valid density, and then save it if it is
                        if (value < 0.0) {
                            std::cout << "ERROR: density is negative: " << value << std::endl;
                        }

                        RUN({
                        region_ics(ic_id).den = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["initial_conditions"][ic_id]["initial_condition"]["density"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., initial_conditions::uniform 
                            switch(scalar_ics_type_map[type]){

                                case initial_conditions::uniform:
                                    std::cout << "Setting density initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_ics(ic_id).den_field = initial_conditions::uniform;
                                    });
                                    break;

                                case initial_conditions::tgVortexScalar:
                                    std::cout << "Setting density initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_ics(ic_id).den_field = initial_conditions::tgVortexScalar;
                                    });
                                    break;

                                case initial_conditions::noICsScalar:
                                    std::cout << "Setting density initial conditions type to no density" << std::endl;
                                    RUN({ 
                                        region_ics(ic_id).den_field = initial_conditions::noICsScalar;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_ics(ic_id).den_field = initial_conditions::noICsScalar;
                                    });

                                    std::cout << "ERROR: No valid density intial conditions type input " << std::endl;
                                    std::cout << "Valid IC types are: " << std::endl;
                                    
                                    for (const auto& pair : scalar_ics_type_map) {
                                        std::cout << pair.second << std::endl;
                                    }

                                    throw std::runtime_error("**** density Initial Conditions Type Not Understood ****");
                                    break;
                            } // end switch

                        }
                        else{
                            std::cout << "ERROR: invalid input: " << type << std::endl;
                            throw std::runtime_error("**** density IC Not Understood ****");
                        } // end if on density type
                        
                    } // end if on density type
                    else {
                        std::cout << "ERROR: invalid input: " << a_subfield_word << std::endl;
                        std::cout << "Valid options are: " << std::endl;
                        for (const auto& element : str_ics_den_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** ICS density Inputs Not Understood ****");
                    } // end if on all subfields under density

                } // end for loop over text

            } // den
            else if (a_word.compare("specific_internal_energy") == 0) {
                // specific internal energy

                // check to see if specific_internal_energy enum was saved
                bool store = true;
                for (auto field : fill_gauss_states){
                    if (field == fill_gauss_state::specific_internal_energy){store = false;}
                }
                // store specific_internal_energy name if it has not been stored already
                if(store){
                    fill_gauss_states.push_back(fill_gauss_state::specific_internal_energy);
                }

                // -----
                // loop over the sub fields under sie
                // -----
                Yaml::Node& inps_subfields_yaml = root["initial_conditions"][ic_id]["initial_condition"]["specific_internal_energy"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_ics_sie_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_ics_sie_inputs, str_ics_sie_inps, ics_sie_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_ics_sie_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        // specific_internal_energy value
                        double value = root["initial_conditions"][ic_id]["initial_condition"]["specific_internal_energy"]["value"].As<double>();
                        RUN({
                        region_ics(ic_id).sie = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["initial_conditions"][ic_id]["initial_condition"]["specific_internal_energy"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., initial_conditions::uniform 
                            switch(scalar_ics_type_map[type]){

                                case initial_conditions::uniform:
                                    std::cout << "Setting specific_internal_energy initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_ics(ic_id).sie_field = initial_conditions::uniform;
                                    });
                                    break;

                                case initial_conditions::tgVortexScalar:
                                    std::cout << "Setting specific_internal_energy initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_ics(ic_id).sie_field = initial_conditions::tgVortexScalar;
                                    });
                                    break;

                                case initial_conditions::noICsScalar:
                                    std::cout << "Setting specific_internal_energy initial conditions type to no specific_internal_energy" << std::endl;
                                    RUN({ 
                                        region_ics(ic_id).sie_field = initial_conditions::noICsScalar;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_ics(ic_id).sie_field = initial_conditions::noICsScalar;
                                    });

                                    std::cout << "ERROR: No valid specific_internal_energy intial conditions type input " << std::endl;
                                    std::cout << "Valid IC types are: " << std::endl;
                                    
                                    for (const auto& pair : scalar_ics_type_map) {
                                        std::cout << pair.second << std::endl;
                                    }

                                    throw std::runtime_error("**** specific_internal_energy Initial Conditions Type Not Understood ****");
                                    break;
                            } // end switch
                        }
                        else{
                            std::cout << "ERROR: invalid input: " << type << std::endl;
                            throw std::runtime_error("**** specific_internal_energy IC Not Understood ****");
                        } // end if on specific_internal_energy type
                        
                    } // end if on specific_internal_energy type
                    else {
                        std::cout << "ERROR: invalid input: " << a_subfield_word << std::endl;
                        std::cout << "Valid options are: " << std::endl;
                        for (const auto& element : str_ics_sie_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** ICS specific_internal_energy Inputs Not Understood ****");
                    } // end if on all subfields under specific_internal_energy

                } // end for loop over text

            } // sie
            else if (a_word.compare("internal_energy") == 0) {
                // extensive internal energy

                // check to see if internal_energy enum was saved
                bool store = true;
                for (auto field : fill_gauss_states){
                    if (field == fill_gauss_state::internal_energy){store = false;}
                }
                // store internal_energy name if it has not been stored already
                if(store){
                    fill_gauss_states.push_back(fill_gauss_state::internal_energy);
                }

                // -----
                // loop over the sub fields under internal_energy
                // -----
                Yaml::Node& inps_subfields_yaml = root["initial_conditions"][ic_id]["initial_condition"]["internal_energy"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_ics_ie_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_ics_ie_inputs, str_ics_ie_inps, ics_ie_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_ics_ie_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        // extensive internal_energy
                        double value = root["initial_conditions"][ic_id]["initial_condition"]["internal_energy"]["value"].As<double>();
                        RUN({
                        region_ics(ic_id).ie = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["initial_conditions"][ic_id]["initial_condition"]["internal_energy"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., initial_conditions::uniform 
                            switch(scalar_ics_type_map[type]){

                                case initial_conditions::uniform:
                                    std::cout << "Setting extensive internal energy initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_ics(ic_id).ie_field = initial_conditions::uniform;
                                    });
                                    break;

                                case initial_conditions::tgVortexScalar:
                                    std::cout << "Setting extensive  internal energy initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_ics(ic_id).ie_field = initial_conditions::tgVortexScalar;
                                    });
                                    break;

                                case initial_conditions::noICsScalar:
                                    std::cout << "Setting extensive internal energy initial conditions type to no internal energy" << std::endl;
                                    RUN({ 
                                        region_ics(ic_id).ie_field = initial_conditions::noICsScalar;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_ics(ic_id).ie_field = initial_conditions::noICsScalar;
                                    });

                                    std::cout << "ERROR: No valid extensive internal energy intial conditions type input " << std::endl;
                                    std::cout << "Valid IC types are: " << std::endl;
                                    
                                    for (const auto& pair : scalar_ics_type_map) {
                                        std::cout << pair.second << std::endl;
                                    }

                                    throw std::runtime_error("**** internal energy Initial Conditions Type Not Understood ****");
                                    break;
                            } // end switch
                        }
                        else{
                            std::cout << "ERROR: invalid input: " << type << std::endl;
                            throw std::runtime_error("**** internal energy IC Not Understood ****");
                        } // end if on internal energy type
                        
                    } // end if on internal energy type
                    else {
                        std::cout << "ERROR: invalid input: " << a_subfield_word << std::endl;
                        std::cout << "Valid options are: " << std::endl;
                        for (const auto& element : str_ics_ie_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** ICS internal energy Inputs Not Understood ****");
                    } // end if on all subfields under internal energy

                } // end for loop over text

            } // ie            
            else if (a_word.compare("specific_heat") == 0) {

                // check to see if specific_heat enum was saved
                bool store = true;
                for (auto field : fill_gauss_states){
                    if (field == fill_gauss_state::specific_heat){store = false;}
                }
                // store specific_heat name if it has not been stored already
                if(store){
                    fill_gauss_states.push_back(fill_gauss_state::specific_heat);
                }

                // -----
                // loop over the sub fields under specific_heat
                // -----
                Yaml::Node& inps_subfields_yaml = root["initial_conditions"][ic_id]["initial_condition"]["specific_heat"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_ics_specific_heat_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_ics_specific_heat_inputs, str_ics_specific_heat_inps, ics_specific_heat_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_ics_specific_heat_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        // x-component of specific_heat
                        double value = root["initial_conditions"][ic_id]["initial_condition"]["specific_heat"]["value"].As<double>();
                        RUN({
                        region_ics(ic_id).specific_heat = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["initial_conditions"][ic_id]["initial_condition"]["specific_heat"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., initial_conditions::uniform 
                            switch(scalar_ics_type_map[type]){

                                case initial_conditions::uniform:
                                    std::cout << "Setting specific_heat initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_ics(ic_id).specific_heat_field = initial_conditions::uniform;
                                    });
                                    break;

                                case initial_conditions::tgVortexScalar:
                                    std::cout << "Setting specific_heat initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_ics(ic_id).specific_heat_field = initial_conditions::tgVortexScalar;
                                    });
                                    break;

                                case initial_conditions::noICsScalar:
                                    std::cout << "Setting specific_heat initial conditions type to no specific_heat" << std::endl;
                                    RUN({ 
                                        region_ics(ic_id).specific_heat_field = initial_conditions::noICsScalar;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_ics(ic_id).specific_heat_field = initial_conditions::noICsScalar;
                                    });

                                    std::cout << "ERROR: No valid specific_heat intial conditions type input " << std::endl;
                                    std::cout << "Valid IC types are: " << std::endl;
                                    
                                    for (const auto& pair : scalar_ics_type_map) {
                                        std::cout << pair.second << std::endl;
                                    }

                                    throw std::runtime_error("**** specific_heat Initial Conditions Type Not Understood ****");
                                    break;
                            } // end switch

                        }
                        else{
                            std::cout << "ERROR: invalid input: " << type << std::endl;
                            throw std::runtime_error("**** specific_heat IC Not Understood ****");
                        } // end if on specific_heat type
                        
                    } // end if on specific_heat type
                    else {
                        std::cout << "ERROR: invalid input: " << a_subfield_word << std::endl;
                        std::cout << "Valid options are: " << std::endl;
                        for (const auto& element : str_ics_specific_heat_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** ICS specific_heat Inputs Not Understood ****");
                    } // end if on all subfields under specific_heat

                } // end for loop over text

            } // specific_heat
            else if (a_word.compare("thermal_conductivity") == 0) {

                // check to see if thermal_conductivity enum was saved
                bool store = true;
                for (auto field : fill_gauss_states){
                    if (field == fill_gauss_state::thermal_conductivity){store = false;}
                }
                // store thermal_conductivity name if it has not been stored already
                if(store){
                    fill_gauss_states.push_back(fill_gauss_state::thermal_conductivity);
                }

                // -----
                // loop over the sub fields under thermal_conductivity
                // -----
                Yaml::Node& inps_subfields_yaml = root["initial_conditions"][ic_id]["initial_condition"]["thermal_conductivity"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_ics_thermal_conductivity_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_ics_thermal_conductivity_inputs, str_ics_thermal_conductivity_inps, ics_thermal_conductivity_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_ics_thermal_conductivity_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        // thermal_conductivity
                        double value = root["initial_conditions"][ic_id]["initial_condition"]["thermal_conductivity"]["value"].As<double>();

                        RUN({
                        region_ics(ic_id).thermal_conductivity = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["initial_conditions"][ic_id]["initial_condition"]["thermal_conductivity"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., initial_conditions::uniform 
                            switch(scalar_ics_type_map[type]){

                                case initial_conditions::uniform:
                                    std::cout << "Setting thermal_conductivity initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_ics(ic_id).thermal_conductivity_field = initial_conditions::uniform;
                                    });
                                    break;

                                case initial_conditions::tgVortexScalar:
                                    std::cout << "Setting thermal_conductivity initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_ics(ic_id).thermal_conductivity_field = initial_conditions::tgVortexScalar;
                                    });
                                    break;

                                case initial_conditions::noICsScalar:
                                    std::cout << "Setting thermal_conductivity initial conditions type to no thermal_conductivity" << std::endl;
                                    RUN({ 
                                        region_ics(ic_id).thermal_conductivity_field = initial_conditions::noICsScalar;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_ics(ic_id).thermal_conductivity_field = initial_conditions::noICsScalar;
                                    });

                                    std::cout << "ERROR: No valid thermal_conductivity intial conditions type input " << std::endl;
                                    std::cout << "Valid IC types are: " << std::endl;
                                    
                                    for (const auto& pair : scalar_ics_type_map) {
                                        std::cout << pair.second << std::endl;
                                    }

                                    throw std::runtime_error("**** thermal_conductivity Initial Conditions Type Not Understood ****");
                                    break;
                            } // end switch

                        }
                        else{
                            std::cout << "ERROR: invalid input: " << type << std::endl;
                            throw std::runtime_error("**** thermal_conductivity IC Not Understood ****");
                        } // end if on thermal_conductivity type
                        
                    } // end if on thermal_conductivity type
                    else {
                        std::cout << "ERROR: invalid input: " << a_subfield_word << std::endl;
                        std::cout << "Valid options are: " << std::endl;
                        for (const auto& element : str_ics_thermal_conductivity_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** ICS thermal_conductivity Inputs Not Understood ****");
                    } // end if on all subfields under thermal_conductivity

                } // end for loop over text

            } // thermal_conductivity
            else if (a_word.compare("material_volume_fraction") == 0){

                // always build, so no need to add this varname to guass_point

                // -----
                // loop over the sub fields under volfrac
                // -----
                Yaml::Node& inps_subfields_yaml = root["initial_conditions"][ic_id]["initial_condition"]["material_volume_fraction"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_ics_mat_volfrac_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_ics_mat_volfrac_inputs, str_ics_mat_volfrac_inps, ics_volfrac_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_ics_mat_volfrac_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        // volfrac value or the intercept if linear variation
                        double value = root["initial_conditions"][ic_id]["initial_condition"]["material_volume_fraction"]["value"].As<double>();
      
                        RUN({
                            region_ics(ic_id).volfrac = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("slope") == 0) {
                        // volfrac slope
                        double slope = root["initial_conditions"][ic_id]["initial_condition"]["material_volume_fraction"]["slope"].As<double>();
      
                        RUN({
                            region_ics(ic_id).volfrac_slope = slope;
                        });
                    } // slope
                    else if (a_subfield_word.compare("origin") == 0) {
                        std::string origin = root["initial_conditions"][ic_id]["initial_condition"]["material_volume_fraction"]["origin"].As<std::string>();

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
                            region_ics(ic_id).volfrac_origin[0] = x1;
                            region_ics(ic_id).volfrac_origin[1] = y1;
                            region_ics(ic_id).volfrac_origin[2] = z1;
                        });
                    } // origin
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["initial_conditions"][ic_id]["initial_condition"]["material_volume_fraction"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., initial_conditions::uniform 
                            switch(scalar_ics_type_map[type]){

                                case initial_conditions::uniform:
                                    std::cout << "Setting volfrac initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_ics(ic_id).volfrac_field = initial_conditions::uniform;
                                    });
                                    break;

                                case initial_conditions::radialScalar:
                                    std::cout << "Setting volfrac initial conditions type to radial scalar " << std::endl;
                                    RUN({
                                        region_ics(ic_id).volfrac_field = initial_conditions::radialScalar;
                                    });
                                    break;

                                case initial_conditions::sphericalScalar:
                                    std::cout << "Setting volfrac initial conditions type to spherical scalar " << std::endl;
                                    RUN({
                                        region_ics(ic_id).volfrac_field = initial_conditions::sphericalScalar;
                                    });
                                    break;

                                case initial_conditions::xlinearScalar:
                                    std::cout << "Setting volfrac initial conditions type to xlinearScalar " << std::endl;
                                    RUN({
                                        region_ics(ic_id).volfrac_field = initial_conditions::xlinearScalar;
                                    });
                                    break;
                                
                                case initial_conditions::ylinearScalar:
                                    std::cout << "Setting volfrac initial conditions type to ylinearScalar " << std::endl;
                                    RUN({
                                        region_ics(ic_id).volfrac_field = initial_conditions::ylinearScalar;
                                    });
                                    break;
                                
                                case initial_conditions::zlinearScalar:
                                    std::cout << "Setting volfrac initial conditions type to zlinearScalar " << std::endl;
                                    RUN({
                                        region_ics(ic_id).volfrac_field = initial_conditions::zlinearScalar;
                                    });
                                    break;

                                case initial_conditions::tgVortexScalar:
                                    std::cout << "Setting volfrac initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_ics(ic_id).volfrac_field = initial_conditions::tgVortexScalar;
                                    });
                                    break;

                                case initial_conditions::noICsScalar:
                                    std::cout << "Default Volume Fraction Used:" << std::endl;
                                    std::cout << "Setting volume fraction to uniform field with a value equal to 1" << std::endl;
                                    RUN({ 
                                        region_ics(ic_id).volfrac_field = initial_conditions::uniform;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_ics(ic_id).volfrac_field = initial_conditions::noICsScalar;
                                    });

                                    std::cout << "ERROR: No valid material volume fraction intial conditions type input " << std::endl;
                                    std::cout << "Valid IC types are: " << std::endl;
                                    
                                    for (const auto& pair : scalar_ics_type_map) {
                                        std::cout << pair.second << std::endl;
                                    }

                                    throw std::runtime_error("**** Material Volume Fraction Initial Conditions Type Not Understood ****");
                                    break;
                            } // end switch

                        }
                        else{
                            std::cout << "ERROR: invalid input: " << type << std::endl;
                            throw std::runtime_error("**** Material Volume Fraction IC Not Understood ****");
                        } // end if on Volume Fraction type
                        
                    } // end if on Volume Fraction type
                    else {
                        std::cout << "ERROR: invalid input: " << a_subfield_word << std::endl;
                        std::cout << "Valid options are: " << std::endl;
                        for (const auto& element : str_ics_mat_volfrac_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** ICS Material Volume Fraction Inputs Not Understood ****");
                    } // end if on all subfields under Volume Fraction

                } // end for loop over text
            }
            // ----------  nodal variables ----------
            else if (a_word.compare("temperature") == 0) {

                // check to see if temperature enum was saved
                bool store = true;
                for (auto field : fill_node_states){
                    if (field == fill_node_state::temperature){store = false;}
                }
                // store temperature name if it has not been stored already
                if(store){
                    fill_node_states.push_back(fill_node_state::temperature);
                }

                // -----
                // loop over the sub fields under temperature
                // -----
                Yaml::Node& inps_subfields_yaml = root["initial_conditions"][ic_id]["initial_condition"]["temperature"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_ics_temperature_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_ics_temperature_inputs, 
                                str_ics_temperature_inps, ics_temperature_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_ics_temperature_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        //temperature
                        double value = root["initial_conditions"][ic_id]["initial_condition"]["temperature"]["value"].As<double>();
                        RUN({
                        region_ics(ic_id).temperature = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["initial_conditions"][ic_id]["initial_condition"]["temperature"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., initial_conditions::uniform 
                            switch(scalar_ics_type_map[type]){

                                case initial_conditions::uniform:
                                    std::cout << "Setting temperature initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_ics(ic_id).temperature_field = initial_conditions::uniform;
                                    });
                                    break;

                                case initial_conditions::tgVortexScalar:
                                    std::cout << "Setting temperature initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_ics(ic_id).temperature_field = initial_conditions::tgVortexScalar;
                                    });
                                    break;

                                case initial_conditions::noICsScalar:
                                    std::cout << "Setting temperature initial conditions type to no temperature" << std::endl;
                                    RUN({ 
                                        region_ics(ic_id).temperature_field = initial_conditions::noICsScalar;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_ics(ic_id).temperature_field = initial_conditions::noICsScalar;
                                    });

                                    std::cout << "ERROR: No valid temperature intial conditions type input " << std::endl;
                                    std::cout << "Valid IC types are: " << std::endl;
                                    
                                    for (const auto& pair : scalar_ics_type_map) {
                                        std::cout << pair.second << std::endl;
                                    }

                                    throw std::runtime_error("**** temperature Initial Conditions Type Not Understood ****");
                                    break;
                            } // end switch

                        }
                        else{
                            std::cout << "ERROR: invalid input: " << type << std::endl;
                            throw std::runtime_error("**** temperature IC Not Understood ****");
                        } // end if on temperature type
                        
                    } // end if on temperature type
                    else {
                        std::cout << "ERROR: invalid input: " << a_subfield_word << std::endl;
                        std::cout << "Valid options are: " << std::endl;
                        for (const auto& element : str_ics_temperature_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** ICS temperature Inputs Not Understood ****");
                    } // end if on all subfields under temperature

                } // end for loop over text

            } // temperature
            else if (a_word.compare("level_set") == 0) {

                // check to see if level set enum was saved
                bool store = true;
                for (auto field : fill_gauss_states){
                    if (field == fill_gauss_state::level_set){store = false;}
                }
                // store level set name if it has not been stored already
                if(store){
                    fill_gauss_states.push_back(fill_gauss_state::level_set);
                }
                
                // -----
                // loop over the sub fields under level set
                // -----
                Yaml::Node& inps_subfields_yaml = root["initial_conditions"][ic_id]["initial_condition"]["level_set"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_ics_level_set_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_ics_level_set_inputs, str_ics_level_set_inps, ics_level_set_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_ics_level_set_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        // level set
                        double value = root["initial_conditions"][ic_id]["initial_condition"]["level_set"]["value"].As<double>();

                        RUN({
                            region_ics(ic_id).level_set = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("slope") == 0) {
                        // volfrac slope
                        double slope = root["initial_conditions"][ic_id]["initial_condition"]["level_set"]["slope"].As<double>();
        
                        RUN({
                            region_ics(ic_id).level_set_slope = slope;
                        });
                    } // slope
                    else if (a_subfield_word.compare("origin") == 0) {
                        std::string origin = root["initial_conditions"][ic_id]["initial_condition"]["level_set"]["origin"].As<std::string>();

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
                            region_ics(ic_id).level_set_origin[0] = x1;
                            region_ics(ic_id).level_set_origin[1] = y1;
                            region_ics(ic_id).level_set_origin[2] = z1;
                        });
                    } // origin
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["initial_conditions"][ic_id]["initial_condition"]["level_set"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., initial_conditions::uniform 
                            switch(scalar_ics_type_map[type]){

                                case initial_conditions::uniform:
                                    std::cout << "Setting level set initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_ics(ic_id).level_set_field = initial_conditions::uniform;
                                    });
                                    break;
                                case initial_conditions::radialScalar:
                                    std::cout << "Setting level set initial conditions type to radial scalar " << std::endl;
                                    RUN({
                                        region_ics(ic_id).level_set_field = initial_conditions::radialScalar;
                                    });
                                    break;
                                case initial_conditions::sphericalScalar:
                                    std::cout << "Setting level set initial conditions type to spherical scalar " << std::endl;
                                    RUN({
                                        region_ics(ic_id).level_set_field = initial_conditions::sphericalScalar;
                                    });
                                    break;
                                case initial_conditions::xlinearScalar:
                                    std::cout << "Setting volfrac initial conditions type to xlinearScalar " << std::endl;
                                    RUN({
                                        region_ics(ic_id).level_set_field = initial_conditions::xlinearScalar;
                                    });
                                    break;
                                
                                case initial_conditions::ylinearScalar:
                                    std::cout << "Setting volfrac initial conditions type to ylinearScalar " << std::endl;
                                    RUN({
                                        region_ics(ic_id).level_set_field = initial_conditions::ylinearScalar;
                                    });
                                    break;
                                
                                case initial_conditions::zlinearScalar:
                                    std::cout << "Setting volfrac initial conditions type to zlinearScalar " << std::endl;
                                    RUN({
                                        region_ics(ic_id).level_set_field = initial_conditions::zlinearScalar;
                                    });
                                    break;
                                case initial_conditions::tgVortexScalar:
                                    std::cout << "Setting level set initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_ics(ic_id).level_set_field = initial_conditions::tgVortexScalar;
                                    });
                                    break;

                                case initial_conditions::noICsScalar:
                                    std::cout << "Setting level set initial conditions type to no level set" << std::endl;
                                    RUN({ 
                                        region_ics(ic_id).level_set_field = initial_conditions::noICsScalar;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_ics(ic_id).level_set_field = initial_conditions::noICsScalar;
                                    });

                                    std::cout << "ERROR: No valid level set intial conditions type input " << std::endl;
                                    std::cout << "Valid IC types are: " << std::endl;
                                    
                                    for (const auto& pair : scalar_ics_type_map) {
                                        std::cout << pair.second << std::endl;
                                    }

                                    throw std::runtime_error("**** level set Initial Conditions Type Not Understood ****");
                                    break;
                            } // end switch

                        }
                        else{
                            std::cout << "ERROR: invalid input: " << type << std::endl;
                            throw std::runtime_error("**** level set IC Not Understood ****");
                        } // end if on level set type
                        
                    } // end if on level set type
                    else {
                        std::cout << "ERROR: invalid input: " << a_subfield_word << std::endl;
                        std::cout << "Valid options are: " << std::endl;
                        for (const auto& element : str_ics_level_set_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** ICS level set Inputs Not Understood ****");
                    } // end if on all subfields under level set

                } // end for loop over text

            } // level set
            else if (a_word.compare("velocity") == 0) {

                // check to see if velocity enum was saved
                bool store = true;
                for (auto field : fill_node_states){
                    if (field == fill_node_state::velocity){store = false;}
                }
                // store velocity name if it has not been stored already
                if(store){
                    fill_node_states.push_back(fill_node_state::velocity);
                }

                // -----
                // loop over the sub fields under velocity
                // -----
                Yaml::Node& inps_subfields_yaml = root["initial_conditions"][ic_id]["initial_condition"]["velocity"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_ics_vel_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_ics_vel_inputs, str_ics_vel_inps, ics_vel_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_ics_vel_inputs){ 

                    if (a_subfield_word.compare("u") == 0) {
                        // x-component of velocity
                        double u = root["initial_conditions"][ic_id]["initial_condition"]["velocity"]["u"].As<double>();

                        RUN({
                        region_ics(ic_id).u = u;
                        });
                    } // u
                    else if (a_subfield_word.compare("v") == 0) {
                        // y-component of velocity
                        double v = root["initial_conditions"][ic_id]["initial_condition"]["velocity"]["v"].As<double>();

                        RUN({
                            region_ics(ic_id).v = v;
                        });
                    } // v
                    else if (a_subfield_word.compare("w") == 0) {
                        // z-component of velocity

                        double w = root["initial_conditions"][ic_id]["initial_condition"]["velocity"]["w"].As<double>();

                        RUN({
                            region_ics(ic_id).w = w;
                        });
                    } // w
                    else if (a_subfield_word.compare("speed") == 0) {
                        double speed = root["initial_conditions"][ic_id]["initial_condition"]["velocity"]["speed"].As<double>();

                        RUN({
                            region_ics(ic_id).speed = speed;
                        });
                    } // speed
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["initial_conditions"][ic_id]["initial_condition"]["velocity"]["type"].As<std::string>();

                        // set the volume tag type
                        if (vector_ics_type_map.find(type) != vector_ics_type_map.end()) {
                        
                            // vector_ics_type_map[type] returns enum value, e.g., initial_conditions::velocity 
                            switch(vector_ics_type_map[type]){

                                case initial_conditions::stationary:
                                    std::cout << "Setting velocity initial conditions type to static " << std::endl;
                                    RUN({
                                        region_ics(ic_id).vel_field = initial_conditions::stationary;
                                    });
                                    break;

                                case initial_conditions::cartesian:
                                    std::cout << "Setting velocity initial conditions type to cartesian " << std::endl;
                                    RUN({
                                        region_ics(ic_id).vel_field = initial_conditions::cartesian;
                                    });
                                    break;

                                case initial_conditions::radialVec:
                                    std::cout << "Setting velocity initial conditions type to radial " << std::endl;
                                    RUN({
                                        region_ics(ic_id).vel_field = initial_conditions::radialVec;
                                    });
                                    break;

                                case initial_conditions::sphericalVec:
                                    std::cout << "Setting velocity initial conditions type to spherical " << std::endl;
                                    RUN({
                                        region_ics(ic_id).vel_field = initial_conditions::sphericalVec;
                                    });
                                    break;

                                case initial_conditions::radialLinearVec:
                                    std::cout << "Setting velocity initial conditions type to radial_linear " << std::endl;
                                    RUN({
                                        region_ics(ic_id).vel_field = initial_conditions::radialLinearVec;
                                    });
                                    break;

                                case initial_conditions::sphericalLinearVec:
                                    std::cout << "Setting velocity initial conditions type to spherical_linear " << std::endl;
                                    RUN({
                                        region_ics(ic_id).vel_field = initial_conditions::sphericalLinearVec;
                                    });
                                    break;

                                case initial_conditions::tgVortexVec:
                                    std::cout << "Setting velocity initial conditions type to tg_vortex " << std::endl;
                                    RUN({
                                        region_ics(ic_id).vel_field = initial_conditions::tgVortexVec;
                                    });
                                    break;

                                case initial_conditions::noICsVec:
                                    std::cout << "Setting velocity initial conditions type to no velocity" << std::endl;
                                    RUN({ 
                                        region_ics(ic_id).vel_field = initial_conditions::noICsVec;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_ics(ic_id).vel_field = initial_conditions::noICsVec;
                                    });

                                    std::cout << "ERROR: No valid velocity intial conditions type input " << std::endl;
                                    std::cout << "Valid IC types are: " << std::endl;
                                    
                                    for (const auto& pair : vector_ics_type_map) {
                                        std::cout << pair.second << std::endl;
                                    }

                                    throw std::runtime_error("**** Velocity Initial Conditions Type Not Understood ****");
                                    break;
                            } // end switch

                        }
                        else{
                            std::cout << "ERROR: invalid input: " << type << std::endl;
                            throw std::runtime_error("**** Velocity IC Not Understood ****");
                        } // end if on velocity type
                        
                    } // end if on velocity type
                    else {
                        std::cout << "ERROR: invalid input: " << a_subfield_word << std::endl;
                        std::cout << "Valid options are: " << std::endl;
                        for (const auto& element : str_ics_vel_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** ICS Velocity Inputs Not Understood ****");
                    } // end if on all subfields under velocity

                } // end for loop over text
            } // end if on velocity
            //
            //
            else {
                std::cout << "ERROR: invalid input: " << a_word << std::endl;
                std::cout << "Valid options are: " << std::endl;
                for (const auto& element : str_ics_inps) {
                    std::cout << element << std::endl;
                }
                throw std::runtime_error("**** ICS Not Understood ****");
            }
        } // end for words in fill region


        // -----------------------------------------------
        // check for consistency in input settings

        if(fill_gauss_states.size() == 0){
            throw std::runtime_error("**** No Initial Conditions Set, Please Specify Fields in the Region Fills ****");
        }

        // -----------------------------------------------


    } // end loop over regions
} // end of function to parse region
