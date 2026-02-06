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

#include "Yaml.hpp"
#include "matar.h"
#include "parse_tools.hpp"
#include "simulation_parameters.h"
#include "parse_regions.hpp"

#include "state.h"




// =================================================================================
//    Parse Fill regions
// =================================================================================
void parse_regions(Yaml::Node& root, 
                   DCArrayKokkos<size_t>& reg_fills_in_solver,  
                   DCArrayKokkos<size_t>& num_reg_fills_in_solver,
                   CArrayKokkos<RegionFill_t>& region_fills, 
                   CArray<RegionFill_host_t>&  region_fills_host,
                   std::vector <fill_gauss_state>& fill_gauss_states,
                   std::vector <fill_node_state>& fill_node_states,
                   const size_t num_solvers)

{
    // allocate memory
    num_reg_fills_in_solver = DCArrayKokkos<size_t>(num_solvers, "sim_param.region_setup.num_reg_fills_in_solver");
    num_reg_fills_in_solver.set_values(0);
    Kokkos::fence();
    num_reg_fills_in_solver.update_host(); // initiallizing host side to 0


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



    Yaml::Node& region_yaml = root["regions"];

    size_t num_regions = region_yaml.Size();

    // reg_fills_in_solver(solver_id, fill_lid) = fill_id
    reg_fills_in_solver = DCArrayKokkos<size_t>(num_solvers, num_regions, "sim_param.region_setup.reg_fills_in_solver");
    region_fills = CArrayKokkos<RegionFill_t>(num_regions , "sim_param.region_setup.region_fills");
    region_fills_host = CArray<RegionFill_host_t>(num_regions); 

    // loop over the fill regions specified
    for (int reg_id = 0; reg_id < num_regions; reg_id++) {
        // read the variables names
        Yaml::Node& inps_yaml = root["regions"][reg_id]["region"];

        // get the material variables names set by the user
        std::vector<std::string> user_str_region_inps;

        // extract words from the input file and validate they are correct
        validate_inputs(inps_yaml, user_str_region_inps, str_region_inps, region_required_inps);

        // loop over the words in the material input definition
        for (auto& a_word : user_str_region_inps) {


            Yaml::Node& material_inps_yaml = root["regions"][reg_id]["region"][a_word];

            // set the values
            if (a_word.compare("solver_id") == 0) {
                int solver_id = root["regions"][reg_id]["region"][a_word].As<int>();
                
                // get the local id for filling this region
                size_t fill_lid = num_reg_fills_in_solver.host(solver_id);

                // save the fill_id, which is the reg_id
                reg_fills_in_solver.host(solver_id, fill_lid) = reg_id;
                num_reg_fills_in_solver.host(solver_id) ++;

                RUN({
                    region_fills(reg_id).solver_id = solver_id;
                });
            } // mat_id
            else if (a_word.compare("material_id") == 0) {
                int id = root["regions"][reg_id]["region"][a_word].As<int>();

                RUN({
                    region_fills(reg_id).material_id = id;
                });
            } // mat_id
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
                Yaml::Node& inps_subfields_yaml = root["regions"][reg_id]["region"]["density"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_region_den_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_region_den_inputs, str_region_den_inps, region_den_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_region_den_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        // density
                        double value = root["regions"][reg_id]["region"]["density"]["value"].As<double>();

                        // check for a valid density, and then save it if it is
                        if (value < 0.0) {
                            std::cout << "ERROR: density is negative: " << value << std::endl;
                        }

                        RUN({
                        region_fills(reg_id).den = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["regions"][reg_id]["region"]["density"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., init_conds::uniform 
                            switch(scalar_ics_type_map[type]){

                                case init_conds::uniform:
                                    std::cout << "Setting density initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_fills(reg_id).den_field = init_conds::uniform;
                                    });
                                    break;

                                case init_conds::tgVortexScalar:
                                    std::cout << "Setting density initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_fills(reg_id).den_field = init_conds::tgVortexScalar;
                                    });
                                    break;

                                case init_conds::noICsScalar:
                                    std::cout << "Setting density initial conditions type to no density" << std::endl;
                                    RUN({ 
                                        region_fills(reg_id).den_field = init_conds::noICsScalar;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_fills(reg_id).den_field = init_conds::noICsScalar;
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
                        for (const auto& element : str_region_den_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** Region density Inputs Not Understood ****");
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
                Yaml::Node& inps_subfields_yaml = root["regions"][reg_id]["region"]["specific_internal_energy"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_region_sie_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_region_sie_inputs, str_region_sie_inps, region_sie_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_region_sie_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        // specific_internal_energy value
                        double value = root["regions"][reg_id]["region"]["specific_internal_energy"]["value"].As<double>();
                        RUN({
                        region_fills(reg_id).sie = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["regions"][reg_id]["region"]["specific_internal_energy"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., init_conds::uniform 
                            switch(scalar_ics_type_map[type]){

                                case init_conds::uniform:
                                    std::cout << "Setting specific_internal_energy initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_fills(reg_id).sie_field = init_conds::uniform;
                                    });
                                    break;

                                case init_conds::tgVortexScalar:
                                    std::cout << "Setting specific_internal_energy initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_fills(reg_id).sie_field = init_conds::tgVortexScalar;
                                    });
                                    break;

                                case init_conds::noICsScalar:
                                    std::cout << "Setting specific_internal_energy initial conditions type to no specific_internal_energy" << std::endl;
                                    RUN({ 
                                        region_fills(reg_id).sie_field = init_conds::noICsScalar;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_fills(reg_id).sie_field = init_conds::noICsScalar;
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
                        for (const auto& element : str_region_sie_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** Region specific_internal_energy Inputs Not Understood ****");
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
                Yaml::Node& inps_subfields_yaml = root["regions"][reg_id]["region"]["internal_energy"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_region_ie_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_region_ie_inputs, str_region_ie_inps, region_ie_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_region_ie_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        // extensive internal_energy
                        double value = root["regions"][reg_id]["region"]["internal_energy"]["value"].As<double>();
                        RUN({
                        region_fills(reg_id).ie = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["regions"][reg_id]["region"]["internal_energy"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., init_conds::uniform 
                            switch(scalar_ics_type_map[type]){

                                case init_conds::uniform:
                                    std::cout << "Setting extensive internal energy initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_fills(reg_id).ie_field = init_conds::uniform;
                                    });
                                    break;

                                case init_conds::tgVortexScalar:
                                    std::cout << "Setting extensive  internal energy initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_fills(reg_id).ie_field = init_conds::tgVortexScalar;
                                    });
                                    break;

                                case init_conds::noICsScalar:
                                    std::cout << "Setting extensive internal energy initial conditions type to no internal energy" << std::endl;
                                    RUN({ 
                                        region_fills(reg_id).ie_field = init_conds::noICsScalar;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_fills(reg_id).ie_field = init_conds::noICsScalar;
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
                        for (const auto& element : str_region_ie_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** Region internal energy Inputs Not Understood ****");
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
                Yaml::Node& inps_subfields_yaml = root["regions"][reg_id]["region"]["specific_heat"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_region_specific_heat_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_region_specific_heat_inputs, str_region_specific_heat_inps, region_specific_heat_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_region_specific_heat_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        // x-component of specific_heat
                        double value = root["regions"][reg_id]["region"]["specific_heat"]["value"].As<double>();
                        RUN({
                        region_fills(reg_id).specific_heat = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["regions"][reg_id]["region"]["specific_heat"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., init_conds::uniform 
                            switch(scalar_ics_type_map[type]){

                                case init_conds::uniform:
                                    std::cout << "Setting specific_heat initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_fills(reg_id).specific_heat_field = init_conds::uniform;
                                    });
                                    break;

                                case init_conds::tgVortexScalar:
                                    std::cout << "Setting specific_heat initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_fills(reg_id).specific_heat_field = init_conds::tgVortexScalar;
                                    });
                                    break;

                                case init_conds::noICsScalar:
                                    std::cout << "Setting specific_heat initial conditions type to no specific_heat" << std::endl;
                                    RUN({ 
                                        region_fills(reg_id).specific_heat_field = init_conds::noICsScalar;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_fills(reg_id).specific_heat_field = init_conds::noICsScalar;
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
                        for (const auto& element : str_region_specific_heat_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** Region specific_heat Inputs Not Understood ****");
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
                Yaml::Node& inps_subfields_yaml = root["regions"][reg_id]["region"]["thermal_conductivity"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_region_thermal_conductivity_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_region_thermal_conductivity_inputs, str_region_thermal_conductivity_inps, region_thermal_conductivity_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_region_thermal_conductivity_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        // thermal_conductivity
                        double value = root["regions"][reg_id]["region"]["thermal_conductivity"]["value"].As<double>();

                        RUN({
                        region_fills(reg_id).thermal_conductivity = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["regions"][reg_id]["region"]["thermal_conductivity"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., init_conds::uniform 
                            switch(scalar_ics_type_map[type]){

                                case init_conds::uniform:
                                    std::cout << "Setting thermal_conductivity initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_fills(reg_id).thermal_conductivity_field = init_conds::uniform;
                                    });
                                    break;

                                case init_conds::tgVortexScalar:
                                    std::cout << "Setting thermal_conductivity initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_fills(reg_id).thermal_conductivity_field = init_conds::tgVortexScalar;
                                    });
                                    break;

                                case init_conds::noICsScalar:
                                    std::cout << "Setting thermal_conductivity initial conditions type to no thermal_conductivity" << std::endl;
                                    RUN({ 
                                        region_fills(reg_id).thermal_conductivity_field = init_conds::noICsScalar;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_fills(reg_id).thermal_conductivity_field = init_conds::noICsScalar;
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
                        for (const auto& element : str_region_thermal_conductivity_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** Region thermal_conductivity Inputs Not Understood ****");
                    } // end if on all subfields under thermal_conductivity

                } // end for loop over text

            } // thermal_conductivity
            else if (a_word.compare("volume_fraction") == 0){

                // always built, so no need to add this varname to guass_point

                // -----
                // loop over the sub fields under volfrac
                // -----
                Yaml::Node& inps_subfields_yaml = root["regions"][reg_id]["region"]["volume_fraction"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_region_volfrac_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_region_volfrac_inputs, str_region_volfrac_inps, region_volfrac_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_region_volfrac_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        // volfrac value or the intercept if linear variation
                        double value = root["regions"][reg_id]["region"]["volume_fraction"]["value"].As<double>();
      
                        RUN({
                            region_fills(reg_id).volfrac = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("slope") == 0) {
                        // volfrac slope
                        double slope = root["regions"][reg_id]["region"]["volume_fraction"]["slope"].As<double>();
      
                        RUN({
                            region_fills(reg_id).volfrac_slope = slope;
                        });
                    } // slope
                    else if (a_subfield_word.compare("origin") == 0) {
                        std::string origin = root["regions"][reg_id]["region"]["volume_fraction"]["origin"].As<std::string>();

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
                            region_fills(reg_id).volfrac_origin[0] = x1;
                            region_fills(reg_id).volfrac_origin[1] = y1;
                            region_fills(reg_id).volfrac_origin[2] = z1;
                        });
                    } // origin
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["regions"][reg_id]["region"]["volume_fraction"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., init_conds::uniform 
                            switch(scalar_ics_type_map[type]){

                                case init_conds::uniform:
                                    std::cout << "Setting volfrac initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_fills(reg_id).volfrac_field = init_conds::uniform;
                                    });
                                    break;

                                case init_conds::radialScalar:
                                    std::cout << "Setting volfrac initial conditions type to radial scalar " << std::endl;
                                    RUN({
                                        region_fills(reg_id).volfrac_field = init_conds::radialScalar;
                                    });
                                    break;

                                case init_conds::sphericalScalar:
                                    std::cout << "Setting volfrac initial conditions type to spherical scalar " << std::endl;
                                    RUN({
                                        region_fills(reg_id).volfrac_field = init_conds::sphericalScalar;
                                    });
                                    break;

                                case init_conds::xlinearScalar:
                                    std::cout << "Setting volfrac initial conditions type to xlinearScalar " << std::endl;
                                    RUN({
                                        region_fills(reg_id).volfrac_field = init_conds::xlinearScalar;
                                    });
                                    break;
                                
                                case init_conds::ylinearScalar:
                                    std::cout << "Setting volfrac initial conditions type to ylinearScalar " << std::endl;
                                    RUN({
                                        region_fills(reg_id).volfrac_field = init_conds::ylinearScalar;
                                    });
                                    break;
                                
                                case init_conds::zlinearScalar:
                                    std::cout << "Setting volfrac initial conditions type to zlinearScalar " << std::endl;
                                    RUN({
                                        region_fills(reg_id).volfrac_field = init_conds::zlinearScalar;
                                    });
                                    break;

                                case init_conds::tgVortexScalar:
                                    std::cout << "Setting volfrac initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_fills(reg_id).volfrac_field = init_conds::tgVortexScalar;
                                    });
                                    break;

                                case init_conds::noICsScalar:
                                    std::cout << "Default Volume Fraction Used:" << std::endl;
                                    std::cout << "Setting volume fraction to uniform field with a value equal to 1" << std::endl;
                                    RUN({ 
                                        region_fills(reg_id).volfrac_field = init_conds::uniform;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_fills(reg_id).volfrac_field = init_conds::noICsScalar;
                                    });

                                    std::cout << "ERROR: No valid volume fraction intial conditions type input " << std::endl;
                                    std::cout << "Valid IC types are: " << std::endl;
                                    
                                    for (const auto& pair : scalar_ics_type_map) {
                                        std::cout << pair.second << std::endl;
                                    }

                                    throw std::runtime_error("**** Volume Fraction Initial Conditions Type Not Understood ****");
                                    break;
                            } // end switch

                        }
                        else{
                            std::cout << "ERROR: invalid input: " << type << std::endl;
                            throw std::runtime_error("**** Volume Fraction IC Not Understood ****");
                        } // end if on Volume Fraction type
                        
                    } // end if on Volume Fraction type
                    else {
                        std::cout << "ERROR: invalid input: " << a_subfield_word << std::endl;
                        std::cout << "Valid options are: " << std::endl;
                        for (const auto& element : str_region_volfrac_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** Region Volume Fraction Inputs Not Understood ****");
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
                Yaml::Node& inps_subfields_yaml = root["regions"][reg_id]["region"]["temperature"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_region_temperature_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_region_temperature_inputs, 
                                str_region_temperature_inps, region_temperature_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_region_temperature_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        //temperature
                        double value = root["regions"][reg_id]["region"]["temperature"]["value"].As<double>();
                        RUN({
                        region_fills(reg_id).temperature = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["regions"][reg_id]["region"]["temperature"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., init_conds::uniform 
                            switch(scalar_ics_type_map[type]){

                                case init_conds::uniform:
                                    std::cout << "Setting temperature initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_fills(reg_id).temperature_field = init_conds::uniform;
                                    });
                                    break;

                                case init_conds::tgVortexScalar:
                                    std::cout << "Setting temperature initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_fills(reg_id).temperature_field = init_conds::tgVortexScalar;
                                    });
                                    break;

                                case init_conds::noICsScalar:
                                    std::cout << "Setting temperature initial conditions type to no temperature" << std::endl;
                                    RUN({ 
                                        region_fills(reg_id).temperature_field = init_conds::noICsScalar;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_fills(reg_id).temperature_field = init_conds::noICsScalar;
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
                        for (const auto& element : str_region_temperature_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** Region temperature Inputs Not Understood ****");
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
                Yaml::Node& inps_subfields_yaml = root["regions"][reg_id]["region"]["level_set"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_region_level_set_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_region_level_set_inputs, str_region_level_set_inps, region_level_set_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_region_level_set_inputs){ 

                    if (a_subfield_word.compare("value") == 0) {
                        // level set
                        double value = root["regions"][reg_id]["region"]["level_set"]["value"].As<double>();

                        RUN({
                            region_fills(reg_id).level_set = value;
                        });
                    } // value
                    else if (a_subfield_word.compare("slope") == 0) {
                        // volfrac slope
                        double slope = root["regions"][reg_id]["region"]["level_set"]["slope"].As<double>();
        
                        RUN({
                            region_fills(reg_id).level_set_slope = slope;
                        });
                    } // slope
                    else if (a_subfield_word.compare("origin") == 0) {
                        std::string origin = root["regions"][reg_id]["region"]["level_set"]["origin"].As<std::string>();

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
                            region_fills(reg_id).level_set_origin[0] = x1;
                            region_fills(reg_id).level_set_origin[1] = y1;
                            region_fills(reg_id).level_set_origin[2] = z1;
                        });
                    } // origin
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["regions"][reg_id]["region"]["level_set"]["type"].As<std::string>();

                        // set the IC tag type
                        if (scalar_ics_type_map.find(type) != scalar_ics_type_map.end()) {
                        
                            // scalar_ics_type_map[type] returns enum value, e.g., init_conds::uniform 
                            switch(scalar_ics_type_map[type]){

                                case init_conds::uniform:
                                    std::cout << "Setting level set initial conditions type to uniform " << std::endl;
                                    RUN({
                                        region_fills(reg_id).level_set_field = init_conds::uniform;
                                    });
                                    break;
                                case init_conds::radialScalar:
                                    std::cout << "Setting level set initial conditions type to radial scalar " << std::endl;
                                    RUN({
                                        region_fills(reg_id).level_set_field = init_conds::radialScalar;
                                    });
                                    break;
                                case init_conds::sphericalScalar:
                                    std::cout << "Setting level set initial conditions type to spherical scalar " << std::endl;
                                    RUN({
                                        region_fills(reg_id).level_set_field = init_conds::sphericalScalar;
                                    });
                                    break;
                                case init_conds::xlinearScalar:
                                    std::cout << "Setting volfrac initial conditions type to xlinearScalar " << std::endl;
                                    RUN({
                                        region_fills(reg_id).level_set_field = init_conds::xlinearScalar;
                                    });
                                    break;
                                
                                case init_conds::ylinearScalar:
                                    std::cout << "Setting volfrac initial conditions type to ylinearScalar " << std::endl;
                                    RUN({
                                        region_fills(reg_id).level_set_field = init_conds::ylinearScalar;
                                    });
                                    break;
                                
                                case init_conds::zlinearScalar:
                                    std::cout << "Setting volfrac initial conditions type to zlinearScalar " << std::endl;
                                    RUN({
                                        region_fills(reg_id).level_set_field = init_conds::zlinearScalar;
                                    });
                                    break;
                                case init_conds::tgVortexScalar:
                                    std::cout << "Setting level set initial conditions type to TG Vortex " << std::endl;
                                    RUN({
                                        region_fills(reg_id).level_set_field = init_conds::tgVortexScalar;
                                    });
                                    break;

                                case init_conds::noICsScalar:
                                    std::cout << "Setting level set initial conditions type to no level set" << std::endl;
                                    RUN({ 
                                        region_fills(reg_id).level_set_field = init_conds::noICsScalar;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_fills(reg_id).level_set_field = init_conds::noICsScalar;
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
                        for (const auto& element : str_region_level_set_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** Region level set Inputs Not Understood ****");
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
                Yaml::Node& inps_subfields_yaml = root["regions"][reg_id]["region"]["velocity"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_region_vel_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_region_vel_inputs, str_region_vel_inps, region_vel_required_inps);

                // loop over the subfield words
                for(auto& a_subfield_word : user_region_vel_inputs){ 

                    if (a_subfield_word.compare("u") == 0) {
                        // x-component of velocity
                        double u = root["regions"][reg_id]["region"]["velocity"]["u"].As<double>();

                        RUN({
                        region_fills(reg_id).u = u;
                        });
                    } // u
                    else if (a_subfield_word.compare("v") == 0) {
                        // y-component of velocity
                        double v = root["regions"][reg_id]["region"]["velocity"]["v"].As<double>();

                        RUN({
                            region_fills(reg_id).v = v;
                        });
                    } // v
                    else if (a_subfield_word.compare("w") == 0) {
                        // z-component of velocity

                        double w = root["regions"][reg_id]["region"]["velocity"]["w"].As<double>();

                        RUN({
                            region_fills(reg_id).w = w;
                        });
                    } // w
                    else if (a_subfield_word.compare("speed") == 0) {
                        double speed = root["regions"][reg_id]["region"]["velocity"]["speed"].As<double>();

                        RUN({
                            region_fills(reg_id).speed = speed;
                        });
                    } // speed
                    else if (a_subfield_word.compare("type") == 0){

                        std::string type = root["regions"][reg_id]["region"]["velocity"]["type"].As<std::string>();

                        // set the volume tag type
                        if (vector_ics_type_map.find(type) != vector_ics_type_map.end()) {
                        
                            // vector_ics_type_map[type] returns enum value, e.g., init_conds::velocity 
                            switch(vector_ics_type_map[type]){

                                case init_conds::stationary:
                                    std::cout << "Setting velocity initial conditions type to static " << std::endl;
                                    RUN({
                                        region_fills(reg_id).vel_field = init_conds::stationary;
                                    });
                                    break;

                                case init_conds::cartesian:
                                    std::cout << "Setting velocity initial conditions type to cartesian " << std::endl;
                                    RUN({
                                        region_fills(reg_id).vel_field = init_conds::cartesian;
                                    });
                                    break;

                                case init_conds::radialVec:
                                    std::cout << "Setting velocity initial conditions type to radial " << std::endl;
                                    RUN({
                                        region_fills(reg_id).vel_field = init_conds::radialVec;
                                    });
                                    break;

                                case init_conds::sphericalVec:
                                    std::cout << "Setting velocity initial conditions type to spherical " << std::endl;
                                    RUN({
                                        region_fills(reg_id).vel_field = init_conds::sphericalVec;
                                    });
                                    break;

                                case init_conds::radialLinearVec:
                                    std::cout << "Setting velocity initial conditions type to radial_linear " << std::endl;
                                    RUN({
                                        region_fills(reg_id).vel_field = init_conds::radialLinearVec;
                                    });
                                    break;

                                case init_conds::sphericalLinearVec:
                                    std::cout << "Setting velocity initial conditions type to spherical_linear " << std::endl;
                                    RUN({
                                        region_fills(reg_id).vel_field = init_conds::sphericalLinearVec;
                                    });
                                    break;

                                case init_conds::tgVortexVec:
                                    std::cout << "Setting velocity initial conditions type to tg_vortex " << std::endl;
                                    RUN({
                                        region_fills(reg_id).vel_field = init_conds::tgVortexVec;
                                    });
                                    break;

                                case init_conds::noICsVec:
                                    std::cout << "Setting velocity initial conditions type to no velocity" << std::endl;
                                    RUN({ 
                                        region_fills(reg_id).vel_field = init_conds::noICsVec;
                                    });
                                    break;

                                default:

                                    RUN({ 
                                        region_fills(reg_id).vel_field = init_conds::noICsVec;
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
                        for (const auto& element : str_region_vel_inps) {
                            std::cout << element << std::endl;
                        }
                        throw std::runtime_error("**** Region Velocity Inputs Not Understood ****");
                    } // end if on all subfields under velocity

                } // end for loop over text
            } // end if on velocity
            else if (a_word.compare("volume") == 0) {

                // -----
                // loop over the sub fields under volume
                // -----
                Yaml::Node& inps_subfields_yaml = root["regions"][reg_id]["region"]["volume"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_region_volume_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_region_volume_inputs, str_region_volume_inps, region_volume_required_inps);


                // loop over the subfield words
                for(auto& a_subfield_word : user_region_volume_inputs){ 

                    if (a_subfield_word.compare("radius1") == 0) {
                        // inner radius of sphere/cylinder

                        double radius1 = root["regions"][reg_id]["region"]["volume"]["radius1"].As<double>();

                        RUN({
                            region_fills(reg_id).radius1 = radius1;
                        });
                    } // radius1
                    else if (a_subfield_word.compare("radius2") == 0) {
                        // outer radius of sphere/cylinder

                        double radius2 = root["regions"][reg_id]["region"]["volume"]["radius2"].As<double>();

                        RUN({
                            region_fills(reg_id).radius2 = radius2;
                        });
                    } // radius2
                    else if (a_subfield_word.compare("x1") == 0) {
                        // inner plane

                        double x1 = root["regions"][reg_id]["region"]["volume"]["x1"].As<double>();

                        RUN({
                            region_fills(reg_id).x1 = x1;
                        });
                    } // x1
                    else if (a_subfield_word.compare("x2") == 0) {
                        // outer plane

                        double x2 = root["regions"][reg_id]["region"]["volume"]["x2"].As<double>();

                        RUN({
                            region_fills(reg_id).x2 = x2;
                        });
                    } // x2
                    else if (a_subfield_word.compare("y1") == 0) {
                        // inner plane

                        double y1 = root["regions"][reg_id]["region"]["volume"]["y1"].As<double>();

                        RUN({
                            region_fills(reg_id).y1 = y1;
                        });
                    } // y1
                    else if (a_subfield_word.compare("y2") == 0) {
                        // outer plane

                        double y2 = root["regions"][reg_id]["region"]["volume"]["y2"].As<double>();

                        RUN({
                            region_fills(reg_id).y2 = y2;
                        });
                    } // y2
                    else if (a_subfield_word.compare("z1") == 0) {
                        // inner plane

                        double z1 = root["regions"][reg_id]["region"]["volume"]["z1"].As<double>();

                        RUN({
                            region_fills(reg_id).z1 = z1;
                        });
                    } // z1
                    else if (a_subfield_word.compare("z2") == 0) {
                        // outer plane

                        double z2 = root["regions"][reg_id]["region"]["volume"]["z2"].As<double>();

                        RUN({
                            region_fills(reg_id).z2 = z2;
                        });
                    } // z2
                    else if (a_subfield_word.compare("scale_x") == 0) {
                        // outer plane

                        double scale_x = root["regions"][reg_id]["region"]["volume"]["scale_x"].As<double>();

                        // on the host side because it relates to reading a mesh file
                        region_fills_host(reg_id).scale_x = scale_x;

                    } // scale_x
                    else if (a_subfield_word.compare("scale_y") == 0) {
                        // outer plane

                        double scale_y = root["regions"][reg_id]["region"]["volume"]["scale_y"].As<double>();

                        // on the host side because it relates to reading a mesh file
                        region_fills_host(reg_id).scale_y = scale_y;

                    } // scale_y
                    else if (a_subfield_word.compare("scale_z") == 0) {
                        // outer plane

                        double scale_z = root["regions"][reg_id]["region"]["volume"]["scale_z"].As<double>();

                        // on the host side because it relates to reading a mesh file
                        region_fills_host(reg_id).scale_z = scale_z;

                    } // scale_z
                    else if (a_subfield_word.compare("part_id") == 0) {
                        // part_id in 

                        int part_id = root["regions"][reg_id]["region"]["volume"]["part_id"].As<int>();

                        RUN({
                            region_fills(reg_id).part_id = part_id;
                        });

                    } // scale_z
                    //
                    else if (a_subfield_word.compare("type") == 0) {

                        // region volume fill type
                        std::string type = root["regions"][reg_id]["region"]["volume"]["type"].As<std::string>();

                        // set the velocity tag type
                        if (region_type_map.find(type) != region_type_map.end()) {
                        
                            // region_type_map[type] returns enum value, e.g., init_conds::velocity 
                            switch(region_type_map[type]){

                                case region::global:
                                    std::cout << "Setting volume fill type to global " << std::endl;
                                    RUN({
                                        region_fills(reg_id).volume = region::global;
                                    });
                                    break;

                                case region::box:
                                    std::cout << "Setting volume fill type to box " << std::endl;
                                    RUN({
                                        region_fills(reg_id).volume = region::box;
                                    });
                                    break;

                                case region::cylinder:
                                    std::cout << "Setting volume fill type to cylinder " << std::endl;
                                    RUN({
                                        region_fills(reg_id).volume = region::cylinder;
                                    });
                                    break;

                                case region::sphere:
                                    std::cout << "Setting volume fill type to sphere " << std::endl;
                                    RUN({
                                        region_fills(reg_id).volume = region::sphere;
                                    });
                                    break;

                                case region::readVoxelFile:
                                    std::cout << "Setting volume fill type to readVoxelFile " << std::endl;
                                    RUN({
                                        region_fills(reg_id).volume = region::readVoxelFile;
                                    });
                                    break;
                                case region::readVTUFile:
                                    std::cout << "Setting volume fill type to readVTUFile " << std::endl;
                                    RUN({
                                        region_fills(reg_id).volume = region::readVTUFile;
                                    });
                                    break;
                                case region::no_volume:
                                    std::cout << "Setting volume fill type to none " << std::endl;
                                    RUN({
                                        region_fills(reg_id).volume = region::no_volume;
                                    });
                                    break;
                                default:

                                    RUN({ 
                                        region_fills(reg_id).volume = region::no_volume;
                                    });

                                    std::cout << "ERROR: No valid region volume fill type input " << std::endl;
                                    std::cout << "Valid IC volume fill types are: " << std::endl;
                                    
                                    for (const auto& pair : region_type_map) {
                                        std::cout << pair.second << std::endl;
                                    }

                                    throw std::runtime_error("**** Region Volume Fill Type Not Understood ****");
                                    break;
                            } // end switch
                        } // end if on setting volume tag

                    }  // end if on volume type
                    // Get mesh file path
                    else if (a_subfield_word.compare("file_path") == 0) {
                        // region volume fill type
                        std::string path = root["regions"][reg_id]["region"]["volume"]["file_path"].As<std::string>();

                        // absolute path to file or local to the director where exe is run
                        region_fills_host(reg_id).file_path = path;   // saving the absolute file path
        

                    } // end file path
                    //
                    else if (a_subfield_word.compare("origin") == 0) {
                        std::string origin = root["regions"][reg_id]["region"]["volume"]["origin"].As<std::string>();

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
                            region_fills(reg_id).origin[0] = x1;
                            region_fills(reg_id).origin[1] = y1;
                            region_fills(reg_id).origin[2] = z1;
                        });
                    } // origin
                    else{
                        std::cout << "ERROR: invalid input: " << a_subfield_word << std::endl;
                        throw std::runtime_error("**** Volume Fill Not Understood ****");
                    } // end if

                } // end for loop over subfields under volume
            } // end volume fill type
            else {
                std::cout << "ERROR: invalid input: " << a_word << std::endl;
                std::cout << "Valid options are: " << std::endl;
                for (const auto& element : str_region_inps) {
                    std::cout << element << std::endl;
                }
                throw std::runtime_error("**** Region Not Understood ****");
            }
        } // end for words in fill region

        // update the device
        reg_fills_in_solver.update_device();    
        num_reg_fills_in_solver.update_device();

        // -----------------------------------------------
        // check for consistency in input settings

        if(fill_gauss_states.size() == 0){
            throw std::runtime_error("**** No Initial Conditions Set, Please Specify Fields in the Region Fills ****");
        }

        // NOTE:
        // Each solver checks to in the intiialize function to see if the fills specified 
        // are sufficient for setting up the fields in the simulation


        // check to see if a file path is empty
        if(region_fills_host(reg_id).file_path.empty()){

            RUN({
                // if the following is true, stop simulation; must add all mesh read options
                if (region_fills(reg_id).volume == region::readVoxelFile) {
                    Kokkos::abort("\n********************************************************************************************\n"
                                    "ERROR: \n"
                                    "When using a file to initialize a region, a file_path must be set to point to the mesh file\n"
                                    "********************************************************************************************\n");
                }
            });
        } // end if check

        // check to see if a file path was set
        if(region_fills_host(reg_id).file_path.size()>0){
            RUN({
                if (region_fills(reg_id).volume != region::readVoxelFile){  
                    // this means it is a geometric definition of the region
                    Kokkos::abort("\n********************************************************************************************\n"
                                    "ERROR: \n"
                                    "When a geometric entity defines the region, a mesh file cannot be passed to set the region\n"
                                    "********************************************************************************************\n");
                }
            });

        }  // end if        

        // -----------------------------------------------


    } // end loop over regions
} // end of function to parse region
