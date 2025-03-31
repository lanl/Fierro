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
#include "parse_material_inputs.hpp"


// eos files
#include "gamma_law_eos.h"
#include "no_eos.h"
#include "user_defined_eos.h"
#include "void_eos.h"
#include "host_user_defined_eos.h"

// ----
#if __has_include("analytic_defined_eos.h")
#include "analytic_defined_eos.h"
#endif
// ----

// strength
#include "no_strength.h"
#include "user_defined_strength.h"
#include "host_user_defined_strength.h"
#include "host_ann_strength.h"

// ----
#if __has_include("decoupled_strength.h")
#include "decoupled_strength.h"
#endif
// ----

// erosion files
#include "basic_erosion.h"
#include "no_erosion.h"

// dissipation files
#include "mars.h"
#include "no_dissipation.h"

// fracture files
#include "user_defined_fracture.h"



// ==============================================================================
//   Function Definitions
// ==============================================================================




// =================================================================================
//    Parse Material Definitions
// =================================================================================
void parse_materials(Yaml::Node& root, Material_t& Materials, const size_t num_dims)
{
    Yaml::Node& material_yaml = root["materials"];

    size_t num_materials = material_yaml.Size();
    std::cout << "Number of materials =  "<< num_materials << std::endl;

    Materials.num_mats = num_materials;

    // --- allocate memory for arrays inside material struct ---

    // setup
    Materials.MaterialSetup = DCArrayKokkos<MaterialSetup_t>(num_materials, "material_setup");

    // function pointers to material models
    Materials.MaterialFunctions = DCArrayKokkos<MaterialFunctions_t>(num_materials, "material_functions");

    // enums
    Materials.MaterialEnums = DCArrayKokkos<MaterialEnums_t>(num_materials, "material_enums");

    // these are temp arrays to store global variables given in the yaml input file for each material, 100 vars is the max allowable
    DCArrayKokkos<double> tempGlobalEOSVars(num_materials, 100, "temp_array_eos_vars");
    DCArrayKokkos<double> tempGlobalStrengthVars(num_materials, 100, "temp_array_strength_vars");
    DCArrayKokkos<double> tempGlobalDissipationVars(num_materials, 10, "temp_array_dissipation_vars");

    Materials.num_eos_global_vars      =  CArrayKokkos <size_t> (num_materials, "num_eos_global_vars");
    Materials.num_strength_global_vars =  CArrayKokkos <size_t> (num_materials, "num_strength_global_vars");
    Materials.num_dissipation_global_vars = CArrayKokkos <size_t> (num_materials, "num_dissipations_vars");

    // initialize the num of global vars to 0 for all models
    FOR_ALL(mat_id, 0, num_materials, {

        Materials.num_eos_global_vars(mat_id) = 0;
        Materials.num_strength_global_vars(mat_id) = 0;
        Materials.num_dissipation_global_vars(mat_id) = 0;   // a minimum of 6 inputs  
        
    }); // end parallel for

    // a check on material_id not being specified more than once or not at all
    CArray <bool> check_mat_ids(num_materials);
    check_mat_ids.set_values(false);

    // loop over the materials specified in the input file
    for (int m_id = 0; m_id < num_materials; m_id++) {

        // Important: m_id corresponds to the order of the materials entered in the input file

        // read the variables names
        Yaml::Node& inps_yaml = root["materials"][m_id]["material"];

        size_t num_vars_set = inps_yaml.Size();

        std::cout << "Number of vars set =  "<< num_vars_set << std::endl;

        // get the material variables names set by the user
        std::vector<std::string> user_str_material_inps;

        // extract words from the input file and validate they are correct
        validate_inputs(inps_yaml, user_str_material_inps, str_material_inps, material_hydrodynamics_required_inps);

        // loop over the words in the material input definition and find the material id
        int mat_id = -1;
        for (auto& a_word : user_str_material_inps) {

            Yaml::Node& material_inps_yaml = root["materials"][m_id]["material"][a_word];

            if (a_word.compare("id") == 0) {
                mat_id = root["materials"][m_id]["material"]["id"].As<int>();

                if (mat_id<0 || mat_id>=num_materials){
                    std::cout << "ERROR: invalid material_id specified in the material definition " << std::endl;
            
                    throw std::runtime_error("**** Material_id is out of bounds ****");
                } // end check on m_id range

                if (check_mat_ids(mat_id) == true){
                    std::cout << "ERROR: material_id = " << mat_id << " was already specified "<< std::endl;
                    throw std::runtime_error("**** Multiple materials used the same material_id ****");
                }
                else {
                    check_mat_ids(mat_id) = true;
                } // end check on mat_id

            } // end id
        } // end loop over all material inputs

        if (mat_id<0){
            std::cout << "ERROR: material_id must be specified in the material definition " << std::endl;
            
            throw std::runtime_error("**** Material_id is missing ****");
        } // end check on m_id being specified



        // loop over the words in the material input definition again
        for (auto& a_word : user_str_material_inps) {

            Yaml::Node& material_inps_yaml = root["materials"][m_id]["material"][a_word];

            
            //extract eos model
            if (a_word.compare("id") == 0) {
                // do nothing
                // this id was read in an earlier loop
            }
            else if (a_word.compare("eos_model_type") == 0) {
                std::string type = root["materials"][m_id]["material"]["eos_model_type"].As<std::string>();

                // set the eos type
                if (eos_type_map.find(type) != eos_type_map.end()) {

                    // eos_type_map[type] returns enum value, e.g., model::decoupled
                    switch(eos_type_map[type]){
                        case model::decoupledEOSType:
                            std::cout << "Setting EOS type to decoupled " << std::endl;
                            RUN({
                                Materials.MaterialEnums(mat_id).EOSType = model::decoupledEOSType;
                            });

                            Materials.MaterialEnums.host(mat_id).EOSType = model::decoupledEOSType;
                            break;

                        case model::coupledEOSType:
                            std::cout << "Setting EOS type to coupled " << std::endl;
                            RUN({
                                Materials.MaterialEnums(mat_id).EOSType = model::coupledEOSType;
                            });
                            Materials.MaterialEnums.host(mat_id).EOSType = model::coupledEOSType;
                            break;

                        default:
                            RUN({ 
                                Materials.MaterialEnums(mat_id).EOSType = model::noEOSType;
                            });
                            Materials.MaterialEnums.host(mat_id).EOSType = model::noEOSType;
                            std::cout << "ERROR: No valid EOS type input " << std::endl;
                            std::cout << "Valid EOS types are: " << std::endl;
                            
                            for (const auto& pair : eos_type_map) {
                                std::cout << pair.second << std::endl;
                            }

                            throw std::runtime_error("**** EOS Type Not Understood ****");
                            break;
                    } // end switch
                } 
                else{
                    std::cout << "ERROR: invalid eos type input: " << type << std::endl;
                } // end if
            }
            //
            // set the eos_model
            else if (a_word.compare("eos_model") == 0) {
                std::string eos = root["materials"][m_id]["material"]["eos_model"].As<std::string>();

                // set the EOS
                if (eos_models_map.find(eos) != eos_models_map.end()) {
                    
                    switch(eos_models_map[eos]){

                        case model::noEOS:

                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_pressure    = &NoEOSModel::calc_pressure;
                                Materials.MaterialFunctions(mat_id).calc_sound_speed = &NoEOSModel::calc_sound_speed;

                                Materials.MaterialEnums(mat_id).EOSModels = model::noEOS;
                            });
                            Materials.MaterialEnums.host(mat_id).EOSModels = model::noEOS;
                            break;

                        case model::gammaLawGasEOS:
                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_pressure    = &GammaLawGasEOSModel::calc_pressure;
                                Materials.MaterialFunctions(mat_id).calc_sound_speed = &GammaLawGasEOSModel::calc_sound_speed;

                                Materials.MaterialEnums(mat_id).EOSModels = model::gammaLawGasEOS;
                            });
                            Materials.MaterialEnums.host(mat_id).EOSModels = model::gammaLawGasEOS;
                            break;

                        case model::voidEOS:
                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_pressure    = &VoidEOSModel::calc_pressure;
                                Materials.MaterialFunctions(mat_id).calc_sound_speed = &VoidEOSModel::calc_sound_speed;

                                Materials.MaterialEnums(mat_id).EOSModels = model::voidEOS;
                            });
                            Materials.MaterialEnums.host(mat_id).EOSModels = model::voidEOS;
                            break;

                        case model::userDefinedEOS:
                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_pressure    = &UserDefinedEOSModel::calc_pressure;
                                Materials.MaterialFunctions(mat_id).calc_sound_speed = &UserDefinedEOSModel::calc_sound_speed;

                                Materials.MaterialEnums(mat_id).EOSModels = model::userDefinedEOS;
                            });
                            Materials.MaterialEnums.host(mat_id).EOSModels = model::userDefinedEOS;
                            break;
                        // --------
                        // adding host run EOSs
                        case model::hostUserDefinedEOS:
                            Materials.MaterialFunctions.host(mat_id).calc_pressure    = &HostUserDefinedEOSModel::calc_pressure;
                            Materials.MaterialFunctions.host(mat_id).calc_sound_speed = &HostUserDefinedEOSModel::calc_sound_speed;
                            Materials.MaterialEnums.host(mat_id).EOSModels = model::hostUserDefinedEOS;

                            RUN({
                                Materials.MaterialEnums(mat_id).EOSRunLocation = model::host;
                            });
                            Materials.MaterialEnums.host(mat_id).EOSRunLocation = model::host;
                            break;
#ifdef ANALYTIC_DEFINED_EOS_H
                        // call Gruneisen
                        case model::mieGruneisenEOS:
                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_pressure    = &MieGruneisenEOSModel::calc_pressure;
                                Materials.MaterialFunctions(mat_id).calc_sound_speed = &MieGruneisenEOSModel::calc_sound_speed;

                                Materials.MaterialEnums(mat_id).EOSModels = model::mieGruneisenEOS;
                            });
                            Materials.MaterialEnums.host(mat_id).EOSModels = model::mieGruneisenEOS;
                            break;  

                        // add other analytic EOS models here, e.g., Johnson-Cook etc.
                        // ....

#endif
                        default:
                            std::cout << "ERROR: invalid input: " << eos << std::endl;
                            throw std::runtime_error("**** EOS Not Understood ****");
                            break;
                    } // end switch on EOS type
                }
                else{
                    std::cout << "ERROR: invalid EOS input: " << eos << std::endl;
                    throw std::runtime_error("**** EOS Not Understood ****");
                } // end if
            } // EOS model

            // Type of strength model
            else if (a_word.compare("strength_model_type") == 0) {
                std::string strength_model_type = root["materials"][m_id]["material"]["strength_model_type"].As<std::string>();

                // set the EOS
                if (strength_type_map.find(strength_model_type) != strength_type_map.end()) {
                    
                    switch(strength_type_map[strength_model_type]){

                        case model::noStrengthType:
                            RUN({
                                Materials.MaterialEnums(mat_id).StrengthType = model::noStrengthType;
                            });
                            Materials.MaterialEnums.host(mat_id).StrengthType = model::noStrengthType;
                            break;

                        case model::incrementBased:
                            RUN({
                                Materials.MaterialEnums(mat_id).StrengthType = model::incrementBased;
                            });
                            Materials.MaterialEnums.host(mat_id).StrengthType = model::incrementBased;
                            
                            break;
                        case model::stateBased:
                            RUN({
                                Materials.MaterialEnums(mat_id).StrengthType = model::stateBased;
                            });
                            Materials.MaterialEnums.host(mat_id).StrengthType = model::stateBased;
                            
                            std::cout << "ERROR: state_based models not yet defined: " << std::endl;
                            throw std::runtime_error("**** ERROR: state_based models not yet defined ****");
                            break;
                        default:
                            std::cout << "ERROR: invalid strength type input: " << strength_model_type << std::endl;
                            throw std::runtime_error("**** Strength Model Type Not Understood ****");
                            break;
                    } // end switch on EOS type
                }
                else{
                    std::cout << "ERROR: Invalid strength model type input: " << strength_model_type << std::endl;
                    throw std::runtime_error("**** Strength model type not understood ****");
                } // end if
            } // Strength model type
            
            // Set specific strength model
            else if (a_word.compare("strength_model") == 0) {
                std::string strength_model = root["materials"][m_id]["material"]["strength_model"].As<std::string>();

                // set the strength
                if (strength_models_map.find(strength_model) != strength_models_map.end()) {

                    std::cout << "strength model = \n" << strength_models_map[strength_model] << std::endl;
                    
                    switch(strength_models_map[strength_model]){

                        case model::noStrengthModel:
                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_stress = &NoStrengthModel::calc_stress;
                            });
                            // note: default run location for strength is device

                            Materials.MaterialFunctions.host(mat_id).init_strength_state_vars = &NoStrengthModel::init_strength_state_vars;
                            
                            // the default run location for the initialization function is host, but below here shows how to set it
                            RUN({
                                Materials.MaterialEnums(mat_id).StrengthSetupLocation = model::host;
                            });
                            Materials.MaterialEnums.host(mat_id).StrengthSetupLocation = model::host;

                            break;

                        case model::userDefinedStrength:
                            
                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_stress = &UserDefinedStrengthModel::calc_stress;
                            });
                            // note: default run location for strength is device

                            Materials.MaterialFunctions.host(mat_id).init_strength_state_vars = &UserDefinedStrengthModel::init_strength_state_vars;
                            // note: default run location for initialization is always host

                            break;

                        case model::hostANNStrength:
                            
                            // set the stress function
                            Materials.MaterialFunctions.host(mat_id).calc_stress = &HostANNStrengthModel::calc_stress;

                            // set the run location for strength
                            Materials.MaterialEnums(mat_id).StrengthRunLocation = model::host;
                            Materials.MaterialEnums.host(mat_id).StrengthRunLocation = model::host;

                            // set the strength initialization function
                            Materials.MaterialFunctions.host(mat_id).init_strength_state_vars = &HostANNStrengthModel::init_strength_state_vars;
                            // note: default run location for initialization is always host

                            break;
#ifdef DECOUPLED_STRENGTH_H
                        // call elastic plastic model
                        case model::hypoElasticPlasticStrength:

                            if(num_dims == 2){
                                std::cout << "ERROR: specified 2D but this is a 3D strength model: " << strength_model << std::endl;
                                throw std::runtime_error("**** Strength model is not valid in 2D ****");
                            }

                            // set the stress function
                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_stress = &HypoElasticPlasticModel::calc_stress;
                            });
                            // note: default run location for strength is device

                            // set the strength initialization function
                            Materials.MaterialFunctions.host(mat_id).init_strength_state_vars = &HypoElasticPlasticModel::init_strength_state_vars;
                            // note: default run location for initialization is always host

                            break;  

                        
                        case model::hypoElasticPlasticStrengthRZ:

                            if(num_dims == 3){
                                std::cout << "ERROR: specified 3D but this is a 2D-RZ strength model: " << strength_model << std::endl;
                                throw std::runtime_error("**** Strength model is not valid in 3D ****");
                            }

                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_stress = &HypoElasticPlasticRZModel::calc_stress;
                            });
                            // note: default run location for strength is device

                            // set the strength initialization function
                            Materials.MaterialFunctions.host(mat_id).init_strength_state_vars = &HypoElasticPlasticRZModel::init_strength_state_vars;
                            // note: default run location for initialization is always host

                            break;  

                        // add other elastic plastic models here, e.g., Johnson-Cook strength etc.
                        // ....
                        
#endif
                        default:
                            std::cout << "ERROR: invalid strength input: " << strength_model << std::endl;
                            throw std::runtime_error("**** Strength model Not Understood ****");
                            break;
                    } // end switch on strength model name
                }
                else{
                    std::cout << "ERROR: invalid Strength model input: " << strength_model << std::endl;
                    throw std::runtime_error("**** Strength model Not Understood ****");
                } // end if
            } // Strength model
            
            //extract erosion model
            else if (a_word.compare("erosion_model") == 0) {
                std::string erosion_model = root["materials"][m_id]["material"]["erosion_model"].As<std::string>();

                // set the erosion model
                if (erosion_model_map.find(erosion_model) != erosion_model_map.end()) {

                    // erosion_model_map[erosion_model] returns enum value, e.g., model::erosion
                    switch(erosion_model_map[erosion_model]){
                        case model::basicErosion:
                            Materials.MaterialEnums.host(mat_id).ErosionModels = model::basicErosion;
                            RUN({
                                Materials.MaterialEnums(mat_id).ErosionModels = model::basicErosion;
                                Materials.MaterialFunctions(mat_id).erode = &BasicErosionModel::erode;
                            });
                            break;
                        case model::noErosion:
                            Materials.MaterialEnums.host(mat_id).ErosionModels = model::noErosion;
                            RUN({
                                Materials.MaterialEnums(mat_id).ErosionModels = model::noErosion;
                                Materials.MaterialFunctions(mat_id).erode = &NoErosionModel::erode;
                            });
                            break;
                        default:
                            std::cout << "ERROR: invalid erosion input: " << erosion_model << std::endl;
                            throw std::runtime_error("**** Erosion model Not Understood ****");
                            break;
                    } // end switch
                } 
                else{
                    std::cout << "ERROR: invalid erosion type input: " << erosion_model<< std::endl;
                    throw std::runtime_error("**** Erosion model Not Understood ****");
                    break;
                } // end if

            } // erosion model variables
            //extract dissipation (artificial viscosity) model
            else if (a_word.compare("dissipation_model") == 0) {
                std::string dissipation_model = root["materials"][m_id]["material"]["dissipation_model"].As<std::string>();

                // set the erosion model
                if (dissipation_model_map.find(dissipation_model) != dissipation_model_map.end()) {

                    // dissipation_model_map[dissipation_model] returns enum value, e.g., model::dissipation
                    switch(dissipation_model_map[dissipation_model]){
                        case model::MARS:
                            
                            if(num_dims == 2){
                                std::cout << "ERROR: specified 2D but this is a 3D MARS model: " << dissipation_model << std::endl;
                                throw std::runtime_error("**** Dissipation model is not valid in 2D ****");
                            }

                            Materials.MaterialEnums.host(mat_id).DissipationModels = model::MARS;
                            RUN({
                                Materials.MaterialEnums(mat_id).DissipationModels = model::MARS;
                                Materials.MaterialFunctions(mat_id).calc_dissipation = &MARSDissipationModel::calc_dissipation;
                            });
                            break;
                        //
                        case model::MARSRZ:
                            
                            if(num_dims == 3){
                                std::cout << "ERROR: specified 3D but this is a 2D-RZ MARS model: " << dissipation_model << std::endl;
                                throw std::runtime_error("**** Dissipation model is not valid in 3D ****");
                            }

                            Materials.MaterialEnums.host(mat_id).DissipationModels = model::MARSRZ;
                            RUN({
                                Materials.MaterialEnums(mat_id).DissipationModels = model::MARSRZ;
                                Materials.MaterialFunctions(mat_id).calc_dissipation = &MARSRZDissipationModel::calc_dissipation;
                            });
                            break;
                        //
                        case model::directionalMARS:
                            
                            if(num_dims == 2){
                                std::cout << "ERROR: specified 2D but this is a 3D MARS model: " << dissipation_model << std::endl;
                                throw std::runtime_error("**** Dissipation model is not valid in 2D ****");
                            }

                            Materials.MaterialEnums.host(mat_id).DissipationModels = model::directionalMARS;
                            RUN({
                                Materials.MaterialEnums(mat_id).DissipationModels = model::directionalMARS;
                                Materials.MaterialFunctions(mat_id).calc_dissipation = &MARSDissipationModel::calc_dissipation;
                            });
                            break;
                        //
                        case model::directionalMARSRZ:
                            
                            if(num_dims == 3){
                                std::cout << "ERROR: specified 3D but this is a 2D-RZ MARS model: " << dissipation_model << std::endl;
                                throw std::runtime_error("**** Dissipation model is not valid in 3D ****");
                            }

                            Materials.MaterialEnums.host(mat_id).DissipationModels = model::directionalMARSRZ;
                            RUN({
                                Materials.MaterialEnums(mat_id).DissipationModels = model::directionalMARSRZ;
                                Materials.MaterialFunctions(mat_id).calc_dissipation = &MARSRZDissipationModel::calc_dissipation;
                            });
                            break;                        
                        default:
                            std::cout << "ERROR: invalid dissipation input: " << dissipation_model << std::endl;
                            throw std::runtime_error("**** Dissipation model Not Understood ****");
                            break;
                    } // end switch

                } 
                else{
                    std::cout << "ERROR: invalid disspation type input: " << dissipation_model << std::endl;
                    throw std::runtime_error("**** Dissipation model Not Understood ****");
                    break;
                } // end if

            } // erosion model variables
            //
            else if (a_word.compare("erode_tension_val") == 0) {
                double erode_tension_val = root["materials"][m_id]["material"]["erode_tension_val"].As<double>();

                RUN({
                    Materials.MaterialFunctions(mat_id).erode_tension_val = erode_tension_val;
                });
            } // erode_tension_val
            else if (a_word.compare("erode_density_val") == 0) {
                double erode_density_val = root["materials"][m_id]["material"]["erode_density_val"].As<double>();

                RUN({
                    Materials.MaterialFunctions(mat_id).erode_density_val = erode_density_val;
                });
            } // erode_density_val
            
            // exact the eos_global_vars
            else if (a_word.compare("eos_global_vars") == 0) {
                Yaml::Node & mat_global_vars_yaml = root["materials"][m_id]["material"][a_word];

                size_t num_global_vars = mat_global_vars_yaml.Size();
                
                RUN({ 
                    Materials.num_eos_global_vars(mat_id) = num_global_vars;
                });

                if(num_global_vars>100){
                    throw std::runtime_error("**** Per material, the code only supports up to 100 eos global vars in the input file ****");
                } // end check on num_global_vars

                // store the global eos model parameters
                for (int global_var_id = 0; global_var_id < num_global_vars; global_var_id++) {
                    double eos_var = root["materials"][m_id]["material"]["eos_global_vars"][global_var_id].As<double>();
                    

                    RUN({
                        tempGlobalEOSVars(mat_id, global_var_id) = eos_var;
                    });

                } // end loop over global vars
            } // "eos_global_vars"
            
            // exact the strength_global_vars
            else if (a_word.compare("strength_global_vars") == 0) {
                Yaml::Node & mat_global_vars_yaml = root["materials"][m_id]["material"][a_word];

                size_t num_global_vars = mat_global_vars_yaml.Size();
                
                RUN({ 
                    Materials.num_strength_global_vars(mat_id) = num_global_vars;
                });

                if(num_global_vars>100){
                    throw std::runtime_error("**** Per material, the code only supports up to 100 strength global vars in the input file ****");
                } // end check on num_global_vars

                // store the global strength model parameters
                for (int global_var_id = 0; global_var_id < num_global_vars; global_var_id++) {
                    double strength_var = root["materials"][m_id]["material"]["strength_global_vars"][global_var_id].As<double>();
                    
                    RUN({
                        tempGlobalStrengthVars(mat_id,global_var_id) = strength_var;
                    });

                } // end loop over global vars
            } // "strength_global_vars"
            else if (a_word.compare("dissipation_global_vars") == 0) {
                Yaml::Node & mat_global_vars_yaml = root["materials"][m_id]["material"][a_word];

                size_t num_global_vars = mat_global_vars_yaml.Size();

                RUN({ 
                    Materials.num_dissipation_global_vars(mat_id) = num_global_vars;
                });

                if(num_global_vars<6){
                    throw std::runtime_error("**** Per material, must specify 6 dissipation global vars in the input file ****");
                } // end check on num_global_vars

                if(num_global_vars>10){
                    throw std::runtime_error("**** Per material, the code only supports up to 10 dissipation global vars in the input file ****");
                } // end check on num_global_vars

                // store the global eos model parameters
                for (int global_var_id = 0; global_var_id < num_global_vars; global_var_id++) {
                    double dissipation_var = root["materials"][m_id]["material"]["dissipation_global_vars"][global_var_id].As<double>();
                    
                    RUN({
                        tempGlobalDissipationVars(mat_id, global_var_id) = dissipation_var;
                    });

                } // end loop over global vars

            } // end else if
            //
            // print and error because text is unknown
            else {
                std::cout << "ERROR: invalid input: " << a_word << std::endl;
                std::cout << "Valid options are: " << std::endl;
            
                for (const auto& element : str_material_inps) {
                    std::cout << element << std::endl;
                }
                throw std::runtime_error("**** Material Input Not Understood ****");
            }
        } // end for words in material
    } // end loop over materials

    // allocate ragged right memory to hold the model global variables
    Materials.eos_global_vars = RaggedRightArrayKokkos <double> (Materials.num_eos_global_vars, "Materials.eos_global_vars");
    Materials.strength_global_vars = RaggedRightArrayKokkos <double> (Materials.num_strength_global_vars, "Materials.strength_global_vars");
    Materials.dissipation_global_vars = RaggedRightArrayKokkos <double> (Materials.num_dissipation_global_vars, "Materials.dissipation_global_vars");


    // save the global variables
    FOR_ALL(mat_id, 0, num_materials, {
        
        for (size_t var_lid=0; var_lid<Materials.num_eos_global_vars(mat_id); var_lid++){
            Materials.eos_global_vars(mat_id, var_lid) = tempGlobalEOSVars(mat_id, var_lid);
        } // end for eos var_lid

        for (size_t var_lid=0; var_lid<Materials.num_strength_global_vars(mat_id); var_lid++){
            Materials.strength_global_vars(mat_id, var_lid) = tempGlobalStrengthVars(mat_id, var_lid);
        } // end for strength var_lid


        for (size_t var_lid=0; var_lid<Materials.num_dissipation_global_vars(mat_id); var_lid++){
            Materials.dissipation_global_vars(mat_id, var_lid) = tempGlobalDissipationVars(mat_id, var_lid);
        } // end for strength var_lid

    }); // end for loop over materials



    // set defaults, which are no models
    FOR_ALL(mat_id, 0, num_materials, {

        // default dissipation model is no dissipation
        if (Materials.MaterialEnums(mat_id).DissipationModels == model::noDissipation){

            // set the fcn pointer
            Materials.MaterialFunctions(mat_id).calc_dissipation = &NoDissipationModel::calc_dissipation;

        } // end if


        // if the following is true, stop simulation, parameters are models must be specified

        // --- EOS checks ---
        if (Materials.num_eos_global_vars(mat_id)>0 && Materials.MaterialFunctions(mat_id).calc_pressure == NULL) {
            Kokkos::abort("\n********************************************************************************************\n"
                            "ERROR: \n"
                            "EOS parameters were given, but no EOS model was not specified.  Both must be specified. \n"
                            "********************************************************************************************\n");
        }
        if (Materials.num_eos_global_vars(mat_id)>0 && Materials.MaterialEnums(mat_id).EOSModels == model::noEOS) {
            Kokkos::abort("\n********************************************************************************************\n"
                            "ERROR: \n"
                            "EOS parameters were given, but the 'no_eos' model was specified.  An EOS must be specified. \n"
                            "********************************************************************************************\n");
        }
        if (Materials.num_eos_global_vars(mat_id)==0 && Materials.MaterialFunctions(mat_id).calc_pressure != NULL) {

            if(Materials.MaterialEnums(mat_id).EOSModels != model::noEOS){ 
            Kokkos::abort("\n********************************************************************************************\n"
                            "ERROR: \n"
                            "EOS model was specified, but no EOS parameters were given. Both must be specified. \n"
                            "********************************************************************************************\n");
            }

        }


        // --- Dissipation checks ---
        if (Materials.num_dissipation_global_vars(mat_id)>0 && Materials.MaterialFunctions(mat_id).calc_dissipation == NULL) {
            Kokkos::abort("\n********************************************************************************************\n"
                            "ERROR: \n"
                            "Dissipation parameters were given, but no dissipation model was not specified.  \n"
                            "Both must be specified. \n"
                            "********************************************************************************************\n");
        }
        if (Materials.num_dissipation_global_vars(mat_id)>0 && Materials.MaterialEnums(mat_id).DissipationModels == model::noDissipation) {
            Kokkos::abort("\n********************************************************************************************\n"
                            "ERROR: \n"
                            "Dissipation parameters were given, but the 'no_dissipation' model was specified.\n " 
                            "A dissipation model must be specified. \n"
                            "********************************************************************************************\n");
        }
        if (Materials.num_dissipation_global_vars(mat_id)==0 && Materials.MaterialFunctions(mat_id).calc_dissipation != NULL) {

            if(Materials.MaterialEnums(mat_id).DissipationModels != model::noDissipation){            
            Kokkos::abort("\n********************************************************************************************\n"
                            "ERROR: \n"
                            "A dissipation model was specified, but no dissipation parameters were given. \n"
                            "Both must be specified. \n"
                            "********************************************************************************************\n");
            }
        }

    }); // end if check

} // end of function to parse material information
