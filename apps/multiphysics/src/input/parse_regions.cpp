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
#include "parse_regions.hpp"

#include "state.hpp"




// =================================================================================
//    Parse Fill regions
// =================================================================================
void parse_regions(Yaml::Node& root, 
                   DCArrayKokkos<size_t>& reg_fills_in_solver,  
                   DCArrayKokkos<size_t>& num_reg_fills_in_solver,
                   CArrayKokkos<RegionFill_t>& region_fills, 
                   CArray<RegionFill_host_t>&  region_fills_host,
                   const size_t num_solvers)

{
    bool verbose = false; 
    // allocate memory
    num_reg_fills_in_solver = DCArrayKokkos<size_t>(num_solvers, "sim_param.region_setup.num_reg_fills_in_solver");
    num_reg_fills_in_solver.set_values(0);
    Kokkos::fence();
    num_reg_fills_in_solver.update_host(); // initiallizing host side to 0



    Yaml::Node& region_yaml = root["regions"];

    size_t num_regions = region_yaml.Size();

    printf("Number of regions = %zu \n", num_regions);

    // reg_fills_in_solver(solver_id, fill_lid) = fill_id
    reg_fills_in_solver = DCArrayKokkos<size_t>(num_solvers, num_regions, "sim_param.region_setup.reg_fills_in_solver");
    region_fills = CArrayKokkos<RegionFill_t>(num_regions , "sim_param.region_setup.region_fills");
    region_fills_host = CArray<RegionFill_host_t>(num_regions); 


    // a check on region_id not being specified more than once or not at all
    CArray <bool> check_reg_ids(num_regions);
    check_reg_ids.set_values(false);


    // loop over the fill regions specified in yaml file
    for (int r_id = 0; r_id < num_regions; r_id++) {

        // read the variables names
        Yaml::Node& inps_yaml = root["regions"][r_id]["region"];

        // get the region variables names set by the user
        std::vector<std::string> user_str_region_inps;

        // extract words from the input file and validate they are correct
        validate_inputs(inps_yaml, user_str_region_inps, str_region_inps, region_required_inps);


        // loop over the words in r_id region input definition and find the region id
        int reg_id = -1;
        for (auto& a_word : user_str_region_inps) {

            if (a_word.compare("id") == 0) {
                reg_id = root["regions"][r_id]["region"]["id"].As<int>(); // the region id

                if (reg_id<0 || reg_id>=num_regions){
                    std::cout << "ERROR: invalid region_id specified in the region definition " << std::endl;
            
                    throw std::runtime_error("**** region_id is out of bounds ****");
                } // end check on m_id range

                if (check_reg_ids(reg_id) == true){
                    std::cout << "ERROR: region_id = " << reg_id << " was already specified "<< std::endl;
                    throw std::runtime_error("**** Multiple regions used the same region_id ****");
                }
                else {
                    check_reg_ids(reg_id) = true;
                } // end check on reg_id

            } // end if id
        } // end loop over all region inputs in r_id block

        if (reg_id<0){
            std::cout << "ERROR: region_id must be specified in each region definition " << std::endl;
            
            throw std::runtime_error("**** region_id is missing ****");
        } // end check on reg_id being specified



        // loop over the words in the material input definition
        for (auto& a_word : user_str_region_inps) {
            
            
            if (a_word.compare("id") == 0) {
                // do nothing
                // this id was read in an earlier loop
            }
            //extract solver id, currently not used by Fierro's solvers
            else if (a_word.compare("solver_id") == 0) {
                int solver_id = root["regions"][r_id]["region"][a_word].As<int>();
                
                // get the local id for filling this region
                size_t fill_lid = num_reg_fills_in_solver.host(solver_id);

                // save the fill_id, which is the reg_id
                reg_fills_in_solver.host(solver_id, fill_lid) = reg_id; 
                num_reg_fills_in_solver.host(solver_id) ++;

                RUN({
                    region_fills(reg_id).solver_id = solver_id;
                });
            } // solver_id
            //
            else if (a_word.compare("volume") == 0) {

                // -----
                // loop over the sub fields under volume
                // -----
                Yaml::Node& inps_subfields_yaml = root["regions"][r_id]["region"]["volume"];

                // get the bc_geometery variables names set by the user
                std::vector<std::string> user_region_volume_inputs;
                
                // extract words from the input file and validate they are correct
                validate_inputs(inps_subfields_yaml, user_region_volume_inputs, str_region_volume_inps, region_volume_required_inps);


                // loop over the subfield words
                for(auto& a_subfield_word : user_region_volume_inputs){ 

                    if (a_subfield_word.compare("radius1") == 0) {
                        // inner radius of sphere/cylinder

                        double radius1 = root["regions"][r_id]["region"]["volume"]["radius1"].As<double>();

                        RUN({
                            region_fills(reg_id).radius1 = radius1;
                        });
                    } // radius1
                    else if (a_subfield_word.compare("radius2") == 0) {
                        // outer radius of sphere/cylinder

                        double radius2 = root["regions"][r_id]["region"]["volume"]["radius2"].As<double>();

                        RUN({
                            region_fills(reg_id).radius2 = radius2;
                        });
                    } // radius2
                    else if (a_subfield_word.compare("x1") == 0) {
                        // inner plane

                        double x1 = root["regions"][r_id]["region"]["volume"]["x1"].As<double>();

                        RUN({
                            region_fills(reg_id).x1 = x1;
                        });
                    } // x1
                    else if (a_subfield_word.compare("x2") == 0) {
                        // outer plane

                        double x2 = root["regions"][r_id]["region"]["volume"]["x2"].As<double>();

                        RUN({
                            region_fills(reg_id).x2 = x2;
                        });
                    } // x2
                    else if (a_subfield_word.compare("y1") == 0) {
                        // inner plane

                        double y1 = root["regions"][r_id]["region"]["volume"]["y1"].As<double>();

                        RUN({
                            region_fills(reg_id).y1 = y1;
                        });
                    } // y1
                    else if (a_subfield_word.compare("y2") == 0) {
                        // outer plane

                        double y2 = root["regions"][r_id]["region"]["volume"]["y2"].As<double>();

                        RUN({
                            region_fills(reg_id).y2 = y2;
                        });
                    } // y2
                    else if (a_subfield_word.compare("z1") == 0) {
                        // inner plane

                        double z1 = root["regions"][r_id]["region"]["volume"]["z1"].As<double>();

                        RUN({
                            region_fills(reg_id).z1 = z1;
                        });
                    } // z1
                    else if (a_subfield_word.compare("z2") == 0) {
                        // outer plane

                        double z2 = root["regions"][r_id]["region"]["volume"]["z2"].As<double>();

                        RUN({
                            region_fills(reg_id).z2 = z2;
                        });
                    } // z2
                    else if (a_subfield_word.compare("half_angle") == 0) {
                        // half angle

                        double half_angle = root["regions"][r_id]["region"]["volume"]["half_angle"].As<double>();

                        RUN({
                            region_fills(reg_id).half_angle = half_angle;
                        });
                    } // z2
                    else if (a_subfield_word.compare("scale_x") == 0) {
                        // outer plane

                        double scale_x = root["regions"][r_id]["region"]["volume"]["scale_x"].As<double>();

                        // on the host side because it relates to reading a mesh file
                        region_fills_host(reg_id).scale_x = scale_x;

                    } // scale_x
                    else if (a_subfield_word.compare("scale_y") == 0) {
                        // outer plane

                        double scale_y = root["regions"][r_id]["region"]["volume"]["scale_y"].As<double>();

                        // on the host side because it relates to reading a mesh file
                        region_fills_host(reg_id).scale_y = scale_y;

                    } // scale_y
                    else if (a_subfield_word.compare("scale_z") == 0) {
                        // outer plane

                        double scale_z = root["regions"][r_id]["region"]["volume"]["scale_z"].As<double>();

                        // on the host side because it relates to reading a mesh file
                        region_fills_host(reg_id).scale_z = scale_z;

                    } // scale_z
                    else if (a_subfield_word.compare("part_id") == 0) {
                        // part_id in 

                        int part_id = root["regions"][r_id]["region"]["volume"]["part_id"].As<int>();

                        RUN({
                            region_fills(reg_id).part_id = part_id;
                        });

                    } // scale_z
                    //
                    else if (a_subfield_word.compare("type") == 0) {

                        // region volume fill type
                        std::string type = root["regions"][r_id]["region"]["volume"]["type"].As<std::string>();

                        // set the velocity tag type
                        if (region_type_map.find(type) != region_type_map.end()) {
                        
                            // region_type_map[type] returns enum value, e.g., init_conds::velocity 
                            switch(region_type_map[type]){

                                case region::global:
                                    if (verbose) std::cout << "Setting volume fill type to global " << std::endl;
                                    region_fills_host(reg_id).volume = region::global;
                                    RUN({
                                        region_fills(reg_id).volume = region::global;
                                    });
                                    break;

                                case region::box:
                                    if (verbose) std::cout << "Setting volume fill type to box " << std::endl;
                                    region_fills_host(reg_id).volume = region::box;
                                    RUN({
                                        region_fills(reg_id).volume = region::box;
                                    });
                                    break;

                                case region::cylinder:
                                    if (verbose) std::cout << "Setting volume fill type to cylinder " << std::endl;
                                    region_fills_host(reg_id).volume = region::cylinder;
                                    RUN({
                                        region_fills(reg_id).volume = region::cylinder;
                                    });
                                    break;

                                case region::sphere:
                                    if (verbose) std::cout << "Setting volume fill type to sphere " << std::endl;
                                    region_fills_host(reg_id).volume = region::sphere;
                                    RUN({
                                        region_fills(reg_id).volume = region::sphere;
                                    });
                                    break;

                                case region::readVoxelFile:
                                    if (verbose) std::cout << "Setting volume fill type to readVoxelFile " << std::endl;
                                    region_fills_host(reg_id).volume = region::readVoxelFile;
                                    RUN({
                                        region_fills(reg_id).volume = region::readVoxelFile;
                                    });
                                    break;
                                
                                case region::readSTLFile:
                                    std::cout << "Setting volume fill type to readSTLFile " << std::endl;
                                    region_fills_host(reg_id).volume = region::readSTLFile;
                                    RUN({
                                        region_fills(reg_id).volume = region::readSTLFile;
                                    });
                                    break;

                                case region::readVTUFile:
                                    if (verbose) std::cout << "Setting volume fill type to readVTUFile " << std::endl;
                                    region_fills_host(reg_id).volume = region::readVTUFile;
                                    RUN({
                                        region_fills(reg_id).volume = region::readVTUFile;
                                    });
                                    break;
                                case region::no_volume:
                                    if (verbose) std::cout << "Setting volume fill type to none " << std::endl;
                                    region_fills_host(reg_id).volume = region::no_volume;
                                    RUN({
                                        region_fills(reg_id).volume = region::no_volume;
                                    });
                                    break;
                                default:
                                    region_fills_host(reg_id).volume = region::no_volume;
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
                        std::string path = root["regions"][r_id]["region"]["volume"]["file_path"].As<std::string>();

                        // absolute path to file or local to the director where exe is run
                        region_fills_host(reg_id).file_path = path;   // saving the absolute file path
        

                    } // end file path
                    //
                    else if (a_subfield_word.compare("unit_vector") == 0) {
                        std::string unit_vector = root["regions"][r_id]["region"]["volume"]["unit_vector"].As<std::string>();

                        // get the origin numbers, values are words
                        std::vector<std::string> numbers = exact_array_values(unit_vector, ",");

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
                        RUN({
                            region_fills(reg_id).unit_vector[0] = x1;
                            region_fills(reg_id).unit_vector[1] = y1;
                            region_fills(reg_id).unit_vector[2] = z1;
                        });
                    } // unit vector
                    else if (a_subfield_word.compare("origin") == 0) {
                        std::string origin = root["regions"][r_id]["region"]["volume"]["origin"].As<std::string>();

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
                        region_fills_host(reg_id).origin[0] = x1;
                        region_fills_host(reg_id).origin[1] = y1;
                        region_fills_host(reg_id).origin[2] = z1;
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
        // NOTE:
        // Each solver checks within the intialize function to see if the fills specified 
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
                
                // if the following is true, stop simulation; must add all mesh read options
                if (region_fills(reg_id).volume == region::readSTLFile) {
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
                if (region_fills(reg_id).volume != region::readVoxelFile &&
                    region_fills(reg_id).volume != region::readSTLFile){  
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
