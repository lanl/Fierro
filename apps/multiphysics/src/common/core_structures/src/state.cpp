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

#include "state.hpp"
#include <iostream>
#include <stdexcept>

void fillGaussState_t::initialize(size_t num_gauss_points, 
                    size_t max_mat_storage_in_elem, 
                    size_t num_dims, 
                    std::vector<fill_gauss_state> fill_gauss_states)
{           

    max_mats_in_elem = max_mat_storage_in_elem;

    for (auto field : fill_gauss_states){
        switch(field){
            case fill_gauss_state::density:
                if (den.size() == 0) this->den = DCArrayKokkos<double>(num_gauss_points, max_mats_in_elem, "fill_gauss_point_density");
                break;
            case fill_gauss_state::stress:
                if (stress.size() == 0) this->stress = DCArrayKokkos<double>(num_gauss_points, max_mats_in_elem, num_dims, num_dims, "fill_gauss_point_stress");
                break;
            case fill_gauss_state::strain:
                if (strain.size() == 0) this->strain = DCArrayKokkos<double>(num_gauss_points, max_mats_in_elem, num_dims, num_dims, "fill_gauss_point_strain");
                break;
            case fill_gauss_state::elastic_modulii:
                if (elastic_modulii.size() == 0) this->elastic_modulii = DCArrayKokkos<double>(num_gauss_points, max_mats_in_elem, 3, "fill_gauss_point_elastic_modulii");
                break;
            case fill_gauss_state::shear_modulii:
                if (shear_modulii.size() == 0) this->shear_modulii = DCArrayKokkos<double>(num_gauss_points, max_mats_in_elem, 3, "fill_gauss_point_shear_modulii");
                break;
            case fill_gauss_state::poisson_ratios:
                if (poisson_ratios.size() == 0) this->poisson_ratios = DCArrayKokkos<double>(num_gauss_points, max_mats_in_elem, 3, "fill_gauss_point_poisson_ratios");
                break;
            case fill_gauss_state::specific_internal_energy:
                if (sie.size() == 0) this->sie = DCArrayKokkos<double>(num_gauss_points, max_mats_in_elem, "fill_gauss_point_sie");
                if (use_sie.size() == 0){ 
                    this->use_sie = DCArrayKokkos<bool>(num_gauss_points, max_mats_in_elem, "fill_gauss_point_use_sie");
                    use_sie.set_values(false);
                }
                break;
            case fill_gauss_state::internal_energy:
                if (sie.size() == 0) this->ie = DCArrayKokkos<double>(num_gauss_points, max_mats_in_elem, "fill_gauss_point_ie");
                if (use_sie.size() == 0){ 
                    this->use_sie = DCArrayKokkos<bool>(num_gauss_points, max_mats_in_elem, "fill_gauss_point_use_sie");
                    use_sie.set_values(false);
                }
                break;
            case fill_gauss_state::thermal_conductivity:
                if (thermal_conductivity.size() == 0) this->thermal_conductivity = DCArrayKokkos<double>(num_gauss_points, max_mats_in_elem, "fill_gauss_point_thermal_conductivity");
                break;
            case fill_gauss_state::specific_heat:
                if (specific_heat.size() == 0) this->specific_heat = DCArrayKokkos<double>(num_gauss_points, max_mats_in_elem, "fill_gauss_point_specific_heat");
                break;
            case fill_gauss_state::level_set:
                if (level_set.size() == 0) this->level_set = DCArrayKokkos<double>(num_gauss_points,max_mats_in_elem, "fill_gauss_level_set");
                break;
            default:
                std::cout<<"Desired Gauss point fill state not understood in initialize"<<std::endl;
                throw std::runtime_error("**** Error in State Fill Field Name ****");
        } // end switch
    }
} // end method

void fillElemState_t::initialize(size_t num_elems, 
                    size_t max_mat_storage_in_elem,
                    size_t num_mats)
{
    this-> max_mats_in_elem = max_mat_storage_in_elem;

    if (mat_volfrac.size() == 0){
        this->mat_volfrac = DCArrayKokkos<double>(num_elems, max_mats_in_elem, "elem_mat_volfrac");
    }

    if (geo_volfrac.size() == 0){
        this->geo_volfrac = DCArrayKokkos<double>(num_elems, max_mats_in_elem, "elem_geo_volfrac");
    }

    // voxel_elem_mat_id is allocated in the voxel file read
    
} // end method

bool check_fill_node_states(
        const std::vector<fill_node_state>& required_fill_node_states,
        const std::vector<fill_node_state>& fill_node_states)
{
    
    bool states_filled = false;
    size_t count = 0;

    // check if all required fields were set
    for (auto required_field : required_fill_node_states){
        
        for (auto field : fill_node_states){
            if (field == required_field){
                count++;  // it has the required field
            }
        } // end user specified fields loop

    } // end for required field
    
    if (count == required_fill_node_states.size()){
        states_filled = true;
    }

    return states_filled;
} // end check_fill_states

bool check_fill_mat_states(
        const std::vector<fill_gauss_state>& required_fill_matpt_states,
        const std::vector<fill_gauss_state>& fill_gauss_states)
{
    
    bool states_filled = false;
    size_t count = 0;

    // check if all required fields were set
    for (auto required_field : required_fill_matpt_states){
        
        for (auto field : fill_gauss_states){
            if (field == required_field){
                count++;  // it has the required field
            }
        } // end user specified fields loop

    } // end for required field
    
    if (count == required_fill_matpt_states.size()){
        states_filled = true;
    }

    return states_filled;
} // end check_fill_states

void node_t::initialize(size_t num_nodes, size_t num_dims, std::vector<node_state> node_states)
{
    for (auto field : node_states){
        switch(field){
            case node_state::coords:
                if (coords.size() == 0) this->coords = MPICArrayKokkos<double>(num_nodes, num_dims, "node_coordinates");
                if (coords_n0.size() == 0) this->coords_n0 = MPICArrayKokkos<double>(num_nodes, num_dims, "node_coordinates_n0");
                break;
            case node_state::coords_t0:
                if (coords_t0.size() == 0) this->coords_t0 = MPICArrayKokkos<double>(num_nodes, num_dims, "node_reference_coordinates");
                break;
            case node_state::velocity:
                if (vel.size() == 0) this->vel = MPICArrayKokkos<double>(num_nodes, num_dims, "node_velocity");
                if (vel_n0.size() == 0) this->vel_n0 = MPICArrayKokkos<double>(num_nodes, num_dims, "node_velocity_n0");
                break;
            case node_state::force:
                if (force.size() == 0) this->force = DCArrayKokkos<double>(num_nodes, num_dims, "node_force");
                break;
            case node_state::mass:
                if (mass.size() == 0) this->mass = DCArrayKokkos<double>(num_nodes, "node_mass");
                break;
            case node_state::temp:
                if (temp.size() == 0) this->temp = MPICArrayKokkos<double>(num_nodes, "node_temp");
                if (temp_n0.size() == 0) this->temp_n0 = MPICArrayKokkos<double>(num_nodes, "node_temp_n0");
                break;
            case node_state::heat_transfer:
                if (q_transfer.size() == 0) this->q_transfer = DCArrayKokkos<double>(num_nodes, "node_q_transfer");
                break;
            case node_state::gradient_level_set:
                if (gradient_level_set.size() == 0) this->gradient_level_set = DCArrayKokkos<double>(num_nodes, num_dims, "node_grad_levelset");
                break;
            case node_state::displacement:
                if (displacement.size() == 0) this->displacement = MPICArrayKokkos<double>(num_nodes, num_dims, "node_force_displacement");
                break;
            default:
                std::cout<<"Desired node state not understood in node_t initialize"<<std::endl;
                throw std::runtime_error("**** Error in State Field Name ****");
        }
    }
} // end method

void node_t::initialize(size_t num_nodes, size_t num_dims, std::vector<node_state> node_states, CommunicationPlan& comm_plan)
{
    for (auto field : node_states){
        switch(field){
            case node_state::coords:
                if (coords.size() == 0){
                    this->coords = MPICArrayKokkos<double>(num_nodes, num_dims, "node_coordinates");
                    this->coords.initialize_comm_plan(comm_plan);
                }
                if (coords_n0.size() == 0){
                    this->coords_n0 = MPICArrayKokkos<double>(num_nodes, num_dims, "node_coordinates_n0");
                    this->coords_n0.initialize_comm_plan(comm_plan);
                }
                break;

            case node_state::velocity:
                if (vel.size() == 0){
                    this->vel = MPICArrayKokkos<double>(num_nodes, num_dims, "node_velocity");
                    this->vel.initialize_comm_plan(comm_plan);
                } 
                if (vel_n0.size() == 0){
                    this->vel_n0 = MPICArrayKokkos<double>(num_nodes, num_dims, "node_velocity_n0");
                    this->vel_n0.initialize_comm_plan(comm_plan);
                }
                break;

            case node_state::temp:
                if (temp.size() == 0){
                    this->temp = MPICArrayKokkos<double>(num_nodes, "node_temp");
                    this->temp.initialize_comm_plan(comm_plan);
                }
                if (temp_n0.size() == 0){
                    this->temp_n0 = MPICArrayKokkos<double>(num_nodes, "node_temp_n0");
                    this->temp_n0.initialize_comm_plan(comm_plan);
                }

            default:
                std::cout<<"Desired node state not understood in node_t initialize with communication plan"<<std::endl;
                throw std::runtime_error("**** Error in State Field Name ****");
        }
    }
} // end method

void GaussPoint_t::initialize(size_t num_gauss_pnts, size_t num_dims, std::vector<gauss_pt_state> gauss_pt_states)
{

    for (auto field : gauss_pt_states){
        switch(field){
            case gauss_pt_state::volume:
                if (vol.size() == 0) this->vol = DCArrayKokkos<double>(num_gauss_pnts, "gauss_point_volume");
                break;
            case gauss_pt_state::divergence_velocity:
                if (div.size() == 0) this->div = DCArrayKokkos<double>(num_gauss_pnts, "gauss_point_div");
                break;
            case gauss_pt_state::gradient_velocity:
                // Note: all current solvers treat velocity gradient as a 3d tensor, even the rz solver.
                if (vel_grad.size() == 0) this->vel_grad = DCArrayKokkos<double>(num_gauss_pnts, 3, 3, "gauss_point_vel_grad");
                break;
            case gauss_pt_state::level_set:
                if (level_set.size() == 0) this->level_set = MPICArrayKokkos<double>(num_gauss_pnts, "gauss_point_level_set");
                if (level_set_n0.size() == 0) this->level_set_n0 = MPICArrayKokkos<double>(num_gauss_pnts, "gauss_point_level_set_n0");
                break;
            case gauss_pt_state::shock_detector:
                if (shock_detector.size() == 0) this->shock_detector = MPICArrayKokkos<double>(num_gauss_pnts, "gauss_point_shock_detector");
                break;
            default:
                std::cout<<"Desired gauss point state not understood in GaussPoint_t initialize"<<std::endl;
                throw std::runtime_error("**** Error in State Field Name ****");
        }
    }
} // end method

void GaussPoint_t::initialize(size_t num_gauss_pnts, size_t num_dims, std::vector<gauss_pt_state> gauss_pt_states, CommunicationPlan& comm_plan)
{
    for (auto field : gauss_pt_states){
        switch(field){
            case gauss_pt_state::shock_detector:
                if (shock_detector.size() == 0){
                    this->shock_detector = MPICArrayKokkos<double>(num_gauss_pnts, "gauss_point_shock_detector");
                    this->shock_detector.initialize_comm_plan(comm_plan);
                }
                break;
            case gauss_pt_state::level_set:

            default:
                std::cout<<"Desired MPI distributed Gauss point state not understood in GaussPoint_t initialize with communication plan"<<std::endl;
                throw std::runtime_error("**** Error in State Field Name ****");
        }
    }
} // end method

void MeshtoMaterialMap_t::initialize(size_t num_elem_max, size_t num_mats_per_elem_max)
{
    if (num_mats_in_elem.size() == 0){
        this->num_mats_in_elem = DCArrayKokkos<size_t>(num_elem_max, "num_mats_in_elem");
        this->num_mats_in_elem.set_values(0); // initialize all elems to storing 0 materials
        this->num_mats_in_elem.update_host(); // copy from GPU to CPU
    }
    if (mats_in_elem.size() == 0){
        this->mats_in_elem = DCArrayKokkos<size_t>(num_elem_max, num_mats_per_elem_max, "mat_id_in_elem");
    }
    if (mat_elems_in_elem.size() == 0){
        this->mat_elems_in_elem = DCArrayKokkos<size_t>(num_elem_max, num_mats_per_elem_max, "mat_elem_in_elem");
    }
    
} // end method

void MaterialToMeshMap_t::initialize()
{
    if (elem_in_mat_elem.size() == 0){ 
        this->elem_in_mat_elem = DRaggedRightArrayKokkos<size_t>(this->num_mat_elems_buffer, "material_space_to_elem");
    }

} // end method

void MaterialToMeshMap_t::initialize_num_mats(size_t num_mats)
{
    // Note: num_mat_elems is allocated in problem setup
    if (num_mat_elems.size() == 0){
        this->num_mat_elems = DCArrayKokkos <size_t> (num_mats, "num_mat_elems"); 
    }

    // Note: num_mat_elems_buffer is allocated in problem setup, the values are set in region_fill.cpp routine
    if (num_mat_elems_buffer.size() == 0){
        this->num_mat_elems_buffer = DCArrayKokkos <size_t> (num_mats, "num_mat_elems_with_buffer"); 
    }

} // end method

void MaterialPoint_t::initialize_num_mats(size_t num_mats)
{
    // Note: num_material_points is allocated in problem setup
    if (num_material_points.size() == 0){
        this->num_material_points = DCArrayKokkos <size_t> (num_mats, "num_material_points"); 
    }

    // Note: num_material_points_buffer is allocated in problem setup, the values are set in region_fill.cpp routine
    if (num_material_points_buffer.size() == 0){
        this->num_material_points_buffer = DCArrayKokkos <size_t> (num_mats, "num_material_points_with_buffer"); 
    }

} // end method

void MaterialPoint_t::initialize(size_t num_dims, std::vector<material_pt_state> material_pt_states)
{
    // WARNING: this should be moved to for loop over mat_pt_states
    if (temp_grad.size() == 0){
        this->temp_grad = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, 3, "material_point_temperature_gradient"); 
    }

    for (auto field : material_pt_states){
        switch(field){
            case material_pt_state::density:
                if (den.size() == 0) this->den = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_density");
                break;
            case material_pt_state::pressure:
                if (pres.size() == 0) this->pres = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_pressure");
                break;
            case material_pt_state::stress:
                // Note: all current solvers treat stress as a 3d tensor, even the rz solver.
                if (stress.size() == 0) this->stress = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, 3, 3, "material_point_stress");  
                if (stress_n0.size() == 0) this->stress_n0 = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, 3, 3, "material_point_stress_n0"); 
                break;
            case material_pt_state::strain:
                if (strain.size() == 0) this->strain = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, num_dims, num_dims, "material_point_strain");  
                break;
            case material_pt_state::elastic_modulii:
                if (elastic_modulii.size() == 0) this->elastic_modulii = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, 3, "material_elastic_modulii");
                break;
            case material_pt_state::shear_modulii:
                if (shear_modulii.size() == 0) this->shear_modulii = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, 3, "material_shear_modulii");
                break;
            case material_pt_state::poisson_ratios:
                if (poisson_ratios.size() == 0) this->poisson_ratios = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, 3, "material_poisson_ratios");
                break;
            case material_pt_state::sound_speed:
                if (sspd.size() == 0) this->sspd = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_sspd");
                break;
            case material_pt_state::mass:
                if (mass.size() == 0) this->mass = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_mass");
                break;
            case material_pt_state::volume_fraction:
                if (mat_volfrac.size() == 0) this->mat_volfrac = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_mat_volfrac");
                if (geo_volfrac.size() == 0) this->geo_volfrac = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_geo_volfrac");
                // changes in volume fraction
                if (delta_mat_volfrac.size() == 0) this->delta_mat_volfrac = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_mat_volfrac_delta");
                if (delta_geo_volfrac.size() == 0) this->delta_geo_volfrac = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_geo_volfrac_delta");
                break;
            case material_pt_state::specific_internal_energy:
                if (sie.size() == 0)  this->sie = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_sie");
                if (sie_n0.size() == 0) this->sie_n0 = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_sie_n0");
                break;
            case material_pt_state::heat_flux:
                if (q_flux.size() == 0) this->q_flux = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, num_dims, "material_point_q_flux");
                break;
            case material_pt_state::thermal_conductivity:
                if (conductivity.size() == 0) this->conductivity = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_thermal_conductivity");
                break;
            case material_pt_state::specific_heat:
                if (specific_heat.size() == 0) this->specific_heat = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_specific_heat");
                break;
            case material_pt_state::eroded_flag:
                if (eroded.size() == 0) this->eroded = DRaggedRightArrayKokkos<bool>(this->num_material_points_buffer, "material_point_eroded");
                break;
            default:
                std::cout<<"Desired material point state not understood in MaterialPoint_t initialize"<<std::endl;
                throw std::runtime_error("**** Error in State Field Name ****");
        }
    }
} // end method

void MaterialZone_t::initialize_num_mats(size_t num_mats)
{
    // Note: num_material_zones is allocated in problem setup
    if (num_material_zones.size() == 0){
        this->num_material_zones = DCArrayKokkos <size_t> (num_mats, "num_material_zones"); 
    }

    // Note: num_material_zones_buffer is allocated in problem setup, the values are set in region_fill.cpp routine
    if (num_material_zones_buffer.size() == 0){
        this->num_material_zones_buffer = DCArrayKokkos <size_t> (num_mats, "num_material_zones_with_buffer"); 
    }

} // end method

void MaterialZone_t::initialize_Pn(std::vector<material_zone_state> material_zone_states)
{
    for (auto field : material_zone_states){
        switch(field){
            case material_zone_state::specific_internal_energy:
                if (sie.size() == 0) this->sie = DRaggedRightArrayKokkos<double>(this->num_material_zones_buffer, "material_zone_sie");
                if (sie_n0.size() == 0) this->sie_n0 = DRaggedRightArrayKokkos<double>(this->num_material_zones_buffer, "material_zone_sie_n0");
                break;
            default:
                std::cout<<"Desired material zone state not understood in MaterialZone_t initialize"<<std::endl;
                throw std::runtime_error("**** Error in State Field Name ****");
        }
    }
    
} // end method

void MaterialCorner_t::initialize_num_mats(size_t num_mats)
{
    // Note: num_material_corners is allocated in problem setup
    if (num_material_corners.size() == 0){
        this->num_material_corners = DCArrayKokkos <size_t> (num_mats, "num_material_corners"); 
    }

    // Note: num_material_corners_buffer is allocated in problem setup, the values are set in region_fill.cpp routine
    if (num_material_corners_buffer.size() == 0){
        this->num_material_corners_buffer = DCArrayKokkos <size_t> (num_mats, "num_material_corners_with_buffer"); 
    }

} // end method

void MaterialCorner_t::initialize(size_t num_dims, std::vector<material_corner_state> material_corner_states)
{
    for (auto field : material_corner_states){
        switch(field){
            case material_corner_state::force:
                if (force.size() == 0) this->force = DRaggedRightArrayKokkos<double>(this->num_material_corners_buffer, num_dims, "material_corner_force");
                break;
            case material_corner_state::heat_transfer:
                if (q_transfer.size() == 0) this->q_transfer =DRaggedRightArrayKokkos<double>(this->num_material_corners_buffer, "material_corner_q_transfer"); 
                break;
            case material_corner_state::density:
                if (density.size() == 0) this->density = DRaggedRightArrayKokkos<double>(this->num_material_corners_buffer, "material_corner_density");
                break;
            case material_corner_state::specific_internal_energy:
                if (specific_internal_energy.size() == 0) this->specific_internal_energy = DRaggedRightArrayKokkos<double>(this->num_material_corners_buffer, "material_corner_specific_internal_energy");
                break;
            case material_corner_state::kinetic_energy:
                if (kinetic_energy.size() == 0) this->kinetic_energy = DRaggedRightArrayKokkos<double>(this->num_material_corners_buffer, "material_corner_kinetic_energy");
                break;
            case material_corner_state::velocity:
                if (velocity.size() == 0) this->velocity = DRaggedRightArrayKokkos<double>(this->num_material_corners_buffer, num_dims, "material_corner_velocity");
                break;
            default:
                std::cout<<"Desired material corner state not understood in MaterialCorner_t initialize"<<std::endl;
                throw std::runtime_error("**** Error in State Field Name ****");
        }
    }
} // end method

void corner_t::initialize(size_t num_corners, size_t num_dims, std::vector<corner_state> corner_states)
{

    for (auto field : corner_states){
        switch(field){
            case corner_state::force:
                if (force.size() == 0) this->force = DCArrayKokkos<double>(num_corners, num_dims, "corner_force");
                break;
            case corner_state::mass:
                if (mass.size() == 0) this->mass  = DCArrayKokkos<double>(num_corners, "corner_mass");
                break;
            case corner_state::heat_transfer:
                if (q_transfer.size() == 0) this->q_transfer = DCArrayKokkos<double>(num_corners, "corner_q_transfer"); 
                break;
            case corner_state::normal:
                if (normal.size() == 0) this->normal = DCArrayKokkos<double>(num_corners, num_dims, "corner_normal");
                break;
            case corner_state::volume:
                if (volume.size() == 0) this->volume  = DCArrayKokkos<double>(num_corners, "corner_volume");
                break;
            default:
                std::cout<<"Desired corner state not understood in corner_t initialize"<<std::endl;
                throw std::runtime_error("**** Error in State Field Name ****");
        }
    }

    
    
    
} // end method
