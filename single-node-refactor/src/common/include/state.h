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
#ifndef STATE_H
#define STATE_H

#include "matar.h"

using namespace mtr;


enum class fill_node_state
{
    velocity,
    temperature
};

enum class fill_gauss_state
{
    density,
    stress,
    specific_internal_energy,
    internal_energy,
    elastic_modulii,
    shear_modulii,
    poisson_ratios,
    thermal_conductivity,
    specific_heat,
    level_set
};

//distributed vector type in use
using DistributedMap = TpetraPartitionMap<>;
template <typename T>
using DistributedDFArray = TpetraDFArray<T>;
template <typename T>
using DistributedDCArray = TpetraDCArray<T>;
template <typename T>
using CommPlan = TpetraLRCommunicationPlan<T>;


template <typename T>
void super_vector_initialization(DistributedDCArray<T> super_vector, DistributedDCArray<T> sub_vector, size_t extent);

template <typename T>
void super_vector_initialization(DistributedDFArray<T> super_vector, DistributedDFArray<T> sub_vector, size_t extent);

/////////////////////////////////////////////////////////////////////////////
///
/// \struct fillGaussState_t
///
/// \brief Stores state to setup a problem
///
/////////////////////////////////////////////////////////////////////////////
// Possible states, used to initialize fillState_t
struct fillGaussState_t
{
    size_t max_mats_in_elem;      ///< the max number of materials possible per element

    DCArrayKokkos<double> den;    ///< Gauss Point density
    DCArrayKokkos<double> sie;    ///< Gauss Point specific internal energy
    DCArrayKokkos<double> ie;     ///< Gauss Point extensive internal energy
    DCArrayKokkos<bool> use_sie;  ///< use sie to set sie, else use ie

    DCArrayKokkos<double> stress; ///< Gauss Point stress

    DCArrayKokkos<double> thermal_conductivity; ///< Thermal conductivity
    DCArrayKokkos<double> specific_heat;        ///< Specific Heat

    DCArrayKokkos<double> elastic_modulii;  ///<  Gauss Point elastic modulii Exx, Eyy, Ezz
    DCArrayKokkos<double> shear_modulii;    ///<  Gauss Point shear modulii Gxy, Gxz, Gyz
    DCArrayKokkos<double> poisson_ratios;   ///<  Gauss Point poisson ratios nu_xy, nu_xz, nu_yz

    DCArrayKokkos<double> level_set;        ///< level set


    // initialization method 
    void initialize(size_t num_gauss_points, 
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
}; // end Gauss fill states


/////////////////////////////////////////////////////////////////////////////
///
/// \struct fillElemState_t
///
/// \brief Stores state to setup of a problem
///
/////////////////////////////////////////////////////////////////////////////
// Possible states, used to initialize fillState_t
struct fillElemState_t
{
    size_t max_mats_in_elem;    ///< the max number of materials possible per element

    DCArrayKokkos<double> volfrac;     ///< element volume fraction
    DCArrayKokkos<double> geo_volfrac;  ///< element geometric (the part) volume fraction

    // arrays for building material index space:
    //    mat_id                     material ids in the element (num_elems, num_mats_saved)
    //    num_mats_saved_in_elem     material ids in the element (num_elems, num_mats_saved)
    //    num_elems_saved_for_mat    the number of elements the material resides in, (num_mats)
    // are in the MeshToMaterialMap struct
   
    
    // initialization method 
    void initialize(size_t num_elems, 
                    size_t max_mat_storage_in_elem,
                    size_t num_mats)
    {
        this-> max_mats_in_elem = max_mat_storage_in_elem;

        if (volfrac.size() == 0){
            this->volfrac = DCArrayKokkos<double>(num_elems, max_mats_in_elem, "elem_volfrac");
        }

        if (geo_volfrac.size() == 0){
            this->geo_volfrac = DCArrayKokkos<double>(num_elems, max_mats_in_elem, "elem_geo_volfrac");
        }

        // voxel_elem_mat_id is allocated in the voxel file read
        
    } // end method

}; // end fill elem states

// a function to verify that the required nodal state has initial conditions (i.e., fill values)
static bool check_fill_node_states(
        const std::vector<fill_node_state>& required_fill_node_states,
        const std::vector<fill_node_state>& fill_node_states){
    
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


// a function to verify that the required matpt state has initial conditions (i.e., fill values)
// fills occure on gauss points, where matpts live on top of those gauss points
static bool check_fill_mat_states(
        const std::vector<fill_gauss_state>& required_fill_matpt_states,
        const std::vector<fill_gauss_state>& fill_gauss_states){
    
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




// Possible node states, used to initialize node_t
enum class node_state
{
    coords,
    velocity,
    mass,
    temp,
    heat_transfer,
    force,
    gradient_level_set
};


/////////////////////////////////////////////////////////////////////////////
///
/// \struct node_t
///
/// \brief Stores state information associated with a node
///
/////////////////////////////////////////////////////////////////////////////
struct node_t
{
    DistributedDCArray<double> coords, local_coords, ghost_coords;    ///< Nodal coordinates
    DistributedDCArray<double> coords_n0, local_coords_n0, ghost_coords_n0; ///< Nodal coordinates at tn=0 of time integration
    DistributedDCArray<double> vel, local_vel, ghost_vel;       ///< Nodal velocity
    DistributedDCArray<double> vel_n0, local_vel_n0, ghost_vel_n0;    ///< Nodal velocity at tn=0 of time integration
    DistributedDCArray<double> mass, local_mass, ghost_mass;      ///< Nodal mass
    DistributedDCArray<double> force, local_force, ghost_force;     ///< Nodal force
    DistributedDCArray<double> temp, local_temp, ghost_temp;      ///< Nodal temperature
    DistributedDCArray<double> temp_n0, local_temp_n0, ghost_temp_n0;   ///< Nodal temperature at tn=0 of time integration
    DistributedDCArray<double> q_transfer, local_q_transfer, ghost_q_transfer; ///< Nodal heat flux
    DistributedDCArray<double> gradient_level_set, local_gradient_level_set, ghost_gradient_level_set;   ///< Nodal gradient of the level set function

    // initialization method (num_nodes, num_dims, state to allocate)
    void initialize(size_t num_nodes, size_t num_dims, std::vector<node_state> node_states)
    {
        for (auto field : node_states){
            switch(field){
                case node_state::coords:
                    if (coords.size() == 0) this->coords = DistributedDCArray<double>(num_nodes, num_dims, "node_coordinates");
                    if (coords_n0.size() == 0) this->coords_n0 = DistributedDCArray<double>(num_nodes, num_dims, "node_coordinates_n0");
                    break;
                case node_state::velocity:
                    if (vel.size() == 0) this->vel = DistributedDCArray<double>(num_nodes, num_dims, "node_velocity");
                    if (vel_n0.size() == 0) this->vel_n0 = DistributedDCArray<double>(num_nodes, num_dims, "node_velocity_n0");
                    break;
                case node_state::force:
                    if (force.size() == 0) this->force = DistributedDCArray<double>(num_nodes, num_dims, "node_force");
                    break;
                case node_state::mass:
                    if (mass.size() == 0) this->mass = DistributedDCArray<double>(num_nodes, "node_mass");
                    break;
                case node_state::temp:
                    if (temp.size() == 0) this->temp = DistributedDCArray<double>(num_nodes, "node_temp");
                    if (temp_n0.size() == 0) this->temp_n0 = DistributedDCArray<double>(num_nodes, "node_temp_n0");
                    break;
                case node_state::heat_transfer:
                    if (q_transfer.size() == 0) this->q_transfer = DistributedDCArray<double>(num_nodes, "node_q_transfer");
                    break;
                case node_state::gradient_level_set:
                    if (gradient_level_set.size() == 0) this->gradient_level_set = DistributedDCArray<double>(num_nodes, num_dims, "node_grad_levelset");
                    break;
                default:
                    std::cout<<"Desired node state not understood in node_t initialize"<<std::endl;
                    throw std::runtime_error("**** Error in State Field Name ****");
            }
        }
    }; // end method

    /*initialize overload with a partitioned map; also has the option of rebuilding the same vector with a larger map to encompass its previous contents
      last option is used for num_all = num_local + num_ghost contiguous storage */
    void initialize(DistributedMap partitioned_map, size_t num_dims, std::vector<node_state> node_states, DistributedMap subview_map = DistributedMap(), size_t offset=0)
    {   
        if(subview_map.size() == 0){
            for (auto field : node_states){
                switch(field){
                    case node_state::coords:
                        if (coords.size() == 0) this->coords = DistributedDCArray<double>(partitioned_map, num_dims, "node_coordinates");
                        if (coords_n0.size() == 0) this->coords_n0 = DistributedDCArray<double>(partitioned_map, num_dims, "node_coordinates_n0");
                        break;
                    case node_state::velocity:
                        if (vel.size() == 0) this->vel = DistributedDCArray<double>(partitioned_map, num_dims, "node_velocity");
                        if (vel_n0.size() == 0) this->vel_n0 = DistributedDCArray<double>(partitioned_map, num_dims, "node_velocity_n0");
                        break;
                    case node_state::force:
                        if (force.size() == 0) this->force = DistributedDCArray<double>(partitioned_map, num_dims, "node_force");
                        break;
                    case node_state::mass:
                        if (mass.size() == 0) this->mass = DistributedDCArray<double>(partitioned_map, "node_mass");
                        break;
                    case node_state::temp:
                        if (temp.size() == 0) this->temp = DistributedDCArray<double>(partitioned_map, "node_temp");
                        if (temp_n0.size() == 0) this->temp_n0 = DistributedDCArray<double>(partitioned_map, "node_temp_n0");
                        break;
                    case node_state::heat_transfer:
                        if (q_transfer.size() == 0) this->q_transfer = DistributedDCArray<double>(partitioned_map, "node_q_transfer");
                        break;
                    case node_state::gradient_level_set:
                        if (gradient_level_set.size() == 0) this->gradient_level_set = DistributedDCArray<double>(partitioned_map, num_dims, "node_grad_levelset");
                        break;
                    default:
                        std::cout<<"Desired node state not understood in node_t initialize"<<std::endl;
                        throw std::runtime_error("**** Error in State Field Name ****");
                }
            }
        }
        else{
            //first assign already partitioned local array to the local variable (if allocated) since up until now there was no local vs all distinction
            //then create array storing all = local + ghost array using the corresponding local array as a subview to avoid duplicate storage
            for (auto field : node_states){
                switch(field){
                    case node_state::coords:
                        //store local data with existing managed view made by mesh read for now
                        if(this->local_coords.size()==0&&this->coords.size()!=0){
                            this->local_coords = this->coords;
                        }
                        if(this->local_coords_n0.size()==0&&this->coords_n0.size()!=0){
                            this->local_coords_n0 = this->coords_n0;
                        }
                        //storage for nlocal+nghost
                        this->coords = DistributedDCArray<double>(partitioned_map, num_dims, "node_coordinates");
                        this->coords_n0 = DistributedDCArray<double>(partitioned_map, num_dims, "node_coordinates_n0");
                        //assign local data to new storage if local data was allocated
                        if(this->local_coords.size()!=0){
                            super_vector_initialization(this->coords, this->local_coords, subview_map.size());
                        }
                        if(this->local_coords_n0.size()!=0){
                            super_vector_initialization(this->coords_n0, this->local_coords_n0, subview_map.size());
                        }
                        //replace local data storage with subview of nlocal+nghost; previous managed view should self-destruct here
                        this->local_coords = DistributedDCArray<double>(this->coords,subview_map);
                        this->local_coords_n0 = DistributedDCArray<double>(this->coords_n0,subview_map);
                        break;
                    case node_state::velocity:
                        //store local data with existing managed view made by mesh read for now
                        if(this->local_vel.size()==0&&this->vel.size()!=0){
                            this->local_vel = this->vel;
                        }
                        if(this->local_vel_n0.size()==0&&this->vel_n0.size()!=0){
                            this->local_vel_n0 = this->vel_n0;
                        }
                        //storage for nlocal+nghost
                        this->vel = DistributedDCArray<double>(partitioned_map, num_dims, "node_velocity");
                        this->vel_n0 = DistributedDCArray<double>(partitioned_map, num_dims, "node_velocity_n0");
                        //assign local data to new storage
                        if(this->local_vel.size()!=0){
                            super_vector_initialization(this->vel, this->local_vel, subview_map.size());
                        }
                        if(this->local_vel_n0.size()!=0){
                            super_vector_initialization(this->vel_n0, this->local_vel_n0, subview_map.size());
                        }
                        //replace local data storage with subview of nlocal+nghost; previous managed view should self-destruct here
                        this->local_vel = DistributedDCArray<double>(this->vel,subview_map);
                        this->local_vel_n0 = DistributedDCArray<double>(this->vel_n0,subview_map);
                        break;
                    case node_state::force:
                        //store local data with existing managed view made by mesh read for now
                        if(this->local_force.size()==0&&this->force.size()!=0){
                            this->local_force = this->force;
                        }
                        //storage for nlocal+nghost
                        this->force = DistributedDCArray<double>(partitioned_map, num_dims, "node_force");
                        //assign local data to new storage
                        if(this->local_force.size()!=0){
                            super_vector_initialization(this->force, this->local_force, subview_map.size());
                        }
                        //replace local data storage with subview of nlocal+nghost; previous managed view should self-destruct here
                        this->local_force = DistributedDCArray<double>(this->force,subview_map);
                        break;
                    case node_state::mass:
                        //store local data with existing managed view made by mesh read for now
                        if(this->local_mass.size()==0&&this->mass.size()!=0){
                            this->local_mass = this->mass;
                        }
                        //storage for nlocal+nghost
                        this->mass = DistributedDCArray<double>(partitioned_map, "node_mass");
                        //assign local data to new storage
                        if(this->local_mass.size()!=0){
                            super_vector_initialization(this->mass, this->local_mass, subview_map.size());
                        }
                        //replace local data storage with subview of nlocal+nghost; previous managed view should self-destruct here
                        this->local_mass = DistributedDCArray<double>(this->mass,subview_map);
                        break;
                    case node_state::temp:
                        //store local data with existing managed view made by mesh read for now
                        if(this->local_temp.size()==0&&this->temp.size()!=0){
                            this->local_temp = this->temp;
                        }
                        if(this->local_temp_n0.size()==0&&this->temp_n0.size()!=0){
                            this->local_temp_n0 = this->temp_n0;
                        }
                        //storage for nlocal+nghost
                        this->temp = DistributedDCArray<double>(partitioned_map, "node_temp");
                        this->temp_n0 = DistributedDCArray<double>(partitioned_map, "node_temp_n0");
                        //assign local data to new storage
                        if(this->local_temp.size()!=0){
                            super_vector_initialization(this->temp, this->local_temp, subview_map.size());
                        }
                        if(this->local_temp_n0.size()!=0){
                            super_vector_initialization(this->temp_n0, this->local_temp_n0, subview_map.size());
                        }
                        //replace local data storage with subview of nlocal+nghost; previous managed view should self-destruct here
                        this->local_temp = DistributedDCArray<double>(this->temp,subview_map);
                        this->local_temp_n0 = DistributedDCArray<double>(this->temp_n0,subview_map);
                        break;
                    case node_state::heat_transfer:
                        //store local data with existing managed view made by mesh read for now
                        if(this->local_q_transfer.size()==0&&this->q_transfer.size()!=0){
                            this->local_q_transfer = this->q_transfer;
                        }
                        //storage for nlocal+nghost
                        this->q_transfer = DistributedDCArray<double>(partitioned_map, "node_q_transfer");
                        
                        //assign local data to new storage
                        if(this->local_q_transfer.size()!=0){
                            super_vector_initialization(this->q_transfer, this->local_q_transfer, subview_map.size());
                        }
                        //replace local data storage with subview of nlocal+nghost; previous managed view should self-destruct here
                        this->local_q_transfer = DistributedDCArray<double>(this->q_transfer,subview_map);
                        break;
                    case node_state::gradient_level_set:
                        //store local data with existing managed view made by mesh read for now
                        if(this->local_gradient_level_set.size()==0&&this->gradient_level_set.size()!=0){
                            this->local_gradient_level_set = this->gradient_level_set;
                        }
                        //storage for nlocal+nghost
                        this->gradient_level_set = DistributedDCArray<double>(partitioned_map, num_dims, "node_grad_levelset");
                        //assign local data to new storage
                        if(this->gradient_level_set.size()!=0){
                            super_vector_initialization(this->gradient_level_set, this->local_gradient_level_set, subview_map.size());
                        }
                        //replace local data storage with subview of nlocal+nghost; previous managed view should self-destruct here
                        this->local_gradient_level_set = DistributedDCArray<double>(this->gradient_level_set,subview_map);
                        break;
                    default:
                        std::cout<<"Desired node state not understood in node_t initialize"<<std::endl;
                        throw std::runtime_error("**** Error in State Field Name ****");
                }
            }
        }
    }; // end method

}; // end node_t


// Possible gauss point states, used to initialize GaussPoint_t
enum class gauss_pt_state
{
    volume,
    divergence_velocity,
    gradient_velocity,
    level_set
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct GaussPoint_t
///
/// \brief Stores state information associated with the Gauss point
///
/////////////////////////////////////////////////////////////////////////////
struct GaussPoint_t
{

    DCArrayKokkos<double> vol;  ///< GaussPoint volume
    DCArrayKokkos<double> div;  ///< GaussPoint divergence of velocity
    DCArrayKokkos<double> vel_grad;  ///< GaussPoint velocity gradient tensor

    DCArrayKokkos<double> level_set;  ///< GaussPoint level set field
    DCArrayKokkos<double> level_set_n0;  ///< GaussPoint level set field

    // initialization method (num_cells, num_dims)
    void initialize(size_t num_gauss_pnts, size_t num_dims, std::vector<gauss_pt_state> gauss_pt_states)
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
                    if (vel_grad.size() == 0) this->vel_grad = DCArrayKokkos<double>(num_gauss_pnts, num_dims, num_dims, "gauss_point_vel_grad");
                    break;
                case gauss_pt_state::level_set:
                    if (level_set.size() == 0) this->level_set = DCArrayKokkos<double>(num_gauss_pnts, "gauss_point_level_set");
                    if (level_set_n0.size() == 0) this->level_set_n0 = DCArrayKokkos<double>(num_gauss_pnts, "gauss_point_level_set_n0");
                    break;
                default:
                    std::cout<<"Desired gauss point state not understood in GaussPoint_t initialize"<<std::endl;
                    throw std::runtime_error("**** Error in State Field Name ****");
            }
        }
    }; // end method
};  // end GaussPoint_t

/////////////////////////////////////////////////////////////////////////////
///
/// \struct MeshtoMaterialMap_t
///
/// \brief Stores state information associated with maps from material to mesh maps
///
/////////////////////////////////////////////////////////////////////////////
struct MeshtoMaterialMap_t
{
    DCArrayKokkos<size_t> num_mats_in_elem; ///< returns the exact number of materials in elem
    DCArrayKokkos<size_t> mat_id;           ///< returns the mat_id 
    DCArrayKokkos<size_t> mat_storage_lid;  ///< returns the material storage local index

    // initialization method for FE-SGH and MPM methods (max number of elems needed)
    void initialize(size_t num_elem_max, size_t num_mats_per_elem_max)
    {
        if (num_mats_in_elem.size() == 0){
            this->num_mats_in_elem = DCArrayKokkos<size_t>(num_elem_max, "num_mats_in_elem");
            this->num_mats_in_elem.set_values(0); // initialize all elems to storing 0 materials
            this->num_mats_in_elem.update_host(); // copy from GPU to CPU
        }
        if (mat_id.size() == 0){
            this->mat_id = DCArrayKokkos<size_t>(num_elem_max, num_mats_per_elem_max, "mat_id_in_elem");
        }
        if (mat_storage_lid.size() == 0){
            this->mat_storage_lid = DCArrayKokkos<size_t>(num_elem_max, num_mats_per_elem_max, "mat_storage_lid_in_elem");
        }
        
    }; // end method
}; // end MeshtoMaterialMaps_t


/////////////////////////////////////////////////////////////////////////////
///
/// \struct MaterialtoMeshMap_t
///
/// \brief Stores state information associated with maps from material to mesh maps
///
/////////////////////////////////////////////////////////////////////////////
struct MaterialToMeshMap_t
{
    DCArrayKokkos <size_t> num_material_elems;        ///< returns the exact number of matpts
    DCArrayKokkos <size_t> num_material_elems_buffer; ///< returns the number of matpts plus buffer

    DRaggedRightArrayKokkos<size_t> elem;             ///< returns the elem for this material

    // initialization method for FE-SGH and MPM methods (max number of elems needed)
    void initialize()
    {
        if (elem.size() == 0){ 
            this->elem = DRaggedRightArrayKokkos<size_t>(this->num_material_elems_buffer, "material_space_to_elem");
        }

    }; // end method

    void initialize_num_mats(size_t num_mats)
    {
        // Note: num_material_elems is allocated in problem setup
        if (num_material_elems.size() == 0){
            this->num_material_elems = DCArrayKokkos <size_t> (num_mats, "num_material_elems"); 
        }

        // Note: num_material_elems_buffer is allocated in problem setup, the values are set in region_fill.cpp routine
        if (num_material_elems_buffer.size() == 0){
            this->num_material_elems_buffer = DCArrayKokkos <size_t> (num_mats, "num_material_elems_with_buffer"); 
        }

    }; // end method

}; // end MaterialtoMeshMaps_t



// Possible material point states, used to initialize MaterialPoint_t
enum class material_pt_state
{
    density,
    pressure,
    stress,
    sound_speed,
    specific_internal_energy,
    mass,
    volume_fraction,
    heat_flux,
    eroded_flag,
    elastic_modulii,
    shear_modulii,
    poisson_ratios,
    thermal_conductivity,
    specific_heat
};
/////////////////////////////////////////////////////////////////////////////
///
/// \struct MaterialPoint_t
///
/// \brief Stores state information associated with the material
///
/////////////////////////////////////////////////////////////////////////////
struct MaterialPoint_t
{
    DCArrayKokkos <size_t> num_material_points;    ///< the actual number of material points, omitting the buffer
    DCArrayKokkos <size_t> num_material_points_buffer; ///< number of material points plus a buffer

    DRaggedRightArrayKokkos<double> den;    ///< MaterialPoint density
    DRaggedRightArrayKokkos<double> pres;   ///< MaterialPoint pressure

    DRaggedRightArrayKokkos<double> stress;    ///< MaterialPoint stress
    DRaggedRightArrayKokkos<double> stress_n0; ///< MaterialPoint stress at t=n0 of time integration

    DRaggedRightArrayKokkos<double> sie;    ///< coefficients for the sie in strong form, only used in some methods e.g., FE-SGH and MPM
    DRaggedRightArrayKokkos<double> sie_n0; ///< coefficients for the sie in strong form at t=n0 of time integration

    DRaggedRightArrayKokkos<double> sspd;   ///< MaterialPoint sound speed
    DRaggedRightArrayKokkos<double> mass;   ///< MaterialPoint mass

    DRaggedRightArrayKokkos<double> q_flux; ///< Divergence of heat flux
    DRaggedRightArrayKokkos<double> conductivity;  ///< Thermal conductivity
    DRaggedRightArrayKokkos<double> specific_heat; ///< Specific Heat

    DRaggedRightArrayKokkos<double> elastic_modulii;  ///<  MaterialPoint elastic modulii Exx, Eyy, Ezz
    DRaggedRightArrayKokkos<double> shear_modulii;    ///<  MaterialPoint shear modulii Gxy, Gxz, Gyz
    DRaggedRightArrayKokkos<double> poisson_ratios;   ///<  MaterialPoint poisson ratios nu_xy, nu_xz, nu_yz
    

    // Material Models are stored on Material points
    DRaggedRightArrayKokkos<double> eos_state_vars;        ///< Array of state variables for the EOS, accessed as (mat_id, elem_lid, var_lid)
    DRaggedRightArrayKokkos<double> strength_state_vars;   ///< Array of state variables for the strength

    DRaggedRightArrayKokkos<double> temp_grad;     ///< Temperature gradient
    DRaggedRightArrayKokkos<double> volfrac;       ///< MaterialPoint volume fraction
    DRaggedRightArrayKokkos<double> delta_volfrac; ///< change in MaterialPoint volume fraction
    DRaggedRightArrayKokkos<double> geo_volfrac;   ///< change in MaterialPoint geometric (part) volume fraction (interface reconstruction)
    DRaggedRightArrayKokkos<double> delta_geo_volfrac; ///< change in MaterialPoint geometric (part) volume fraction (interface reconstruction)
    DRaggedRightArrayKokkos<bool> eroded;              ///< MaterialPoint eroded or not flag


    void initialize_num_mats(size_t num_mats)
    {
        // Note: num_material_points is allocated in problem setup
        if (num_material_points.size() == 0){
            this->num_material_points = DCArrayKokkos <size_t> (num_mats, "num_material_points"); 
        }

        // Note: num_material_points_buffer is allocated in problem setup, the values are set in region_fill.cpp routine
        if (num_material_points_buffer.size() == 0){
            this->num_material_points_buffer = DCArrayKokkos <size_t> (num_mats, "num_material_points_with_buffer"); 
        }

    }; // end method

    // initialization method (num_dims)
    void initialize(size_t num_dims, std::vector<material_pt_state> material_pt_states)
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
                    if (stress.size() == 0) this->stress = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, num_dims, num_dims, "material_point_stress");  
                    if (stress_n0.size() == 0) this->stress_n0 = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, num_dims, num_dims, "material_point_stress_n0"); 
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
                    if (volfrac.size() == 0) this->volfrac = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_volfrac");
                    if (geo_volfrac.size() == 0) this->geo_volfrac = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_geo_volfrac");
                    // changes in volume fraction
                    if (delta_volfrac.size() == 0) this->delta_volfrac = DRaggedRightArrayKokkos<double>(this->num_material_points_buffer, "material_point_volfrac_delta");
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
    }; // end method

}; // end MaterialPoint



// Possible material zone states, used to initialize MaterialZone_t
enum class material_zone_state
{
    specific_internal_energy
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct MaterialZone_t
///
/// \brief Stores state information associated with zone index space
///
/////////////////////////////////////////////////////////////////////////////
struct MaterialZone_t
{
    DCArrayKokkos <size_t> num_material_zones; ///< the actual number of material zones, omitting the buffer
    DCArrayKokkos <size_t> num_material_zones_buffer; ///< the number of material zones plus a buffer

    DRaggedRightArrayKokkos<double> sie;     ///< coefficients for the sie polynomial field
    DRaggedRightArrayKokkos<double> sie_n0;  ///< coefficients for the sie polynomial field at t=n0 of time integration


    void initialize_num_mats(size_t num_mats)
    {
        // Note: num_material_zones is allocated in problem setup
        if (num_material_zones.size() == 0){
            this->num_material_zones = DCArrayKokkos <size_t> (num_mats, "num_material_zones"); 
        }

        // Note: num_material_zones_buffer is allocated in problem setup, the values are set in region_fill.cpp routine
        if (num_material_zones_buffer.size() == 0){
            this->num_material_zones_buffer = DCArrayKokkos <size_t> (num_mats, "num_material_zones_with_buffer"); 
        }

    }; // end method

    // initialization method for arbitrary-order FE (num_zones)
    void initialize_Pn(std::vector<material_zone_state> material_zone_states)
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
        
    }; // end method

}; // end MaterialZone_t



// Possible material corner states, used to initialize MaterialCorner_t
enum class material_corner_state
{
    force, 
    heat_transfer
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct MaterialCorner_t
///
/// \brief Stores state information associated with a material in an element corner
///
/////////////////////////////////////////////////////////////////////////////
struct MaterialCorner_t
{
    DCArrayKokkos<size_t>num_material_corners;   ///< the actual number of material corners, omitting the buffer
    DCArrayKokkos<size_t>num_material_corners_buffer;   ///< the number of material corners plus a buffer

    DRaggedRightArrayKokkos<double> force;   ///< Corner force for the material

    DRaggedRightArrayKokkos<double> q_transfer;  ///< Corner heat tranfer per material

    
    void initialize_num_mats(size_t num_mats)
    {
        // Note: num_material_corners is allocated in problem setup
        if (num_material_corners.size() == 0){
            this->num_material_corners = DCArrayKokkos <size_t> (num_mats, "num_material_corners"); 
        }

        // Note: num_material_corners_buffer is allocated in problem setup, the values are set in region_fill.cpp routine
        if (num_material_corners_buffer.size() == 0){
            this->num_material_corners_buffer = DCArrayKokkos <size_t> (num_mats, "num_material_corners_with_buffer"); 
        }

    }; // end method

    // initialization method (num_corners, num_dims)
    void initialize(size_t num_dims, std::vector<material_corner_state> material_corner_states)
    {
        for (auto field : material_corner_states){
            switch(field){
                case material_corner_state::force:
                    if (force.size() == 0) this->force = DRaggedRightArrayKokkos<double>(this->num_material_corners_buffer, num_dims, "material_corner_force");
                    break;
                case material_corner_state::heat_transfer:
                    if (q_transfer.size() == 0) this->q_transfer =DRaggedRightArrayKokkos<double>(this->num_material_corners_buffer, "material_corner_q_transfer"); 
                    break;
                default:
                    std::cout<<"Desired material corner state not understood in MaterialCorner_t initialize"<<std::endl;
                    throw std::runtime_error("**** Error in State Field Name ****");
            }
        }
    }; // end method
}; // end material corner



// Possible material corner states, used to initialize MaterialCorner_t
enum class corner_state
{
    force, 
    mass,
    heat_transfer,
    normal,
    volume
};
/////////////////////////////////////////////////////////////////////////////
///
/// \struct corner_t
///
/// \brief Stores state information associated with a corner
///
/////////////////////////////////////////////////////////////////////////////
struct corner_t
{
    DCArrayKokkos<double> force;  ///< Corner force
    DCArrayKokkos<double> mass;   ///< Corner mass
    DCArrayKokkos<double> q_transfer;  ///< Corner heat transfer
    DCArrayKokkos<double> normal; ///< Corner normal
    DCArrayKokkos<double> volume; ///< Corner volume

    // initialization method (num_corners, num_dims)
    void initialize(size_t num_corners, size_t num_dims, std::vector<corner_state> corner_states)
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

        
        
        
    }; // end method
}; // end corner_t

/////////////////////////////////////////////////////////////////////////////
///
/// \struct map for getting the corners in the material index
///
/// \brief Stores state information associated with material corner index space
///
/////////////////////////////////////////////////////////////////////////////
struct corners_in_mat_t
{
    private:
        size_t num_corners_in_elem_;
    public:
        corners_in_mat_t() {
        };

        corners_in_mat_t(const size_t num_corners_in_elem_inp) {
            this->num_corners_in_elem_ = num_corners_in_elem_inp;
        };

        // return global corner index for given local corner index in a material storage
        size_t  host(const size_t mat_storage_lid, const size_t corner_lid) const
        {
            return mat_storage_lid * num_corners_in_elem_ + corner_lid;
        };

        // Return the global corner ID given a material storage gloabl ID and a local corner ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t mat_storage_lid, const size_t corner_lid) const
        {
            return mat_storage_lid * num_corners_in_elem_ + corner_lid;
        };
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct maps for high-order FE methods
///
/// \brief Stores state information associated with other mesh index spaces
///
/////////////////////////////////////////////////////////////////////////////

struct zones_in_mat_t
{
    private:
        size_t num_zones_in_elem_;
    public:
        zones_in_mat_t() {
        };

        zones_in_mat_t(const size_t num_zones_in_elem_inp) {
            this->num_zones_in_elem_ = num_zones_in_elem_inp;
        };

        // return global zone index for given local zone index in a material storage
        size_t  host(const size_t mat_storage_lid, const size_t zone_lid) const
        {
            return mat_storage_lid * num_zones_in_elem_ + zone_lid;
        };

        // Return the global zone ID given a material storage gloabl ID and a local zone ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t mat_storage_lid, const size_t zone_lid) const
        {
            return mat_storage_lid * num_zones_in_elem_ + zone_lid;
        };
};

// if material points are defined strictly internal to the element.
struct legendre_in_mat_t
{
    private:
        size_t num_leg_gauss_in_elem_;
    public:
        legendre_in_mat_t() {
        };

        legendre_in_mat_t(const size_t num_leg_gauss_in_elem_inp) {
                this->num_leg_gauss_in_elem_ = num_leg_gauss_in_elem_inp;
        };

        // return global gauss index for given local gauss index in a material storage
        size_t  host(const size_t mat_storage_lid, const size_t leg_gauss_lid) const
        {
            return mat_storage_lid * num_leg_gauss_in_elem_ + leg_gauss_lid;
        };

        // Return the global gauss ID given a material storage gloabl ID and a local gauss ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t mat_storage_lid, const size_t leg_gauss_lid) const
        {
            return mat_storage_lid * num_leg_gauss_in_elem_ + leg_gauss_lid;
        };
};

/// if material points are defined at element interfaces, e.g., for nodal DG
struct lobatto_in_mat_t
{
    private:
        size_t num_lob_gauss_in_elem_;
    public:
        lobatto_in_mat_t() {
        };

        lobatto_in_mat_t(const size_t num_lob_gauss_in_elem_inp) {
                this->num_lob_gauss_in_elem_ = num_lob_gauss_in_elem_inp;
        };

        // return global gauss index for given local gauss index in a material storage
        size_t  host(const size_t mat_storage_lid, const size_t lob_gauss_lid) const
        {
            return mat_storage_lid * num_lob_gauss_in_elem_ + lob_gauss_lid;
        };

        // Return the global gauss ID given a material storage ID and a local gauss ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t mat_storage_lid, const size_t lob_gauss_lid) const
        {
            return mat_storage_lid * num_lob_gauss_in_elem_ + lob_gauss_lid;
        };
};

// the local id for material points in elem
struct points_in_mat_t
{
    private:
        size_t num_points_in_elem_;
    public:
        points_in_mat_t() {
        };

        points_in_mat_t(const size_t num_points_in_elem_inp) {
                this->num_points_in_elem_ = num_points_in_elem_inp;
        };

        // return global gauss index for given local gauss index in a material storage
        size_t  host(const size_t mat_storage_lid, const size_t points_lid) const
        {
            return mat_storage_lid * num_points_in_elem_ + points_lid;
        };

        // Return the global gauss ID given a material storage gloabl ID and a local gauss ID
        KOKKOS_INLINE_FUNCTION
        size_t operator()(const size_t mat_storage_lid, const size_t points_lid) const
        {
            return mat_storage_lid * num_points_in_elem_ + points_lid;
        };
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct state_t
///
/// \brief Stores all state
///
/////////////////////////////////////////////////////////////////////////////
struct State_t
{
    // ---------------------------------------------------------------------
    //    state data on mesh declarations
    // ---------------------------------------------------------------------
    node_t node;              ///< access as node.coords(node_gid,dim)
    GaussPoint_t GaussPoints; ///< access as GaussPoints.vol(gauss_pt_gid)
    corner_t corner;          ///< access as corner.force(corner_gid,dim)

    // ---------------------------------------------------------------------
    //    material to mesh maps and mesh to material maps
    // ---------------------------------------------------------------------
    MaterialToMeshMap_t MaterialToMeshMaps; ///< access as MaterialToMeshMaps.elem(mat_id, mat_storage_lid)
    MeshtoMaterialMap_t MeshtoMaterialMaps; ///< acces as MeshtoMaterialMaps.mat_id(elem, mat_lid)

    // ---------------------------------------------------------------------
    //    material to material maps
    // ---------------------------------------------------------------------
    corners_in_mat_t corners_in_mat_elem; ///< access the corner mat lid using (mat_elem_lid, corn_lid)
    points_in_mat_t points_in_mat_elem;   ///< for accessing e.g., gauss points mat lid with arbitrary-order FE (mat_elem_lid, gauss_lid)
    zones_in_mat_t zones_in_mat_elem;     ///< for accessing sub-zones mat lid with arbitrary-order FE

    // ---------------------------------------------------------------------
    //    material state, compressed, and sequentially accessed
    // ---------------------------------------------------------------------
    MaterialPoint_t MaterialPoints;   ///< access as MaterialPoints.var(mat_id, mat_pt)
    MaterialCorner_t MaterialCorners; ///< access as MaterialCorners.var(mat_id, mat_corner)
    MaterialZone_t MaterialZones;     ///< access as MaterialZones.var(mat_id, mat_zone), only used with arbitrary-order FE
}; // end state_t

template <typename T>
void super_vector_initialization(DistributedDCArray<T> super_vector, DistributedDCArray<T> sub_vector, size_t extent){
    FOR_ALL(i, 0, extent, {
        for (size_t dim = 0; dim < super_vector.component_length(); dim++) {
            super_vector(i, dim) = sub_vector(i, dim);
        }
    }); // end parallel for corners
}

template <typename T>
void super_vector_initialization(DistributedDFArray<T> super_vector, DistributedDFArray<T> sub_vector, size_t extent){
    FOR_ALL(i, 0, extent, {
        for (size_t dim = 0; dim < super_vector.component_length(); dim++) {
            super_vector(i, dim) = sub_vector(i, dim);
        }
    }); // end parallel for corners
}



#endif