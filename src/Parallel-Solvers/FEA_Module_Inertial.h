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

#ifndef FEA_MODULE_INERTIAL_H
#define FEA_MODULE_INERTIAL_H

#include "FEA_Module.h"

class Inertial_Parameters;

class FEA_Module_Inertial : public FEA_Module
{
public:
    FEA_Module_Inertial(Inertial_Parameters& params, Solver* Solver_Pointer, const int my_fea_module_index = 0);
    ~FEA_Module_Inertial();

    void comm_variables(Teuchos::RCP<const MV> zp);

    void compute_element_volumes();

    void compute_element_masses(const_host_vec_array design_variables, bool max_flag, bool use_initial_coords = false);
    
    void compute_element_masses_TO(const_host_vec_array design_densities, bool max_flag, bool use_initial_coords = false);
    
    void compute_element_masses_SO(const_host_vec_array design_coords, bool max_flag, bool use_initial_coords = false);

    void compute_element_moments(const_host_vec_array design_densities, bool max_flag, int moment_component, bool use_initial_coords = false);

    void compute_element_moments_of_inertia(const_host_vec_array design_densities, bool max_flag, int inertia_component, bool use_initial_coords = false);

    void compute_nodal_gradients(const_host_vec_array design_variables, host_vec_array gradients, bool use_initial_coords = false);

    void compute_TO_gradients(const_host_vec_array design_densities, host_vec_array gradients, bool use_initial_coords = false);

    void compute_shape_gradients(const_host_vec_array design_coords, host_vec_array gradients, bool use_initial_coords = false);

    void compute_moment_gradients(const_host_vec_array design_densities, host_vec_array gradients, int moment_component, bool use_initial_coords = false);

    void compute_moment_of_inertia_gradients(const_host_vec_array design_densities, host_vec_array gradients, int intertia_component, bool use_initial_coords = false);

    // forward declare
    Inertial_Parameters* module_params;

    // Global FEA data
    Teuchos::RCP<MV> mass_gradients_distributed;
    Teuchos::RCP<MV> center_of_mass_gradients_distributed;
    Teuchos::RCP<MV> Global_Element_Volumes;
    Teuchos::RCP<MV> Global_Element_Masses;
    Teuchos::RCP<MV> Global_Element_Moments_x;
    Teuchos::RCP<MV> Global_Element_Moments_y;
    Teuchos::RCP<MV> Global_Element_Moments_z;
    Teuchos::RCP<MV> Global_Element_Moments_of_Inertia_xx;
    Teuchos::RCP<MV> Global_Element_Moments_of_Inertia_yy;
    Teuchos::RCP<MV> Global_Element_Moments_of_Inertia_zz;
    Teuchos::RCP<MV> Global_Element_Moments_of_Inertia_xy;
    Teuchos::RCP<MV> Global_Element_Moments_of_Inertia_xz;
    Teuchos::RCP<MV> Global_Element_Moments_of_Inertia_yz;

    // inertial properties
    real_t mass, center_of_mass[3], moments_of_inertia[6];

    bool use_initial_density; // if density variable is from initial configuration then jacobian is not needed

    // runtime flags
    bool mass_init, com_init[3];

    // update counters (first attempt at reducing redundant calls through ROL for Moments of Inertia and Center of Mass)
    int mass_update, com_update[3];
    int mass_gradient_update, com_gradient_update[3];
};

#endif // end HEADER_H
