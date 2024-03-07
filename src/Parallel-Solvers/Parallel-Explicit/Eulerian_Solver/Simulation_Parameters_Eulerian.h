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

#ifndef SIMULATION_PARAMETERS_EULERIAN_H
#define SIMULATION_PARAMETERS_EULERIAN_H

#include "matar.h"
#include <cmath>
#include <stdexcept>
#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>

#include "mesh.h"
#include "state.h"
#include "yaml-serializable.h"
#include "Simulation_Parameters.h"
#include "Simulation_Parameters_Base.h"

using namespace mtr;

struct Simulation_Parameters_Eulerian : Simulation_Parameters
{
    Time_Variables time_variables;
    std::vector<MaterialFill> region_options;
    std::vector<Material> material_options;
    std::vector<Boundary> boundary_conditions;
    Graphics_Options graphics_options;

    bool gravity_flag   = false;
    bool report_runtime = true;

    size_t rk_num_stages = 2;
    int NB   = 6; // number of boundaries
    int NBSF = 4; // number of surface density force conditions
    int NBV  = 2; // number of surface sets used to specify a fixed displacement on nodes belonging to respective surfaces

    // Non-serialized fields
    int num_gauss_points = 2;
    size_t max_num_eos_global_vars;
    size_t max_num_strength_global_vars;
    size_t rk_num_bins;
    double time_value = 0.0;
    DCArrayKokkos<double> eos_global_vars;
    DCArrayKokkos<double> strength_global_vars;

    DCArrayKokkos<mat_fill_t> mat_fill;
    DCArrayKokkos<material_t> material;
    DCArrayKokkos<boundary_t> boundary;
    std::vector<double> gravity_vector { 9.81, 0., 0. };

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn init_material_variable_arrays
    ///
    /// \brief Initialize data for EOS and material strength models
    ///
    /// \param Number of global EOS variables
    /// \param Number of global strength variables
    ///
    /////////////////////////////////////////////////////////////////////////////
    void init_material_variable_arrays(size_t num_eos_global_vars, size_t num_strength_global_vars)
    {
        eos_global_vars = DCArrayKokkos<double>(material_options.size(), num_eos_global_vars);
        strength_global_vars = DCArrayKokkos<double>(material_options.size(), num_strength_global_vars);

        for (size_t i = 0; i < material_options.size(); i++) {
            auto mat = material_options[i];

            for (size_t j = 0; j < mat.eos_global_vars.size(); j++) {
                eos_global_vars.host(i, j) = mat.eos_global_vars[j];
            }

            for (size_t j = 0; j < mat.strength_global_vars.size(); j++) {
                strength_global_vars.host(i, j) = mat.strength_global_vars[j];
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn from_vector
    ///
    /// \brief Copies data from a vector on the host device
    ///
    /// \param Array data is coped into
    /// \param Vector on the host device
    ///
    /////////////////////////////////////////////////////////////////////////////
    template<typename T, typename K> void from_vector(DCArrayKokkos<T>& array, const std::vector<K>& vec)
    {
        array = DCArrayKokkos<T>(vec.size());
        for (size_t i = 0; i < vec.size(); i++) {
            array.host(i) = *(T*)&vec[i];
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn derive_kokkos_arrays
    ///
    /// \brief Updates arrays on the device
    ///
    /////////////////////////////////////////////////////////////////////////////
    void derive_kokkos_arrays()
    {
        max_num_eos_global_vars = 0;
        max_num_strength_global_vars = 0;
        for (auto mo : material_options) {
            max_num_eos_global_vars = std::max(max_num_eos_global_vars, mo.eos_global_vars.size());
            max_num_strength_global_vars = std::max(max_num_strength_global_vars, mo.strength_global_vars.size());
        }

        init_material_variable_arrays(max_num_eos_global_vars, max_num_strength_global_vars);

        from_vector(mat_fill, region_options);
        from_vector(material, material_options);
        from_vector(boundary, boundary_conditions);

        // Send to device.
        mat_fill.update_device();
        boundary.update_device();
        material.update_device();
        eos_global_vars.update_device();
        strength_global_vars.update_device();
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn derive
    ///
    /// \brief Calls derive_kokkos_arrays and sets the number or Runge Kutta integration
    ///
    /////////////////////////////////////////////////////////////////////////////
    void derive()
    {
        derive_kokkos_arrays();
        rk_num_bins = rk_num_stages;
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn validate
    ///
    /// \brief Validates the specified FEA module
    ///
    /////////////////////////////////////////////////////////////////////////////
    void validate()
    {
        validate_module_is_specified(FEA_MODULE_TYPE::Eulerian);
    }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Simulation_Parameters_Eulerian, Simulation_Parameters,
  time_variables, material_options, region_options,
  boundary_conditions, gravity_flag, report_runtime, rk_num_stages,
  NB, NBSF, NBV,
  graphics_options
    )
#endif // end HEADER_H
