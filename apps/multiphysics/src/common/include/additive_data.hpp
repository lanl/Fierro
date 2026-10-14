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
#ifndef FIERRO_ADDITIVE_DATA_H
#define FIERRO_ADDITIVE_DATA_H

#include <stdio.h>
#include "matar.h"
#include "table.hpp"


/////////////////////////////////////////////////////////////////////////////
///
/// \class ToolPathInfo
///
/// \brief Stores and manages additive manufacturing tool paths with time-parameterized 3D points
///        and provides utilities for obtaining the current tool position.
///
/// ToolPathInfo uses MATAR data types (e.g., CArray) and Table_t to represent
/// a sequence of tool path points, each defined by time, x, y, z, and power.
/// Provides methods for setting data points, updating device views, and interpolating
/// tool position and power at arbitrary times along the path.
///
/////////////////////////////////////////////////////////////////////////////
class ToolPathInfo {
public:
    
    enum Fields
    {
        time = 0,   ///<  time
        x = 1,      ///<  x position
        y = 2,      ///<  y position
        z = 3,      ///<  z position
        power = 4,  ///<  power
    };
    Table_t tool_path_table;

    size_t num_columns = 5;

    ToolPathInfo() = default;

    // CSV rows are time,x,y,z,power. Equal times are permitted for instantaneous
    // path transitions, matching the existing table interpolation convention.
    static ToolPathInfo load_laser_path(const char* filename);

    // Default constructor that takes the number of data points
    ToolPathInfo(size_t npoints)
    {
        tool_path_table = Table_t(npoints, num_columns, "ToolPathInfo_table");
    }

    // Set a data point (time, x, y, z, power) at index i
    void set_data_point(size_t i, double time, double x, double y, double z, double power) {
        tool_path_table.set_value(i, Fields::time, time);
        tool_path_table.set_value(i, Fields::x, x);
        tool_path_table.set_value(i, Fields::y, y);
        tool_path_table.set_value(i, Fields::z, z);
        tool_path_table.set_value(i, Fields::power, power);
    }
    
    // Update the device views (copy to the GPU from the CPU)
    void update_device() {
        tool_path_table.update_device();
    }

    // Compute current position of tool at time t, assuming linear motion between path points.
    // Returns a std::array<double,3> {x, y, z}
    KOKKOS_INLINE_FUNCTION
    void get_position(const double& t, double& x, double& y, double& z) const {

        // get the x position at time t
        x = tool_path_table.linear_interpolation(t, Fields::x, Fields::time);
        // get the y position at time t
        y = tool_path_table.linear_interpolation(t, Fields::y, Fields::time);
        // get the z position at time t
        z = tool_path_table.linear_interpolation(t, Fields::z, Fields::time);
    } // end function


    // Earliest time at or after start when tool elevation reaches a cell center.
    // A false return means the path never reaches this height. Duplicate-time
    // upward jumps activate at that timestamp; subsequent descent cannot undo it.
    KOKKOS_INLINE_FUNCTION
    bool first_height_crossing(double height, double start, double& event_time) const {
        const auto& table = tool_path_table;
        if (table.num_rows == 0) return false;
        if (table.linear_interpolation(start, Fields::z) >= height) {
            event_time = start;
            return true;
        }
        for (size_t row = 0; row < table.num_rows; ++row) {
            const double time = table.data(row, Fields::time);
            const double elevation = table.data(row, Fields::z);
            if (time < start || elevation < height) continue;
            event_time = time;
            if (row > 0) {
                const double previous_time = table.data(row - 1, Fields::time);
                const double previous_elevation = table.data(row - 1, Fields::z);
                if (time > previous_time && previous_elevation < height)
                    event_time = previous_time + (time - previous_time) *
                        ((height - previous_elevation) / (elevation - previous_elevation));
            }
            event_time = Kokkos::max(start, event_time);
            return true;
        }
        return false;
    }

    // Compute current power of tool at time t, assuming linear interpolation between path points.
    // Returns the power at the time t
    KOKKOS_INLINE_FUNCTION
    double get_power(const double& t) const {
        return tool_path_table.linear_interpolation(t, Fields::power, Fields::time);
    } // end function
};

enum class LaserModelType { SPHERICAL, GOLDAK };

struct SphericalParams_t {
    double radius = 8.0;  // Legacy moving-source radius.
};

struct GoldakParams_t {
    double major_front = 0.0;
    double major_rear = 0.0;
    double minor = 0.0;
    double depth = 0.0;
    double absorptivity = 0.0;
};

// Laser configuration contains numeric parameters and MATAR tool-path views.
struct Laser_t {
    bool configured = false;
    ToolPathInfo tool_path_info;
    LaserModelType model = LaserModelType::SPHERICAL;
    SphericalParams_t spherical;
    GoldakParams_t goldak;
};

// Numeric source data sampled on the host once per Runge-Kutta stage. This
// structure contains no host pointers and can be captured by a device kernel.
struct LaserSource_t {
    LaserModelType model = LaserModelType::SPHERICAL;
    SphericalParams_t spherical;
    GoldakParams_t goldak;
    double x = 0.0, y = 0.0, z = 0.0;
    double direction_x = 1.0, direction_y = 0.0;
    double power = 0.0;
};

LaserSource_t sample_laser_source(const Laser_t& laser, double time);

KOKKOS_INLINE_FUNCTION
double laser_power_density(const LaserSource_t& source, double x, double y, double z)
{
    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double sqrt_three = 1.732050807568877293527446341505872367;
    constexpr double sqrt_pi = 1.772453850905516027298167483341145183;
    if (source.power == 0.0) return 0.0;
    const double dx = x - source.x;
    const double dy = y - source.y;
    const double dz = z - source.z;
    if (source.model == LaserModelType::SPHERICAL) {
        const double radius = source.spherical.radius;
        return dx*dx + dy*dy + dz*dz <= radius*radius
            ? source.power / ((4.0/3.0) * pi * radius*radius*radius) : 0.0;
    }

    // Goldak's 6*sqrt(3) normalization deposits eta*P in the lower half-space.
    // Without this surface restriction the same formula deposits twice eta*P.
    if (dz > 0.0) return 0.0;
    const double xi = dx * source.direction_x + dy * source.direction_y;
    const double eta = -dx * source.direction_y + dy * source.direction_x;
    const auto& p = source.goldak;
    const double a = xi >= 0.0 ? p.major_front : p.major_rear;
    const double fraction = 2.0 * a / (p.major_front + p.major_rear);
    const double exponent = -3.0 * (xi*xi/(a*a) + eta*eta/(p.minor*p.minor) +
                                     dz*dz/(p.depth*p.depth));
    return 6.0 * sqrt_three * fraction * p.absorptivity * source.power /
           (a * p.minor * p.depth * pi * sqrt_pi) * Kokkos::exp(exponent);
}

// ----------------------------------
// valid inputs for additive data
// ----------------------------------
static std::vector<std::string> str_additive_data_inps
{
    "time_initial",
    "time_final",
    "dt_min",
    "dt_max",
    "dt_start",
    "dt_cfl",
    "cycle_stop",
    "fuzz",
    "tiny",
    "small",
    "rk_num_stages"
};

// ----------------------------------
// required inputs for dynamic options
// ----------------------------------
static std::vector<std::string> additive_data_required_inps
{
    "time_final",
    "cycle_stop",
};

#endif // end Header Guard
