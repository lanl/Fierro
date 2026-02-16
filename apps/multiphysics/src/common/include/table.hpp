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

#ifndef TABLE_H
#define TABLE_H

#include "matar.h"
#include <cstring>

// -----------------------------------------------------------------------------
//  Table_t: A structure to hold tabular data and perform linear interpolation
//  using MATAR data types for portability (CPU/GPU).
// ------------------------------------------------------------------------------
struct Table_t{
    
    // Store the table name as a fixed-size char array for host/device printing
    static constexpr int max_name_length = 64;
    char name[max_name_length];

    size_t num_rows = 0;
    size_t num_columns = 0;

    // Data stored as a 2D array (num_rows x num_columns) using MATAR's DCArrayKokkos.
    // - Rows correspond to data points.
    // - Columns correspond to variables (e.g., [time, x, y, z, power] or [temp, property]).
    DCArrayKokkos<double> data;


    // Default Constructor
    Table_t() : num_rows(0), num_columns(0) {
        name[0] = '\0';
    }

    // Constructor  
    Table_t(size_t num_rows_, size_t num_columns_, const char* name_) 
        : num_rows(num_rows_), num_columns(num_columns_)
    {
        // Initialize data with dimensions and name
        data = DCArrayKokkos<double>(num_rows_, num_columns_, name_);
        
        // Hardcoded NaN. Here, we use a bit pattern for portability in device code. This is to verify that data has been initialized correctly. 
        union { uint64_t i; double d; } nan_union;
        nan_union.i = 0x7ff8000000000000ULL;
        data.set_values(nan_union.d);

        // Copy name safely
        size_t i = 0;
        for (; i < max_name_length - 1 && name_[i] != '\0'; ++i) {
            name[i] = name_[i];
        }
        name[i] = '\0';
    }

    // Set a value in the table on the host
    void set_value(size_t row, size_t col, double val) {
        assert(row < num_rows && col < num_columns && "Row and column indices are out of bounds in Table_t::set_value");
        data.host(row, col) = val;
    }
    
    // Update the device view with data from the host
    void update_device() {
        data.update_device();
    }

    // Update the host view with data from the device
    void update_host() {
        data.update_host();
    }

    // Verify that the table is valid and the data in at least one column is sorted in ascending order.
    void verify_data(size_t sorted_column_idx) {
        assert(data.host.extent(0) == num_rows && data.host.extent(1) == num_columns && "Data dimensions do not match the table dimensions in Table_t::verify_data");
        for (size_t row = 0; row < num_rows - 1; row++) {
            assert(data.host(row, sorted_column_idx) <= data.host(row + 1, sorted_column_idx) && "Data in column " + std::to_string(sorted_column_idx) + " is not sorted in ascending order in Table_t::verify_data");
        }

        // Check for NaN in each element; throw assertion if found
        for (size_t row = 0; row < num_rows; ++row) {
            for (size_t col = 0; col < num_columns; ++col) {
                assert(!std::isnan(data.host(row, col)) &&
                    "NaN detected in Table_t::verify_data at row " + std::to_string(row) + ", column " + std::to_string(col));
            }
        }
    }

    // Linear interpolation function for host side access
    // Returns the interpolated value for the dependent variable in column `dependent_idx`
    // based on the independent variable value `x` in column `independent_idx`.
    // Assumes that the data in `independent_idx` is sorted in ascending order.
    KOKKOS_INLINE_FUNCTION
    double linear_interpolation_host(double independent_value, size_t dependent_idx, size_t independent_idx = 0) const {
        if (num_rows == 0) return 0.0;
        
        // Handle cases where independent_value is outside the range (clamping/extrapolation)
        // If independent_value is before the first point, return the first point's value
        if (independent_value <= data.host(0, independent_idx)) {
            return data.host(0, dependent_idx);
        }
        
        // If independent_value is after the last point, return the last point's value
        if (independent_value >= data.host(num_rows - 1, independent_idx)) {
            return data.host(num_rows - 1, dependent_idx);
        }

        // Find the segment [x0, x1] containing independent_value using linear search
        // (Note: For very large tables, binary search would be more efficient)
        for (size_t i = 0; i < num_rows - 1; ++i) {
            double x0 = data.host(i, independent_idx);
            double x1 = data.host(i + 1, independent_idx);

            double x = independent_value;
            
            if (x >= x0 && x <= x1) {
                // Avoid division by zero
                if (x1 == x0) return data.host(i, dependent_idx);

                double alpha = (x - x0) / (x1 - x0); // Interpolation factor
                double y0 = data.host(i, dependent_idx);
                double y1 = data.host(i + 1, dependent_idx);
                
                return (1.0 - alpha) * y0 + alpha * y1;
            }
        }
        
        // Should not be reached given the boundary checks above
        return 0.0;
    }

    // Linear interpolation function for device side access
    // Returns the interpolated value for the dependent variable in column `dependent_idx`
    // based on the independent variable value `x` in column `independent_idx`.
    // Assumes that the data in `independent_idx` is sorted in ascending order.
    KOKKOS_INLINE_FUNCTION
    double linear_interpolation(double independent_value, size_t dependent_idx, size_t independent_idx = 0) const {
        if (num_rows == 0) return 0.0;
        
        // Handle cases where independent_value is outside the range (clamping/extrapolation)
        // If independent_value is before the first point, return the first point's value
        if (independent_value <= data(0, independent_idx)) {
            return data(0, dependent_idx);
        }
        
        // If independent_value is after the last point, return the last point's value
        if (independent_value >= data(num_rows - 1, independent_idx)) {
            return data(num_rows - 1, dependent_idx);
        }

        // Find the segment [x0, x1] containing independent_value using linear search
        // (Note: For very large tables, binary search would be more efficient)
        for (size_t i = 0; i < num_rows - 1; ++i) {
            double x0 = data(i, independent_idx);
            double x1 = data(i + 1, independent_idx);

            double x = independent_value;
            
            if (x >= x0 && x <= x1) {
                // Avoid division by zero
                if (x1 == x0) return data(i, dependent_idx);

                double alpha = (x - x0) / (x1 - x0); // Interpolation factor
                double y0 = data(i, dependent_idx);
                double y1 = data(i + 1, dependent_idx);
                
                return (1.0 - alpha) * y0 + alpha * y1;
            }
        }
        
        // Should not be reached given the boundary checks above
        return 0.0;
    }

    // Print the table to the console
    void print_table() const {
        for (size_t i = 0; i < num_rows; ++i) {
            for (size_t j = 0; j < num_columns; ++j) {
                std::cout << data.host(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }

};



#endif // end Header Guard
