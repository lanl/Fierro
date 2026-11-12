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

#include "table.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>

Table_t::Table_t() : num_rows(0), num_columns(0)
{
    name[0] = '\0';
}

Table_t::Table_t(size_t num_rows_, size_t num_columns_, const char* name_) : num_rows(num_rows_), num_columns(num_columns_)
{
    // Initialize data with dimensions and name
    data = DCArrayKokkos<real_t>(num_rows_, num_columns_, name_);
    
    // Hardcoded NaN. Here, we use a bit pattern for portability in device code. This is to verify that data has been initialized correctly. 
    union { uint64_t i; real_t d; } nan_union;
    nan_union.i = 0x7ff8000000000000ULL;
    data.set_values(nan_union.d);

    // Copy name safely
    size_t i = 0;
    for (; i < max_name_length - 1 && name_[i] != '\0'; ++i) {
        name[i] = name_[i];
    }
    name[i] = '\0';
}

void Table_t::set_value(size_t row, size_t col, real_t val)
{
    assert(row < num_rows && col < num_columns && "Row and column indices are out of bounds in Table_t::set_value");
    data.host(row, col) = val;
}

void Table_t::update_device()
{
    data.update_device();
}

void Table_t::update_host()
{
    data.update_host();
}

void Table_t::verify_data(size_t sorted_column_idx)
{
    assert(sorted_column_idx < num_columns && "Sorted column index is out of bounds in Table_t::verify_data");
    assert(data.host.dims(0) == num_rows && data.host.dims(1) == num_columns && "Data dimensions do not match the table dimensions in Table_t::verify_data");
    for (size_t row = 0; row + 1 < num_rows; row++) {
        assert(data.host(row, sorted_column_idx) <= data.host(row + 1, sorted_column_idx) &&
               "Data in sorted column is not ascending in Table_t::verify_data");
    }

    // Check for NaN in each element; throw assertion if found
    for (size_t row = 0; row < num_rows; ++row) {
        for (size_t col = 0; col < num_columns; ++col) {
            assert(!std::isnan(data.host(row, col)) &&
                   "NaN detected in Table_t::verify_data");
        }
    }
}

void Table_t::print_table() const
{
    for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < num_columns; ++j) {
            printf("%g ", data.host(i, j));
        }
        printf("\n");
    }
}
