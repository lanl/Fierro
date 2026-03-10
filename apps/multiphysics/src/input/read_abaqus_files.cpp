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


#include "read_abaqus_files.hpp"

// ==============================================================================
//   Function Definitions
// ==============================================================================

void AbaqusReader::read_tabular_jmatpro_file(const std::string& file_path, const std::string& table_name, Table_t& table)
{
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open Abaqus JMatPro file at path: " + file_path);
    }

    std::string search_keyword;
    std::string lower_name = table_name;
    std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);

    // Map table_name to Abaqus keyword
    if (lower_name == "density") {
        search_keyword = "*Density";
    } else if (lower_name == "thermal_conductivity") {
        search_keyword = "*CONDUCTIVITY";
    } else if (lower_name == "specific_heat" || lower_name == "specific heat") {
        search_keyword = "*SPECIFIC HEAT";
    } else if (lower_name == "elastic") {
        search_keyword = "*ELASTIC";
    } else {
        // Fallback: assume the user provided the keyword name (without *)
        search_keyword = "*" + table_name;
    }

    // Ensure search_keyword is uppercase for case-insensitive comparison
    std::transform(search_keyword.begin(), search_keyword.end(), search_keyword.begin(), ::toupper);

    std::vector<std::vector<double>> data_buffer;
    std::string line;
    bool in_block = false;

    while (std::getline(file, line)) {
        std::string trimmed_line = trim(line);
        if (trimmed_line.empty()) continue;

        // Check if line is a keyword (starts with * but not **)
        if (trimmed_line.rfind("*", 0) == 0 && trimmed_line.rfind("**", 0) != 0) {
            std::string line_upper = trimmed_line;
            std::transform(line_upper.begin(), line_upper.end(), line_upper.begin(), ::toupper);

            // Check if this is the keyword we are looking for
            // We check if the line starts with the keyword.
            // Abaqus keywords can have parameters, e.g., *ELASTIC, TYPE=ISOTROPIC
            if (line_upper.find(search_keyword) == 0) {
                in_block = true;
                continue;
            } else {
                if (in_block) {
                    // We hit a new keyword, so the block for our table has ended
                    break;
                }
            }
        }

        if (in_block) {
            // Skip comments
            if (trimmed_line.rfind("**", 0) == 0) continue;

            // Parse data line
            std::vector<std::string> tokens = split(trimmed_line, ",");
            std::vector<double> row_data;
            bool parse_error = false;
            
            for (const auto& token : tokens) {
                std::string trimmed_token = trim(token);
                if (trimmed_token.empty()) continue;
                try {
                    row_data.push_back(std::stod(trimmed_token));
                } catch (...) {
                    parse_error = true;
                    break; 
                }
            }

            if (!parse_error && !row_data.empty()) {
                data_buffer.push_back(row_data);
            }
        }
    }

    if (data_buffer.empty()) {
        throw std::runtime_error("Could not find data for table: " + table_name + " in file: " + file_path);
    }

    size_t num_rows = data_buffer.size();
    size_t num_cols = data_buffer[0].size();

    // Verify column consistency
    for (const auto& row : data_buffer) {
        if (row.size() != num_cols) {
            throw std::runtime_error("Inconsistent column count in table: " + table_name + " in file: " + file_path);
        }
    }

    // Initialize the table
    // Table_t constructor: (num_rows, num_columns, name)
    table = Table_t(num_rows, num_cols, table_name.c_str());

    // Populate data
    for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < num_cols; ++j) {
            // Data needs to be stores as (temperature, value), they come in from the file as (value, temperature)
            table.set_value(i, j, data_buffer[i][num_cols - j - 1]);
        }
    }

    table.verify_data(0);

    // Sync to device
    table.update_device();
}