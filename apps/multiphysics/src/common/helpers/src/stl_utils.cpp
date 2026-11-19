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

#include "stl_utils.hpp"
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>

std::tuple<
    CArray<real_t>,   // normal
    CArray<real_t>, CArray<real_t>, CArray<real_t>,   // v0X, v0Y, v0Z 
    CArray<real_t>, CArray<real_t>, CArray<real_t>,   // v1X, v1Y, v1Z
    CArray<real_t>, CArray<real_t>, CArray<real_t>,   // v2X, v2Y, v2Z
    size_t // n_facets
> binary_stl_reader(const std::string& path)
{
    
    // Open file in binary mode and position at end
    std::ifstream in(path, std::ios::binary | std::ios::ate);

    // Check file opened successfully
    if (!in) {
        throw std::runtime_error(
            "ERROR: Failed to open STL file: " + path + "\n"
        );
    }

    // Get file size
    const std::streamoff filesize = in.tellg();

    // Basic sanity check
    if (filesize < 100) {
        throw std::runtime_error(
            "ERROR: STL file too small to be valid: " + path + "\n"
        );
    }

    // Move back to beginning of file
    in.seekg(0);

    std::cout << "Successfully opened STL file.\n";
    std::cout << "File size = " << filesize << " bytes\n";

    // Continue reading/parsing STL file here...


    // ---- check if ASCII -------------------------------------------------
    char magic[6] = { 0 };
    in.read(magic, 5);          // read first 5 chars
    in.seekg(0);               // rewind
    if (std::strncmp(magic, "solid", 5) == 0) {
        std::cerr
            << "ERROR: \"" << path
            << "\" looks like an **ASCII** STL (starts with \"solid\").\n"
            << "Re‑export it as *binary* or implement an ASCII parser.\n";
        std::exit(EXIT_FAILURE);        // or call ascii_stl_reader();
    }

    // ---- read 80‑byte header + nominal facet count ----------------------
    char header[80];                in.read(header, 80);
    size_t n_facets_nominal;  in.read(reinterpret_cast<char*>(&n_facets_nominal), 4);

    // ---- compute expected count from file size to sanity‑check ----------
    // binary facet record = 50 bytes (12×4 + 12×4 + 12×4 + 2)
    const size_t n_facets_from_size =
        static_cast<size_t>((filesize - 84) / 50);

    uint32_t n_facets = n_facets_nominal;
    if (n_facets_nominal != n_facets_from_size) {
        std::cout << "WARNING: facet count in header (" << n_facets_nominal
            << ") disagrees with file size (" << n_facets_from_size
            << ").  Using size‑derived value.\n";
        n_facets = n_facets_from_size;
    }
    std::cout << "STL facets: " << n_facets << '\n';

    // ---- allocate MATAR arrays -----------------------------------------
    CArray<real_t> normal(n_facets, 3);
    CArray<real_t> v0X(n_facets), v0Y(n_facets), v0Z(n_facets);
    CArray<real_t> v1X(n_facets), v1Y(n_facets), v1Z(n_facets);
    CArray<real_t> v2X(n_facets), v2Y(n_facets), v2Z(n_facets);

    // ---- read facet records --------------------------------------------
    float nrm[3], v0[3], v1[3], v2[3];
    for (size_t i = 0; i < n_facets; ++i) {
        in.read(reinterpret_cast<char*>(nrm), 12);
        in.read(reinterpret_cast<char*>(v0), 12);
        in.read(reinterpret_cast<char*>(v1), 12);
        in.read(reinterpret_cast<char*>(v2), 12);
        in.ignore(2);                        // attribute byte count

        for (int dim = 0; dim < 3; ++dim){ 
             normal(i, dim) = (double)nrm[dim]; 
        } // end for


        v0X(i) = (double)v0[0]; 
        v0Y(i) = (double)v0[1]; 
        v0Z(i) = (double)v0[2];

        v1X(i) = (double)v1[0]; 
        v1Y(i) = (double)v1[1]; 
        v1Z(i) = (double)v1[2];

        v2X(i) = (double)v2[0]; 
        v2Y(i) = (double)v2[1]; 
        v2Z(i) = (double)v2[2];
        
    } // end for facet

    return { normal,v0X,v0Y,v0Z,v1X,v1Y,v1Z,v2X,v2Y,v2Z,n_facets };

} // end of function to read STL file
