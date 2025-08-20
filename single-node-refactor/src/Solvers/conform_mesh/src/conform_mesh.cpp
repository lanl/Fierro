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

#include "conform_mesh.h"
#include "state.h"
#include "mesh.h"
#include "region_fill.h"
#include "material.h"
#include "boundary_conditions.h"
#include "state.h"
#include "simulation_parameters.h"
#include "geometry_new.h"



void ConformMesh::initialize(SimulationParameters_t& SimulationParamaters, 
                	   Material_t& Materials, 
                	   Mesh_t& mesh, 
                	   BoundaryCondition_t& Boundary,
                	   State_t& State) const
{
	const size_t num_nodes = mesh.num_nodes;
    const size_t num_gauss_pts = mesh.num_elems;
    const size_t num_corners = mesh.num_corners;
    const size_t num_dims = mesh.num_dims;
} // end solver initialization


/////////////////////////////////////////////////////////////////////////////
///
/// \fn setup the ConformMesh method
///
/// \brief Allocate state, setup models, and fill mesh regions per the YAML input
///
/////////////////////////////////////////////////////////////////////////////
void ConformMesh::setup(SimulationParameters_t& SimulationParamaters, 
                Material_t& Materials, 
                Mesh_t& mesh, 
                BoundaryCondition_t& Boundary,
                State_t& State)
{
    
} // end ConformMesh setup

void ConformMesh::execute(SimulationParameters_t& SimulationParamaters, 
                    Material_t& Materials, 
                    BoundaryCondition_t& BoundaryConditions, 
                    Mesh_t& mesh, 
                    State_t& State)
{

    double fuzz  = SimulationParamaters.dynamic_options.fuzz;  // 1.e-16
    double tiny  = SimulationParamaters.dynamic_options.tiny;  // 1.e-12
    double small = SimulationParamaters.dynamic_options.small; // 1.e-8

    double graphics_dt_ival  = SimulationParamaters.output_options.graphics_time_step;
    int    graphics_cyc_ival = SimulationParamaters.output_options.graphics_iteration_step;

    // double time_initial = SimulationParamaters.dynamic_options.time_initial;
    double time_final   = this->time_end; //SimulationParamaters.dynamic_options.time_final;
    double dt_min   = SimulationParamaters.dynamic_options.dt_min;
    double dt_max   = SimulationParamaters.dynamic_options.dt_max;
    double dt_start = SimulationParamaters.dynamic_options.dt_start;
    double dt_cfl   = SimulationParamaters.dynamic_options.dt_cfl;


} // end of SGH execute


// ==============================================================
// ------------------- VOXELIZATION FUNCTIONS -------------------
// ==============================================================

// BINARY STL READER - (Note: it can ONLY read binary stl files)
std::tuple<CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, CArray<float>, unsigned int> ConformMesh::binary_stl_reader(std::string stl_file_path){
    
    // Open the binary file
    std::string filename = stl_file_path;
    auto input = std::ifstream{filename, std::ifstream::in | std::ifstream::binary};
    
    // Check that the file is actually open
    if (input.is_open()) {
        std::cout << "Opening .stl file... \n";
    }
    else {
        std::cout << "WARNING: .stl file is not open \n";
    }
    
    // Initialize variables
    char heading[81];
    unsigned int n_facets;
    
    // Read the file's heading
    char* ptr1 = heading;
    input.read(ptr1, 80);
    heading[80] = '\0'; // to ensure a proper C string
    std::cout << "File heading: " << heading << "\n";
    
    // Read the number of facets in the file
    unsigned int* ptr2 = &n_facets;
    input.read(reinterpret_cast<char*>(ptr2), sizeof(unsigned int));
    std::cout << "Voxelizing " << n_facets << " facets..." << "\n";
    
    // Initialize storage variables
    CArray<float> normal(n_facets,3);
    CArray<float> v1X(n_facets);
    CArray<float> v1Y(n_facets);
    CArray<float> v1Z(n_facets);
    CArray<float> v2X(n_facets);
    CArray<float> v2Y(n_facets);
    CArray<float> v2Z(n_facets);
    CArray<float> v3X(n_facets);
    CArray<float> v3Y(n_facets);
    CArray<float> v3Z(n_facets);
    float normalp[3];
    float v1p[3];
    float v2p[3];
    float v3p[3];
    
    // Pull data from file
    float* ptr3 = normalp;
    float* ptr4 = v1p;
    float* ptr5 = v2p;
    float* ptr6 = v3p;
    for (size_t i = 0; i < n_facets; ++i) {
        size_t n_bytes = 3 * sizeof(float);
        input.read(reinterpret_cast<char*>(ptr3), n_bytes);
        input.read(reinterpret_cast<char*>(ptr4), n_bytes);
        input.read(reinterpret_cast<char*>(ptr5), n_bytes);
        input.read(reinterpret_cast<char*>(ptr6), n_bytes);
        input.seekg(2,input.cur);
        for (int j=0; j<3; j++) {
            normal(i,j) = normalp[j];
        }
        v1X(i) = v1p[0];
        v1Y(i) = v1p[1];
        v1Z(i) = v1p[2];
        v2X(i) = v2p[0];
        v2Y(i) = v2p[1];
        v2Z(i) = v2p[2];
        v3X(i) = v3p[0];
        v3Y(i) = v3p[1];
        v3Z(i) = v3p[2];
//        v1X.host(i) = v1p[0];
//        v1Y.host(i) = v1p[1];
//        v1Z.host(i) = v1p[2];
//        v2X.host(i) = v2p[0];
//        v2Y.host(i) = v2p[1];
//        v2Z.host(i) = v2p[2];
//        v3X.host(i) = v3p[0];
//        v3Y.host(i) = v3p[1];
//        v3Z.host(i) = v3p[2];
    }
    input.close();
    return {normal, v1X, v1Y, v1Z, v2X, v2Y, v2Z, v3X, v3Y, v3Z, n_facets};
}