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
#ifndef FIERRO_IO_H
#define FIERRO_IO_H

#include "matar.h"
#include "mesh.h"
#include "state.h"
#include "simulation_parameters.h"
#include "region.h"
#include "string_utils.h"

#include <map>
#include <memory>
#include <cstring>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_id
///
/// \brief This gives the index value of the point or the elem
///
/// Assumes that the grid has an i,j,k structure
/// the elem = i + (j)*(num_points_i-1) + (k)*(num_points_i-1)*(num_points_j-1)
/// the point = i + (j)*num_points_i + (k)*num_points_i*num_points_j
///
/// \param i index
/// \param j index
/// \param k index
/// \param Number of i indices
/// \param Number of j indices
///
/////////////////////////////////////////////////////////////////////////////
inline int get_id(int i, int j, int k, int num_i, int num_j)
{
    return i + j * num_i + k * num_i * num_j;
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_id_device
///
/// \brief This gives the index value of the point or the elem
///
/// Assumes that the grid has an i,j,k structure
/// the elem = i + (j)*(num_points_i-1) + (k)*(num_points_i-1)*(num_points_j-1)
/// the point = i + (j)*num_points_i + (k)*num_points_i*num_points_j
///
/// \param i index
/// \param j index
/// \param k index
/// \param Number of i indices
/// \param Number of j indices
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_INLINE_FUNCTION
int get_id_device(int i, int j, int k, int num_i, int num_j)
{
    return i + j * num_i + k * num_i * num_j;
}


/////////////////////////////////////////////////////////////////////////////
///
/// \class MeshReader
///
/// \brief Class for simplifying reading meshes
///
/// This class contains the requisite functions required to read different
/// mesh formats. The idea is to set the mesh file name, and parse the
/// extension to decide which reader to use. Currently, only ensight .geo
/// files are supported.
///
/////////////////////////////////////////////////////////////////////////////
class MeshReader
{
private:
    // Handy structs for parsing input meshes
    struct Node {
        int id;
        double x, y, z;
    };

    struct Element {
        int id;
        std::vector<int> connectivity; 
    };

public:

    char* mesh_file_ = NULL;

    MeshReader() {} // Simulation_Parameters& _simparam);

    ~MeshReader() = default;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn set_mesh_file
    ///
    /// \brief Sets the mesh file path for reading in a mesh
    ///
    /// \param Path to mesh file
    ///
    /////////////////////////////////////////////////////////////////////////////
    void set_mesh_file(char* MESH)
    {
        mesh_file_ = MESH;
    }

    // Reads and initializes the mesh and geometric state entities
    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn read_mesh
    ///
    /// \brief Read mesh from file
    ///
    /// \param Simulation mesh
    /// \param Simulation state
    /// \param Number of dimensions
    /// \param Number of RK bins
    ///
    ///
    /////////////////////////////////////////////////////////////////////////////
    void read_mesh(Mesh_t& mesh,
                   State_t& State,
                   mesh_input_t& mesh_inps,
                   int      num_dims,
                   int      rk_num_bins)
    {
        if (mesh_file_ == NULL) {
            throw std::runtime_error("**** No mesh path given for read_mesh ****");
        }

        std::ifstream file(mesh_file_);
        if (file.is_open()) {
            std::cout << "The file exists." << std::endl;
            file.close();
        } else {
            throw std::runtime_error("**** Mesh path given does not exists ****");
        }

        // Check mesh file extension
        // and read based on extension
        std::string filePathStr(mesh_file_);
        std::string extension;

        size_t pos = filePathStr.rfind('.');
        if (pos != std::string::npos) {
            extension = filePathStr.substr(pos + 1);
        } else {
            extension =  "";
        }

        std::cout << "File extension is: " << extension << std::endl;

        if(extension == "geo"){ // Ensight meshfile extension
            read_ensight_mesh(mesh, State.GaussPoints, State.node, State.corner, mesh_inps, num_dims, rk_num_bins);
        }
        else if(extension == "inp"){ // Abaqus meshfile extension
            read_Abaqus_mesh(mesh, State, num_dims, rk_num_bins);
        }
        else if(extension == "vtk"){ // vtk file format

            read_vtk_mesh(mesh, State.GaussPoints, State.node, State.corner, mesh_inps, num_dims, rk_num_bins);

            //throw std::runtime_error("**** VTK mesh reader not yet implemented ****");
            // read_VTK_mesh(mesh, State, num_dims, rk_num_bins);
        }
        else{
            throw std::runtime_error("**** Mesh file extension not understood ****");
        }

    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn read_ensight_mesh
    ///
    /// \brief Read .geo mesh file
    ///
    /// \param Simulation mesh
    /// \param Element state struct
    /// \param Node state struct
    /// \param Corner state struct
    /// \param Number of dimensions
    /// \param Number of RK bins
    ///
    /////////////////////////////////////////////////////////////////////////////
    void read_ensight_mesh(Mesh_t& mesh,
                           GaussPoint_t& GaussPoints,
                           node_t&   node,
                           corner_t& corner,
                           mesh_input_t& mesh_inps,
                           int num_dims,
                           int rk_num_bins)
    {
        FILE* in;
        char  ch;

        size_t num_nodes_in_elem = 1;
        for (int dim = 0; dim < num_dims; dim++) {
            num_nodes_in_elem *= 2;
        }

        // read the mesh    WARNING: assumes a .geo file
        in = fopen(mesh_file_, "r");

        // skip 8 lines
        for (int j = 1; j <= 8; j++) {
            int i = 0;
            while ((ch = (char)fgetc(in)) != '\n') {
                i++;
            }
        }

        // --- Read in the nodes in the mesh ---

        size_t num_nodes = 0;

        fscanf(in, "%lu", &num_nodes);
        printf("Number if nodes read in %lu\n", num_nodes);

        
        mesh.initialize_nodes(num_nodes);

        // initialize node state variables, for now, we just need coordinates, the rest will be initialize by the respective solvers
        std::vector<node_state> required_node_state = { node_state::coords };
        node.initialize(rk_num_bins, num_nodes, num_dims, required_node_state);

        // read the initial mesh coordinates
        // x-coords
        for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
            fscanf(in, "%le", &node.coords.host(0, node_id, 0));
            node.coords.host(0, node_id, 0)*= mesh_inps.scale_x;
        }

        // y-coords
        for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
            fscanf(in, "%le", &node.coords.host(0, node_id, 1));
            node.coords.host(0, node_id, 1)*= mesh_inps.scale_y;
        }

        // z-coords
        for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
            if (num_dims == 3) {
                fscanf(in, "%le", &node.coords.host(0, node_id, 2));
                node.coords.host(0, node_id, 2)*= mesh_inps.scale_z;
            }
            else{
                double dummy;
                fscanf(in, "%le", &dummy);
            }
        } // end for

        // save the node coords to the current RK value
        for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
            for (int rk = 1; rk < rk_num_bins; rk++) {
                for (int dim = 0; dim < num_dims; dim++) {
                    node.coords.host(rk, node_gid, dim) = node.coords.host(0, node_gid, dim);
                } // end for dim
            } // end for rk
        } // end parallel for

        // Update device nodal positions
        node.coords.update_device();

        ch = (char)fgetc(in);

        // skip 1 line
        for (int j = 1; j <= 1; j++) {
            int i = 0;
            while ((ch = (char)fgetc(in)) != '\n') {
                i++;
            }
        }

        // --- read in the elements in the mesh ---
        size_t num_elem = 0;

        fscanf(in, "%lu", &num_elem);
        printf("Number of elements read in %lu\n", num_elem);

        // initialize elem variables
        mesh.initialize_elems(num_elem, num_dims);
        // GaussPoints.initialize(rk_num_bins, num_elem, 3); // always 3D here, even for 2D

        
        // for each cell read the list of associated nodes
        for (int elem_gid = 0; elem_gid < num_elem; elem_gid++) {
            for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                fscanf(in, "%lu", &mesh.nodes_in_elem.host(elem_gid, node_lid));  // %d vs zu

                // shift to start node index space at 0
                mesh.nodes_in_elem.host(elem_gid, node_lid) -= 1;
            }
        }

        // Convert from ensight to IJK mesh
        int convert_ensight_to_ijk[8];
        convert_ensight_to_ijk[0] = 0;
        convert_ensight_to_ijk[1] = 1;
        convert_ensight_to_ijk[2] = 3;
        convert_ensight_to_ijk[3] = 2;
        convert_ensight_to_ijk[4] = 4;
        convert_ensight_to_ijk[5] = 5;
        convert_ensight_to_ijk[6] = 7;
        convert_ensight_to_ijk[7] = 6;

        int tmp_ijk_indx[8];

        for (int elem_gid = 0; elem_gid < num_elem; elem_gid++) {
            for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                tmp_ijk_indx[node_lid] = mesh.nodes_in_elem.host(elem_gid, convert_ensight_to_ijk[node_lid]);
            }

            for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
                mesh.nodes_in_elem.host(elem_gid, node_lid) = tmp_ijk_indx[node_lid];
            }
        }
        // update device side
        mesh.nodes_in_elem.update_device();

        // initialize corner variables
        int num_corners = num_elem * mesh.num_nodes_in_elem;
        mesh.initialize_corners(num_corners);
        // corner.initialize(num_corners, num_dims);

        // Close mesh input file
        fclose(in);

        // Build connectivity
        mesh.build_connectivity();

        return;
    } // end read ensight mesh

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn read_Abaqus_mesh
    ///
    /// \brief Read .inp mesh file
    ///
    /// \param Simulation mesh
    /// \param Simulation state
    /// \param Node state struct
    /// \param Number of dimensions
    /// \param Number of RK bins
    ///
    /////////////////////////////////////////////////////////////////////////////
    void read_Abaqus_mesh(Mesh_t& mesh,
                          State_t& State,
                          int      num_dims,
                          int      rk_num_bins)
    {

        std::cout<<"Reading abaqus input file for mesh"<<std::endl;
        std::ifstream inputFile(mesh_file_);
        if (!inputFile.is_open()) {
            std::cerr << "Failed to open the file." << std::endl;

        }

        std::vector<Node> nodes;
        std::vector<Element> elements;

        std::string line;
        bool readingNodes = false;
        bool readingElements = false;

        while (std::getline(inputFile, line)) {
            if (line.find("*Node") != std::string::npos) {
                readingNodes = true;
                std::cout<<"Found *Node"<<std::endl;

            } 
            else if (readingNodes && !line.find("*") ) { // End of nodes
                readingNodes = false;
            } 
            else if (readingNodes) {
                // std::cout<<"Reading Nodes"<<std::endl;
                std::istringstream iss(line);
                std::ws(iss); // Skip leading whitespace
                std::string token;
                Node node;

                if (!(iss >> node.id && std::getline(iss, token, ',') && iss >> node.x &&
                    std::getline(iss, token, ',') && iss >> node.y &&
                    std::getline(iss, token, ',') && iss >> node.z)) {
                    std::cerr << "Failed to parse line: " << line << std::endl;
                    continue; // Skip this line if parsing failed
                }
                nodes.push_back(node);
            }

            if (line.find("*Element") != std::string::npos) {
                readingElements = true;
                std::cout<<"Found *Element*"<<std::endl;
            } 
            else if (readingElements &&  !line.find("*") ) { // End of elements
                readingElements = false;
            } 
            else if (readingElements ) {
                std::istringstream iss(line);
                Element element;
                std::string token;

                if (!(iss >> element.id)){
                    std::cout << "Failed to parse line: " << line << std::endl;
                    continue; // Skip this line if parsing failed
                } 

                while ((std::getline(iss, token, ','))) { 
                    // Now extract the integer, ignoring any trailing whitespace
                    int val;
                    iss >> val;
                    element.connectivity.push_back(val);
                }

                // Convert from abaqus to IJK mesh
                int convert_abq_to_ijk[8];
                convert_abq_to_ijk[0] = 0;
                convert_abq_to_ijk[1] = 1;
                convert_abq_to_ijk[2] = 3;
                convert_abq_to_ijk[3] = 2;
                convert_abq_to_ijk[4] = 4;
                convert_abq_to_ijk[5] = 5;
                convert_abq_to_ijk[6] = 7;
                convert_abq_to_ijk[7] = 6;

                int tmp_ijk_indx[8];

                for (int node_lid = 0; node_lid < 8; node_lid++) {
                    tmp_ijk_indx[node_lid] = element.connectivity[convert_abq_to_ijk[node_lid]];
                }

                for (int node_lid = 0; node_lid < 8; node_lid++){
                    element.connectivity[node_lid] = tmp_ijk_indx[node_lid];
                }

                elements.push_back(element);
            }
        }

        inputFile.close();

        size_t num_nodes = nodes.size();

        printf("Number if nodes read in %lu\n", num_nodes);

        // initialize node variables
        mesh.initialize_nodes(num_nodes);

        // initialize node state, for now, we just need coordinates, the rest will be initialize by the respective solvers
        std::vector<node_state> required_node_state = { node_state::coords };

        State.node.initialize(rk_num_bins, num_nodes, num_dims, required_node_state);


        // Copy nodes to mesh
        for(int node_gid = 0; node_gid < num_nodes; node_gid++){
            State.node.coords.host(0, node_gid, 0) = nodes[node_gid].x;
            State.node.coords.host(0, node_gid, 1) = nodes[node_gid].y;
            State.node.coords.host(0, node_gid, 2) = nodes[node_gid].z;
        }

        // save the node coords to the current RK value
        for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
            for (int rk = 1; rk < rk_num_bins; rk++) {
                for (int dim = 0; dim < num_dims; dim++) {
                    State.node.coords.host(rk, node_gid, dim) = State.node.coords.host(0, node_gid, dim);
                } // end for dim
            } // end for rk
        } // end parallel for

        // Update device nodal positions
        State.node.coords.update_device();


        // --- read in the elements in the mesh ---
        size_t num_elem = elements.size();
        printf("Number of elements read in %lu\n", num_elem);

        // initialize elem variables
        mesh.initialize_elems(num_elem, num_dims);
        // State.GaussPoints.initialize(rk_num_bins, num_elem, 3); // always 3D here, even for 2D


        // for each cell read the list of associated nodes
        for (int elem_gid = 0; elem_gid < num_elem; elem_gid++) {
            for (int node_lid = 0; node_lid < 8; node_lid++) {
                mesh.nodes_in_elem.host(elem_gid, node_lid) = elements[elem_gid].connectivity[node_lid];

                // shift to start node index space at 0
                mesh.nodes_in_elem.host(elem_gid, node_lid) -= 1;
            }
        }

        // update device side
        mesh.nodes_in_elem.update_device();

        // initialize corner variables
        int num_corners = num_elem * mesh.num_nodes_in_elem;
        mesh.initialize_corners(num_corners);
        // State.corner.initialize(num_corners, num_dims);

        // Build connectivity
        mesh.build_connectivity();
    } // end read abaqus mesh


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn read_vtk_mesh
    ///
    /// \brief Read ASCII .vtk mesh file
    ///
    /// \param Simulation mesh
    /// \param Simulation state
    /// \param Node state struct
    /// \param Number of dimensions
    /// \param Number of RK bins
    ///
    /////////////////////////////////////////////////////////////////////////////
    void read_vtk_mesh(Mesh_t& mesh,
                    GaussPoint_t& GaussPoints,
                    node_t&   node,
                    corner_t& corner,
                    mesh_input_t& mesh_inps,
                    int num_dims,
                    int rk_num_bins)
    {

        std::cout<<"Reading VTK mesh"<<std::endl;
    
        int i;           // used for writing information to file
        int node_gid;    // the global id for the point
        int elem_gid;     // the global id for the elem

        size_t num_nodes_in_elem = 1;
        for (int dim = 0; dim < num_dims; dim++) {
            num_nodes_in_elem *= 2;
        }
        

        std::string token;
        
        bool found = false;
        
        std::ifstream in;  // FILE *in;
        in.open(mesh_file_);
        

        // look for POINTS
        i = 0;
        while (found==false) {
            std::string str;
            std::string delimiter = " ";
            std::getline(in, str);
            std::vector<std::string> v = split (str, delimiter);
            
            // looking for the following text:
            //      POINTS %d float
            if(v[0] == "POINTS"){
                size_t num_nodes = std::stoi(v[1]);
                printf("Num nodes read in %zu\n", num_nodes);
                mesh.initialize_nodes(num_nodes);

                std::vector<node_state> required_node_state = { node_state::coords };
                node.initialize(rk_num_bins, num_nodes, num_dims, required_node_state);
                
                found=true;
            } // end if
            
            
            if (i>1000){
                std::cerr << "ERROR: Failed to find POINTS in file" << std::endl;
                break;
            } // end if
            
            i++;
        } // end while
        
        // read the node coordinates
        for (node_gid=0; node_gid<mesh.num_nodes; node_gid++){
            
            std::string str;
            std::getline(in, str);
            
            std::string delimiter = " ";
            std::vector<std::string> v = split (str, delimiter);
            
            // save the nodal coordinates
            node.coords.host(0, node_gid, 0) = mesh_inps.scale_x*std::stod(v[0]); // double
            node.coords.host(0, node_gid, 1) = mesh_inps.scale_y*std::stod(v[1]); // double
            if(num_dims==3){
                node.coords.host(0, node_gid, 2) = mesh_inps.scale_z*std::stod(v[2]); // double
            }

            
        } // end for nodes

        // save the node coords to the current RK value
        for (size_t node_gid = 0; node_gid < mesh.num_nodes; node_gid++) {
            for (int rk = 1; rk < rk_num_bins; rk++) {
                for (int dim = 0; dim < num_dims; dim++) {
                    node.coords.host(rk, node_gid, dim) = node.coords.host(0, node_gid, dim);
                } // end for dim
            } // end for rk
        } // end parallel for

        // Update device nodal positions
        node.coords.update_device();
        

        found=false;

        // look for CELLS
        i = 0;
        size_t num_elem = 0;
        while (found==false) {
            std::string str;
            std::getline(in, str);
            
            std::string delimiter = " ";
            std::vector<std::string> v = split (str, delimiter);
            std::cout << v[0] << std::endl; // printing
            
            // looking for the following text:
            //      CELLS num_elem size
            if(v[0] == "CELLS"){
                num_elem = std::stoi(v[1]);
                printf("Num elements read in %zu\n", num_elem);

                // initialize elem variables
                mesh.initialize_elems(num_elem, num_dims);
                
                found=true;
            } // end if
            
            
            if (i>1000){
                printf("ERROR: Failed to find CELLS \n");
                break;
            } // end if
            
            i++;
        } // end while
        
        
        // read the node ids in the element
        for (elem_gid=0; elem_gid<num_elem; elem_gid++) {
            
            std::string str;
            std::getline(in, str);
            
            std::string delimiter = " ";
            std::vector<std::string> v = split (str, delimiter);
            num_nodes_in_elem = std::stoi(v[0]);
            
            for (size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++){
                mesh.nodes_in_elem.host(elem_gid, node_lid) = std::stod(v[node_lid+1]);
                //printf(" %zu ", elem_point_list(elem_gid,node_lid) ); // printing
            }
            //printf("\n"); // printing
            
        } // end for

        // Convert from ensight to IJK mesh
        size_t convert_ensight_to_ijk[8];
        convert_ensight_to_ijk[0] = 0;
        convert_ensight_to_ijk[1] = 1;
        convert_ensight_to_ijk[2] = 3;
        convert_ensight_to_ijk[3] = 2;
        convert_ensight_to_ijk[4] = 4;
        convert_ensight_to_ijk[5] = 5;
        convert_ensight_to_ijk[6] = 7;
        convert_ensight_to_ijk[7] = 6;

        size_t tmp_ijk_indx[8];

        for (size_t elem_gid = 0; elem_gid < num_elem; elem_gid++) {
            for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                tmp_ijk_indx[node_lid] = mesh.nodes_in_elem.host(elem_gid, convert_ensight_to_ijk[node_lid]);
            }

            for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
                mesh.nodes_in_elem.host(elem_gid, node_lid) = tmp_ijk_indx[node_lid];
            }
        }
        // update device side
        mesh.nodes_in_elem.update_device();


        // initialize corner variables
        size_t num_corners = num_elem * num_nodes_in_elem;
        mesh.initialize_corners(num_corners);


        // Build connectivity
        mesh.build_connectivity();


        found=false;

        printf("\n");
        
        
        // look for CELL_TYPE
        i = 0;
        size_t elem_type = 0;
        while (found==false) {
            std::string str;
            std::string delimiter = " ";
            std::getline(in, str);
            std::vector<std::string> v = split (str, delimiter);
            
            // looking for the following text:
            //      CELLS num_elem size
            if(v[0] == "CELL_TYPES"){

                std::getline(in, str);
                elem_type = std::stoi(str);
                
                found=true;
            } // end if
            
            
            if (i>1000){
                printf("ERROR: Failed to find elem_TYPE \n");
                break;
            } // end if
            
            i++;
        } // end while
        printf("elem type = %zu \n", elem_type);
        // elem types:
        // linear hex = 12, linear quad = 9
        found=false;
        
        
        if(num_nodes_in_elem==8 & elem_type != 12) {
            printf("wrong elem type of %zu \n", elem_type);
            std::cerr << "ERROR: incorrect element type in VTK file" << std::endl;
        }
        
        in.close();
        
    } // end of VTKread function

}; // end of Mesh reader class

/////////////////////////////////////////////////////////////////////////////
///
/// \class MeshBuilder
///
/// \brief Class for building simple meshes
///
/// This class contains the requisite functions required to build simple
/// 2D and 3D Box meshes as well as 2D polar meshes. It uses the parsed
/// simulation parameters to decide what type of mesh to build.
///
/////////////////////////////////////////////////////////////////////////////
class MeshBuilder
{
public:

    MeshBuilder() {}

    ~MeshBuilder()
    {
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_mesh
    ///
    /// \brief Build a mesh for Fierro based on the input instructions
    ///
    /// \param Simulation mesh that is built
    /// \param Element state data
    /// \param Node state data
    /// \param Corner state data
    /// \param Simulation parameters
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_mesh(Mesh_t& mesh,
        GaussPoint_t& GaussPoints,
        node_t&   node,
        corner_t& corner,
        SimulationParameters_t& SimulationParamaters)
    {
        if (SimulationParamaters.mesh_input.num_dims == 2) {
            if (SimulationParamaters.mesh_input.type == mesh_input::Polar) {
                build_2d_polar(mesh, GaussPoints, node, corner, SimulationParamaters);
            }
            else if (SimulationParamaters.mesh_input.type == mesh_input::Box) {
                build_2d_box(mesh, GaussPoints, node, corner, SimulationParamaters);
            }
            else{
                std::cout << "**** 2D MESH TYPE NOT SUPPORTED **** " << std::endl;
                std::cout << "Valid options are: " << std::endl;
                auto map = mesh_input_type_map;
                for (const auto& pair : map) {
                    std::cout << "\t" << pair.first << std::endl;
                }
                throw std::runtime_error("**** 2D MESH TYPE NOT SUPPORTED ****");
            }
        }
        else if (SimulationParamaters.mesh_input.num_dims == 3) {
            build_3d_box(mesh, GaussPoints, node, corner, SimulationParamaters);
        }
        else{
            throw std::runtime_error("**** ONLY 2D RZ OR 3D MESHES ARE SUPPORTED ****");
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_2d_box
    ///
    /// \brief Builds an unstructured 2D rectilinear mesh
    ///
    /// \param Simulation mesh that is built
    /// \param Element state data
    /// \param Node state data
    /// \param Corner state data
    /// \param Simulation parameters
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_2d_box(Mesh_t& mesh,
        GaussPoint_t& GaussPoints,
        node_t&   node,
        corner_t& corner,
        SimulationParameters_t& SimulationParamaters) const
    {
        printf(" Creating a 2D box mesh \n");

        const int num_dim = 2;

        const double lx = SimulationParamaters.mesh_input.length[0];
        const double ly = SimulationParamaters.mesh_input.length[1];

        const int num_elems_i = SimulationParamaters.mesh_input.num_elems[0];
        const int num_elems_j = SimulationParamaters.mesh_input.num_elems[1];

        const int num_points_i = num_elems_i + 1; // num points in x
        const int num_points_j = num_elems_j + 1; // num points in y

        const int num_nodes = num_points_i * num_points_j;

        const double dx = lx / ((double)num_elems_i);  // len/(num_elems_i)
        const double dy = ly / ((double)num_elems_j);  // len/(num_elems_j)

        const int num_elems = num_elems_i * num_elems_j;

        std::vector<double> origin(num_dim);
        // SimulationParamaters.mesh_input.origin.update_host();
        for (int i = 0; i < num_dim; i++) { origin[i] = SimulationParamaters.mesh_input.origin[i]; }

        // --- 2D parameters ---
        // const int num_faces_in_elem  = 4;  // number of faces in elem
        // const int num_points_in_elem = 4;  // number of points in elem
        // const int num_points_in_face = 2;  // number of points in a face
        // const int num_edges_in_elem  = 4;  // number of edges in a elem

        // --- mesh node ordering ---
        // Convert ijk index system to the finite element numbering convention
        // for vertices in elem
        auto convert_point_number_in_quad = CArray<int>(4);
        convert_point_number_in_quad(0) = 0;
        convert_point_number_in_quad(1) = 1;
        convert_point_number_in_quad(2) = 3;
        convert_point_number_in_quad(3) = 2;

        int rk_num_bins = SimulationParamaters.dynamic_options.rk_num_bins;

        // intialize node variables
        mesh.initialize_nodes(num_nodes);

        // initialize node state, for now, we just need coordinates, the rest will be initialize by the respective solvers
        std::vector<node_state> required_node_state = { node_state::coords };
        node.initialize(rk_num_bins, num_nodes, num_dim, required_node_state);

        // --- Build nodes ---

        // populate the point data structures
        for (int j = 0; j < num_points_j; j++) {
            for (int i = 0; i < num_points_i; i++) {
                // global id for the point
                int node_gid = get_id(i, j, 0, num_points_i, num_points_j);

                // store the point coordinates
                node.coords.host(0, node_gid, 0) = origin[0] + (double)i * dx;
                node.coords.host(0, node_gid, 1) = origin[1] + (double)j * dy;
            } // end for i
        } // end for j

        for (int rk_level = 1; rk_level < rk_num_bins; rk_level++) {
            for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
                node.coords.host(rk_level, node_gid, 0) = node.coords.host(0, node_gid, 0);
                node.coords.host(rk_level, node_gid, 1) = node.coords.host(0, node_gid, 1);
            }
        }
        node.coords.update_device();

        // initialize elem variables
        mesh.initialize_elems(num_elems, num_dim);

        // populate the elem center data structures
        for (int j = 0; j < num_elems_j; j++) {
            for (int i = 0; i < num_elems_i; i++) {
                // global id for the elem
                int elem_gid = get_id(i, j, 0, num_elems_i, num_elems_j);

                // store the point IDs for this elem where the range is
                // (i:i+1, j:j+1 for a linear quad
                int this_point = 0;

                for (int jcount = j; jcount <= j + 1; jcount++) {
                    for (int icount = i; icount <= i + 1; icount++) {
                        // global id for the points
                        int node_gid = get_id(icount, jcount, 0, num_points_i, num_points_j);

                        // convert this_point index to the FE index convention
                        int this_index = convert_point_number_in_quad(this_point);

                        // store the points in this elem according the the finite
                        // element numbering convention
                        mesh.nodes_in_elem.host(elem_gid, this_index) = node_gid;

                        // increment the point counting index
                        this_point = this_point + 1;
                    } // end for icount
                } // end for jcount
            } // end for i
        } // end for j

        // update device side
        mesh.nodes_in_elem.update_device();

        // intialize corner variables
        int num_corners = num_elems * mesh.num_nodes_in_elem;
        mesh.initialize_corners(num_corners);
        // corner.initialize(num_corners, num_dim);

        // Build connectivity
        mesh.build_connectivity();
    } // end build_2d_box

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_2d_polar
    ///
    /// \brief Builds an unstructured 2D polar mesh
    ///
    /// \param Simulation mesh that is built
    /// \param Element state data
    /// \param Node state data
    /// \param Corner state data
    /// \param Simulation parameters
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_2d_polar(Mesh_t& mesh,
        GaussPoint_t& GaussPoints,
        node_t&   node,
        corner_t& corner,
        SimulationParameters_t& SimulationParamaters) const
    {
        printf(" Creating a 2D polar mesh \n");

        int num_dim     = 2;
        int rk_num_bins = SimulationParamaters.dynamic_options.rk_num_bins;

        const double inner_radius = SimulationParamaters.mesh_input.inner_radius;
        const double outer_radius = SimulationParamaters.mesh_input.outer_radius;

        const double start_angle = PI / 180.0 * SimulationParamaters.mesh_input.starting_angle;
        const double end_angle   = PI / 180.0 * SimulationParamaters.mesh_input.ending_angle;

        const int num_elems_i = SimulationParamaters.mesh_input.num_radial_elems;
        const int num_elems_j = SimulationParamaters.mesh_input.num_angular_elems;

        const int num_points_i = num_elems_i + 1; // num points in x
        const int num_points_j = num_elems_j + 1; // num points in y

        const int num_nodes = num_points_i * num_points_j;

        const double dx = (outer_radius - inner_radius) / ((double)num_elems_i);  // len/(elems)
        const double dy = (end_angle - start_angle) / ((double)num_elems_j);  // len/(elems)

        const int num_elems = num_elems_i * num_elems_j;

        std::vector<double> origin(num_dim);

        for (int i = 0; i < num_dim; i++) { origin[i] = SimulationParamaters.mesh_input.origin[i]; }

        // --- 2D parameters ---
        // const int num_faces_in_elem  = 4;  // number of faces in elem
        // const int num_points_in_elem = 4;  // number of points in elem
        // const int num_points_in_face = 2;  // number of points in a face
        // const int num_edges_in_elem  = 4;  // number of edges in a elem

        // --- mesh node ordering ---
        // Convert ijk index system to the finite element numbering convention
        // for vertices in elem
        auto convert_point_number_in_quad = CArray<int>(4);
        convert_point_number_in_quad(0) = 0;
        convert_point_number_in_quad(1) = 1;
        convert_point_number_in_quad(2) = 3;
        convert_point_number_in_quad(3) = 2;

        // intialize node variables
        mesh.initialize_nodes(num_nodes);

        // initialize node state, for now, we just need coordinates, the rest will be initialize by the respective solvers
        std::vector<node_state> required_node_state = { node_state::coords };
        node.initialize(rk_num_bins, num_nodes, num_dim, required_node_state);

        // populate the point data structures
        for (int j = 0; j < num_points_j; j++) {
            for (int i = 0; i < num_points_i; i++) {
                // global id for the point
                int node_gid = get_id(i, j, 0, num_points_i, num_points_j);

                double r_i     = inner_radius + (double)i * dx;
                double theta_j = start_angle + (double)j * dy;

                // store the point coordinates
                node.coords.host(0, node_gid, 0) = origin[0] + r_i * cos(theta_j);
                node.coords.host(0, node_gid, 1) = origin[1] + r_i * sin(theta_j);

                if(node.coords.host(0, node_gid, 0) < 0.0){
                    throw std::runtime_error("**** NODE RADIUS FOR RZ MESH MUST BE POSITIVE ****");
                }

            } // end for i
        } // end for j

        for (int rk_level = 1; rk_level < rk_num_bins; rk_level++) {
            for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
                node.coords.host(rk_level, node_gid, 0) = node.coords.host(0, node_gid, 0);
                node.coords.host(rk_level, node_gid, 1) = node.coords.host(0, node_gid, 1);
            }
        }
        node.coords.update_device();

        // initialize elem variables
        mesh.initialize_elems(num_elems, num_dim);

        // populate the elem center data structures
        for (int j = 0; j < num_elems_j; j++) {
            for (int i = 0; i < num_elems_i; i++) {
                // global id for the elem
                int elem_gid = get_id(i, j, 0, num_elems_i, num_elems_j);

                // store the point IDs for this elem where the range is
                // (i:i+1, j:j+1 for a linear quad
                int this_point = 0;

                for (int jcount = j; jcount <= j + 1; jcount++) {
                    for (int icount = i; icount <= i + 1; icount++) {
                        // global id for the points
                        int node_gid = get_id(icount, jcount, 0, num_points_i, num_points_j);

                        // convert this_point index to the FE index convention
                        int this_index = convert_point_number_in_quad(this_point);

                        // store the points in this elem according the the finite
                        // element numbering convention
                        mesh.nodes_in_elem.host(elem_gid, this_index) = node_gid;

                        // increment the point counting index
                        this_point = this_point + 1;
                    } // end for icount
                } // end for jcount
            } // end for i
        } // end for j

        // update device side
        mesh.nodes_in_elem.update_device();

        // intialize corner variables
        int num_corners = num_elems * mesh.num_nodes_in_elem;
        mesh.initialize_corners(num_corners);
        // corner.initialize(num_corners, num_dim);

        // Build connectivity
        mesh.build_connectivity();
    } // end build_2d_box

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_3d_box
    ///
    /// \brief Builds an unstructured 3D rectilinear mesh
    ///
    /// \param Simulation mesh that is built
    /// \param Element state data
    /// \param Node state data
    /// \param Corner state data
    /// \param Simulation parameters
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_3d_box(Mesh_t& mesh,
        GaussPoint_t& GaussPoints,
        node_t&   node,
        corner_t& corner,
        SimulationParameters_t& SimulationParamaters) const
    {
        printf(" Creating a 3D box mesh \n");

        const int num_dim = 3;

        // SimulationParamaters.mesh_input.length.update_host();
        const double lx = SimulationParamaters.mesh_input.length[0];
        const double ly = SimulationParamaters.mesh_input.length[1];
        const double lz = SimulationParamaters.mesh_input.length[2];

        // SimulationParamaters.mesh_input.num_elems.update_host();
        const int num_elems_i = SimulationParamaters.mesh_input.num_elems[0];
        const int num_elems_j = SimulationParamaters.mesh_input.num_elems[1];
        const int num_elems_k = SimulationParamaters.mesh_input.num_elems[2];

        const int num_points_i = num_elems_i + 1; // num points in x
        const int num_points_j = num_elems_j + 1; // num points in y
        const int num_points_k = num_elems_k + 1; // num points in y

        const int num_nodes = num_points_i * num_points_j * num_points_k;

        const double dx = lx / ((double)num_elems_i);  // len/(num_elems_i)
        const double dy = ly / ((double)num_elems_j);  // len/(num_elems_j)
        const double dz = lz / ((double)num_elems_k);  // len/(num_elems_k)

        const int num_elems = num_elems_i * num_elems_j * num_elems_k;

        std::vector<double> origin(num_dim);
        // SimulationParamaters.mesh_input.origin.update_host();
        for (int i = 0; i < num_dim; i++) { origin[i] = SimulationParamaters.mesh_input.origin[i]; }

        // --- 3D parameters ---
        // const int num_faces_in_elem  = 6;  // number of faces in elem
        // const int num_points_in_elem = 8;  // number of points in elem
        // const int num_points_in_face = 4;  // number of points in a face
        // const int num_edges_in_elem  = 12; // number of edges in a elem


        int rk_num_bins = SimulationParamaters.dynamic_options.rk_num_bins;

        // initialize mesh node variables
        mesh.initialize_nodes(num_nodes);

         // initialize node state variables, for now, we just need coordinates, the rest will be initialize by the respective solvers
        std::vector<node_state> required_node_state = { node_state::coords };
        node.initialize(rk_num_bins, num_nodes, num_dim, required_node_state);

        // --- Build nodes ---

        // populate the point data structures
        for (int k = 0; k < num_points_k; k++) {
            for (int j = 0; j < num_points_j; j++) {
                for (int i = 0; i < num_points_i; i++) {
                    // global id for the point
                    int node_gid = get_id(i, j, k, num_points_i, num_points_j);

                    // store the point coordinates
                    node.coords.host(0, node_gid, 0) = origin[0] + (double)i * dx;
                    node.coords.host(0, node_gid, 1) = origin[1] + (double)j * dy;
                    node.coords.host(0, node_gid, 2) = origin[2] + (double)k * dz;
                } // end for i
            } // end for j
        } // end for k

        for (int rk_level = 1; rk_level < rk_num_bins; rk_level++) {
            for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
                node.coords.host(rk_level, node_gid, 0) = node.coords.host(0, node_gid, 0);
                node.coords.host(rk_level, node_gid, 1) = node.coords.host(0, node_gid, 1);
                node.coords.host(rk_level, node_gid, 2) = node.coords.host(0, node_gid, 2);
            }
        }
        node.coords.update_device();

        // initialize elem variables
        mesh.initialize_elems(num_elems, num_dim);

        // --- Build elems  ---

        // populate the elem center data structures
        for (int k = 0; k < num_elems_k; k++) {
            for (int j = 0; j < num_elems_j; j++) {
                for (int i = 0; i < num_elems_i; i++) {
                    // global id for the elem
                    int elem_gid = get_id(i, j, k, num_elems_i, num_elems_j);

                    // store the point IDs for this elem where the range is
                    // (i:i+1, j:j+1, k:k+1) for a linear hexahedron
                    int this_point = 0;
                    for (int kcount = k; kcount <= k + 1; kcount++) {
                        for (int jcount = j; jcount <= j + 1; jcount++) {
                            for (int icount = i; icount <= i + 1; icount++) {
                                // global id for the points
                                int node_gid = get_id(icount, jcount, kcount,
                                                  num_points_i, num_points_j);

                                // convert this_point index to the FE index convention
                                int this_index = this_point; //convert_point_number_in_Hex(this_point);

                                // store the points in this elem according the the finite
                                // element numbering convention
                                mesh.nodes_in_elem.host(elem_gid, this_index) = node_gid;

                                // increment the point counting index
                                this_point = this_point + 1;
                            } // end for icount
                        } // end for jcount
                    }  // end for kcount
                } // end for i
            } // end for j
        } // end for k

        // update device side
        mesh.nodes_in_elem.update_device();

        // initialize corner variables
        int num_corners = num_elems * mesh.num_nodes_in_elem;
        mesh.initialize_corners(num_corners);
        // corner.initialize(num_corners, num_dim);

        // Build connectivity
        mesh.build_connectivity();
    } // end build_3d_box

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_3d_HexN_box
    ///
    /// \brief Builds an unstructured high order 3D rectilinear mesh
    ///
    /// \param Simulation mesh that is built
    /// \param Element state data
    /// \param Node state data
    /// \param Corner state data
    /// \param Simulation parameters
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_3d_HexN_box(Mesh_t& mesh,
        GaussPoint_t& GaussPoints,
        node_t&   node,
        corner_t& corner,
        SimulationParameters_t& SimulationParamaters) const
    {
        printf(" ***** WARNING::  build_3d_HexN_box not yet implemented\n");
        const int num_dim = 3;

        const int rk_num_bins = SimulationParamaters.dynamic_options.rk_num_bins;

        // SimulationParamaters.mesh_input.length.update_host();
        const double lx = SimulationParamaters.mesh_input.length[0];
        const double ly = SimulationParamaters.mesh_input.length[1];
        const double lz = SimulationParamaters.mesh_input.length[2];

        // SimulationParamaters.mesh_input.num_elems.update_host();
        const int num_elems_i = SimulationParamaters.mesh_input.num_elems[0];
        const int num_elems_j = SimulationParamaters.mesh_input.num_elems[1];
        const int num_elems_k = SimulationParamaters.mesh_input.num_elems[2];

        // creating zones for the Pn order
        const int Pn_order = SimulationParamaters.mesh_input.p_order;
        
        if (Pn_order > 19) {
            printf(" Fierro DG and RD solvers are only valid for elements up to Pn = 19 \n");
            return;
        }

        const int num_zones_i = Pn_order*num_elems_i;
        const int num_zones_j = Pn_order*num_elems_j;
        const int num_zones_k = Pn_order*num_elems_k;
        
        const int num_points_i = num_zones_i+1; // num points in x accounting for Pn
        const int num_points_j = num_zones_j+1; // num points in y accounting for Pn
        const int num_points_k = num_zones_k+1; // num points in y accounting for Pn
        
        
        const double dx = lx/((double)num_zones_i);  // len/(num_zones_i)
        const double dy = ly/((double)num_zones_j);  // len/(num_zones_j)
        const double dz = lz/((double)num_zones_k);  // len/(num_zones_k)
        
        const int num_elems = num_elems_i*num_elems_j*num_elems_k;
        // const int num_zones = num_zones_i*num_zones_j*num_zones_k; // accounts for Pn

        std::vector<double> origin(num_dim);
        for (int i = 0; i < num_dim; i++) { origin[i] = SimulationParamaters.mesh_input.origin[i]; }

        // --- 3D parameters ---
        // const int num_faces_in_zone = 6;   // number of faces in zone
        // const int num_points_in_zone = 8;  // number of points in zone
        // const int num_points_in_face = 4;  // number of points in a face
        
        // p_order   = 1, 2, 3, 4, 5
        // num_nodes = 2, 3, 4, 5, 6
        const int num_1D_points = Pn_order+1;
        const int num_points_in_elem = num_1D_points*num_1D_points*num_1D_points;
           
        
        // --- elem ---
        auto elem_coords = CArray <double> (num_elems, num_dim);
        auto elem_point_list = CArray <int> (num_elems, num_points_in_elem);
        
        
        // --- point ---
        int num_points = num_points_i * num_points_j * num_points_k;
        auto pt_coords = CArray <double> (num_points, num_dim);


        // --- Build nodes ---
        
        // initialize node variables
        mesh.initialize_nodes(num_points);

        // 
        std::vector<node_state> required_node_state = { node_state::coords };
        node.initialize(rk_num_bins, num_points, num_dim, required_node_state);
        // populate the point data structures
        for (int k = 0; k < num_points_k; k++){
            for (int j = 0; j < num_points_j; j++){
                for (int i = 0; i < num_points_i; i++){

                
                    // global id for the point
                    int node_gid = get_id(i, j, k, num_points_i, num_points_j);

                    // store the point coordinates
                    node.coords.host(0, node_gid, 0) = origin[0] + (double)i * dx;
                    node.coords.host(0, node_gid, 1) = origin[1] + (double)j * dy;
                    node.coords.host(0, node_gid, 2) = origin[2] + (double)k * dz;
                    
                } // end for k
            } // end for i
        } // end for j

        for (int rk_level = 1; rk_level < rk_num_bins; rk_level++) {
            for (int node_gid = 0; node_gid < num_points; node_gid++) {
                node.coords.host(rk_level, node_gid, 0) = node.coords.host(0, node_gid, 0);
                node.coords.host(rk_level, node_gid, 1) = node.coords.host(0, node_gid, 1);
                node.coords.host(rk_level, node_gid, 2) = node.coords.host(0, node_gid, 2);
            }
        }

        node.coords.update_device();


        // initialize elem variables
        mesh.initialize_elems(num_elems, num_dim);

        // --- Build elems  ---
        
        // populate the elem center data structures accounting for Pn
        for (int k=0; k<num_elems_k; k++){
            for (int j=0; j<num_elems_j; j++){
                for (int i=0; i<num_elems_i; i++){
                  
                    // global id for the elem
                    size_t elem_gid = get_id(i, j, k, num_elems_i, num_elems_j);
                    
                    // store the point IDs for this elem where the range is
                    // (i:i+1, j:j+1, k:k+1) for a linear hexahedron
                    // (i:(i+1)*Pn_order, j:(j+1)*Pn_order, k:(k+1)*Pn_order) for a Pn hexahedron
                    int node_lid = 0;
                    
                    int k_local = 0;
                    for (int kcount=k*Pn_order; kcount<=(k+1)*Pn_order; kcount++){
                        
                        int j_local = 0;
                        for (int jcount=j*Pn_order; jcount<=(j+1)*Pn_order; jcount++){
                            
                            int i_local = 0;
                            for (int icount=i*Pn_order; icount<=(i+1)*Pn_order; icount++){
                                
                                // global id for the points
                                size_t node_gid = get_id(icount, jcount, kcount,
                                                  num_points_i, num_points_j);

                                // Saved using i,j,k indexing
                                mesh.nodes_in_elem.host(elem_gid, node_lid) = node_gid;
                                
                                // increment the point counting index
                                node_lid = node_lid + 1;
                                
                                i_local++;
                            } // end for icount
                            
                            j_local++;
                        } // end for jcount
                        
                        k_local ++;
                    }  // end for kcount
                } // end for i
            } // end for j
        } // end for k

        // update device side
        mesh.nodes_in_elem.update_device();

        // initialize corner variables
        int num_corners = num_elems * mesh.num_nodes_in_elem;
        mesh.initialize_corners(num_corners);
        // corner.initialize(num_corners, num_dim);

        // Build connectivity
        mesh.build_connectivity();

    }
};

/////////////////////////////////////////////////////////////////////////////
///
/// \class MeshWriter
///
/// \brief Class for writing out a mesh with its associated state from Fierro
///
/// This class contains the requisite functions required to write out a mesh
/// with its associated state data from solvers in Fierro.
///
/////////////////////////////////////////////////////////////////////////////
class MeshWriter
{
private:
    int graphics_id = 0;

public:

    MeshWriter() {}

    ~MeshWriter()
    {
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn writes mesh with the format given in the input.yaml file
    ///
    /// \param Simulation mesh
    /// \param Element related state
    /// \param Node related state
    /// \param Corner related state
    /// \param Simulation input parameters
    ///
    /////////////////////////////////////////////////////////////////////////////
    void write_mesh(Mesh_t& mesh,
        State_t& State,
        SimulationParameters_t& SimulationParamaters,
        double dt,
        double time_value,
        CArray<double> graphics_times,
        std::vector<node_state> node_states,
        std::vector<gauss_pt_state> gauss_pt_states,
        std::vector<material_pt_state> material_pt_states)
    {
        if (SimulationParamaters.output_options.format == output_options::vtk) {
            write_vtk(mesh,
                      State,
                      SimulationParamaters,
                      time_value,
                      graphics_times,
                      node_states,
                      gauss_pt_states,
                      material_pt_states);
        }
        else if (SimulationParamaters.output_options.format == output_options::ensight) {
            write_ensight(mesh,
                          State,
                          SimulationParamaters,
                          dt,
                          time_value,
                          graphics_times,
                          node_states,
                          gauss_pt_states,
                          material_pt_states);
        }
        else if (SimulationParamaters.output_options.format == output_options::state) {
            write_material_point_state(mesh,
                                      State,
                                      SimulationParamaters,
                                      time_value,
                                      graphics_times,
                                      node_states,
                                      gauss_pt_states,
                                      material_pt_states);
        }
        else{
            std::cout << "**** MESH OUTPUT TYPE NOT SUPPORTED **** " << std::endl;
            std::cout << "Valid options are: " << std::endl;
            auto map = output_format_map;
            for (const auto& pair : map) {
                std::cout << "\t" << pair.first << std::endl;
            }
            throw std::runtime_error("**** MESH OUTPUT TYPE NOT SUPPORTED ****");
            // return;
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn write_ensight
    ///
    /// \brief Writes an ensight output file
    ///
    /// \param Simulation mesh
    /// \param State data
    /// \param Simulation parameters
    /// \param current time value
    /// \param Vector of all graphics output times
    ///
    /////////////////////////////////////////////////////////////////////////////
    void write_ensight(Mesh_t& mesh,
        State_t& State,
        SimulationParameters_t& SimulationParamaters,
        double dt,
        double time_value,
        CArray<double> graphics_times,
        std::vector<node_state> node_states,
        std::vector<gauss_pt_state> gauss_pt_states,
        std::vector<material_pt_state> material_pt_states)
    {
        size_t num_mats = State.MaterialPoints.size();

        // ---- Update host data ----

        // material point values
        for (int mat_id = 0; mat_id < num_mats; mat_id++) {
            State.MaterialPoints(mat_id).den.update_host();
            State.MaterialPoints(mat_id).pres.update_host();
            State.MaterialPoints(mat_id).stress.update_host();
            State.MaterialPoints(mat_id).sspd.update_host();
            State.MaterialPoints(mat_id).sie.update_host();
            State.MaterialPoints(mat_id).mass.update_host();
            State.MaterialPoints(mat_id).eroded.update_host();
        } // end for mat_id

        // gauss point values
        State.GaussPoints.vol.update_host();

        // nodal values
        State.node.coords.update_host();
        State.node.vel.update_host();
        State.node.mass.update_host();

        Kokkos::fence();

        // --------------------------

        const int num_scalar_vars = 10;
        const int num_vec_vars    = 3;

        std::string name_tmp;
        name_tmp = "Outputs_SGH";

        char* name = new char [name_tmp.length() + 1];
        std::strcpy(name, name_tmp.c_str());

        const char scalar_var_names[num_scalar_vars][15] = {
            "den", "pres", "sie", "vol", "mass", "sspd", "speed", "mat_id", "elem_switch", "eroded"
        };

        const char vec_var_names[num_vec_vars][15] = {
            "pos", "vel", "accel"
        };

        // short hand
        const size_t num_nodes = mesh.num_nodes;
        const size_t num_elems = mesh.num_elems;
        const size_t num_dims  = mesh.num_dims;

        // save the cell state to an array for exporting to graphics files
        auto elem_fields = CArray<double>(num_elems, num_scalar_vars);
        int  elem_switch = 1;


        DCArrayKokkos<double> speed(num_elems);
        FOR_ALL(elem_gid, 0, num_elems, {
            double elem_vel[3]; // note:initialization with a list won't work
            elem_vel[0] = 0.0;
            elem_vel[1] = 0.0;
            elem_vel[2] = 0.0;
            // get the coordinates of the element center
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                elem_vel[0] += State.node.vel(1, mesh.nodes_in_elem(elem_gid, node_lid), 0);
                elem_vel[1] += State.node.vel(1, mesh.nodes_in_elem(elem_gid, node_lid), 1);
                if (mesh.num_dims == 3) {
                    elem_vel[2] += State.node.vel(1, mesh.nodes_in_elem(elem_gid, node_lid), 2);
                }
                else{
                    elem_vel[2] = 0.0;
                }
            } // end loop over nodes in element
            elem_vel[0] = elem_vel[0] / mesh.num_nodes_in_elem;
            elem_vel[1] = elem_vel[1] / mesh.num_nodes_in_elem;
            elem_vel[2] = elem_vel[2] / mesh.num_nodes_in_elem;

            double speed_sqrd = 0.0;
            for (int dim = 0; dim < num_dims; dim++) {
                speed_sqrd += elem_vel[dim] * elem_vel[dim];
            }
            speed(elem_gid) = sqrt(speed_sqrd);
        }); // end parallel for
        speed.update_host();
        Kokkos::fence();

        // save the output scale fields to a single 2D array

        // export material centeric data to the elements
        for (int mat_id = 0; mat_id < num_mats; mat_id++) {
            size_t num_mat_elems = State.MaterialToMeshMaps(mat_id).num_material_elems;

            for (size_t mat_elem_lid = 0; mat_elem_lid < num_mat_elems; mat_elem_lid++) {
                // 1 material per element

                // get elem gid
                size_t elem_gid = State.MaterialToMeshMaps(mat_id).elem.host(mat_elem_lid);

                // save outputs
                elem_fields(elem_gid, 0) = State.MaterialPoints(mat_id).den.host(mat_elem_lid);
                elem_fields(elem_gid, 1) = State.MaterialPoints(mat_id).pres.host(mat_elem_lid);
                elem_fields(elem_gid, 2) = State.MaterialPoints(mat_id).sie.host(1, mat_elem_lid);
                // 3 is guass point vol
                elem_fields(elem_gid, 4) = State.MaterialPoints(mat_id).mass.host(mat_elem_lid);
                elem_fields(elem_gid, 5) = State.MaterialPoints(mat_id).sspd.host(mat_elem_lid);
                // 6 is elem speed
                elem_fields(elem_gid, 7) = (double)mat_id;
                // 8 is the e_switch
                elem_fields(elem_gid, 9) = (double)State.MaterialPoints(mat_id).eroded.host(mat_elem_lid);
            } // end for mat elems storage
        } // end parallel loop over materials

        // export element centric data
        double e_switch = 1;
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            elem_fields(elem_gid, 3) = State.GaussPoints.vol.host(elem_gid);
            elem_fields(elem_gid, 6) = speed.host(elem_gid);
            elem_fields(elem_gid, 8) = e_switch;
            elem_switch *= -1;
        } // end for elem_gid

        // save the vertex vector fields to an array for exporting to graphics files
        CArray<double> vec_fields(num_nodes, num_vec_vars, 3);

        for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
            // position, var 0
            vec_fields(node_gid, 0, 0) = State.node.coords.host(1, node_gid, 0);
            vec_fields(node_gid, 0, 1) = State.node.coords.host(1, node_gid, 1);
            if (num_dims == 2) {
                vec_fields(node_gid, 0, 2) = 0.0;
            }
            else{
                vec_fields(node_gid, 0, 2) = State.node.coords.host(1, node_gid, 2);
            }

            // velocity, var 1
            vec_fields(node_gid, 1, 0) = State.node.vel.host(1, node_gid, 0);
            vec_fields(node_gid, 1, 1) = State.node.vel.host(1, node_gid, 1);
            if (num_dims == 2) {
                vec_fields(node_gid, 1, 2) = 0.0;
            }
            else{
                vec_fields(node_gid, 1, 2) = State.node.vel.host(1, node_gid, 2);
            }

            // accelleration, var 2
            vec_fields(node_gid, 2, 0) = (State.node.vel.host(1, node_gid, 0) - State.node.vel.host(0, node_gid, 0))/dt;
            vec_fields(node_gid, 2, 1) = (State.node.vel.host(1, node_gid, 1) - State.node.vel.host(0, node_gid, 1))/dt;
            if (num_dims == 2) {
                vec_fields(node_gid, 2, 2) = 0.0;
            }
            else{
                vec_fields(node_gid, 2, 2) = (State.node.vel.host(1, node_gid, 2) - State.node.vel.host(0, node_gid, 2))/dt;
            }


        } // end for loop over vertices

        //  ---------------------------------------------------------------------------
        //  Setup of file and directoring for exporting
        //  ---------------------------------------------------------------------------
        FILE* out[20];   // the output files that are written to
        char  filename[128];
        int   max_len = sizeof filename;
        int   str_output_len;

        struct stat st;

        if (stat("ensight", &st) != 0) {
            system("mkdir ensight");
        }

        if (stat("ensight/data", &st) != 0) {
            system("mkdir ensight/data");
        }

        //  ---------------------------------------------------------------------------
        //  Write the Geometry file
        //  ---------------------------------------------------------------------------
        // sprintf(filename, "ensight/data/%s.%05d.geo", name, graphics_id);
        str_output_len = snprintf(filename, max_len, "ensight/data/%s.%05d.geo", name, graphics_id);
        // filename has the full string
        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

        out[0] = fopen(filename, "w");

        fprintf(out[0], "A graphics dump by Fierro \n");

        fprintf(out[0], "%s", "EnSight Gold geometry\n");
        fprintf(out[0], "%s", "node id assign\n");
        fprintf(out[0], "%s", "element id assign\n");

        fprintf(out[0], "part\n");
        fprintf(out[0], "%10d\n", 1);
        fprintf(out[0], "Mesh\n");

        // --- vertices ---
        fprintf(out[0], "coordinates\n");
        fprintf(out[0], "%10lu\n", num_nodes);

        // write all components of the point coordinates
        for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
            fprintf(out[0], "%12.5e\n", State.node.coords.host(1, node_gid, 0));
        }

        for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
            fprintf(out[0], "%12.5e\n", State.node.coords.host(1, node_gid, 1));
        }

        for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
            if (num_dims == 3) {
                fprintf(out[0], "%12.5e\n", State.node.coords.host(1, node_gid, 2));
            }
            else{
                fprintf(out[0], "%12.5e\n", 0.0);
            }
        }

        // --- elements ---
        if (num_dims == 3) {
            fprintf(out[0], "hexa8\n");
        }
        else{
            fprintf(out[0], "quad4\n");
        }
        fprintf(out[0], "%10lu\n", num_elems);


        int convert_ijk_to_ensight[8];
        if(mesh.num_dims==3){
            convert_ijk_to_ensight[0] = 0;
            convert_ijk_to_ensight[1] = 1;
            convert_ijk_to_ensight[2] = 3;
            convert_ijk_to_ensight[3] = 2;
            convert_ijk_to_ensight[4] = 4;
            convert_ijk_to_ensight[5] = 5;
            convert_ijk_to_ensight[6] = 7;
            convert_ijk_to_ensight[7] = 6;
        }
        else{
        
            convert_ijk_to_ensight[0] = 0;
            convert_ijk_to_ensight[1] = 1;
            convert_ijk_to_ensight[2] = 2;
            convert_ijk_to_ensight[3] = 3;
            convert_ijk_to_ensight[4] = 4;
            convert_ijk_to_ensight[5] = 5;
            convert_ijk_to_ensight[6] = 6;
            convert_ijk_to_ensight[7] = 7;
        } // end if


        // write all global point numbers for this cell
        for (int elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                fprintf(out[0], "%10lu\t", mesh.nodes_in_elem.host(elem_gid, convert_ijk_to_ensight[node_lid]) + 1); // note: node_gid starts at 1
            }
            fprintf(out[0], "\n");
        }

        fclose(out[0]);

        // ---------------------------------------------------------------------------
        // Write the Scalar variable files
        // ---------------------------------------------------------------------------

        // ensight_vars = (den, pres,...)
        for (int var = 0; var < num_scalar_vars; var++) {
            // write a scalar value
            // sprintf(filename, "ensight/data/%s.%05d.%s", name, graphics_id, scalar_var_names[var]);
            str_output_len = snprintf(filename, max_len, "ensight/data/%s.%05d.%s", name, graphics_id, scalar_var_names[var]);
            if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

            out[0] = fopen(filename, "w");

            fprintf(out[0], "Per_elem scalar values\n");
            fprintf(out[0], "part\n");
            fprintf(out[0], "%10d\n", 1);
            if (num_dims == 3) {
                fprintf(out[0], "hexa8\n");
            }
            else{
                fprintf(out[0], "quad4\n");
            }

            for (int elem_id = 0; elem_id < num_elems; elem_id++) {
                fprintf(out[0], "%12.5e\n", elem_fields(elem_id, var));
            }

            fclose(out[0]);
        } // end for var

        //  ---------------------------------------------------------------------------
        //  Write the Vector variable files
        //  ---------------------------------------------------------------------------

        // ensight vector vars = (position, velocity, force)
        for (int var = 0; var < num_vec_vars; var++) {
            // sprintf(filename, "ensight/data/%s.%05d.%s", name, graphics_id, vec_var_names[var]);
            str_output_len = snprintf(filename, max_len, "ensight/data/%s.%05d.%s", name, graphics_id, vec_var_names[var]);
            if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

            out[0] = fopen(filename, "w");
            // fprintf(out[0],"Per_node vector values\n");
            // fprintf(out[0],"part\n");
            // fprintf(out[0],"%10d \n",1);
            // fprintf(out[0],"hexa8\n"); // WARNING, maybe bug here?

            fprintf(out[0], "Per_node vector values\n");
            fprintf(out[0], "part\n");
            fprintf(out[0], "%10d\n", 1);
            fprintf(out[0], "block\n");

            for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
                fprintf(out[0], "%12.5e\n", vec_fields(node_gid, var, 0));
            }

            for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
                fprintf(out[0], "%12.5e\n", vec_fields(node_gid, var, 1));
            }

            for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
                fprintf(out[0], "%12.5e\n", vec_fields(node_gid, var, 2));
            }

            fclose(out[0]);
        } // end for var

        // ---------------------------------------------------------------------------
        // Write the case file
        // ---------------------------------------------------------------------------

        // sprintf(filename, "ensight/%s.case", name);
        str_output_len = snprintf(filename, max_len, "ensight/%s.case", name);
        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

        out[0] = fopen(filename, "w");

        fprintf(out[0], "FORMAT\n");
        fprintf(out[0], "type: ensight gold\n");
        fprintf(out[0], "GEOMETRY\n");

        // sprintf(filename, "model: data/%s.*****.geo\n", name);
        str_output_len = snprintf(filename, max_len, "model: data/%s.*****.geo\n", name);
        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

        fprintf(out[0], "%s", filename);
        fprintf(out[0], "VARIABLE\n");

        for (int var = 0; var < num_scalar_vars; var++) {
            // sprintf(filename, "scalar per element: %s data/%s.*****.%s\n", scalar_var_names[var], name, scalar_var_names[var]);
            str_output_len = snprintf(filename, max_len, "scalar per element: %s data/%s.*****.%s\n", scalar_var_names[var], name, scalar_var_names[var]);
            if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }

            fprintf(out[0], "%s", filename);
        }

        for (int var = 0; var < num_vec_vars; var++) {
            // sprintf(filename, "vector per node: %s data/%s.*****.%s\n", vec_var_names[var], name, vec_var_names[var]);
            str_output_len = snprintf(filename, max_len, "vector per node: %s data/%s.*****.%s\n", vec_var_names[var], name, vec_var_names[var]);
            if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }
            fprintf(out[0], "%s", filename);
        }

        fprintf(out[0], "TIME\n");
        fprintf(out[0], "time set: 1\n");
        fprintf(out[0], "number of steps: %4d\n", graphics_id + 1);
        fprintf(out[0], "filename start number: 0\n");
        fprintf(out[0], "filename increment: 1\n");
        fprintf(out[0], "time values: \n");

        graphics_times(graphics_id) = time_value;

        for (int i = 0; i <= graphics_id; i++) {
            fprintf(out[0], "%12.5e\n", graphics_times(i));
        }
        fclose(out[0]);

        // ---------------------------------------------------------------------------
        // Done writing the graphics dump
        // ---------------------------------------------------------------------------

        // increment graphics id counter
        graphics_id++;

        delete[] name;

        return;
    }

    /**\brief Given (i,j,k) coordinates within the Lagrange hex, return an offset into the local connectivity (PointIds) array.
    *"den", "pres", "sie", "vol", "mass", "sspd", "speed", "mat_id", "elem_switch", "eroded"
    * The \a order parameter must point to an array of 3 integers specifying the order
    * along each axis of the hexahedron.
    */
    int PointIndexFromIJK(int i, int j, int k, const int* order)
    {
        bool ibdy = (i == 0 || i == order[0]);
        bool jbdy = (j == 0 || j == order[1]);
        bool kbdy = (k == 0 || k == order[2]);
        // How many boundaries do we lie on at once?
        int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0);

        if (nbdy == 3) { // Vertex DOF
            // ijk is a corner node. Return the proper index (somewhere in [0,7]):
            return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
        }

        int offset = 8;
        if (nbdy == 2) { // Edge DOF
            if (!ibdy) { // On i axis
                return (i - 1) + (j ? order[0] - 1 + order[1] - 1 : 0) + (k ? 2 * (order[0] - 1 + order[1] - 1) : 0) + offset;
            }
            if (!jbdy) { // On j axis
                return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) + (k ? 2 * (order[0] - 1 + order[1] - 1) : 0) + offset;
            }
            // !kbdy, On k axis
            offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
            return (k - 1) + (order[2] - 1) * (i ? (j ? 3 : 1) : (j ? 2 : 0)) + offset;
        }

        offset += 4 * (order[0] - 1 + order[1] - 1 + order[2] - 1);
        if (nbdy == 1) { // Face DOF
            if (ibdy) { // On i-normal face
                return (j - 1) + ((order[1] - 1) * (k - 1)) + (i ? (order[1] - 1) * (order[2] - 1) : 0) + offset;
            }
            offset += 2 * (order[1] - 1) * (order[2] - 1);
            if (jbdy) { // On j-normal face
                return (i - 1) + ((order[0] - 1) * (k - 1)) + (j ? (order[2] - 1) * (order[0] - 1) : 0) + offset;
            }
            offset += 2 * (order[2] - 1) * (order[0] - 1);
            // kbdy, On k-normal face
            return (i - 1) + ((order[0] - 1) * (j - 1)) + (k ? (order[0] - 1) * (order[1] - 1) : 0) + offset;
        }

        // nbdy == 0: Body DOF
        offset += 2 * ( (order[1] - 1) * (order[2] - 1) + (order[2] - 1) * (order[0] - 1) + (order[0] - 1) * (order[1] - 1));
        return offset + (i - 1) + (order[0] - 1) * ( (j - 1) + (order[1] - 1) * ( (k - 1)));
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn write_vtk
    ///
    /// \brief Writes a vtk output file
    ///
    /// \param Simulation mesh
    /// \param State data
    /// \param Simulation parameters
    /// \param current time value
    /// \param Vector of all graphics output times
    ///
    /////////////////////////////////////////////////////////////////////////////
    void write_vtk(Mesh_t& mesh,
        State_t& State,
        SimulationParameters_t& SimulationParamaters,
        double time_value,
        CArray<double> graphics_times,
        std::vector<node_state> node_states,
        std::vector<gauss_pt_state> gauss_pt_states,
        std::vector<material_pt_state> material_pt_states)
    {

        size_t num_mats = State.MaterialPoints.size();

        // ---- Update host data ----

        // material point values
        for (int mat_id = 0; mat_id < num_mats; mat_id++) {
            State.MaterialPoints(mat_id).den.update_host();
            State.MaterialPoints(mat_id).pres.update_host();
            State.MaterialPoints(mat_id).stress.update_host();
            State.MaterialPoints(mat_id).sspd.update_host();
            State.MaterialPoints(mat_id).sie.update_host();
            State.MaterialPoints(mat_id).mass.update_host();
            State.MaterialPoints(mat_id).eroded.update_host();
        } // end for mat_id

        // gauss point values
        State.GaussPoints.vol.update_host();

        // nodal values
        State.node.coords.update_host();
        State.node.vel.update_host();
        State.node.mass.update_host();
        State.node.temp.update_host();

        Kokkos::fence();


        const int num_cell_scalar_vars = 10;
        const int num_cell_vec_vars    = 0;

        const int num_point_scalar_vars = 1;
        const int num_point_vec_vars = 2;


        // Scalar values associated with a cell
        const char cell_scalar_var_names[num_cell_scalar_vars][15] = {
            "den", "pres", "sie", "vol", "mass", "sspd", "speed", "mat_id", "elem_switch", "eroded"
        };
        
        const char cell_vec_var_names[num_cell_vec_vars][15] = {
            
        };

        const char point_scalar_var_names[num_point_scalar_vars][15] = {
            "temp"
        };

        const char point_vec_var_names[num_point_vec_vars][15] = {
            "pos", "vel" 
        };

        // short hand
        const size_t num_nodes = mesh.num_nodes;
        const size_t num_elems = mesh.num_elems;
        const size_t num_dims  = mesh.num_dims;

        // save the cell state to an array for exporting to graphics files
        auto elem_fields = CArray<double>(num_elems, num_cell_scalar_vars);
        int  elem_switch = 1;

        DCArrayKokkos<double> speed(num_elems);
        FOR_ALL(elem_gid, 0, num_elems, {
            double elem_vel[3]; // note:initialization with a list won't work
            elem_vel[0] = 0.0;
            elem_vel[1] = 0.0;
            elem_vel[2] = 0.0;
            // get the coordinates of the element center
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                elem_vel[0] += State.node.vel(1, mesh.nodes_in_elem(elem_gid, node_lid), 0);
                elem_vel[1] += State.node.vel(1, mesh.nodes_in_elem(elem_gid, node_lid), 1);
                if (mesh.num_dims == 3) {
                    elem_vel[2] += State.node.vel(1, mesh.nodes_in_elem(elem_gid, node_lid), 2);
                }
                else{
                    elem_vel[2] = 0.0;
                }
            } // end loop over nodes in element
            elem_vel[0] = elem_vel[0] / mesh.num_nodes_in_elem;
            elem_vel[1] = elem_vel[1] / mesh.num_nodes_in_elem;
            elem_vel[2] = elem_vel[2] / mesh.num_nodes_in_elem;

            double speed_sqrd = 0.0;
            for (int dim = 0; dim < num_dims; dim++) {
                speed_sqrd += elem_vel[dim] * elem_vel[dim];
            }
            speed(elem_gid) = sqrt(speed_sqrd);
        }); // end parallel for
        speed.update_host();
        Kokkos::fence();

        // save the output scale fields to a single 2D array

        // export material centeric data to the elements
        for (int mat_id = 0; mat_id < num_mats; mat_id++) {
            size_t num_mat_elems = State.MaterialToMeshMaps(mat_id).num_material_elems;

            for (size_t mat_elem_lid = 0; mat_elem_lid < num_mat_elems; mat_elem_lid++) {
                // 1 material per element

                // get elem gid
                size_t elem_gid = State.MaterialToMeshMaps(mat_id).elem.host(mat_elem_lid);

                // save outputs
                elem_fields(elem_gid, 0) = State.MaterialPoints(mat_id).den.host(mat_elem_lid);
                elem_fields(elem_gid, 1) = State.MaterialPoints(mat_id).pres.host(mat_elem_lid);
                elem_fields(elem_gid, 2) = State.MaterialPoints(mat_id).sie.host(1, mat_elem_lid);
                // 3 is guass point vol
                elem_fields(elem_gid, 4) = State.MaterialPoints(mat_id).mass.host(mat_elem_lid);
                elem_fields(elem_gid, 5) = State.MaterialPoints(mat_id).sspd.host(mat_elem_lid);
                // 6 is elem speed
                elem_fields(elem_gid, 7) = (double)mat_id;
                // 8 is the e_switch
                elem_fields(elem_gid, 9) = (double)State.MaterialPoints(mat_id).eroded.host(mat_elem_lid);
            } // end for mat elems storage
        } // end parallel loop over materials

        // export element centric data
        double e_switch = 1;
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            elem_fields(elem_gid, 3) = State.GaussPoints.vol.host(elem_gid);
            elem_fields(elem_gid, 6) = speed.host(elem_gid);
            elem_fields(elem_gid, 8) = e_switch;
            elem_switch *= -1;
        } // end for elem_gid

        // save the vertex vector fields to an array for exporting to graphics files
        CArray<double> vec_fields(num_nodes, num_point_vec_vars, 3);
        CArray<double> point_scalar_fields(num_nodes, num_point_scalar_vars);

        for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
            // position, var 0
            vec_fields(node_gid, 0, 0) = State.node.coords.host(1, node_gid, 0);
            vec_fields(node_gid, 0, 1) = State.node.coords.host(1, node_gid, 1);
            if (num_dims == 2) {
                vec_fields(node_gid, 0, 2) = 0.0;
            }
            else{
                vec_fields(node_gid, 0, 2) = State.node.coords.host(1, node_gid, 2);
            }

            // position, var 1
            vec_fields(node_gid, 1, 0) = State.node.vel.host(1, node_gid, 0);
            vec_fields(node_gid, 1, 1) = State.node.vel.host(1, node_gid, 1);
            if (num_dims == 2) {
                vec_fields(node_gid, 1, 2) = 0.0;
            }
            else{
                vec_fields(node_gid, 1, 2) = State.node.vel.host(1, node_gid, 2);
            }

            point_scalar_fields(node_gid, 0) = State.node.temp.host(1,node_gid);
        } // end for loop over vertices

        FILE* out[20];   // the output files that are written to
        char  filename[100]; // char string
        int   max_len = sizeof filename;
        int   str_output_len;

        struct stat st;

        if (stat("vtk", &st) != 0) {
            system("mkdir vtk");
        }

        // snprintf(filename, max_len, "ensight/data/%s.%05d.%s", name, graphics_id, vec_var_names[var]);

        //sprintf(filename, "vtk/meshHexPn.%05d.vtk", graphics_id);  // mesh file
        str_output_len = snprintf(filename, max_len, "vtk/meshHexPn.%05d.vtk", graphics_id);
        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }
         // mesh file
        
        out[0] = fopen(filename, "w");

        fprintf(out[0], "# vtk DataFile Version 2.0\n");  // part 2
        fprintf(out[0], "Mesh for Fierro\n");             // part 2
        fprintf(out[0], "ASCII \n");                      // part 3
        fprintf(out[0], "DATASET UNSTRUCTURED_GRID\n\n"); // part 4

        fprintf(out[0], "POINTS %zu float\n", mesh.num_nodes);

        // write all components of the point coordinates
        for (size_t node_gid = 0; node_gid < mesh.num_nodes; node_gid++) {
            fprintf(out[0],
                    "%f %f %f\n",
                    State.node.coords.host(1, node_gid, 0),
                    State.node.coords.host(1, node_gid, 1),
                    State.node.coords.host(1, node_gid, 2));
        } // end for

        /*
        ---------------------------------------------------------------------------
        Write the elems
        ---------------------------------------------------------------------------
        */

        fprintf(out[0], "\n");
        fprintf(out[0], "CELLS %lu %lu\n", mesh.num_elems, mesh.num_elems + mesh.num_elems * mesh.num_nodes_in_elem);  // size=all printed values

        int Pn_order   = mesh.Pn;
        int order[3]   = { Pn_order, Pn_order, Pn_order };

        // const int num_1D_points = Pn_order+1;

        // write all global point numbers for this elem
        for (size_t elem_gid = 0; elem_gid < mesh.num_elems; elem_gid++) {
            fprintf(out[0], "%lu ", mesh.num_nodes_in_elem); // num points in this elem

            for (int k = 0; k <= Pn_order; k++) {
                for (int j = 0; j <= Pn_order; j++) {
                    for (int i = 0; i <= Pn_order; i++) {
                        size_t node_lid = PointIndexFromIJK(i, j, k, order);
                        fprintf(out[0], "%lu ", mesh.nodes_in_elem.host(elem_gid, node_lid));
                    }
                }
            }

            fprintf(out[0], "\n");
        } // end for

        // Write the element types
        fprintf(out[0], "\n");
        fprintf(out[0], "CELL_TYPES %zu \n", mesh.num_elems);
        // VTK_LAGRANGE_HEXAHEDRON: 72,
        // VTK_HIGHER_ORDER_HEXAHEDRON: 67
        // VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33
        // element types: https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
        // element types: https://kitware.github.io/vtk-js/api/Common_DataModel_CellTypes.html
        // vtk format: https://www.kitware.com//modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/
        for (size_t elem_gid = 0; elem_gid < mesh.num_elems; elem_gid++) {
            fprintf(out[0], "%d \n", 72);
        }

        /*
        ---------------------------------------------------------------------------
        Write the nodal vector variables to file
        ---------------------------------------------------------------------------
        */

        fprintf(out[0], "\n");
        fprintf(out[0], "POINT_DATA %zu \n", mesh.num_nodes);

        // vtk vector vars = (position, velocity)
        for (int var = 0; var < num_point_vec_vars; var++) {
            fprintf(out[0], "VECTORS %s float \n", point_vec_var_names[var]);
            for (size_t node_gid = 0; node_gid < mesh.num_nodes; node_gid++) {
                fprintf(out[0], "%f %f %f\n",
                        vec_fields(node_gid, var, 0),
                        vec_fields(node_gid, var, 1),
                        vec_fields(node_gid, var, 2));
            } // end for nodes
        } // end for vec_vars


        // vtk scalar vars = (temp)
        for (int var = 0; var < num_point_scalar_vars; var++) {
            fprintf(out[0], "SCALARS %s float 1\n", point_scalar_var_names[var]);
            fprintf(out[0], "LOOKUP_TABLE default\n");
            for (size_t node_gid = 0; node_gid < mesh.num_nodes; node_gid++) {
                fprintf(out[0], "%f\n",
                        point_scalar_fields(node_gid, 0));
            } // end for nodes
        } // end for vec_vars

        /*
        ---------------------------------------------------------------------------
        Write the scalar elem variable to file
        ---------------------------------------------------------------------------
        */
        fprintf(out[0], "\n");
        fprintf(out[0], "CELL_DATA %zu \n", mesh.num_elems);

        for (int var = 0; var < num_cell_scalar_vars; var++) {
            fprintf(out[0], "SCALARS %s float 1\n", cell_scalar_var_names[var]); // the 1 is number of scalar components [1:4]
            fprintf(out[0], "LOOKUP_TABLE default\n");
            for (size_t elem_gid = 0; elem_gid < mesh.num_elems; elem_gid++) {
                fprintf(out[0], "%f\n", elem_fields(elem_gid, var));
            } // end for elem
        } // end for cell scalar_vars

        fclose(out[0]);

        graphics_times(graphics_id) = time_value;

        // Write time series metadata
        //sprintf(filename, "vtk/meshHexPn.vtk.series", graphics_id);  // mesh file
        str_output_len = snprintf(filename, max_len, "vtk/meshHexPn.vtk.series"); 
        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }
        // mesh file

        out[0] = fopen(filename, "w");

        fprintf(out[0], "{\n");
        fprintf(out[0], "  \"file-series-version\" : \"1.0\",\n");
        fprintf(out[0], "  \"files\" : [\n");

        for (int i = 0; i <= graphics_id; i++) {
            fprintf(out[0], "    { \"name\" : \"meshHexPn.%05d.vtk\", \"time\" : %12.5e },\n", i, graphics_times(i) );
        }

        // fprintf(out[0], "%12.5e\n", graphics_times(i));
        fprintf(out[0], "  ]\n"); // part 4
        fprintf(out[0], "}"); // part 4

        fclose(out[0]);

        // increment graphics id counter
        graphics_id++;
    } // end write vtk

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn write_material_point_state
    ///
    /// \brief Writes a state output file at each material point
    ///
    /// \param Simulation mesh
    /// \param State data
    /// \param Simulation parameters
    /// \param current time value
    /// \param Vector of all graphics output times
    ///
    /////////////////////////////////////////////////////////////////////////////
    void write_material_point_state(Mesh_t& mesh,
        State_t& State,
        SimulationParameters_t& SimulationParamaters,
        double time_value,
        CArray<double> graphics_times,
        std::vector<node_state> node_states,
        std::vector<gauss_pt_state> gauss_pt_states,
        std::vector<material_pt_state> material_pt_states)
    {
        // WARNING WARNING WARNING:
        // This currently assumes the gauss and material point IDs are the same as the element ID
        // This will need to be updated for high order methods

        // Update host data
        // ---- Update host data ----
        size_t num_mats = State.MaterialPoints.size();

        for (int mat_id = 0; mat_id < num_mats; mat_id++) {
            State.MaterialPoints(mat_id).den.update_host();
            State.MaterialPoints(mat_id).pres.update_host();
            State.MaterialPoints(mat_id).stress.update_host();
            State.MaterialPoints(mat_id).sspd.update_host();
            State.MaterialPoints(mat_id).sie.update_host();
            State.MaterialPoints(mat_id).mass.update_host();
        } // end for mat_id

        State.GaussPoints.vol.update_host();

        State.node.coords.update_host();
        State.node.vel.update_host();
        State.node.mass.update_host();

        Kokkos::fence();

        struct stat st;

        if (stat("state", &st) != 0) {
            system("mkdir state");
        }

        size_t num_dims = mesh.num_dims;

        //  ---------------------------------------------------------------------------
        //  Setup of file and directory for exporting
        //  ---------------------------------------------------------------------------

        // output file
        FILE* out_elem_state;  // element average state
        char  filename[128];

        int max_len = sizeof filename;

        snprintf(filename, max_len, "state/mat_pt_state_t_%6.4e.txt", time_value);

        // output files
        out_elem_state = fopen(filename, "w");

        // write state dump
        fprintf(out_elem_state, "# state dump file\n");
        fprintf(out_elem_state, "# x  y  z  radius_2D  radius_3D  den  pres  sie  sspd  vol  mass \n");

        // write out values for the elem
        for (size_t mat_id = 0; mat_id < num_mats; mat_id++) {
            size_t num_mat_elems = State.MaterialToMeshMaps(mat_id).num_material_elems;
            for (size_t elem_lid = 0; elem_lid < num_mat_elems; elem_lid++)
            {

                size_t elem_gid = State.MaterialToMeshMaps(mat_id).elem.host(elem_lid);


                double elem_coords[3];
                elem_coords[0] = 0.0;
                elem_coords[1] = 0.0;
                elem_coords[2] = 0.0;

                // get the coordinates of the element center
                for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                    elem_coords[0] += State.node.coords.host(1, mesh.nodes_in_elem.host(elem_gid, node_lid), 0);
                    elem_coords[1] += State.node.coords.host(1, mesh.nodes_in_elem.host(elem_gid, node_lid), 1);
                    if (num_dims == 3) {
                        elem_coords[2] += State.node.coords.host(1, mesh.nodes_in_elem.host(elem_gid, node_lid), 2);
                    }
                    else{
                        elem_coords[2] = 0.0;
                    }
                } // end loop over nodes in element

                elem_coords[0] = elem_coords[0] / mesh.num_nodes_in_elem;
                elem_coords[1] = elem_coords[1] / mesh.num_nodes_in_elem;
                elem_coords[2] = elem_coords[2] / mesh.num_nodes_in_elem;

                double rad2 = sqrt(elem_coords[0] * elem_coords[0] +
                                   elem_coords[1] * elem_coords[1]);

                double rad3 = sqrt(elem_coords[0] * elem_coords[0] +
                                   elem_coords[1] * elem_coords[1] +
                                   elem_coords[2] * elem_coords[2]);

                fprintf(out_elem_state, "%4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t \n",
                         elem_coords[0],
                         elem_coords[1],
                         elem_coords[2],
                         rad2,
                         rad3,
                         State.MaterialPoints(mat_id).den.host(elem_gid),
                         State.MaterialPoints(mat_id).pres.host(elem_gid),
                         State.MaterialPoints(mat_id).sie.host(1, elem_gid),
                         State.MaterialPoints(mat_id).sspd.host(elem_gid),
                         State.GaussPoints.vol.host(elem_gid),
                         State.MaterialPoints(mat_id).mass.host(elem_gid) );
            } // end for elements
        } // end for materials
        fclose(out_elem_state);


        // printing nodal state
            
        FILE* out_point_state;  // element average state

        snprintf(filename, max_len, "state/node_state_t_%6.4e.txt", time_value);

        // output files
        out_point_state = fopen(filename, "w");

        // write state dump
        fprintf(out_point_state, "# state node dump file\n");
        fprintf(out_point_state, "# x  y  z  radius_2D  radius_3D  vel_x  vel_y  vel_z  speed  ||err_v_dot_r|| \n");

        // get the coordinates of the node
        for (size_t node_gid = 0; node_gid < mesh.num_nodes; node_gid++) {

            double node_coords[3];

            node_coords[0] = State.node.coords.host(1, node_gid, 0);
            node_coords[1] = State.node.coords.host(1, node_gid, 1);
            if (num_dims == 3) {
                node_coords[2] = State.node.coords.host(1, node_gid, 2);
            }
            else{
                node_coords[2] = 0.0;
            }

            double rad2 = sqrt(node_coords[0] * node_coords[0] +
                               node_coords[1] * node_coords[1]);
            double rad3 = sqrt(node_coords[0] * node_coords[0] +
                               node_coords[1] * node_coords[1] +
                               node_coords[2] * node_coords[2]);

            double node_vel[3];

           node_vel[0] = State.node.vel.host(1, node_gid, 0);
           node_vel[1] = State.node.vel.host(1, node_gid, 1);
            if (num_dims == 3) {
               node_vel[2] = State.node.vel.host(1, node_gid, 2);
            }
            else{
               node_vel[2] = 0.0;
            }

            double speed = sqrt(node_vel[0] * node_vel[0] +
                                node_vel[1] * node_vel[1] +
                                node_vel[2] * node_vel[2]);



            // looking at perfect radial motion
            double unit_r_vec[2];
            unit_r_vec[0] = node_coords[0]/rad2;
            unit_r_vec[1] = node_coords[1]/rad2;

            //the radial motion
            double v_dot_r = node_vel[0] * unit_r_vec[0] +
                             node_vel[1] * unit_r_vec[1];
            

            double err_v_dot_r[3]; 
            err_v_dot_r[0] = node_vel[0]-unit_r_vec[0]*v_dot_r;
            err_v_dot_r[1] = node_vel[1]-unit_r_vec[1]*v_dot_r;

            double mag_err_v_dot_r = sqrt(err_v_dot_r[0]*err_v_dot_r[0] + err_v_dot_r[1]*err_v_dot_r[1]);

            fprintf(out_point_state, "%4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t %4.12e\t  %4.12e\t %4.12e\t \n",
                         node_coords[0],
                         node_coords[1],
                         node_coords[2],
                         rad2,
                         rad3,
                         node_vel[0],
                         node_vel[1],
                         node_vel[2],
                         speed,
                         mag_err_v_dot_r);


        } // end loop over nodes in element


        fclose(out_point_state);


        return;
    } // end of state output
}; // end class

#endif // end Header Guard