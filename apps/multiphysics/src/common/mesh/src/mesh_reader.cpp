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
#include "mesh_reader.hpp"

#include "region.hpp"
#include "string_utils.hpp"

#include <map>
#include <memory>
#include <cstring>
#include <sys/stat.h>
#include <iostream>
#include <regex>
#include <fstream>
#include <sstream>
#include <mpi.h>

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

//-------
// word is the field name e.g., Offsets, connectivity, etc.
// stop is the phrase to stop extracting values
template <typename T>
static bool extract_values_xml(T *values_xml,
                        const std::string& word,
                        const std::string& stop,
                        std::ifstream& in,
                        size_t& size)
{

        bool found = false;

        std::string line;

        size_t i = 0;

        // Read the file line by line looking for specified word
        while (std::getline(in, line)) {

            if (line.find(word) != std::string::npos) { // Check if the portion of the word is in the line
                found = true;
            } 
            if(found) {

                // loop over the lines in the file, extracting the values of the field corresponding to the word
                while (std::getline(in, line)){  
                
                    std::istringstream iss(line);  // Create a stream from the line

                    // extract the individual values from the stream
                    T value;
                    while (iss >> value) {
                        values_xml[i] = value;
                        i++;
                    } // end while

                    if (line.find(stop) != std::string::npos) { // Check if the stop word is in the line
                        break;
                    } // end if

                } // end while

                if(found) break;

            } // end if found

        } // end while

        size = i;

        return found;

} // end function

// find the number of points and number of cells in the mesh
static bool extract_num_points_and_cells_xml(int& numberOfPoints,
                                      int& numberOfCells,
                                      std::ifstream& in)
{
    bool found = false;

    std::string line;

        
    // Read the file line by line looking for NumberOfPoints
    while (std::getline(in, line)) {
        
        std::string word = "NumberOfPoints=";  // A portion of a word

        if (line.find(word) != std::string::npos) { // Check if the portion of the word is in the line
            found = true;
        }
        if(found) {
            // Define regex pattern to match the attributes and capture values
            std::regex pattern(R"(NumberOfPoints=\"(\d+)\" NumberOfCells=\"(\d+)\")");
            std::smatch match;

            if (std::regex_search(line, match, pattern)) {
                //std::cout << "Number of nodes in mesh file: " << match[1] << std::endl;
                //std::cout << "Number of cells in mesh file: " << match[2] << std::endl;

                numberOfPoints = std::stoi(match[1].str());
                numberOfCells = std::stoi(match[2].str());

            } else {
                std::cout << "Error reading the number of points and cells in the mesh!" << std::endl;
            }
            
            break;
        } // end if
        
    } // end while

    return found;

} // end function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn set_mesh_file
///
/// \brief Sets the mesh file path for reading in a mesh
///
/// \param Path to mesh file
///
/////////////////////////////////////////////////////////////////////////////
void MeshReader::set_mesh_file(char* MESH)
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
///
///
/////////////////////////////////////////////////////////////////////////////
void MeshReader::read_mesh(swage::Mesh_t& mesh,
               MPICArrayKokkos<real_t>& node_coords,
               MeshInput_t& mesh_inps,
               int           num_dims,
               bool HexN)
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
        read_ensight_mesh(mesh, node_coords, mesh_inps, num_dims);
    }
    else if(extension == "inp"){ // Abaqus meshfile extension
        read_Abaqus_mesh(mesh, node_coords, mesh_inps, num_dims);
    }
    else if(extension == "vtk"){ // vtk file format
        if (!HexN) {
            read_vtk_mesh(mesh, node_coords, mesh_inps, num_dims);
        }
        else {
            read_vtk_hexN_mesh(mesh, node_coords, mesh_inps, num_dims);
        }
    }
    else if(extension == "vtu"){ // vtu file format
        read_vtu_mesh(mesh, node_coords, mesh_inps, num_dims, HexN);
    }
    else{
        throw std::runtime_error("**** Mesh file extension not understood ****");
    }

    // Initialize/create connectivity

    int num_corners = mesh.num_elems * mesh.num_nodes_in_elem;
    mesh.build_connectivity();
    // corner.initialize(num_corners, num_dims);


    // initialize node state variables, for now, we just need coordinates, the rest will be initialize by the respective solvers
    // std::vector<node_state> required_node_state = { node_state::coords };
    // node.initialize(mesh.num_nodes, num_dims, required_node_state);

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
///
/////////////////////////////////////////////////////////////////////////////
void MeshReader::read_ensight_mesh(swage::Mesh_t& mesh,
                       MPICArrayKokkos<real_t>& node_coords,
                       MeshInput_t& mesh_inps,
                       int num_dims)
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
    printf("Number of nodes read in %lu\n", num_nodes);

    // initialize node variables
    mesh.initialize_nodes(num_nodes);

    node_coords = MPICArrayKokkos<real_t>(num_nodes, num_dims, "Node_coordinates_in_mesh_io");

    // read the initial mesh coordinates
    // x-coords
    for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
        // Read into a double and convert on assignment: fscanf writes
        // sizeof(double) bytes through the pointer, so it must never be
        // handed the address of an array slot whose element type may be
        // narrower than double.
        real_t coord = real_t(0.0);
        fscanf(in, "%le", &coord);
        node_coords.host(node_id, 0) = coord * mesh_inps.scale_x;
    }

    // y-coords
    for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
        real_t coord = real_t(0.0);
        fscanf(in, "%le", &coord);
        node_coords.host(node_id, 1) = coord * mesh_inps.scale_y;
    }

    // z-coords
    for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
        if (num_dims == 3) {
            real_t coord = real_t(0.0);
            fscanf(in, "%le", &coord);
            node_coords.host(node_id, 2) = coord * mesh_inps.scale_z;
        }
        else{
            real_t dummy;
            fscanf(in, "%le", &dummy);
        }
    } // end for


    // Update device nodal positions
    node_coords.update_device();

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
    mesh.initialize_elems(num_elem);
    // GaussPoints.initialize(num_elem, 3); // always 3D here, even for 2D

    
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

    // Close mesh input file
    fclose(in);

    return;
} // end read ensight mesh


/////////////////////////////////////////////////////////////////////////////
///
/// \fn read_Abaqus_mesh
///
/// \brief Read .inp mesh file
///
/// \param Simulation mesh
/// \param Node coordinates
/// \param Number of dimensions
///
/////////////////////////////////////////////////////////////////////////////
void MeshReader::read_Abaqus_mesh(swage::Mesh_t& mesh,
                      MPICArrayKokkos<real_t>& node_coords,
                      MeshInput_t& mesh_inps,
                      int num_dims)
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

    node_coords = MPICArrayKokkos<real_t>(num_nodes, num_dims, "Node_coordinates_in_mesh_io");

    printf("Number of nodes read in %lu\n", num_nodes);

    // initialize node variables
    mesh.initialize_nodes(num_nodes);

    // Copy nodes to mesh
    for(int node_gid = 0; node_gid < num_nodes; node_gid++){
        node_coords.host(node_gid, 0) = nodes[node_gid].x * mesh_inps.scale_x;
        node_coords.host(node_gid, 1) = nodes[node_gid].y * mesh_inps.scale_y;
        node_coords.host(node_gid, 2) = nodes[node_gid].z * mesh_inps.scale_z;
    }

    // Update device nodal positions
    node_coords.update_device();


    // --- read in the elements in the mesh ---
    size_t num_elem = elements.size();
    printf("Number of elements read in %lu\n", num_elem);

    // initialize elem variables
    mesh.initialize_elems(num_elem);


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
///
/////////////////////////////////////////////////////////////////////////////
void MeshReader::read_vtk_mesh(swage::Mesh_t& mesh,
                   MPICArrayKokkos<real_t>& node_coords,
                   MeshInput_t& mesh_inps,
                   int num_dims)
{

    std::cout<<"Reading VTK mesh"<<std::endl;

    int i;           // used for writing information to file
    int node_gid = 0;    // the global id for the point
    int elem_gid = 0;     // the global id for the elem

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
            printf("Number of nodes read in %zu\n", num_nodes);
            mesh.initialize_nodes(num_nodes);

            // std::vector<node_state> required_node_state = { node_state::coords };
            // node.initialize(num_nodes, num_dims, required_node_state);
            
            found=true;
        } // end if
        
        
        if (i>1000){
            std::cerr << "ERROR: Failed to find POINTS in file" << std::endl;
            break;
        } // end if
        
        i++;
    } // end while

    node_coords = MPICArrayKokkos<real_t>(mesh.num_nodes, num_dims, "Node_coordinates_in_mesh_io");
    
    // read the node coordinates
    for (node_gid = 0; node_gid < mesh.num_nodes; node_gid++){
        
        std::string str;
        std::getline(in, str);
        
        std::string delimiter = " ";
        std::vector<std::string> v = split (str, delimiter);
        
        // save the nodal coordinates
        node_coords.host(node_gid, 0) = mesh_inps.scale_x*std::stod(v[0]); // double
        node_coords.host(node_gid, 1) = mesh_inps.scale_y*std::stod(v[1]); // double
        if(num_dims == 3){
            node_coords.host(node_gid, 2) = mesh_inps.scale_z*std::stod(v[2]); // double
        }
        
    } // end for nodes


    // Update device nodal positions
    node_coords.update_device();
    

    found = false;

    // look for CELLS
    i = 0;
    size_t num_elem = 0;
    while (found == false) {
        std::string str;
        std::getline(in, str);
        
        std::string delimiter = " ";
        std::vector<std::string> v = split (str, delimiter);
        std::cout << v[0] << std::endl; // printing
        
        // looking for the following text:
        //      CELLS num_elem size
        if(v[0] == "CELLS"){
            num_elem = std::stoi(v[1]);
            printf("Number of elements read in %zu\n", num_elem);

            // initialize elem variables
            mesh.initialize_elems(num_elem);
            
            found=true;
        } // end if
        
        
        if (i>10000){
            printf("ERROR: Failed to find CELLS \n");
            break;
        } // end if
        
        i++;
    } // end while
    
    
    // read the node ids in the element
    for (elem_gid = 0; elem_gid < num_elem; elem_gid++) {
        
        std::string str;
        std::getline(in, str);
        
        std::string delimiter = " ";
        std::vector<std::string> v = split (str, delimiter);
        num_nodes_in_elem = std::stoi(v[0]);
        
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
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
    printf("Element type = %zu \n", elem_type);
    // elem types:
    // linear hex = 12, linear quad = 9
    found=false;
    
    
    if(num_nodes_in_elem == 8 && elem_type != 12) {
        printf("Wrong element type of %zu \n", elem_type);
        std::cerr << "ERROR: incorrect element type in VTK file" << std::endl;
    }
    
    in.close();
    
} // end of VTKread function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn read_vtk_hexN_mesh
///
/// \brief Read ASCII .vtk mesh file for arbitrary order hex elements (HexN)
///        using Matar CArrays for internal mapping.
///
/// \param Simulation mesh
/// \param Simulation state
/// \param Node state struct
/// \param Number of dimensions
///
/////////////////////////////////////////////////////////////////////////////
void MeshReader::read_vtk_hexN_mesh(swage::Mesh_t& mesh,
                        MPICArrayKokkos<real_t>& node_coords,
                        MeshInput_t& mesh_inps,
                        int num_dims)
{
    std::cout << "Reading VTK HexN mesh" << std::endl;

    if (num_dims != 3) {
        throw std::runtime_error("read_vtk_hexN_mesh: only 3D Lagrange hex meshes are supported");
    }

    std::ifstream in(mesh_file_);
    if (!in.is_open()) {
        throw std::runtime_error("read_vtk_hexN_mesh: could not open mesh file");
    }

    // ---- whitespace-robust tokenizer (handles tabs, \r, repeated spaces) ----
    auto tokenize = [](const std::string& line) {
        std::vector<std::string> tok;
        std::istringstream iss(line);
        std::string t;
        while (iss >> t) tok.push_back(t);
        return tok;
    };

    // read the next NON-EMPTY line; returns false at EOF
    auto next_tokens = [&](std::vector<std::string>& tok) -> bool {
        std::string line;
        while (std::getline(in, line)) {
            tok = tokenize(line);
            if (!tok.empty()) return true;
        }
        return false;
    };

    // ---- 0. Header / version check ----
    {
        std::string header;
        std::getline(in, header);
        const size_t p = header.find("Version");
        if (header.find("vtk DataFile") == std::string::npos || p == std::string::npos) {
            throw std::runtime_error("read_vtk_hexN_mesh: missing '# vtk DataFile Version' header");
        }
        const real_t version = std::stod(header.substr(p + 7));
        if (version >= real_t(5.0)) {
            // 5.x uses OFFSETS/CONNECTIVITY blocks and the newer Lagrange-hex edge ordering
            throw std::runtime_error("read_vtk_hexN_mesh: legacy VTK version 5.x is not supported; "
                                    "save the file as version 4.2 or lower");
        }
    }

    std::vector<std::string> tok;

    // ---- 1. POINTS ----
    bool found = false;
    while (next_tokens(tok)) {
        if (tok[0] == "POINTS") { found = true; break; }
    }
    if (!found || tok.size() < 2) {
        throw std::runtime_error("read_vtk_hexN_mesh: POINTS not found");
    }

    const size_t num_nodes = std::stoul(tok[1]);
    printf("Number of nodes read in %zu\n", num_nodes);
    mesh.initialize_nodes(num_nodes);

    // read 3*num_nodes values regardless of how they are split across lines
    std::vector<double> xyz;
    xyz.reserve(3 * num_nodes);
    while (xyz.size() < 3 * num_nodes && next_tokens(tok)) {
        for (const auto& t : tok) xyz.push_back(std::stod(t));
    }
    if (xyz.size() < 3 * num_nodes) {
        throw std::runtime_error("read_vtk_hexN_mesh: fewer point coordinates than POINTS count");
    }

    node_coords = MPICArrayKokkos<real_t>(num_nodes, num_dims, "Node_coordinates_in_mesh_io");
    for (size_t n = 0; n < num_nodes; n++) {
        node_coords.host(n, 0) = mesh_inps.scale_x * xyz[3 * n + 0];
        node_coords.host(n, 1) = mesh_inps.scale_y * xyz[3 * n + 1];
        node_coords.host(n, 2) = mesh_inps.scale_z * xyz[3 * n + 2];
    }
    node_coords.update_device();

    // ---- 2. CELLS ----
    found = false;
    while (next_tokens(tok)) {
        if (tok[0] == "CELLS") { found = true; break; }
    }
    if (!found || tok.size() < 3) {
        throw std::runtime_error("read_vtk_hexN_mesh: CELLS not found");
    }

    const size_t num_elem  = std::stoul(tok[1]);
    const size_t cell_size = std::stoul(tok[2]);   // total ints in the CELLS block
    printf("Number of elements read in %zu\n", num_elem);

    std::vector<size_t> cell_data;
    cell_data.reserve(cell_size);
    while (cell_data.size() < cell_size && next_tokens(tok)) {
        for (const auto& t : tok) cell_data.push_back(std::stoul(t));
    }
    if (cell_data.size() < cell_size || num_elem == 0) {
        throw std::runtime_error("read_vtk_hexN_mesh: CELLS block is truncated");
    }

    // ---- infer Pn from the file, and check it against the YAML ----
    const size_t n0     = cell_data[0];
    const size_t num_1D = (size_t)std::llround(std::cbrt((double)n0));
    if (num_1D * num_1D * num_1D != n0 || num_1D < 2) {
        throw std::runtime_error("read_vtk_hexN_mesh: cell node count is not (Pn+1)^3");
    }
    const size_t Pn_file = num_1D - 1;

    if ((size_t)mesh_inps.p_order != Pn_file) {
        printf("ERROR: VTK file has Pn = %zu but input file has p_order = %d\n",
            Pn_file, (int)mesh_inps.p_order);
        throw std::runtime_error("read_vtk_hexN_mesh: p_order in YAML does not match the mesh file");
    }

    mesh.initialize_elems_Pn(num_elem, Pn_file, 2 * Pn_file);

    // ---- 3. build the VTK -> ijk map ONCE ----
    const int Pn = (int)Pn_file;
    int order[3] = {Pn, Pn, Pn};
    std::vector<size_t> vtk_lid_of_ijk(mesh.num_nodes_in_elem);
    {
        size_t lid = 0;
        for (int k = 0; k <= Pn; k++)
        for (int j = 0; j <= Pn; j++)
        for (int i = 0; i <= Pn; i++) {
            vtk_lid_of_ijk[lid++] = (size_t)PointIndexFromIJK(i, j, k, order);
        }
    }

    // ---- fill nodes_in_elem with validation ----
    size_t off = 0;
    for (size_t e = 0; e < num_elem; e++) {
        if (off >= cell_data.size()) {
            throw std::runtime_error("read_vtk_hexN_mesh: ran out of CELLS data");
        }
        const size_t n = cell_data[off];
        if (n != mesh.num_nodes_in_elem) {
            throw std::runtime_error("read_vtk_hexN_mesh: mixed element orders are not supported");
        }
        if (off + n >= cell_data.size()) {
            throw std::runtime_error("read_vtk_hexN_mesh: CELLS line is truncated");
        }

        for (size_t lid = 0; lid < n; lid++) {
            const size_t gid = cell_data[off + 1 + vtk_lid_of_ijk[lid]];
            if (gid >= num_nodes) {
                throw std::runtime_error("read_vtk_hexN_mesh: node id out of range in CELLS");
            }
            mesh.nodes_in_elem.host(e, lid) = gid;
        }
        off += n + 1;
    }
    mesh.nodes_in_elem.update_device();

    // ---- 4. CELL_TYPES ----
    found = false;
    while (next_tokens(tok)) {
        if (tok[0] == "CELL_TYPES") { found = true; break; }
    }
    if (!found) {
        throw std::runtime_error("read_vtk_hexN_mesh: CELL_TYPES not found");
    }

    size_t num_types_read = 0;
    while (num_types_read < num_elem && next_tokens(tok)) {
        for (const auto& t : tok) {
            const int type = std::stoi(t);
            // 72 = VTK_LAGRANGE_HEXAHEDRON; 12 = VTK_HEXAHEDRON (same vertex order when Pn = 1)
            if (!(type == 72 || (type == 12 && Pn == 1))) {
                printf("ERROR: unsupported cell type %d for a HexN mesh\n", type);
                throw std::runtime_error("read_vtk_hexN_mesh: unsupported cell type");
            }
            num_types_read++;
        }
    }
    printf("Element type read = 72 (Lagrange hex, Pn = %d)\n", Pn);

    in.close();

    // NOTE: do NOT call mesh.build_connectivity() here; read_mesh() already does it.
} // end of read_vtk_hexN_mesh


/////////////////////////////////////////////////////////////////////////////
///
/// \fn read_vtu_mesh
///
/// \brief Read ASCII .vtu mesh file
///
/// \param mesh Simulation mesh
/// \param node_coords Node coordinates
/// \param mesh_inps Mesh input parameters
/// \param num_dims Number of dimensions
///
/////////////////////////////////////////////////////////////////////////////
void MeshReader::read_vtu_mesh(swage::Mesh_t& mesh,
                   MPICArrayKokkos<real_t>& node_coords,
                   MeshInput_t& mesh_inps,
                   int num_dims,
                   bool HexN)
{

    std::cout<<"Reading VTU file in a multiblock VTK mesh"<<std::endl;

    // int i;           // used for writing information to file
    int node_gid = 0;    // the global id for the point
    int elem_gid = 0;    // the global id for the elem


    //
    int Pn_order = mesh_inps.p_order;
    size_t num_nodes_in_elem = 1;
    for (int dim = 0; dim < num_dims; dim++) {
        num_nodes_in_elem *= (Pn_order + 1);
    }
    
    bool found;
    
    std::ifstream in;  // FILE *in;
    in.open(mesh_file_);
    

    // --- extract the number of points and cells from the XML file ---
    int num_nodes;
    int num_elems;
    found = extract_num_points_and_cells_xml(num_nodes,
                                             num_elems,
                                             in);
    if(found == false){
        throw std::runtime_error("ERROR: number of points and/or cells not found in the XML file!");
        //std::cout << "ERROR: number of points and cells not found in the XML file!" << std::endl;
    }
    std::cout << "Number of nodes in the mesh file: " << num_nodes << std::endl;
    std::cout << "Number of elements in the mesh file: " << num_elems << std::endl;
    
    //------------------------------------
    // allocate mesh class nodes and elems
    mesh.initialize_nodes(num_nodes);
    if(HexN){
        mesh.initialize_elems_Pn(num_elems, Pn_order, 2*Pn_order);
    } else {
        mesh.initialize_elems(num_elems);
    }


    //------------------------------------
    // allocate the elem object id array
    mesh_inps.object_ids = DCArrayKokkos <int> (num_elems, "ObjectIDs");


    // ------------------------
    // Mesh file storage order:
    //     objectId
    //     Points
    //     connectivity
    //     offsets
    //     types
    // ------------------------
    
    // temporary arrays
    DCArrayKokkos<real_t> node_coords_vtu(num_nodes,3, "node_coords_vtu_file"); // always 3 with vtu files
    DCArrayKokkos<int> connectivity(num_elems,num_nodes_in_elem, "connectivity_vtu_file");
    DCArrayKokkos<int> elem_types(num_elems, "elem_types_vtu_file"); // element types

    node_coords = MPICArrayKokkos<real_t>(num_nodes, num_dims, "Node_coordinates_in_mesh_io");


    // for all fields, we stop recording when we get to "<"
    std::string stop = "<";

    // the size of 1D storage from reading the mesh file
    size_t size;

    // ---
    //  Object ids
    // ---

    // the object id in the element
    // array dims are (num_elems)
    found = extract_values_xml(mesh_inps.object_ids.host.pointer(),
                            "\"ObjectId\"",
                            stop,
                            in,
                            size);
    if(found==false){
        // throw std::runtime_error("ERROR: ObjectIDs were not found in the XML file!");
        std::cout << "ERROR: ObjectIDs were not found in the XML file!" << std::endl;
        mesh_inps.object_ids.host.set_values(real_t(0.0));
    }
    mesh_inps.object_ids.update_device();

    // Reset file stream to beginning to search for Points
    in.clear();
    in.seekg(0, std::ios::beg);

    // ---
    //  Nodal coordinates of mesh
    // ---

    // coordinates of the node
    // array dims are (num_nodes,dims)
    // must use the quotes around Points to read the point values
    found = extract_values_xml(node_coords_vtu.host.pointer(),
                            "\"Points\"",
                            stop,
                            in,
                            size);
    if(found==false){
        throw std::runtime_error("**** ERROR: mesh nodes were not found in the XML file! ****");
        //std::cout << "ERROR: mesh nodes were not found in the XML file!" << std::endl;
    }
    if (size!=num_nodes*3){
        throw std::runtime_error("ERROR: failed to read all the mesh nodes!");
        //std::cout << "ERROR: failed to read all the mesh nodes!" << std::endl;
    }
    node_coords_vtu.update_device();

    // dimensional scaling of the mesh
    const real_t scl_x = mesh_inps.scale_x;
    const real_t scl_y = mesh_inps.scale_y;
    const real_t scl_z = mesh_inps.scale_z;

    // save the node coordinates to the state array
    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        
        // save the nodal coordinates
        node_coords(node_gid, 0) = scl_x*node_coords_vtu(node_gid, 0); // double
        node_coords(node_gid, 1) = scl_y*node_coords_vtu(node_gid, 1); // double
        if(num_dims == 3){
            node_coords(node_gid, 2) = scl_z*node_coords_vtu(node_gid, 2); // double
        }

    }); // end for parallel nodes
    node_coords.update_host();


    // ---
    //  Nodes in the element 
    // ---

    // fill temporary nodes in the element array
    // array dims are (num_elems,num_nodes_in_elem)
    found = extract_values_xml(connectivity.host.pointer(),
                            "\"connectivity\"",
                            stop,
                            in,
                            size);
    if(found==false){
        std::cout << "ERROR: mesh connectivity was not found in the XML file!" << std::endl;
    }
    connectivity.update_device();

    // array dims are the (num_elems) 
    //    8  = pixal i,j,k linear quad format
    //    9  = linear quad ensight ordering
    //    12 = linear ensight hex ordering
    //    72 = VTK_LAGRANGE_HEXAHEDRON
    // ....
    found = extract_values_xml(elem_types.host.pointer(),
                            "\"types\"",
                            stop,
                            in,
                            size);
    if(found==false){
        std::cout << "ERROR: element types were not found in the XML file!" << std::endl;
    }
    elem_types.update_device();

    // check that the element type is supported by Fierro
    FOR_ALL (elem_gid, 0, mesh.num_elems, {
        if(elem_types(elem_gid) == element_types::linear_quad || 
           elem_types(elem_gid) == element_types::linear_hex_ijk ||
           elem_types(elem_gid) == element_types::linear_hex ||
           elem_types(elem_gid) == element_types::arbitrary_hex )
        {
            // at least one of them is true
        }
        else 
        {
           // unknown element used
           Kokkos::abort("Unknown element type in the mesh \n");
        }
    });

    // Convert from ensight linear hex to a IJK mesh
    CArrayKokkos <size_t> convert_ensight_to_ijk(8, "convert_ensight_to_ijk");

    // Convert the arbitrary order hex to a IJK mesh
    DCArrayKokkos <size_t> convert_pn_vtk_to_ijk(mesh.num_nodes_in_elem, "convert_pn_vtk_to_ijk");

    //build the connectivity for element type 12
    // elem_types.host(0)
    switch(elem_types.host(0)){

        case element_types::linear_quad:
            // the node order is correct, no changes required

            FOR_ALL (elem_gid, 0, mesh.num_elems, {
                
                for (size_t node_lid=0; node_lid<mesh.num_nodes_in_elem; node_lid++){
                    mesh.nodes_in_elem(elem_gid, node_lid) = connectivity(elem_gid,node_lid);
                }
                
            }); // end for

            break;
            // next case

        case element_types::linear_hex_ijk:

            // read the node ids in the element, no maps required
            FOR_ALL (elem_gid, 0, mesh.num_elems, {
                
                for (size_t node_lid=0; node_lid<mesh.num_nodes_in_elem; node_lid++){
                    mesh.nodes_in_elem(elem_gid, node_lid) = connectivity(elem_gid,node_lid);
                }
                
            }); // end for

            break;
            // next case

        case element_types::linear_hex:

            RUN({
                convert_ensight_to_ijk(0) = 0;
                convert_ensight_to_ijk(1) = 1;
                convert_ensight_to_ijk(2) = 3;
                convert_ensight_to_ijk(3) = 2;
                convert_ensight_to_ijk(4) = 4;
                convert_ensight_to_ijk(5) = 5;
                convert_ensight_to_ijk(6) = 7;
                convert_ensight_to_ijk(7) = 6;
            });

            // read the node ids in the element
            FOR_ALL (elem_gid, 0, mesh.num_elems, {
                
                for (size_t node_lid=0; node_lid<mesh.num_nodes_in_elem; node_lid++){
                    mesh.nodes_in_elem(elem_gid, node_lid) = connectivity(elem_gid,convert_ensight_to_ijk(node_lid));
                }
                
            }); // end for

            break;
            // next case

        case element_types::arbitrary_hex:

            // re-order the nodes to be in i,j,k format for Fierro
            size_t this_node = 0;
            for (int k=0; k<=Pn_order; k++){
                for (int j=0; j<=Pn_order; j++){
                    for (int i=0; i<=Pn_order; i++){
                        
                        // convert this_node index to the FE index convention
                        int order[3] = {Pn_order, Pn_order, Pn_order};
                        int this_index = PointIndexFromIJK(i, j, k, order);
                        
                        // store the points in this elem according the the finite
                        // element numbering convention
                        convert_pn_vtk_to_ijk.host(this_node) = this_index;
                        
                        // increment the point counting index
                        this_node = this_node + 1;
                        
                    } // end for icount
                } // end for jcount
            }  // end for kcount
            convert_pn_vtk_to_ijk.update_device();
            Kokkos::fence();

            // read the node ids in the element
            FOR_ALL (elem_gid, 0, mesh.num_elems, {
                
                for (size_t node_lid=0; node_lid<mesh.num_nodes_in_elem; node_lid++){
                    int vtk_index = convert_pn_vtk_to_ijk(node_lid); 
                    mesh.nodes_in_elem(elem_gid, node_lid) = connectivity(elem_gid,vtk_index);
                }
                
            }); // end for

            break;
            // next case

    } // end switch
    mesh.nodes_in_elem.update_host();

    in.close();
        
} // end of VTMread function
MeshReader::MeshReader()
{} // Simulation_Parameters& _simparam);
