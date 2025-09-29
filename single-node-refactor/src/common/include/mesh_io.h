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
#include <regex>    // for string pattern recoginition
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#define BUFFER_SIZE 100000
#define MAX_WORD 30

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
/// \fn PointIndexFromIJK
///
/// \brief Given (i,j,k) coordinates within the Lagrange hex, return an 
/// offset into the local connectivity (PointIds) array. The order parameter
/// must point to an array of 3 integers specifying the order along each 
/// axis of the hexahedron.
///
/////////////////////////////////////////////////////////////////////////////
inline int PointIndexFromIJK(int i, int j, int k, const int* order)
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


//-------
// word is the field name e.g., Offsets, connectivity, etc.
// stop is the phrase to stop extracting values
template <typename T>
inline bool extract_values_xml(T *values_xml,
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
inline bool extract_num_points_and_cells_xml(size_t& numberOfPoints,
                                      size_t& numberOfCells,
                                      std::ifstream& in)
{
    bool found = false;

    std::string line;

        
    // Read the file line by line looking for NumberOfPoints
    while (std::getline(in, line)) {
        //std::cout << line << std::endl;
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


//    8  = pixal i,j,k linear quad ording
//    9  = linear quad ensight ordering
//    11 = voxel i,j,k linear hex ording
//    12 = linear ensight hex ordering
//    72 = VTK_LAGRANGE_HEXAHEDRON
namespace element_types
{
    enum element_name
    {
        linear_quad_ijk = 8,
        linear_quad = 9,
        linear_hex_ijk = 11,
        linear_hex = 12,
        arbitrary_hex = 72
    };
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
    ///
    ///
    /////////////////////////////////////////////////////////////////////////////
    void read_mesh(Mesh_t& mesh,
                   State_t& State,
                   mesh_input_t& mesh_inps,
                   int      num_dims)
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
            //read_ensight_mesh(mesh, State.GaussPoints, State.node, State.corner, mesh_inps, num_dims);
        }
        else if(extension == "inp"){ // Abaqus meshfile extension
            //read_Abaqus_mesh(mesh, State, num_dims);
        }
        else if(extension == "vtk"){ // vtk file format
            read_vtk_mesh(mesh, State.GaussPoints, State.node, State.corner, mesh_inps, num_dims);
        }
        else if(extension == "vtu"){ // vtu file format
            read_vtu_mesh(mesh, State.GaussPoints, State.node, State.corner, mesh_inps, num_dims);
        }
        else{
            throw std::runtime_error("**** Mesh file extension not understood ****");
        }

    }


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

    void read_vtk_mesh(Mesh_t& mesh,
                            GaussPoint_t& GaussPoints,
                            node_t&   node,
                            corner_t& corner,
                            mesh_input_t& mesh_inps,
                            int num_dims)
    {
        char ch;
        std::string skip_line, read_line, substring;
        std::stringstream line_parse;
        std::ifstream in;  // FILE *in;

        std::cout<<"Reading VTK mesh"<<std::endl;
    
        int myrank, nranks;
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        MPI_Comm_size(MPI_COMM_WORLD,&nranks);

        int local_node_index, current_column_index;
        int buffer_loop, buffer_iteration, buffer_iterations, dof_limit, scan_loop;
        int negative_index_found = 0;
        int global_negative_index_found = 0;

        size_t read_index_start, node_rid, elem_gid;
        long long int node_gid;
        bool zero_index_base = true;
        real_t dof_value;
        real_t unit_scaling = 1;

        CArrayKokkos<char, Kokkos::LayoutRight, Kokkos::HostSpace> read_buffer;


        // read the mesh
        // --- Read the number of nodes in the mesh --- //
        size_t global_num_nodes = 0;
        if (myrank == 0)
        {
            std::cout << " NUM DIM is " << num_dims << std::endl;
            in.open(mesh_file_);
            bool found = false;

            while (found == false&&in.good()) {
                std::getline(in, read_line);
                line_parse.str("");
                line_parse.clear();
                line_parse << read_line;
                line_parse >> substring;

                // looking for the following text:
                //      POINTS %d float
                if (substring == "POINTS")
                {
                    line_parse >> global_num_nodes;
                    std::cout << "declared node count: " << global_num_nodes << std::endl;
                    if (global_num_nodes <= 0)
                    {
                        throw std::runtime_error("ERROR, NO NODES IN MESH");
                    }
                    found = true;
                } // end if
            } // end while

            if (!found){
                throw std::runtime_error("ERROR: Failed to find POINTS");
            } // end if

        } // end if(myrank==0)

        // broadcast number of nodes
        MPI_Bcast(&global_num_nodes, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

        // host_vec_array node_coords = dual_node_coords.view_host();
        // notify that the host view is going to be modified in the file readin
        // dual_node_coords.modify_host();

        // old swage method
        // mesh->init_nodes(local_nrows); // add 1 for index starting at 1

        //std::cout << "Num nodes assigned to task " << myrank << " = " << num_local_nodes << std::endl;

        // read the initial mesh coordinates
        // x-coords
        /*only task 0 reads in nodes and elements from the input file
        stores node data in a buffer and communicates once the buffer cap is reached
        or the data ends*/

        size_t num_nodes_in_elem = 1;
        for (int dim = 0; dim < num_dims; dim++) {
            num_nodes_in_elem *= 2;
        }

        int words_per_line = num_dims;
        int elem_words_per_line = num_nodes_in_elem;

        // allocate read buffer
        read_buffer = CArrayKokkos<char, Kokkos::LayoutRight, Kokkos::HostSpace>(BUFFER_SIZE, words_per_line, MAX_WORD);

        dof_limit = global_num_nodes;
        buffer_iterations = dof_limit / BUFFER_SIZE;
        if (dof_limit % BUFFER_SIZE != 0)
        {
            buffer_iterations++;
        }

        // read coords
        read_index_start = 0;
        size_t num_local_nodes;
        DistributedMap node_map;
        { //scoped so temp FArray data is auto deleted to save memory
            //allocate pre-partition node coords using contiguous decomposition
            //FArray type used since CArray type still doesnt support zoltan2 decomposition
            DistributedDFArray<real_t> node_coords_distributed(global_num_nodes, num_dims);

            // construct contiguous parallel row map now that we know the number of nodes
            node_map = node_coords_distributed.pmap;
            // map->describe(*fos,Teuchos::VERB_EXTREME);

            // set the vertices in the mesh read in
            num_local_nodes = node_map.size();
            for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
            {
                // pack buffer on rank 0
                if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
                {
                    for (buffer_loop = 0; buffer_loop < BUFFER_SIZE; buffer_loop++)
                    {
                        getline(in, read_line);
                        line_parse.clear();
                        line_parse.str(read_line);

                        for (int iword = 0; iword < words_per_line; iword++)
                        {
                            // read portions of the line into the substring variable
                            line_parse >> substring;
                            // debug print
                            // std::cout<<" "<< substring <<std::endl;
                            // assign the substring variable as a word of the read buffer
                            strcpy(&read_buffer(buffer_loop, iword, 0), substring.c_str());
                        }
                    }
                }
                else if (myrank == 0)
                {
                    buffer_loop = 0;
                    while (buffer_iteration * BUFFER_SIZE + buffer_loop < global_num_nodes) {
                        getline(in, read_line);
                        line_parse.clear();
                        line_parse.str(read_line);
                        for (int iword = 0; iword < words_per_line; iword++)
                        {
                            // read portions of the line into the substring variable
                            line_parse >> substring;
                            // debug print
                            // std::cout<<" "<< substring <<std::endl;
                            // assign the substring variable as a word of the read buffer
                            strcpy(&read_buffer(buffer_loop, iword, 0), substring.c_str());
                        }
                        buffer_loop++;
                    }
                }

                // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
                MPI_Bcast(read_buffer.pointer(), BUFFER_SIZE * words_per_line * MAX_WORD, MPI_CHAR, 0, MPI_COMM_WORLD);
                // broadcast how many nodes were read into this buffer iteration
                MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, MPI_COMM_WORLD);

                // debug_print
                // std::cout << "NODE BUFFER LOOP IS: " << buffer_loop << std::endl;
                // for(int iprint=0; iprint < buffer_loop; iprint++)
                // std::cout<<"buffer packing: " << std::string(&read_buffer(iprint,0,0)) << std::endl;
                // return;

                // determine which data to store in the swage mesh members (the local node data)
                // loop through read buffer
                for (scan_loop = 0; scan_loop < buffer_loop; scan_loop++)
                {
                    // set global node id (ensight specific order)
                    node_gid = read_index_start + scan_loop;
                    // let map decide if this node id belongs locally; if yes store data
                    if (node_map.isProcessGlobalIndex(node_gid))
                    {
                        // set local node index in this mpi rank
                        node_rid = node_map.getLocalIndex(node_gid);
                        // extract nodal position from the read buffer
                        // for tecplot format this is the three coords in the same line
                        dof_value = atof(&read_buffer(scan_loop, 0, 0));
                        node_coords_distributed.host(node_rid, 0) = dof_value * unit_scaling;
                        dof_value = atof(&read_buffer(scan_loop, 1, 0));
                        node_coords_distributed.host(node_rid, 1) = dof_value * unit_scaling;
                        if (num_dims == 3)
                        {
                            dof_value = atof(&read_buffer(scan_loop, 2, 0));
                            node_coords_distributed.host(node_rid, 2) = dof_value * unit_scaling;
                        }
                    }
                }
                read_index_start += BUFFER_SIZE;
            }
            // end of coordinate readin
            node_coords_distributed.update_device();
            // repartition node distribution
            node_coords_distributed.repartition_vector();
            //get map from repartitioned Farray and feed it into distributed CArray type; FArray data will be discared after scope
            std::vector<node_state> required_node_state = { node_state::coords };
            node_map = node_coords_distributed.pmap;
            node.initialize(node_map, num_dims, required_node_state);
            //copy coordinate data from repartitioned FArray into CArray
            FOR_ALL(node_id, 0, node_map.size(), {
                for(int idim = 0; idim < num_dims; idim++){
                    node.coords(node_id,idim) = node_coords_distributed(node_id,idim);
                }
            });
        }
        node.coords.update_host();
        //initialize some mesh data
        mesh.initialize_nodes(global_num_nodes);
        num_local_nodes = node_map.size();
        mesh.num_local_nodes = num_local_nodes;
        mesh.node_map = node_map;
        mesh.num_dims = num_dims;
        
        // debug print of nodal data

        // debug print nodal positions and indices
        /*
        std::cout << " ------------NODAL POSITIONS ON TASK " << myrank << " --------------"<<std::endl;
        for (int inode = 0; inode < local_nrows; inode++){
            std::cout << "node: " << map->getGlobalElement(inode) + 1 << " { ";
        for (int istride = 0; istride < num_dims; istride++){
            std::cout << node_coords(inode,istride) << " , ";
        }
        std::cout << " }"<< std::endl;
        }
        */

        // check that local assignments match global total

        // read in element info (ensight file format is organized in element type sections)
        // loop over this later for several element type sections

        size_t global_num_elems  = 0;
        size_t num_elems = 0;
        CArrayKokkos<int, Kokkos::LayoutRight, Kokkos::HostSpace> node_store(elem_words_per_line);

        // --- read the number of cells in the mesh ---
        // --- Read the number of vertices in the mesh --- //
        if (myrank == 0)
        {
            bool found = false;
            while (found == false&&in.good()) {
                std::getline(in, read_line);
                line_parse.str("");
                line_parse.clear();
                line_parse << read_line;
                line_parse >> substring;

                // looking for the following text:
                //      CELLS num_cells size
                if (substring == "CELLS")
                {
                    line_parse >> global_num_elems;
                    std::cout << "declared element count: " << global_num_elems << std::endl;
                    if (global_num_elems <= 0)
                    {
                        throw std::runtime_error("ERROR, NO ELEMENTS IN MESH");
                    }
                    found = true;
                } // end if
            } // end while

            if (!found){
                throw std::runtime_error("ERROR: Failed to find CELLS");
            } // end if
        } // end if(myrank==0)

        // broadcast number of elements
        MPI_Bcast(&global_num_elems, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

        //initialize num elem in mesh struct

        if (myrank == 0)
        {
            std::cout << "before mesh initialization" << std::endl;
        }

        // read in element connectivity
        // we're gonna reallocate for the words per line expected for the element connectivity
        read_buffer = CArrayKokkos<char, Kokkos::LayoutRight, Kokkos::HostSpace>(BUFFER_SIZE, elem_words_per_line, MAX_WORD);

        // calculate buffer iterations to read number of lines
        buffer_iterations = global_num_elems / BUFFER_SIZE;
        int assign_flag;

        // dynamic buffer used to store elements before we know how many this rank needs
        std::vector<size_t> element_temp(BUFFER_SIZE * elem_words_per_line);
        std::vector<size_t> global_indices_temp(BUFFER_SIZE);
        size_t buffer_max = BUFFER_SIZE * elem_words_per_line;
        size_t indices_buffer_max = BUFFER_SIZE;

        if (global_num_elems % BUFFER_SIZE != 0)
        {
            buffer_iterations++;
        }
        read_index_start = 0;
        // std::cout << "ELEMENT BUFFER ITERATIONS: " << buffer_iterations << std::endl;
        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
        {
            // pack buffer on rank 0
            if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
            {
                for (buffer_loop = 0; buffer_loop < BUFFER_SIZE; buffer_loop++)
                {
                    getline(in, read_line);
                    line_parse.clear();
                    line_parse.str(read_line);
                    // disregard node count line since we're using one element type per mesh
                    line_parse >> substring;
                    for (int iword = 0; iword < elem_words_per_line; iword++)
                    {
                        // read portions of the line into the substring variable
                        line_parse >> substring;
                        // debug print
                        // std::cout<<" "<< substring;
                        // assign the substring variable as a word of the read buffer
                        strcpy(&read_buffer(buffer_loop, iword, 0), substring.c_str());
                    }
                    // std::cout <<std::endl;
                }
            }
            else if (myrank == 0)
            {
                buffer_loop = 0;
                while (buffer_iteration * BUFFER_SIZE + buffer_loop < global_num_elems) {
                    getline(in, read_line);
                    line_parse.clear();
                    line_parse.str(read_line);
                    line_parse >> substring;
                    for (int iword = 0; iword < elem_words_per_line; iword++)
                    {
                        // read portions of the line into the substring variable
                        line_parse >> substring;
                        // debug print
                        // std::cout<<" "<< substring;
                        // assign the substring variable as a word of the read buffer
                        strcpy(&read_buffer(buffer_loop, iword, 0), substring.c_str());
                    }
                    // std::cout <<std::endl;
                    buffer_loop++;
                    // std::cout<<" "<< node_coords(node_gid, 0)<<std::endl;
                }
            }

            // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
            MPI_Bcast(read_buffer.pointer(), BUFFER_SIZE * elem_words_per_line * MAX_WORD, MPI_CHAR, 0, MPI_COMM_WORLD);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // store element connectivity that belongs to this rank
            // loop through read buffer
            for (scan_loop = 0; scan_loop < buffer_loop; scan_loop++)
            {
                // set global node id (ensight specific order)
                elem_gid = read_index_start + scan_loop;
                // add this element to the local list if any of its nodes belong to this rank according to the map
                // get list of nodes for each element line and check if they belong to the map
                assign_flag = 0;
                for (int inode = 0; inode < elem_words_per_line; inode++)
                {
                    // as we loop through the nodes belonging to this element we store them
                    // if any of these nodes belongs to this rank this list is used to store the element locally
                    node_gid = atoi(&read_buffer(scan_loop, inode, 0));
                    if (zero_index_base)
                    {
                        node_store(inode) = node_gid; // subtract 1 since file index start is 1 but code expects 0
                    }
                    else
                    {
                        node_store(inode) = node_gid - 1; // subtract 1 since file index start is 1 but code expects 0
                    }
                    if (node_store(inode) < 0)
                    {
                        negative_index_found = 1;
                    }
                    // first we add the elements to a dynamically allocated list
                    if (zero_index_base)
                    {
                        if (node_map.isProcessGlobalIndex(node_gid) && !assign_flag)
                        {
                            assign_flag = 1;
                            num_elems++;
                        }
                    }
                    else
                    {
                        if (node_map.isProcessGlobalIndex(node_gid - 1) && !assign_flag)
                        {
                            assign_flag = 1;
                            num_elems++;
                        }
                    }
                }

                if (assign_flag)
                {
                    for (int inode = 0; inode < elem_words_per_line; inode++)
                    {
                        if ((num_elems - 1) * elem_words_per_line + inode >= buffer_max)
                        {
                            element_temp.resize((num_elems - 1) * elem_words_per_line + inode + BUFFER_SIZE * elem_words_per_line);
                            buffer_max = (num_elems - 1) * elem_words_per_line + inode + BUFFER_SIZE * elem_words_per_line;
                        }
                        element_temp[(num_elems - 1) * elem_words_per_line + inode] = node_store(inode);
                        // std::cout << "VECTOR STORAGE FOR ELEM " << num_elems << " ON TASK " << myrank << " NODE " << inode+1 << " IS " << node_store(inode) + 1 << std::endl;
                    }
                    // assign global element id to temporary list
                    if (num_elems - 1 >= indices_buffer_max)
                    {
                        global_indices_temp.resize(num_elems - 1 + BUFFER_SIZE);
                        indices_buffer_max = num_elems - 1 + BUFFER_SIZE;
                    }
                    global_indices_temp[num_elems - 1] = elem_gid;
                }
            }
            read_index_start += BUFFER_SIZE;
        }

        //set global and local shared element counts
        mesh.global_num_elems = global_num_elems;

        // construct partition mapping for shared elements on each process
        DCArrayKokkos<long long int> All_Element_Global_Indices(num_elems);
        // copy temporary global indices storage to view storage
        for (int ielem = 0; ielem < num_elems; ielem++)
        {
            All_Element_Global_Indices.host(ielem) = global_indices_temp[ielem];
            if (global_indices_temp[ielem] < 0)
            {
                negative_index_found = 1;
            }
        }

        MPI_Allreduce(&negative_index_found, &global_negative_index_found, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (global_negative_index_found)
        {
            if (myrank == 0)
            {
                std::cout << "Node index less than or equal to zero detected; set \"zero_index_base: true\" under \"input_options\" in your yaml file if indices start at 0" << std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Finalize();
            exit(0);
        }

        All_Element_Global_Indices.update_device();

        // construct global map of local and shared elements (since different ranks can own the same elements due to the local node map)
        DistributedMap element_map = DistributedMap(All_Element_Global_Indices);

        //initialize elem data structures
        mesh.initialize_elems(num_elems, num_nodes_in_elem, element_map);

        // copy temporary element storage to distributed storage
        DistributedDCArray<size_t> nodes_in_elem = mesh.nodes_in_elem;

        for (int ielem = 0; ielem < num_elems; ielem++)
        {
            for (int inode = 0; inode < num_nodes_in_elem; inode++)
            {   //assign local indices to element-node connectivity (stores global indices until ghost maps are made later)
                nodes_in_elem.host(ielem, inode) = element_temp[ielem * num_nodes_in_elem + inode];
            }
        }

        // delete temporary element connectivity and index storage
        std::vector<size_t>().swap(element_temp);
        std::vector<size_t>().swap(global_indices_temp);

        // Convert ensight index system to the ijk finite element numbering convention
        // for vertices in cell
        CArrayKokkos<size_t, Kokkos::LayoutRight, Kokkos::HostSpace> convert_ensight_to_ijk(num_nodes_in_elem);
        CArrayKokkos<size_t, Kokkos::LayoutRight, Kokkos::HostSpace> tmp_ijk_indx(num_nodes_in_elem);
        convert_ensight_to_ijk(0) = 0;
        convert_ensight_to_ijk(1) = 1;
        convert_ensight_to_ijk(2) = 3;
        convert_ensight_to_ijk(3) = 2;
        convert_ensight_to_ijk(4) = 4;
        convert_ensight_to_ijk(5) = 5;
        convert_ensight_to_ijk(6) = 7;
        convert_ensight_to_ijk(7) = 6;

        for (int cell_rid = 0; cell_rid < num_elems; cell_rid++)
        {
            for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++)
            {
                tmp_ijk_indx(node_lid) = nodes_in_elem.host(cell_rid, convert_ensight_to_ijk(node_lid));
            }

            for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++)
            {
                nodes_in_elem.host(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
            }
        }
        

        nodes_in_elem.update_device();

        // initialize corner variables
        size_t num_corners = num_elems * num_nodes_in_elem;
        mesh.initialize_corners(num_corners);

        // Close mesh input file
        if (myrank == 0)
        {
            in.close();
        }
    } // end read_mesh


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn read_vtu_mesh
    ///
    /// \brief Read ASCII .vtu mesh file
    ///
    /// \param Simulation mesh
    /// \param Simulation state
    /// \param Node state struct
    /// \param Number of dimensions
    ///
    /////////////////////////////////////////////////////////////////////////////
    void read_vtu_mesh(Mesh_t& mesh,
                    GaussPoint_t& GaussPoints,
                    node_t&   node,
                    corner_t& corner,
                    mesh_input_t& mesh_inps,
                    int num_dims)
    {
        int myrank, nranks;
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        MPI_Comm_size(MPI_COMM_WORLD,&nranks);

        std::ifstream in;  // FILE *in;
        int local_node_index, current_column_index;
        int buffer_loop, buffer_iteration, buffer_iterations, dof_limit, scan_loop;

        size_t read_index_start, node_rid, elem_gid;
        long long int node_gid;
        real_t dof_value;
        real_t unit_scaling = 1;
        bool zero_index_base = true;
        int negative_index_found = 0;
        int global_negative_index_found = 0;

        CArrayKokkos<real_t, Kokkos::LayoutRight, Kokkos::HostSpace> read_buffer(BUFFER_SIZE*num_dims);
        mesh_inps.input_elem_types = true;
        mesh_inps.input_elem_objectids = true;


        // read the mesh

        // ------------------------
        // Mesh file storage order:
        //     objectId
        //     Points
        //     connectivity
        //     offsets
        //     types
        // ------------------------

        // --- Read the number of nodes in the mesh --- //
        size_t global_num_nodes = 0;
        size_t global_num_elems = 0;
        size_t num_elems = 0;
        std::streampos objectid_streampos;

        // dimensional scaling of the mesh
        const double scl_x = mesh_inps.scale_x;
        const double scl_y = mesh_inps.scale_y;
        const double scl_z = mesh_inps.scale_z;


        //
        int Pn_order = mesh_inps.p_order;
        size_t num_nodes_in_elem = 1;
        for (int dim = 0; dim < num_dims; dim++) {
            num_nodes_in_elem *= (Pn_order + 1);
        }

        if(myrank==0){
            std::cout<<"Reading VTU file in a multiblock VTK mesh"<<std::endl;
        
            bool found;
            in.open(mesh_file_);
            

            // --- extract the number of points and cells from the XML file ---
            int num_nodes;
            int num_elems;
            found = extract_num_points_and_cells_xml(global_num_nodes,
                                                    global_num_elems,
                                                    in);
            if(found==false){
                throw std::runtime_error("ERROR: number of points and/or cells not found in the XML file!");
                //std::cout << "ERROR: number of points and cells not found in the XML file!" << std::endl;
            }
            std::cout << "Number of nodes in the mesh file: " << global_num_nodes << std::endl;
            std::cout << "Number of elements in the mesh file: " << global_num_elems << std::endl;
        }
        
        // broadcast number of nodes
        MPI_Bcast(&global_num_nodes, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        
        //------------------------------------
        // allocate mesh class nodes and elems
        mesh.global_num_nodes = global_num_nodes;

        //allocate initial contiguously partitioned node coordinate vector and get map of node indices
        DistributedDFArray<real_t> node_coords_distributed(global_num_nodes, num_dims);

        // construct contiguous parallel row map now that we know the number of nodes
        DistributedMap node_map = node_coords_distributed.pmap;
        // map->describe(*fos,Teuchos::VERB_EXTREME);

        // set the vertices in the mesh read in
        size_t num_local_nodes = node_map.size();


        // for all fields, we stop recording when we get to "<"
        std::string stop = "<";

        // the size of 1D storage from reading the mesh file
        size_t size = 0;

        // ---
        //  Object ids
        // ---

        //save fstream pointer for object ids for later after element map is decided (prevents more comms if object ids is always first)
        if(myrank==0){
            objectid_streampos = in.tellg();
        }

        // ---
        //  Nodal coordinates of mesh
        // ---

        // coordinates of the node
        // array dims are (num_nodes,dims)
        // must use the quotes around Points to read the point values
        dof_limit = global_num_nodes*num_dims;
        buffer_iterations = dof_limit / (BUFFER_SIZE*num_dims);
        size_t remainder_size = dof_limit % (BUFFER_SIZE*num_dims);
        if (remainder_size != 0)
        {
            buffer_iterations++;
        }
        
        read_index_start = 0;

        //first find the block with node coords data
        if(myrank==0){

            bool found = false;

            std::string line;
            const std::string word = "\"Points\"";

            // Read the file line by line looking for specified word
            while (std::getline(in, line)) {
                //std::cout << line << std::endl;
                if (line.find(word) != std::string::npos) { // Check if the portion of the word is in the line
                    found = true;
                } 
                if(found) {

                    break;

                } // end if found

            } // end while

            if(found==false){
                throw std::runtime_error("**** ERROR: mesh nodes were not found in the XML file! ****");
                //std::cout << "ERROR: mesh nodes were not found in the XML file!" << std::endl;
            }
        }


        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
            // pack buffer on rank 0
            size_t buffer_iteration_size;
            if(buffer_iteration < buffer_iterations - 1){
                buffer_iteration_size = BUFFER_SIZE*num_dims;
            }
            else{
                buffer_iteration_size = remainder_size;
            }
            
            if(myrank==0){

                std::string line;

                // loop over the lines in the file until the buffer limit is reached
                for(int idata = 0; idata < buffer_iteration_size; idata++){

                    // extract the individual values from the stream
                    std::string value;
                    in >> value;
                    if (value == stop) { // Check if the stop word is in the line
                        break;
                    } // end if
                    read_buffer(idata) = std::stod(value);
                    size++;
                } // end for
            }
            // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
            MPI_Bcast(read_buffer.pointer(), BUFFER_SIZE*num_dims, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_iteration_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // debug_print
            // std::cout << "NODE BUFFER LOOP IS: " << buffer_loop << std::endl;
            // for(int iprint=0; iprint < buffer_loop; iprint++)
            // std::cout<<"buffer packing: " << std::string(&read_buffer(iprint,0,0)) << std::endl;
            // return;

            // determine which data to store in the swage mesh members (the local node data)
            // loop through read buffer
            for (scan_loop = 0; scan_loop < buffer_iteration_size/num_dims; scan_loop++)
            {
                // set global node id (ensight specific order)
                node_gid = read_index_start + scan_loop;
                // let map decide if this node id belongs locally; if yes store data
                if (node_map.isProcessGlobalIndex(node_gid))
                {
                    // set local node index in this mpi rank
                    node_rid = node_map.getLocalIndex(node_gid);
                    // extract nodal position from the read buffer
                    // for tecplot format this is the three coords in the same line
                    dof_value = read_buffer(scan_loop*num_dims);
                    node_coords_distributed.host(node_rid, 0) = dof_value * scl_x;
                    dof_value = read_buffer(scan_loop*num_dims + 1);
                    node_coords_distributed.host(node_rid, 1) = dof_value * scl_y;
                    if (num_dims == 3)
                    {
                        dof_value = read_buffer(scan_loop*num_dims + 2);
                        node_coords_distributed.host(node_rid, 2) = dof_value * scl_z;
                    }
                }
            }
            read_index_start += buffer_iteration_size/num_dims;
        }

        if(myrank==0){        
            if (size!=global_num_nodes*num_dims){
                throw std::runtime_error("ERROR: failed to read all the mesh nodes!");
                //std::cout << "ERROR: failed to read all the mesh nodes!" << std::endl;
            }
        }

        node_coords_distributed.update_device();
        // repartition node distribution
        node_coords_distributed.repartition_vector();
        //get map from repartitioned Farray and feed it into distributed CArray type; FArray data will be discared after scope
        std::vector<node_state> required_node_state = { node_state::coords };
        node_map = node_coords_distributed.pmap;
        node.initialize(node_map, num_dims, required_node_state);
        //copy coordinate data from repartitioned FArray into CArray
        FOR_ALL(node_id, 0, node_map.size(), {
            for(int idim = 0; idim < num_dims; idim++){
                node.coords(node_id,idim) = node_coords_distributed(node_id,idim);
            }
        });
        node.coords.update_host();
        //initialize some mesh data
        mesh.initialize_nodes(global_num_nodes);
        num_local_nodes = node_map.size();
        mesh.num_local_nodes = num_local_nodes;
        mesh.node_map = node_map;

        /***Element data***/

        // broadcast number of elems
        MPI_Bcast(&global_num_elems, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

        // ---
        //  Nodes in the element 
        // ---

        dof_limit = global_num_elems*num_nodes_in_elem;
        buffer_iterations = dof_limit / (BUFFER_SIZE*num_nodes_in_elem);
        remainder_size = dof_limit % (BUFFER_SIZE*num_nodes_in_elem);
        if (remainder_size != 0)
        {
            buffer_iterations++;
        }

        
        CArrayKokkos<long long int, Kokkos::LayoutRight, Kokkos::HostSpace> elem_read_buffer(BUFFER_SIZE*num_nodes_in_elem);
        CArrayKokkos<long long int, Kokkos::LayoutRight, Kokkos::HostSpace> node_store(num_nodes_in_elem);
        
        std::vector<size_t> element_temp(BUFFER_SIZE * num_nodes_in_elem);
        std::vector<size_t> global_indices_temp(BUFFER_SIZE);
        size_t buffer_max = BUFFER_SIZE * num_nodes_in_elem;
        size_t indices_buffer_max = BUFFER_SIZE;
        int assign_flag;
        read_index_start = 0;
        size = 0;

        //first find the block with element connectivity data data
        if(myrank==0){

            bool found = false;

            std::string line;
            const std::string word = "\"connectivity\"";

            // Read the file line by line looking for specified word
            while (std::getline(in, line)) {

                if (line.find(word) != std::string::npos) { // Check if the portion of the word is in the line
                    found = true;
                } 
                if(found) {

                    break;

                } // end if found

            } // end while

            if(found==false){
                throw std::runtime_error("ERROR: mesh connectivity was not found in the XML file!");
                //std::cout << "ERROR: mesh nodes were not found in the XML file!" << std::endl;
            }
        }


        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
            // pack buffer on rank 0
            size_t buffer_iteration_size;
            if(buffer_iteration < buffer_iterations - 1){
                buffer_iteration_size = BUFFER_SIZE*num_nodes_in_elem;
            }
            else{
                buffer_iteration_size = remainder_size;
            }
            
            if(myrank==0){

                std::string line;

                // loop over the lines in the file until the buffer limit is reached
                for(int idata = 0; idata < buffer_iteration_size; idata++){

                    // extract the individual values from the stream
                    std::string value;
                    in >> value;
                    if (value == stop) { // Check if the stop word is in the line
                        break;
                    } // end if
                    elem_read_buffer(idata) = std::stoll(value);
                    size++;
                } // end for
            }
            // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
            MPI_Bcast(elem_read_buffer.pointer(), BUFFER_SIZE*num_nodes_in_elem, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_iteration_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // debug_print
            // std::cout << "NODE BUFFER LOOP IS: " << buffer_loop << std::endl;
            // for(int iprint=0; iprint < buffer_loop; iprint++)
            // std::cout<<"buffer packing: " << std::string(&read_buffer(iprint,0,0)) << std::endl;
            // return;

            // determine which data to store in the swage mesh members (the local node data)
            // loop through read buffer
            for (scan_loop = 0; scan_loop < buffer_iteration_size/num_nodes_in_elem; scan_loop++)
            {
                // set global node id (ensight specific order)
                elem_gid = read_index_start + scan_loop;
                // add this element to the local list if any of its nodes belong to this rank according to the map
                // get list of nodes for each element line and check if they belong to the map
                assign_flag = 0;
                for (int inode = 0; inode < num_nodes_in_elem; inode++)
                {
                    // as we loop through the nodes belonging to this element we store them
                    // if any of these nodes belongs to this rank this list is used to store the element locally
                    node_gid = elem_read_buffer(scan_loop*num_nodes_in_elem + inode);
                    if (zero_index_base)
                    {
                        node_store(inode) = node_gid; // subtract 1 since file index start is 1 but code expects 0
                    }
                    else
                    {
                        node_store(inode) = node_gid - 1; // subtract 1 since file index start is 1 but code expects 0
                    }
                    if (node_store(inode) < 0)
                    {
                        negative_index_found = 1;
                    }
                    // first we add the elements to a dynamically allocated list
                    if (zero_index_base)
                    {
                        if (node_map.isProcessGlobalIndex(node_gid) && !assign_flag)
                        {
                            assign_flag = 1;
                            num_elems++;
                        }
                    }
                    else
                    {
                        if (node_map.isProcessGlobalIndex(node_gid - 1) && !assign_flag)
                        {
                            assign_flag = 1;
                            num_elems++;
                        }
                    }
                }

                if (assign_flag)
                {
                    for (int inode = 0; inode < num_nodes_in_elem; inode++)
                    {
                        if ((num_elems - 1) * num_nodes_in_elem + inode >= buffer_max)
                        {
                            element_temp.resize((num_elems - 1) * num_nodes_in_elem + inode + BUFFER_SIZE * num_nodes_in_elem);
                            buffer_max = (num_elems - 1) * num_nodes_in_elem + inode + BUFFER_SIZE * num_nodes_in_elem;
                        }
                        element_temp[(num_elems - 1) * num_nodes_in_elem + inode] = node_store(inode);
                        // std::cout << "VECTOR STORAGE FOR ELEM " << num_elems << " ON TASK " << myrank << " NODE " << inode+1 << " IS " << node_store(inode) + 1 << std::endl;
                    }
                    // assign global element id to temporary list
                    if (num_elems - 1 >= indices_buffer_max)
                    {
                        global_indices_temp.resize(num_elems - 1 + BUFFER_SIZE);
                        indices_buffer_max = num_elems - 1 + BUFFER_SIZE;
                    }
                    global_indices_temp[num_elems - 1] = elem_gid;
                }
            }
            read_index_start += buffer_iteration_size/num_nodes_in_elem;
        }

        if(myrank==0){        
            if (size!=global_num_elems*num_nodes_in_elem){
                throw std::runtime_error("ERROR: failed to read all the mesh elements!");
                //std::cout << "ERROR: failed to read all the mesh nodes!" << std::endl;
            }
        }

        //interface with permanent storage of nodes in elem and construct element maps
        DCArrayKokkos<long long int> All_Element_Global_Indices(num_elems);
        // copy temporary global indices storage to view storage
        for (int ielem = 0; ielem < num_elems; ielem++)
        {
            All_Element_Global_Indices.host(ielem) = global_indices_temp[ielem];
            if (global_indices_temp[ielem] < 0)
            {
                negative_index_found = 1;
            }
        }

        MPI_Allreduce(&negative_index_found, &global_negative_index_found, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (global_negative_index_found)
        {
            if (myrank == 0)
            {
                std::cout << "Node index less than or equal to zero detected; set \"zero_index_base: true\" under \"input_options\" in your yaml file if indices start at 0" << std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Finalize();
            exit(0);
        }

        All_Element_Global_Indices.update_device();

        // construct global map of local and shared elements (since different ranks can own the same elements due to the local node map)
        DistributedMap element_map = DistributedMap(All_Element_Global_Indices);

        //initialize elem data structures
        mesh.initialize_elems(num_elems, num_nodes_in_elem, element_map);

        // copy temporary element storage to distributed storage
        DistributedDCArray<size_t> nodes_in_elem = mesh.nodes_in_elem;

        for (int ielem = 0; ielem < num_elems; ielem++)
        {
            for (int inode = 0; inode < num_nodes_in_elem; inode++)
            {   //assign local indices to element-node connectivity (stores global indices until ghost maps are made later)
                nodes_in_elem.host(ielem, inode) = element_temp[ielem * num_nodes_in_elem + inode];
            }
        }
        nodes_in_elem.update_device();

        //debug print
        //for (int inode = 0; inode < num_nodes_in_elem; inode++)
            //std::cout << "Element nodes " << inode << " " << nodes_in_elem(num_elems-1,inode) << std::endl;

        // delete temporary element connectivity and index storage
        std::vector<size_t>().swap(element_temp);
        std::vector<size_t>().swap(global_indices_temp);

        /****read in element types ****/
        DCArrayKokkos<int> elem_types(num_elems, "elem_types_vtu_file"); // element types

        // array dims are the (num_elems) 
        //    8  = pixal i,j,k linear quad format
        //    9  = linear quad ensight ordering
        //    12 = linear ensight hex ordering
        //    72 = VTK_LAGRANGE_HEXAHEDRON
        // ....

        dof_limit = global_num_elems;
        buffer_iterations = dof_limit / (BUFFER_SIZE);
        remainder_size = dof_limit % (BUFFER_SIZE);
        if (remainder_size != 0)
        {
            buffer_iterations++;
        }

        read_index_start = 0;
        size = 0;
        elem_read_buffer = CArrayKokkos<long long int, Kokkos::LayoutRight, Kokkos::HostSpace>(BUFFER_SIZE);

        //allocate types array

        //first find the block with element connectivity data data
        if(myrank==0){

            bool found = false;

            std::string line;
            const std::string word = "\"types\"";

            // Read the file line by line looking for specified word
            while (std::getline(in, line)) {

                if (line.find(word) != std::string::npos) { // Check if the portion of the word is in the line
                    found = true;
                } 
                if(found) {

                    break;

                } // end if found

            } // end while

            if(found==false){
                throw std::runtime_error("ERROR: mesh element types were not found in the XML file!");
                //std::cout << "ERROR: mesh nodes were not found in the XML file!" << std::endl;
            }
        }


        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
            // pack buffer on rank 0
            size_t buffer_iteration_size;
            if(buffer_iteration < buffer_iterations - 1){
                buffer_iteration_size = BUFFER_SIZE;
            }
            else{
                buffer_iteration_size = remainder_size;
            }
            
            if(myrank==0){

                std::string line;

                // loop over the lines in the file until the buffer limit is reached
                for(int idata = 0; idata < buffer_iteration_size; idata++){

                    // extract the individual values from the stream
                    std::string value;
                    in >> value;
                    if (value == stop) { // Check if the stop word is in the line
                        break;
                    } // end if
                    elem_read_buffer(idata) = std::stoi(value);
                    // if(elem_read_buffer(i)!=11){
                    //     std::cout << "FOUND TYPE ERROR " << std::endl;
                    // }
                    size++;
                } // end for
            }
            // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
            MPI_Bcast(elem_read_buffer.pointer(), BUFFER_SIZE, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_iteration_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // debug_print
            // std::cout << "NODE BUFFER LOOP IS: " << buffer_loop << std::endl;
            // for(int iprint=0; iprint < buffer_loop; iprint++)
            // std::cout<<"buffer packing: " << std::string(&read_buffer(iprint,0,0)) << std::endl;
            // return;

            // determine which data to store in the swage mesh members (the local node data)
            // loop through read buffer
            for (scan_loop = 0; scan_loop < buffer_iteration_size; scan_loop++)
            {
                // set global node id (ensight specific order)
                elem_gid = read_index_start + scan_loop;
                //add to the local type array if this elem gid belongs to this rank
                if (element_map.isProcessGlobalIndex(elem_gid)){
                    //std::cout << "elocal id " << element_map.getLocalIndex(elem_gid) << " " << elem_read_buffer(scan_loop) << std::endl;
                    elem_types.host(element_map.getLocalIndex(elem_gid)) = elem_read_buffer(scan_loop);
                }
                

            }
            read_index_start += buffer_iteration_size;
        }

        elem_types.update_device();
        //std::cout << "Type read size " << size << std::endl;

        // check that the element type is supported by Fierro
        FOR_ALL (elem_gid, 0, num_elems, {
            //std::cout << "Element type is " << elem_types(elem_gid) << std::endl;
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
        DCArrayKokkos <size_t> convert_pn_vtk_to_ijk(num_nodes_in_elem, "convert_pn_vtk_to_ijk");

        //build the connectivity for element type 12
        // elem_types.host(0)
        switch(elem_types.host(0)){

            case element_types::linear_quad:
                // the node order is correct, no changes required

                // FOR_ALL (elem_gid, 0, mesh.num_elems, {
                    
                //     for (size_t node_lid=0; node_lid<mesh.num_nodes_in_elem; node_lid++){
                //         mesh.nodes_in_elem(elem_gid, node_lid) = connectivity(elem_gid,node_lid);
                //     }
                    
                // }); // end for

                break;
                // next case

            case element_types::linear_hex_ijk:

                // read the node ids in the element, no maps required
                // FOR_ALL (elem_gid, 0, mesh.num_elems, {
                    
                //     for (size_t node_lid=0; node_lid<mesh.num_nodes_in_elem; node_lid++){
                //         mesh.nodes_in_elem(elem_gid, node_lid) = connectivity(elem_gid,node_lid);
                //     }
                    
                // }); // end for

                break;
                // next case

            case element_types::linear_hex:
                std::cout << "Converting mesh element ordering" << std::endl;
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
                FOR_ALL (elem_id, 0, num_elems, {
                    long long int temp[8];
                    for (size_t node_lid=0; node_lid<8; node_lid++){
                        temp[node_lid] =  nodes_in_elem(elem_id,convert_ensight_to_ijk(node_lid));
                    }
                    for (size_t node_lid=0; node_lid<8; node_lid++){
                        nodes_in_elem(elem_id, node_lid) = temp[node_lid];
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
                            convert_pn_vtk_to_ijk.host(this_index) = this_node;
                            
                            // increment the point counting index
                            this_node = this_node + 1;
                            
                        } // end for icount
                    } // end for jcount
                }  // end for kcount
                convert_pn_vtk_to_ijk.update_device();
                Kokkos::fence();

                // read the node ids in the element
                FOR_ALL (elem_id, 0, mesh.num_elems, {
                    
                    for (size_t node_lid=0; node_lid<mesh.num_nodes_in_elem; node_lid++){
                        mesh.nodes_in_elem(elem_id, node_lid) = mesh.nodes_in_elem(elem_id,convert_pn_vtk_to_ijk(node_lid));
                    }
                    
                }); // end for

                break;
                // next case

        } // end switch
        mesh.nodes_in_elem.update_host();


        // initialize corner variables
        size_t num_corners = mesh.num_elems * mesh.num_nodes_in_elem;
        mesh.initialize_corners(num_corners);
        
        //------------------------------------
        // allocate the elem object id array
        mesh_inps.object_ids = DCArrayKokkos <int> (num_elems, "ObjectIDs");
    
        //reset file stream pointer to area where ObjectId is located before node coordinates
        in.seekg(objectid_streampos);

        dof_limit = global_num_elems;
        buffer_iterations = dof_limit / (BUFFER_SIZE);
        remainder_size = dof_limit % (BUFFER_SIZE);
        if (remainder_size != 0)
        {
            buffer_iterations++;
        }

        read_index_start = 0;
        size = 0;

        //allocate types array

        //first find the block with element connectivity data data
        if(myrank==0){

            bool found = false;

            std::string line;
            const std::string word = "\"ObjectId\"";

            // Read the file line by line looking for specified word
            while (std::getline(in, line)) {

                if (line.find(word) != std::string::npos) { // Check if the portion of the word is in the line
                    found = true;
                } 
                if(found) {

                    break;

                } // end if found

            } // end while

            if(found==false){
                throw std::runtime_error("ERROR: ObjectIDs were not found in the XML file!");
                //std::cout << "ERROR: mesh nodes were not found in the XML file!" << std::endl;
            }
        }


        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
            // pack buffer on rank 0
            size_t buffer_iteration_size;
            if(buffer_iteration < buffer_iterations - 1){
                buffer_iteration_size = BUFFER_SIZE;
            }
            else{
                buffer_iteration_size = remainder_size;
            }
            
            if(myrank==0){

                std::string line;

                // loop over the lines in the file until the buffer limit is reached
                for(int idata = 0; idata < buffer_iteration_size; idata++){

                    // extract the individual values from the stream
                    std::string value;
                    in >> value;
                    if (value == stop) { // Check if the stop word is in the line
                        break;
                    } // end if
                    elem_read_buffer(idata) = std::stoi(value);
                    size++;
                } // end for
            }
            // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
            MPI_Bcast(elem_read_buffer.pointer(), BUFFER_SIZE, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_iteration_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // debug_print
            // std::cout << "NODE BUFFER LOOP IS: " << buffer_loop << std::endl;
            // for(int iprint=0; iprint < buffer_loop; iprint++)
            // std::cout<<"buffer packing: " << std::string(&read_buffer(iprint,0,0)) << std::endl;
            // return;

            // determine which data to store in the swage mesh members (the local node data)
            // loop through read buffer
            for (scan_loop = 0; scan_loop < buffer_iteration_size; scan_loop++)
            {
                // set global node id (ensight specific order)
                elem_gid = read_index_start + scan_loop;
                //add to the local type array if this elem gid belongs to this rank
                if (element_map.isProcessGlobalIndex(elem_gid)){
                    mesh_inps.object_ids.host(element_map.getLocalIndex(elem_gid)) = elem_read_buffer(scan_loop);
                }
                

            }
            read_index_start += buffer_iteration_size;
        }
        
        mesh_inps.object_ids.update_device();

        in.close();
            
    } // end of VTMread function


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
        SimulationParameters_t& SimulationParameters)
    {
        if (SimulationParameters.mesh_input.num_dims == 2) {
            if (SimulationParameters.mesh_input.type == mesh_input::Polar) {
                build_2d_polar(mesh, GaussPoints, node, corner, SimulationParameters);
            }
            else if (SimulationParameters.mesh_input.type == mesh_input::Box) {
                build_2d_box(mesh, GaussPoints, node, corner, SimulationParameters);
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
        else if (SimulationParameters.mesh_input.num_dims == 3) {
            build_3d_box(mesh, GaussPoints, node, corner, SimulationParameters);
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
        SimulationParameters_t& SimulationParameters) const
    {
        int myrank, nranks;
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        MPI_Comm_size(MPI_COMM_WORLD,&nranks);
        /*currently we just build the global mesh data on rank 0 and then broadcast relevant data to each rank
          before the global mesh data on rank 0 falls out of scope*/
        int global_num_nodes, global_num_elems;
        CArrayKokkos<real_t, Kokkos::LayoutRight, Kokkos::HostSpace> read_buffer;
        CArrayKokkos<long long int, Kokkos::LayoutRight, Kokkos::HostSpace> read_buffer_edof;
        CArrayKokkos<double, Kokkos::LayoutLeft, Kokkos::HostSpace> global_coords;
        CArrayKokkos<double, Kokkos::LayoutLeft, Kokkos::HostSpace> global_nodes_in_elem;
        auto convert_point_number_in_quad = CArray<int>(4);
        convert_point_number_in_quad(0) = 0;
        convert_point_number_in_quad(1) = 1;
        convert_point_number_in_quad(2) = 3;
        convert_point_number_in_quad(3) = 2;

        int local_node_index, current_column_index;
        int buffer_loop, buffer_iteration, buffer_iterations, dof_limit, scan_loop;
        int negative_index_found = 0;
        int global_negative_index_found = 0;

        size_t read_index_start, node_rid, elem_gid;
        long long int node_gid;
        real_t dof_value;
        real_t unit_scaling = 1;

        const int num_dims = 2;
        size_t num_nodes_in_elem = 1;
        for (int dim = 0; dim < num_dims; dim++) {
            num_nodes_in_elem *= 2;
        }
        if(myrank==0){
            printf("Creating a 2D box mesh \n");
        }

        // SimulationParameters.mesh_input.length.update_host();
        const double lx = SimulationParameters.mesh_input.length[0];
        const double ly = SimulationParameters.mesh_input.length[1];

        // SimulationParameters.mesh_input.num_elems.update_host();
        const int num_elems_i = SimulationParameters.mesh_input.num_elems[0];
        const int num_elems_j = SimulationParameters.mesh_input.num_elems[1];

        const int num_points_i = num_elems_i + 1; // num points in x
        const int num_points_j = num_elems_j + 1; // num points in y

        global_num_nodes = num_points_i * num_points_j;

        const double dx = lx / ((double)num_elems_i);  // len/(num_elems_i)
        const double dy = ly / ((double)num_elems_j);  // len/(num_elems_j)

        global_num_elems = num_elems_i * num_elems_j;
            
        if(myrank==0){
            std::vector<double> origin(num_dims);
            // SimulationParameters.mesh_input.origin.update_host();
            for (int i = 0; i < num_dims; i++) { origin[i] = SimulationParameters.mesh_input.origin[i]; }

            // --- 3D parameters ---
            // const int num_faces_in_elem  = 6;  // number of faces in elem
            // const int num_points_in_elem = 8;  // number of points in elem
            // const int num_points_in_face = 4;  // number of points in a face
            // const int num_edges_in_elem  = 12; // number of edges in a elem
            
            // node coords data on rank 0 for all global nodes
            global_coords = CArrayKokkos<double, Kokkos::LayoutLeft, Kokkos::HostSpace>(global_num_nodes, num_dims, "global_mesh_build_node_coordinates");

            // --- Build nodes ---

            // populate the point data structures
            for (int j = 0; j < num_points_j; j++) {
                for (int i = 0; i < num_points_i; i++) {
                    // global id for the point
                    int node_gid = get_id(i, j, 0, num_points_i, num_points_j);

                    // store the point coordinates
                    global_coords(node_gid, 0) = origin[0] + (double)i * dx;
                    global_coords(node_gid, 1) = origin[1] + (double)j * dy;
                } // end for i
            } // end for j


            // initialize elem variables
            global_nodes_in_elem = CArrayKokkos<double, Kokkos::LayoutLeft, Kokkos::HostSpace>(global_num_elems, num_nodes_in_elem, "global_mesh_build_nodes_in_elem");

            // --- Build elems  ---

            // populate the elem center data structures
            for (int j = 0; j < num_elems_j; j++) {
                for (int i = 0; i < num_elems_i; i++) {
                    // global id for the elem
                    int elem_gid = get_id(i, j, 0, num_elems_i, num_elems_j);

                    // store the point IDs for this elem where the range is
                    // (i:i+1, j:j+1, k:k+1) for a linear hexahedron
                    int this_point = 0;
                    for (int jcount = j; jcount <= j + 1; jcount++) {
                        for (int icount = i; icount <= i + 1; icount++) {
                            // global id for the points
                            int node_gid = get_id(icount, jcount, 0,
                                            num_points_i, num_points_j);

                            // convert this_point index to the FE index convention
                            int this_index = convert_point_number_in_quad(this_point);

                            // store the points in this elem according the the finite
                            // element numbering convention
                            global_nodes_in_elem(elem_gid, this_index) = node_gid;

                            // increment the point counting index
                            this_point = this_point + 1;
                        } // end for icount
                    } // end for jcount
                } // end for i
            } // end for j
        }

        //distribute partitioned data from the global mesh build data on rank 0
        size_t num_local_nodes;
        DistributedMap node_map;
        // read coords
        read_index_start = 0;
        
        buffer_iterations = global_num_nodes / BUFFER_SIZE;
        if (global_num_nodes % BUFFER_SIZE != 0)
        {
            buffer_iterations++;
        }

        read_buffer = CArrayKokkos<real_t, Kokkos::LayoutRight, Kokkos::HostSpace>(BUFFER_SIZE, num_dims);

        { //scoped so temp FArray data is auto deleted to save memory
            //allocate pre-partition node coords using contiguous decomposition
            //FArray type used since CArray type still doesnt support zoltan2 decomposition
            DistributedDFArray<real_t> node_coords_distributed(global_num_nodes, num_dims);

            // construct contiguous parallel row map now that we know the number of nodes
            node_map = node_coords_distributed.pmap;
            // map->describe(*fos,Teuchos::VERB_EXTREME);

            // set the vertices in the mesh read in
            num_local_nodes = node_map.size();
            for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
            {
                // pack buffer on rank 0
                if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
                {
                    for (buffer_loop = 0; buffer_loop < BUFFER_SIZE; buffer_loop++)
                    {

                        for (int idim = 0; idim < num_dims; idim++)
                        {
                            // debug print
                            // std::cout<<" "<< substring <<std::endl;
                            // assign the substring variable as a word of the read buffer
                            read_buffer(buffer_loop, idim) = global_coords(buffer_iteration * BUFFER_SIZE + buffer_loop, idim);
                        }
                    }
                }
                else if (myrank == 0)
                {
                    buffer_loop = 0;
                    while (buffer_iteration * BUFFER_SIZE + buffer_loop < global_num_nodes) {
                        for (int idim = 0; idim < num_dims; idim++)
                        {
                            // debug print
                            // std::cout<<" "<< substring <<std::endl;
                            // assign the substring variable as a word of the read buffer
                            read_buffer(buffer_loop, idim) = global_coords(buffer_iteration * BUFFER_SIZE + buffer_loop, idim);
                        }
                        buffer_loop++;
                    }
                }

                // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
                MPI_Bcast(read_buffer.pointer(), BUFFER_SIZE * num_dims, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                // broadcast how many nodes were read into this buffer iteration
                MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, MPI_COMM_WORLD);

                // debug_print
                // std::cout << "NODE BUFFER LOOP IS: " << buffer_loop << std::endl;
                // for(int iprint=0; iprint < buffer_loop; iprint++)
                // std::cout<<"buffer packing: " << std::string(&read_buffer(iprint,0,0)) << std::endl;
                // return;

                // determine which data to store in the swage mesh members (the local node data)
                // loop through read buffer
                for (scan_loop = 0; scan_loop < buffer_loop; scan_loop++)
                {
                    // set global node id (ensight specific order)
                    node_gid = read_index_start + scan_loop;
                    // let map decide if this node id belongs locally; if yes store data
                    if (node_map.isProcessGlobalIndex(node_gid))
                    {
                        // set local node index in this mpi rank
                        node_rid = node_map.getLocalIndex(node_gid);
                        // extract nodal position from the read buffer
                        // for tecplot format this is the three coords in the same line
                        dof_value = read_buffer(scan_loop,0);
                        node_coords_distributed.host(node_rid, 0) = dof_value * unit_scaling;
                        dof_value = read_buffer(scan_loop,1);
                        node_coords_distributed.host(node_rid, 1) = dof_value * unit_scaling;
                        if (num_dims == 3)
                        {
                            dof_value = read_buffer(scan_loop,2);
                            node_coords_distributed.host(node_rid, 2) = dof_value * unit_scaling;
                        }
                    }
                }
                read_index_start += BUFFER_SIZE;
            }
            // end of coordinate readin
            node_coords_distributed.update_device();
            // repartition node distribution
            node_coords_distributed.repartition_vector();
            //get map from repartitioned Farray and feed it into distributed CArray type; FArray data will be discared after scope
            std::vector<node_state> required_node_state = { node_state::coords };
            node_map = node_coords_distributed.pmap;
            node.initialize(node_map, num_dims, required_node_state);
            //copy coordinate data from repartitioned FArray into CArray
            FOR_ALL(node_id, 0, node_map.size(), {
                for(int idim = 0; idim < num_dims; idim++){
                    node.coords(node_id,idim) = node_coords_distributed(node_id,idim);
                }
            });
        }

        //initialize some mesh data
        mesh.initialize_nodes(global_num_nodes);
        num_local_nodes = node_map.size();
        mesh.num_local_nodes = num_local_nodes;
        mesh.node_map = node_map;
        mesh.num_dims = num_dims;
        //node.coords.print();
        
        // debug print of nodal data

        // debug print nodal positions and indices
        /*
        std::cout << " ------------NODAL POSITIONS ON TASK " << myrank << " --------------"<<std::endl;
        for (int inode = 0; inode < local_nrows; inode++){
            std::cout << "node: " << map->getGlobalElement(inode) + 1 << " { ";
        for (int istride = 0; istride < num_dims; istride++){
            std::cout << node_coords(inode,istride) << " , ";
        }
        std::cout << " }"<< std::endl;
        }
        */

        // check that local assignments match global total

        // read in element info (ensight file format is organized in element type sections)
        // loop over this later for several element type sections
        CArrayKokkos<int, Kokkos::LayoutRight, Kokkos::HostSpace> node_store(num_nodes_in_elem);

        // broadcast number of elements
        MPI_Bcast(&global_num_elems, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

        //initialize num elem in mesh struct

        if (myrank == 0)
        {
            std::cout << "before mesh initialization" << std::endl;
        }

        // read in element connectivity
        // we're gonna reallocate for the words per line expected for the element connectivity
        read_buffer_edof = CArrayKokkos<long long int, Kokkos::LayoutRight, Kokkos::HostSpace>(BUFFER_SIZE, num_nodes_in_elem);

        // calculate buffer iterations to read number of lines
        buffer_iterations = global_num_elems / BUFFER_SIZE;
        int assign_flag;

        // dynamic buffer used to store elements before we know how many this rank needs
        std::vector<size_t> element_temp(BUFFER_SIZE * num_nodes_in_elem);
        std::vector<size_t> global_indices_temp(BUFFER_SIZE);
        size_t buffer_max = BUFFER_SIZE * num_nodes_in_elem;
        size_t indices_buffer_max = BUFFER_SIZE;
        size_t num_elems = 0;

        if (global_num_elems % BUFFER_SIZE != 0)
        {
            buffer_iterations++;
        }
        read_index_start = 0;
        // std::cout << "ELEMENT BUFFER ITERATIONS: " << buffer_iterations << std::endl;
        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
        {
            // pack buffer on rank 0
            if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
            {
                for (buffer_loop = 0; buffer_loop < BUFFER_SIZE; buffer_loop++)
                {
                    for (int inode = 0; inode < num_nodes_in_elem; inode++)
                    {
                        read_buffer_edof(buffer_loop,inode) = global_nodes_in_elem(buffer_iteration * BUFFER_SIZE + buffer_loop, inode);
                    }
                    // std::cout <<std::endl;
                }
            }
            else if (myrank == 0)
            {
                buffer_loop = 0;
                while (buffer_iteration * BUFFER_SIZE + buffer_loop < global_num_elems) {
                    for (int inode = 0; inode < num_nodes_in_elem; inode++)
                    {
                        read_buffer_edof(buffer_loop,inode) = global_nodes_in_elem(buffer_iteration * BUFFER_SIZE + buffer_loop, inode);
                    }
                    // std::cout <<std::endl;
                    buffer_loop++;
                    // std::cout<<" "<< node_coords(node_gid, 0)<<std::endl;
                }
            }

            // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
            MPI_Bcast(read_buffer_edof.pointer(), BUFFER_SIZE * num_nodes_in_elem, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // store element connectivity that belongs to this rank
            // loop through read buffer
            for (scan_loop = 0; scan_loop < buffer_loop; scan_loop++)
            {
                // set global node id (ensight specific order)
                elem_gid = read_index_start + scan_loop;
                // add this element to the local list if any of its nodes belong to this rank according to the map
                // get list of nodes for each element line and check if they belong to the map
                assign_flag = 0;
                for (int inode = 0; inode < num_nodes_in_elem; inode++)
                {
                    // as we loop through the nodes belonging to this element we store them
                    // if any of these nodes belongs to this rank this list is used to store the element locally
                    node_gid = read_buffer_edof(scan_loop, inode);
                    node_store(inode) = node_gid;
                    // first we add the elements to a dynamically allocated list
                    if (node_map.isProcessGlobalIndex(node_gid) && !assign_flag)
                    {
                        assign_flag = 1;
                        num_elems++;
                    }
                }

                if (assign_flag)
                {
                    for (int inode = 0; inode < num_nodes_in_elem; inode++)
                    {
                        if ((num_elems - 1) * num_nodes_in_elem + inode >= buffer_max)
                        {
                            element_temp.resize((num_elems - 1) * num_nodes_in_elem + inode + BUFFER_SIZE * num_nodes_in_elem);
                            buffer_max = (num_elems - 1) * num_nodes_in_elem + inode + BUFFER_SIZE * num_nodes_in_elem;
                        }
                        element_temp[(num_elems - 1) * num_nodes_in_elem + inode] = node_store(inode);
                        // std::cout << "VECTOR STORAGE FOR ELEM " << num_elems << " ON TASK " << myrank << " NODE " << inode+1 << " IS " << node_store(inode) + 1 << std::endl;
                    }
                    // assign global element id to temporary list
                    if (num_elems - 1 >= indices_buffer_max)
                    {
                        global_indices_temp.resize(num_elems - 1 + BUFFER_SIZE);
                        indices_buffer_max = num_elems - 1 + BUFFER_SIZE;
                    }
                    global_indices_temp[num_elems - 1] = elem_gid;
                }
            }
            read_index_start += BUFFER_SIZE;
        }
        //set global and local shared element counts
        mesh.global_num_elems = global_num_elems;

        // construct partition mapping for shared elements on each process
        DCArrayKokkos<long long int> All_Element_Global_Indices(num_elems);
        // copy temporary global indices storage to view storage
        for (int ielem = 0; ielem < num_elems; ielem++)
        {
            All_Element_Global_Indices.host(ielem) = global_indices_temp[ielem];
        }

        All_Element_Global_Indices.update_device();

        // construct global map of local and shared elements (since different ranks can own the same elements due to the local node map)
        DistributedMap element_map = DistributedMap(All_Element_Global_Indices);

        //initialize elem data structures
        mesh.initialize_elems(num_elems, num_nodes_in_elem, element_map);

        // copy temporary element storage to distributed storage
        DistributedDCArray<size_t> nodes_in_elem = mesh.nodes_in_elem;

        for (int ielem = 0; ielem < num_elems; ielem++)
        {
            for (int inode = 0; inode < num_nodes_in_elem; inode++)
            {   //assign local indices to element-node connectivity (stores global indices until ghost maps are made later)
                nodes_in_elem.host(ielem, inode) = element_temp[ielem * num_nodes_in_elem + inode];
            }
        }

        nodes_in_elem.update_device();
        
        // delete temporary element connectivity and index storage
        //std::vector<size_t>().swap(element_temp);
        //std::vector<size_t>().swap(global_indices_temp);

        // initialize corner variables
        size_t num_corners = num_elems * num_nodes_in_elem;
        mesh.initialize_corners(num_corners);

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
        SimulationParameters_t& SimulationParameters) const
    {
        int myrank, nranks;
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        MPI_Comm_size(MPI_COMM_WORLD,&nranks);
        /*currently we just build the global mesh data on rank 0 and then broadcast relevant data to each rank
          before the global mesh data on rank 0 falls out of scope*/
        int global_num_nodes, global_num_elems;
        CArrayKokkos<real_t, Kokkos::LayoutRight, Kokkos::HostSpace> read_buffer;
        CArrayKokkos<long long int, Kokkos::LayoutRight, Kokkos::HostSpace> read_buffer_edof;
        CArrayKokkos<double, Kokkos::LayoutLeft, Kokkos::HostSpace> global_coords;
        CArrayKokkos<double, Kokkos::LayoutLeft, Kokkos::HostSpace> global_nodes_in_elem;
        auto convert_point_number_in_quad = CArray<int>(4);
        convert_point_number_in_quad(0) = 0;
        convert_point_number_in_quad(1) = 1;
        convert_point_number_in_quad(2) = 3;
        convert_point_number_in_quad(3) = 2;

        int local_node_index, current_column_index;
        int buffer_loop, buffer_iteration, buffer_iterations, dof_limit, scan_loop;
        int negative_index_found = 0;
        int global_negative_index_found = 0;

        size_t read_index_start, node_rid, elem_gid;
        long long int node_gid;
        real_t dof_value;
        real_t unit_scaling = 1;

        const int num_dims = 2;
        size_t num_nodes_in_elem = 1;
        for (int dim = 0; dim < num_dims; dim++) {
            num_nodes_in_elem *= 2;
        }
        if(myrank==0){
            printf("Creating a 2D box mesh \n");
        }

        // SimulationParameters.mesh_input.length.update_host();
        const double inner_radius = SimulationParameters.mesh_input.inner_radius;
        const double outer_radius = SimulationParameters.mesh_input.outer_radius;

        const double start_angle = PI / 180.0 * SimulationParameters.mesh_input.starting_angle;
        const double end_angle   = PI / 180.0 * SimulationParameters.mesh_input.ending_angle;

        const int num_elems_i = SimulationParameters.mesh_input.num_radial_elems;
        const int num_elems_j = SimulationParameters.mesh_input.num_angular_elems;

        const int num_points_i = num_elems_i + 1; // num points in x
        const int num_points_j = num_elems_j + 1; // num points in y

        const double dx = (outer_radius - inner_radius) / ((double)num_elems_i);  // len/(elems)
        const double dy = (end_angle - start_angle) / ((double)num_elems_j);  // len/(elems)

        global_num_elems = num_elems_i * num_elems_j;

        global_num_nodes = num_points_i * num_points_j;
            
        if(myrank==0){
            std::vector<double> origin(num_dims);
            // SimulationParameters.mesh_input.origin.update_host();
            for (int i = 0; i < num_dims; i++) { origin[i] = SimulationParameters.mesh_input.origin[i]; }

            // --- 3D parameters ---
            // const int num_faces_in_elem  = 6;  // number of faces in elem
            // const int num_points_in_elem = 8;  // number of points in elem
            // const int num_points_in_face = 4;  // number of points in a face
            // const int num_edges_in_elem  = 12; // number of edges in a elem
            
            // node coords data on rank 0 for all global nodes
            global_coords = CArrayKokkos<double, Kokkos::LayoutLeft, Kokkos::HostSpace>(global_num_nodes, num_dims, "global_mesh_build_node_coordinates");

            // --- Build nodes ---

            // populate the point data structures
            for (int j = 0; j < num_points_j; j++) {
                for (int i = 0; i < num_points_i; i++) {
                    // global id for the point
                    int node_gid = get_id(i, j, 0, num_points_i, num_points_j);

                    double r_i     = inner_radius + (double)i * dx;
                    double theta_j = start_angle + (double)j * dy;

                    // store the point coordinates
                    global_coords(node_gid, 0) = origin[0] + r_i * cos(theta_j);
                    global_coords(node_gid, 1) = origin[1] + r_i * sin(theta_j);

                    if(global_coords(node_gid, 0) < 0.0){
                        throw std::runtime_error("**** NODE RADIUS FOR RZ MESH MUST BE POSITIVE ****");
                    }
                } // end for i
            } // end for j


            // initialize elem variables
            global_nodes_in_elem = CArrayKokkos<double, Kokkos::LayoutLeft, Kokkos::HostSpace>(global_num_elems, num_nodes_in_elem, "global_mesh_build_nodes_in_elem");

            // --- Build elems  ---

            // populate the elem center data structures
            for (int j = 0; j < num_elems_j; j++) {
                for (int i = 0; i < num_elems_i; i++) {
                    // global id for the elem
                    int elem_gid = get_id(i, j, 0, num_elems_i, num_elems_j);

                    // store the point IDs for this elem where the range is
                    // (i:i+1, j:j+1, k:k+1) for a linear hexahedron
                    int this_point = 0;
                    for (int jcount = j; jcount <= j + 1; jcount++) {
                        for (int icount = i; icount <= i + 1; icount++) {
                            // global id for the points
                            int node_gid = get_id(icount, jcount, 0,
                                            num_points_i, num_points_j);

                            // convert this_point index to the FE index convention
                            int this_index = convert_point_number_in_quad(this_point);

                            // store the points in this elem according the the finite
                            // element numbering convention
                            global_nodes_in_elem(elem_gid, this_index) = node_gid;

                            // increment the point counting index
                            this_point = this_point + 1;
                        } // end for icount
                    } // end for jcount
                } // end for i
            } // end for j
        }

        //distribute partitioned data from the global mesh build data on rank 0
        size_t num_local_nodes;
        DistributedMap node_map;
        // read coords
        read_index_start = 0;
        
        buffer_iterations = global_num_nodes / BUFFER_SIZE;
        if (global_num_nodes % BUFFER_SIZE != 0)
        {
            buffer_iterations++;
        }

        read_buffer = CArrayKokkos<real_t, Kokkos::LayoutRight, Kokkos::HostSpace>(BUFFER_SIZE, num_dims);

        { //scoped so temp FArray data is auto deleted to save memory
            //allocate pre-partition node coords using contiguous decomposition
            //FArray type used since CArray type still doesnt support zoltan2 decomposition
            DistributedDFArray<real_t> node_coords_distributed(global_num_nodes, num_dims);

            // construct contiguous parallel row map now that we know the number of nodes
            node_map = node_coords_distributed.pmap;
            // map->describe(*fos,Teuchos::VERB_EXTREME);

            // set the vertices in the mesh read in
            num_local_nodes = node_map.size();
            for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
            {
                // pack buffer on rank 0
                if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
                {
                    for (buffer_loop = 0; buffer_loop < BUFFER_SIZE; buffer_loop++)
                    {

                        for (int idim = 0; idim < num_dims; idim++)
                        {
                            // debug print
                            // std::cout<<" "<< substring <<std::endl;
                            // assign the substring variable as a word of the read buffer
                            read_buffer(buffer_loop, idim) = global_coords(buffer_iteration * BUFFER_SIZE + buffer_loop, idim);
                        }
                    }
                }
                else if (myrank == 0)
                {
                    buffer_loop = 0;
                    while (buffer_iteration * BUFFER_SIZE + buffer_loop < global_num_nodes) {
                        for (int idim = 0; idim < num_dims; idim++)
                        {
                            // debug print
                            // std::cout<<" "<< substring <<std::endl;
                            // assign the substring variable as a word of the read buffer
                            read_buffer(buffer_loop, idim) = global_coords(buffer_iteration * BUFFER_SIZE + buffer_loop, idim);
                        }
                        buffer_loop++;
                    }
                }

                // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
                MPI_Bcast(read_buffer.pointer(), BUFFER_SIZE * num_dims, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                // broadcast how many nodes were read into this buffer iteration
                MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, MPI_COMM_WORLD);

                // debug_print
                // std::cout << "NODE BUFFER LOOP IS: " << buffer_loop << std::endl;
                // for(int iprint=0; iprint < buffer_loop; iprint++)
                // std::cout<<"buffer packing: " << std::string(&read_buffer(iprint,0,0)) << std::endl;
                // return;

                // determine which data to store in the swage mesh members (the local node data)
                // loop through read buffer
                for (scan_loop = 0; scan_loop < buffer_loop; scan_loop++)
                {
                    // set global node id (ensight specific order)
                    node_gid = read_index_start + scan_loop;
                    // let map decide if this node id belongs locally; if yes store data
                    if (node_map.isProcessGlobalIndex(node_gid))
                    {
                        // set local node index in this mpi rank
                        node_rid = node_map.getLocalIndex(node_gid);
                        // extract nodal position from the read buffer
                        // for tecplot format this is the three coords in the same line
                        dof_value = read_buffer(scan_loop,0);
                        node_coords_distributed.host(node_rid, 0) = dof_value * unit_scaling;
                        dof_value = read_buffer(scan_loop,1);
                        node_coords_distributed.host(node_rid, 1) = dof_value * unit_scaling;
                        if (num_dims == 3)
                        {
                            dof_value = read_buffer(scan_loop,2);
                            node_coords_distributed.host(node_rid, 2) = dof_value * unit_scaling;
                        }
                    }
                }
                read_index_start += BUFFER_SIZE;
            }
            // end of coordinate readin
            node_coords_distributed.update_device();
            // repartition node distribution
            node_coords_distributed.repartition_vector();
            //get map from repartitioned Farray and feed it into distributed CArray type; FArray data will be discared after scope
            std::vector<node_state> required_node_state = { node_state::coords };
            node_map = node_coords_distributed.pmap;
            node.initialize(node_map, num_dims, required_node_state);
            //copy coordinate data from repartitioned FArray into CArray
            FOR_ALL(node_id, 0, node_map.size(), {
                for(int idim = 0; idim < num_dims; idim++){
                    node.coords(node_id,idim) = node_coords_distributed(node_id,idim);
                }
            });
        }

        //initialize some mesh data
        mesh.initialize_nodes(global_num_nodes);
        num_local_nodes = node_map.size();
        mesh.num_local_nodes = num_local_nodes;
        mesh.node_map = node_map;
        mesh.num_dims = num_dims;
        //node.coords.print();
        
        // debug print of nodal data

        // debug print nodal positions and indices
        /*
        std::cout << " ------------NODAL POSITIONS ON TASK " << myrank << " --------------"<<std::endl;
        for (int inode = 0; inode < local_nrows; inode++){
            std::cout << "node: " << map->getGlobalElement(inode) + 1 << " { ";
        for (int istride = 0; istride < num_dims; istride++){
            std::cout << node_coords(inode,istride) << " , ";
        }
        std::cout << " }"<< std::endl;
        }
        */

        // check that local assignments match global total

        // read in element info (ensight file format is organized in element type sections)
        // loop over this later for several element type sections
        CArrayKokkos<int, Kokkos::LayoutRight, Kokkos::HostSpace> node_store(num_nodes_in_elem);

        // broadcast number of elements
        MPI_Bcast(&global_num_elems, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

        //initialize num elem in mesh struct

        if (myrank == 0)
        {
            std::cout << "before mesh initialization" << std::endl;
        }

        // read in element connectivity
        // we're gonna reallocate for the words per line expected for the element connectivity
        read_buffer_edof = CArrayKokkos<long long int, Kokkos::LayoutRight, Kokkos::HostSpace>(BUFFER_SIZE, num_nodes_in_elem);

        // calculate buffer iterations to read number of lines
        buffer_iterations = global_num_elems / BUFFER_SIZE;
        int assign_flag;

        // dynamic buffer used to store elements before we know how many this rank needs
        std::vector<size_t> element_temp(BUFFER_SIZE * num_nodes_in_elem);
        std::vector<size_t> global_indices_temp(BUFFER_SIZE);
        size_t buffer_max = BUFFER_SIZE * num_nodes_in_elem;
        size_t indices_buffer_max = BUFFER_SIZE;
        size_t num_elems = 0;

        if (global_num_elems % BUFFER_SIZE != 0)
        {
            buffer_iterations++;
        }
        read_index_start = 0;
        // std::cout << "ELEMENT BUFFER ITERATIONS: " << buffer_iterations << std::endl;
        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
        {
            // pack buffer on rank 0
            if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
            {
                for (buffer_loop = 0; buffer_loop < BUFFER_SIZE; buffer_loop++)
                {
                    for (int inode = 0; inode < num_nodes_in_elem; inode++)
                    {
                        read_buffer_edof(buffer_loop,inode) = global_nodes_in_elem(buffer_iteration * BUFFER_SIZE + buffer_loop, inode);
                    }
                    // std::cout <<std::endl;
                }
            }
            else if (myrank == 0)
            {
                buffer_loop = 0;
                while (buffer_iteration * BUFFER_SIZE + buffer_loop < global_num_elems) {
                    for (int inode = 0; inode < num_nodes_in_elem; inode++)
                    {
                        read_buffer_edof(buffer_loop,inode) = global_nodes_in_elem(buffer_iteration * BUFFER_SIZE + buffer_loop, inode);
                    }
                    // std::cout <<std::endl;
                    buffer_loop++;
                    // std::cout<<" "<< node_coords(node_gid, 0)<<std::endl;
                }
            }

            // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
            MPI_Bcast(read_buffer_edof.pointer(), BUFFER_SIZE * num_nodes_in_elem, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // store element connectivity that belongs to this rank
            // loop through read buffer
            for (scan_loop = 0; scan_loop < buffer_loop; scan_loop++)
            {
                // set global node id (ensight specific order)
                elem_gid = read_index_start + scan_loop;
                // add this element to the local list if any of its nodes belong to this rank according to the map
                // get list of nodes for each element line and check if they belong to the map
                assign_flag = 0;
                for (int inode = 0; inode < num_nodes_in_elem; inode++)
                {
                    // as we loop through the nodes belonging to this element we store them
                    // if any of these nodes belongs to this rank this list is used to store the element locally
                    node_gid = read_buffer_edof(scan_loop, inode);
                    node_store(inode) = node_gid;
                    // first we add the elements to a dynamically allocated list
                    if (node_map.isProcessGlobalIndex(node_gid) && !assign_flag)
                    {
                        assign_flag = 1;
                        num_elems++;
                    }
                }

                if (assign_flag)
                {
                    for (int inode = 0; inode < num_nodes_in_elem; inode++)
                    {
                        if ((num_elems - 1) * num_nodes_in_elem + inode >= buffer_max)
                        {
                            element_temp.resize((num_elems - 1) * num_nodes_in_elem + inode + BUFFER_SIZE * num_nodes_in_elem);
                            buffer_max = (num_elems - 1) * num_nodes_in_elem + inode + BUFFER_SIZE * num_nodes_in_elem;
                        }
                        element_temp[(num_elems - 1) * num_nodes_in_elem + inode] = node_store(inode);
                        // std::cout << "VECTOR STORAGE FOR ELEM " << num_elems << " ON TASK " << myrank << " NODE " << inode+1 << " IS " << node_store(inode) + 1 << std::endl;
                    }
                    // assign global element id to temporary list
                    if (num_elems - 1 >= indices_buffer_max)
                    {
                        global_indices_temp.resize(num_elems - 1 + BUFFER_SIZE);
                        indices_buffer_max = num_elems - 1 + BUFFER_SIZE;
                    }
                    global_indices_temp[num_elems - 1] = elem_gid;
                }
            }
            read_index_start += BUFFER_SIZE;
        }
        //set global and local shared element counts
        mesh.global_num_elems = global_num_elems;

        // construct partition mapping for shared elements on each process
        DCArrayKokkos<long long int> All_Element_Global_Indices(num_elems);
        // copy temporary global indices storage to view storage
        for (int ielem = 0; ielem < num_elems; ielem++)
        {
            All_Element_Global_Indices.host(ielem) = global_indices_temp[ielem];
        }

        All_Element_Global_Indices.update_device();

        // construct global map of local and shared elements (since different ranks can own the same elements due to the local node map)
        DistributedMap element_map = DistributedMap(All_Element_Global_Indices);

        //initialize elem data structures
        mesh.initialize_elems(num_elems, num_nodes_in_elem, element_map);

        // copy temporary element storage to distributed storage
        DistributedDCArray<size_t> nodes_in_elem = mesh.nodes_in_elem;

        for (int ielem = 0; ielem < num_elems; ielem++)
        {
            for (int inode = 0; inode < num_nodes_in_elem; inode++)
            {   //assign local indices to element-node connectivity (stores global indices until ghost maps are made later)
                nodes_in_elem.host(ielem, inode) = element_temp[ielem * num_nodes_in_elem + inode];
            }
        }

        nodes_in_elem.update_device();
        
        // delete temporary element connectivity and index storage
        //std::vector<size_t>().swap(element_temp);
        //std::vector<size_t>().swap(global_indices_temp);

        // initialize corner variables
        size_t num_corners = num_elems * num_nodes_in_elem;
        mesh.initialize_corners(num_corners);

    } // end build_2d_polar

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
        SimulationParameters_t& SimulationParameters) const
    {   
        int myrank, nranks;
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        MPI_Comm_size(MPI_COMM_WORLD,&nranks);
        /*currently we just build the global mesh data on rank 0 and then broadcast relevant data to each rank
          before the global mesh data on rank 0 falls out of scope*/
        int global_num_nodes, global_num_elems;
        CArrayKokkos<real_t, Kokkos::LayoutRight, Kokkos::HostSpace> read_buffer;
        CArrayKokkos<long long int, Kokkos::LayoutRight, Kokkos::HostSpace> read_buffer_edof;
        CArrayKokkos<double, Kokkos::LayoutLeft, Kokkos::HostSpace> global_coords;
        CArrayKokkos<double, Kokkos::LayoutLeft, Kokkos::HostSpace> global_nodes_in_elem;

        int local_node_index, current_column_index;
        int buffer_loop, buffer_iteration, buffer_iterations, dof_limit, scan_loop;
        int negative_index_found = 0;
        int global_negative_index_found = 0;

        size_t read_index_start, node_rid, elem_gid;
        long long int node_gid;
        real_t dof_value;
        real_t unit_scaling = 1;

        const int num_dims = 3;
        size_t num_nodes_in_elem = 1;
        for (int dim = 0; dim < num_dims; dim++) {
            num_nodes_in_elem *= 2;
        }
        if(myrank==0){
            printf("Creating a 3D box mesh \n");
        }

        // SimulationParameters.mesh_input.length.update_host();
        const double lx = SimulationParameters.mesh_input.length[0];
        const double ly = SimulationParameters.mesh_input.length[1];
        const double lz = SimulationParameters.mesh_input.length[2];

        // SimulationParameters.mesh_input.num_elems.update_host();
        const int num_elems_i = SimulationParameters.mesh_input.num_elems[0];
        const int num_elems_j = SimulationParameters.mesh_input.num_elems[1];
        const int num_elems_k = SimulationParameters.mesh_input.num_elems[2];

        const int num_points_i = num_elems_i + 1; // num points in x
        const int num_points_j = num_elems_j + 1; // num points in y
        const int num_points_k = num_elems_k + 1; // num points in y

        global_num_nodes = num_points_i * num_points_j * num_points_k;

        const double dx = lx / ((double)num_elems_i);  // len/(num_elems_i)
        const double dy = ly / ((double)num_elems_j);  // len/(num_elems_j)
        const double dz = lz / ((double)num_elems_k);  // len/(num_elems_k)

        global_num_elems = num_elems_i * num_elems_j * num_elems_k;
            
        if(myrank==0){
            std::vector<double> origin(num_dims);
            // SimulationParameters.mesh_input.origin.update_host();
            for (int i = 0; i < num_dims; i++) { origin[i] = SimulationParameters.mesh_input.origin[i]; }

            // --- 3D parameters ---
            // const int num_faces_in_elem  = 6;  // number of faces in elem
            // const int num_points_in_elem = 8;  // number of points in elem
            // const int num_points_in_face = 4;  // number of points in a face
            // const int num_edges_in_elem  = 12; // number of edges in a elem
            
            // node coords data on rank 0 for all global nodes
            global_coords = CArrayKokkos<double, Kokkos::LayoutLeft, Kokkos::HostSpace>(global_num_nodes, num_dims, "global_mesh_build_node_coordinates");

            // --- Build nodes ---

            // populate the point data structures
            for (int k = 0; k < num_points_k; k++) {
                for (int j = 0; j < num_points_j; j++) {
                    for (int i = 0; i < num_points_i; i++) {
                        // global id for the point
                        int node_gid = get_id(i, j, k, num_points_i, num_points_j);

                        // store the point coordinates
                        global_coords(node_gid, 0) = origin[0] + (double)i * dx;
                        global_coords(node_gid, 1) = origin[1] + (double)j * dy;
                        global_coords(node_gid, 2) = origin[2] + (double)k * dz;
                    } // end for i
                } // end for j
            } // end for k


            // initialize elem variables
            global_nodes_in_elem = CArrayKokkos<double, Kokkos::LayoutLeft, Kokkos::HostSpace>(global_num_elems, num_nodes_in_elem, "global_mesh_build_nodes_in_elem");

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
                                    global_nodes_in_elem(elem_gid, this_index) = node_gid;

                                    // increment the point counting index
                                    this_point = this_point + 1;
                                } // end for icount
                            } // end for jcount
                        }  // end for kcount
                    } // end for i
                } // end for j
            } // end for k
        }

        //distribute partitioned data from the global mesh build data on rank 0
        size_t num_local_nodes;
        DistributedMap node_map;
        // read coords
        read_index_start = 0;
        
        buffer_iterations = global_num_nodes / BUFFER_SIZE;
        if (global_num_nodes % BUFFER_SIZE != 0)
        {
            buffer_iterations++;
        }

        read_buffer = CArrayKokkos<real_t, Kokkos::LayoutRight, Kokkos::HostSpace>(BUFFER_SIZE, num_dims);

        { //scoped so temp FArray data is auto deleted to save memory
            //allocate pre-partition node coords using contiguous decomposition
            //FArray type used since CArray type still doesnt support zoltan2 decomposition
            DistributedDFArray<real_t> node_coords_distributed(global_num_nodes, num_dims);

            // construct contiguous parallel row map now that we know the number of nodes
            node_map = node_coords_distributed.pmap;
            // map->describe(*fos,Teuchos::VERB_EXTREME);

            // set the vertices in the mesh read in
            num_local_nodes = node_map.size();
            for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
            {
                // pack buffer on rank 0
                if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
                {
                    for (buffer_loop = 0; buffer_loop < BUFFER_SIZE; buffer_loop++)
                    {

                        for (int idim = 0; idim < num_dims; idim++)
                        {
                            // debug print
                            // std::cout<<" "<< substring <<std::endl;
                            // assign the substring variable as a word of the read buffer
                            read_buffer(buffer_loop, idim) = global_coords(buffer_iteration * BUFFER_SIZE + buffer_loop, idim);
                        }
                    }
                }
                else if (myrank == 0)
                {
                    buffer_loop = 0;
                    while (buffer_iteration * BUFFER_SIZE + buffer_loop < global_num_nodes) {
                        for (int idim = 0; idim < num_dims; idim++)
                        {
                            // debug print
                            // std::cout<<" "<< substring <<std::endl;
                            // assign the substring variable as a word of the read buffer
                            read_buffer(buffer_loop, idim) = global_coords(buffer_iteration * BUFFER_SIZE + buffer_loop, idim);
                        }
                        buffer_loop++;
                    }
                }

                // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
                MPI_Bcast(read_buffer.pointer(), BUFFER_SIZE * num_dims, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                // broadcast how many nodes were read into this buffer iteration
                MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, MPI_COMM_WORLD);

                // debug_print
                // std::cout << "NODE BUFFER LOOP IS: " << buffer_loop << std::endl;
                // for(int iprint=0; iprint < buffer_loop; iprint++)
                // std::cout<<"buffer packing: " << std::string(&read_buffer(iprint,0,0)) << std::endl;
                // return;

                // determine which data to store in the swage mesh members (the local node data)
                // loop through read buffer
                for (scan_loop = 0; scan_loop < buffer_loop; scan_loop++)
                {
                    // set global node id (ensight specific order)
                    node_gid = read_index_start + scan_loop;
                    // let map decide if this node id belongs locally; if yes store data
                    if (node_map.isProcessGlobalIndex(node_gid))
                    {
                        // set local node index in this mpi rank
                        node_rid = node_map.getLocalIndex(node_gid);
                        // extract nodal position from the read buffer
                        // for tecplot format this is the three coords in the same line
                        dof_value = read_buffer(scan_loop,0);
                        node_coords_distributed.host(node_rid, 0) = dof_value * unit_scaling;
                        dof_value = read_buffer(scan_loop,1);
                        node_coords_distributed.host(node_rid, 1) = dof_value * unit_scaling;
                        if (num_dims == 3)
                        {
                            dof_value = read_buffer(scan_loop,2);
                            node_coords_distributed.host(node_rid, 2) = dof_value * unit_scaling;
                        }
                    }
                }
                read_index_start += BUFFER_SIZE;
            }
            // end of coordinate readin
            node_coords_distributed.update_device();
            // repartition node distribution
            node_coords_distributed.repartition_vector();
            //get map from repartitioned Farray and feed it into distributed CArray type; FArray data will be discared after scope
            std::vector<node_state> required_node_state = { node_state::coords };
            node_map = node_coords_distributed.pmap;
            node.initialize(node_map, num_dims, required_node_state);
            //copy coordinate data from repartitioned FArray into CArray
            FOR_ALL(node_id, 0, node_map.size(), {
                for(int idim = 0; idim < num_dims; idim++){
                    node.coords(node_id,idim) = node_coords_distributed(node_id,idim);
                }
            });
        }

        //initialize some mesh data
        mesh.initialize_nodes(global_num_nodes);
        num_local_nodes = node_map.size();
        mesh.num_local_nodes = num_local_nodes;
        mesh.node_map = node_map;
        //node.coords.print();
        
        // debug print of nodal data

        // debug print nodal positions and indices
        /*
        std::cout << " ------------NODAL POSITIONS ON TASK " << myrank << " --------------"<<std::endl;
        for (int inode = 0; inode < local_nrows; inode++){
            std::cout << "node: " << map->getGlobalElement(inode) + 1 << " { ";
        for (int istride = 0; istride < num_dims; istride++){
            std::cout << node_coords(inode,istride) << " , ";
        }
        std::cout << " }"<< std::endl;
        }
        */

        // check that local assignments match global total

        // read in element info (ensight file format is organized in element type sections)
        // loop over this later for several element type sections
        CArrayKokkos<int, Kokkos::LayoutRight, Kokkos::HostSpace> node_store(num_nodes_in_elem);

        // broadcast number of elements
        MPI_Bcast(&global_num_elems, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

        //initialize num elem in mesh struct

        if (myrank == 0)
        {
            std::cout << "before mesh initialization" << std::endl;
        }

        // read in element connectivity
        // we're gonna reallocate for the words per line expected for the element connectivity
        read_buffer_edof = CArrayKokkos<long long int, Kokkos::LayoutRight, Kokkos::HostSpace>(BUFFER_SIZE, num_nodes_in_elem);

        // calculate buffer iterations to read number of lines
        buffer_iterations = global_num_elems / BUFFER_SIZE;
        int assign_flag;

        // dynamic buffer used to store elements before we know how many this rank needs
        std::vector<size_t> element_temp(BUFFER_SIZE * num_nodes_in_elem);
        std::vector<size_t> global_indices_temp(BUFFER_SIZE);
        size_t buffer_max = BUFFER_SIZE * num_nodes_in_elem;
        size_t indices_buffer_max = BUFFER_SIZE;
        size_t num_elems = 0;

        if (global_num_elems % BUFFER_SIZE != 0)
        {
            buffer_iterations++;
        }
        read_index_start = 0;
        // std::cout << "ELEMENT BUFFER ITERATIONS: " << buffer_iterations << std::endl;
        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
        {
            // pack buffer on rank 0
            if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
            {
                for (buffer_loop = 0; buffer_loop < BUFFER_SIZE; buffer_loop++)
                {
                    for (int inode = 0; inode < num_nodes_in_elem; inode++)
                    {
                        read_buffer_edof(buffer_loop,inode) = global_nodes_in_elem(buffer_iteration * BUFFER_SIZE + buffer_loop, inode);
                    }
                    // std::cout <<std::endl;
                }
            }
            else if (myrank == 0)
            {
                buffer_loop = 0;
                while (buffer_iteration * BUFFER_SIZE + buffer_loop < global_num_elems) {
                    for (int inode = 0; inode < num_nodes_in_elem; inode++)
                    {
                        read_buffer_edof(buffer_loop,inode) = global_nodes_in_elem(buffer_iteration * BUFFER_SIZE + buffer_loop, inode);
                    }
                    // std::cout <<std::endl;
                    buffer_loop++;
                    // std::cout<<" "<< node_coords(node_gid, 0)<<std::endl;
                }
            }

            // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
            MPI_Bcast(read_buffer_edof.pointer(), BUFFER_SIZE * num_nodes_in_elem, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // store element connectivity that belongs to this rank
            // loop through read buffer
            for (scan_loop = 0; scan_loop < buffer_loop; scan_loop++)
            {
                // set global node id (ensight specific order)
                elem_gid = read_index_start + scan_loop;
                // add this element to the local list if any of its nodes belong to this rank according to the map
                // get list of nodes for each element line and check if they belong to the map
                assign_flag = 0;
                for (int inode = 0; inode < num_nodes_in_elem; inode++)
                {
                    // as we loop through the nodes belonging to this element we store them
                    // if any of these nodes belongs to this rank this list is used to store the element locally
                    node_gid = read_buffer_edof(scan_loop, inode);
                    node_store(inode) = node_gid;
                    // first we add the elements to a dynamically allocated list
                    if (node_map.isProcessGlobalIndex(node_gid) && !assign_flag)
                    {
                        assign_flag = 1;
                        num_elems++;
                    }
                }

                if (assign_flag)
                {
                    for (int inode = 0; inode < num_nodes_in_elem; inode++)
                    {
                        if ((num_elems - 1) * num_nodes_in_elem + inode >= buffer_max)
                        {
                            element_temp.resize((num_elems - 1) * num_nodes_in_elem + inode + BUFFER_SIZE * num_nodes_in_elem);
                            buffer_max = (num_elems - 1) * num_nodes_in_elem + inode + BUFFER_SIZE * num_nodes_in_elem;
                        }
                        element_temp[(num_elems - 1) * num_nodes_in_elem + inode] = node_store(inode);
                        // std::cout << "VECTOR STORAGE FOR ELEM " << num_elems << " ON TASK " << myrank << " NODE " << inode+1 << " IS " << node_store(inode) + 1 << std::endl;
                    }
                    // assign global element id to temporary list
                    if (num_elems - 1 >= indices_buffer_max)
                    {
                        global_indices_temp.resize(num_elems - 1 + BUFFER_SIZE);
                        indices_buffer_max = num_elems - 1 + BUFFER_SIZE;
                    }
                    global_indices_temp[num_elems - 1] = elem_gid;
                }
            }
            read_index_start += BUFFER_SIZE;
        }
        //set global and local shared element counts
        mesh.global_num_elems = global_num_elems;

        // construct partition mapping for shared elements on each process
        DCArrayKokkos<long long int> All_Element_Global_Indices(num_elems);
        // copy temporary global indices storage to view storage
        for (int ielem = 0; ielem < num_elems; ielem++)
        {
            All_Element_Global_Indices.host(ielem) = global_indices_temp[ielem];
        }

        All_Element_Global_Indices.update_device();

        // construct global map of local and shared elements (since different ranks can own the same elements due to the local node map)
        DistributedMap element_map = DistributedMap(All_Element_Global_Indices);

        //initialize elem data structures
        mesh.initialize_elems(num_elems, num_nodes_in_elem, element_map);

        // copy temporary element storage to distributed storage
        DistributedDCArray<size_t> nodes_in_elem = mesh.nodes_in_elem;

        for (int ielem = 0; ielem < num_elems; ielem++)
        {
            for (int inode = 0; inode < num_nodes_in_elem; inode++)
            {   //assign local indices to element-node connectivity (stores global indices until ghost maps are made later)
                nodes_in_elem.host(ielem, inode) = element_temp[ielem * num_nodes_in_elem + inode];
            }
        }

        nodes_in_elem.update_device();
        
        // delete temporary element connectivity and index storage
        //std::vector<size_t>().swap(element_temp);
        //std::vector<size_t>().swap(global_indices_temp);

        // initialize corner variables
        size_t num_corners = num_elems * num_nodes_in_elem;
        mesh.initialize_corners(num_corners);
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
    // void build_3d_HexN_box(Mesh_t& mesh,
    //     GaussPoint_t& GaussPoints,
    //     node_t&   node,
    //     corner_t& corner,
    //     SimulationParameters_t& SimulationParameters) const
    // {
    //     printf(" ***** WARNING::  build_3d_HexN_box not yet implemented\n");
    //     const int num_dim = 3;

    //     // SimulationParameters.mesh_input.length.update_host();
    //     const double lx = SimulationParameters.mesh_input.length[0];
    //     const double ly = SimulationParameters.mesh_input.length[1];
    //     const double lz = SimulationParameters.mesh_input.length[2];

    //     // SimulationParameters.mesh_input.num_elems.update_host();
    //     const int num_elems_i = SimulationParameters.mesh_input.num_elems[0];
    //     const int num_elems_j = SimulationParameters.mesh_input.num_elems[1];
    //     const int num_elems_k = SimulationParameters.mesh_input.num_elems[2];

    //     // creating zones for the Pn order
    //     const int Pn_order = SimulationParameters.mesh_input.p_order;
        
    //     if (Pn_order > 19) {
    //         printf("Fierro DG and RD solvers are only valid for elements up to Pn = 19 \n");
    //         return;
    //     }

    //     const int num_zones_i = Pn_order*num_elems_i;
    //     const int num_zones_j = Pn_order*num_elems_j;
    //     const int num_zones_k = Pn_order*num_elems_k;
        
    //     const int num_points_i = num_zones_i+1; // num points in x accounting for Pn
    //     const int num_points_j = num_zones_j+1; // num points in y accounting for Pn
    //     const int num_points_k = num_zones_k+1; // num points in y accounting for Pn
        
        
    //     const double dx = lx/((double)num_zones_i);  // len/(num_zones_i)
    //     const double dy = ly/((double)num_zones_j);  // len/(num_zones_j)
    //     const double dz = lz/((double)num_zones_k);  // len/(num_zones_k)
        
    //     const int num_elems = num_elems_i*num_elems_j*num_elems_k;
    //     // const int num_zones = num_zones_i*num_zones_j*num_zones_k; // accounts for Pn

    //     std::vector<double> origin(num_dim);
    //     for (int i = 0; i < num_dim; i++) { origin[i] = SimulationParameters.mesh_input.origin[i]; }

    //     // --- 3D parameters ---
    //     // const int num_faces_in_zone = 6;   // number of faces in zone
    //     // const int num_points_in_zone = 8;  // number of points in zone
    //     // const int num_points_in_face = 4;  // number of points in a face
        
    //     // p_order   = 1, 2, 3, 4, 5
    //     // num_nodes = 2, 3, 4, 5, 6
    //     const int num_1D_points = Pn_order+1;
    //     const int num_points_in_elem = num_1D_points*num_1D_points*num_1D_points;
           
        
    //     // --- elem ---
    //     auto elem_coords = CArray <double> (num_elems, num_dim);
    //     auto elem_point_list = CArray <int> (num_elems, num_points_in_elem);
        
        
    //     // --- point ---
    //     int num_points = num_points_i * num_points_j * num_points_k;
    //     auto pt_coords = CArray <double> (num_points, num_dim);


    //     // --- Build nodes ---
        
    //     // initialize node variables
    //     mesh.initialize_nodes(num_points);

    //     // 
    //     std::vector<node_state> required_node_state = { node_state::coords };
    //     node.initialize(num_points, num_dim, required_node_state);
    //     // populate the point data structures
    //     for (int k = 0; k < num_points_k; k++){
    //         for (int j = 0; j < num_points_j; j++){
    //             for (int i = 0; i < num_points_i; i++){

                
    //                 // global id for the point
    //                 int node_gid = get_id(i, j, k, num_points_i, num_points_j);

    //                 // store the point coordinates
    //                 node.coords.host(node_gid, 0) = origin[0] + (double)i * dx;
    //                 node.coords.host(node_gid, 1) = origin[1] + (double)j * dy;
    //                 node.coords.host(node_gid, 2) = origin[2] + (double)k * dz;
                    
    //             } // end for k
    //         } // end for i
    //     } // end for j


    //     node.coords.update_device();


    //     // initialize elem variables
    //     mesh.initialize_elems(num_elems, num_dim);

    //     // --- Build elems  ---
        
    //     // populate the elem center data structures accounting for Pn
    //     for (int k=0; k<num_elems_k; k++){
    //         for (int j=0; j<num_elems_j; j++){
    //             for (int i=0; i<num_elems_i; i++){
                  
    //                 // global id for the elem
    //                 size_t elem_gid = get_id(i, j, k, num_elems_i, num_elems_j);
                    
    //                 // store the point IDs for this elem where the range is
    //                 // (i:i+1, j:j+1, k:k+1) for a linear hexahedron
    //                 // (i:(i+1)*Pn_order, j:(j+1)*Pn_order, k:(k+1)*Pn_order) for a Pn hexahedron
    //                 int node_lid = 0;
                    
    //                 int k_local = 0;
    //                 for (int kcount=k*Pn_order; kcount<=(k+1)*Pn_order; kcount++){
                        
    //                     int j_local = 0;
    //                     for (int jcount=j*Pn_order; jcount<=(j+1)*Pn_order; jcount++){
                            
    //                         int i_local = 0;
    //                         for (int icount=i*Pn_order; icount<=(i+1)*Pn_order; icount++){
                                
    //                             // global id for the points
    //                             size_t node_gid = get_id(icount, jcount, kcount,
    //                                               num_points_i, num_points_j);

    //                             // Saved using i,j,k indexing
    //                             mesh.nodes_in_elem.host(elem_gid, node_lid) = node_gid;
                                
    //                             // increment the point counting index
    //                             node_lid = node_lid + 1;
                                
    //                             i_local++;
    //                         } // end for icount
                            
    //                         j_local++;
    //                     } // end for jcount
                        
    //                     k_local ++;
    //                 }  // end for kcount
    //             } // end for i
    //         } // end for j
    //     } // end for k

    //     // update device side
    //     mesh.nodes_in_elem.update_device();

    //     // initialize corner variables
    //     int num_corners = num_elems * mesh.num_nodes_in_elem;
    //     mesh.initialize_corners(num_corners);
    //     // corner.initialize(num_corners, num_dim);

    // }
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
        SimulationParameters_t& SimulationParameters,
        double dt,
        double time_value,
        CArray<double> graphics_times,
        std::vector<node_state> node_states,
        std::vector<gauss_pt_state> gauss_pt_states,
        std::vector<material_pt_state> material_pt_states,
        const size_t solver_id)
    {

        // node_state is an enum for possible fields (e.g., coords, velocity, etc.), see state.h
        // gauss_pt_state is an enum for possible fields (e.g., vol, divergence, etc.)
        // material_pt_state is an enum for possible fields (e.g., den, pres, etc.)


        // *******************
        //  Update host 
        // *******************

        const size_t num_mats = State.MaterialPoints.num_material_points.size();

        // material point values
            
        //  Update host data for mat_pt state
        for (auto field : material_pt_states){
            switch(field){
                // scalar vars to write out
                case material_pt_state::density:
                    State.MaterialPoints.den.update_host();
                    break;
                case material_pt_state::pressure:
                    State.MaterialPoints.pres.update_host();
                    break;
                case material_pt_state::specific_internal_energy:
                    State.MaterialPoints.sie.update_host();
                    break;
                case material_pt_state::sound_speed:
                    State.MaterialPoints.sspd.update_host();
                    break;
                case material_pt_state::mass:
                    State.MaterialPoints.mass.update_host();
                    break;
                case material_pt_state::volume_fraction:
                    State.MaterialPoints.volfrac.update_host();
                    State.MaterialPoints.geo_volfrac.update_host();
                    break;
                case material_pt_state::eroded_flag:
                    State.MaterialPoints.eroded.update_host();
                    break;
                // tensor vars to write out
                case material_pt_state::stress:
                    State.MaterialPoints.stress.update_host();
                    break;
                
                // additional vars for thermal-mechanical solver
                case material_pt_state::thermal_conductivity:
                    State.MaterialPoints.conductivity.update_host();
                    break;
                
                case material_pt_state::specific_heat:
                    State.MaterialPoints.specific_heat.update_host();
                    break;

                // add other variables here
                
                // not used
                case material_pt_state::elastic_modulii:
                    break;
                case material_pt_state::shear_modulii:
                    break;
                case material_pt_state::poisson_ratios:
                    break;
                case material_pt_state::heat_flux:
                    break;
                default:
                    std::cout<<"Desired material point state not understood in outputs"<<std::endl;
            } // end switch
        } // end for over mat_pt_states



        // update gauss point values
        for (auto field : gauss_pt_states){
            switch(field){
                // scalar vars to write out
                case gauss_pt_state::volume:
                    State.GaussPoints.vol.update_host();
                    break;
                case gauss_pt_state::divergence_velocity:
                    State.GaussPoints.div.update_host();
                    break;
                case gauss_pt_state::level_set:
                    State.GaussPoints.level_set.update_host();
                    break;      

                // tensor vars to write out
                case gauss_pt_state::gradient_velocity:
                    State.GaussPoints.vel_grad.update_host();
                    break;
                default:
                    std::cout<<"Desired Gauss point state not understood in vtk outputs"<<std::endl;

            } // end switch
        } // end loop

        // nodal values
        for (auto field : node_states){
            switch(field){
                case node_state::mass:
                    State.node.mass.update_host();
                    break;
                case node_state::temp:
                    State.node.temp.update_host();
                    break;
                case node_state::coords:
                    State.node.coords.update_host();
                    break;
                case node_state::velocity:
                    State.node.vel.update_host();
                    break;
                case node_state::gradient_level_set:
                    State.node.gradient_level_set.update_host();
                    break;  

                case node_state::force:
                    break;

                // heat transer vars
                case node_state::heat_transfer:
                    break;

            } // end switch
        } // end for over 
        Kokkos::fence();


        // ******************************************
        //  Build Material and Element state outputs
        // ******************************************

        size_t num_mat_pt_scalar_vars = 0;
        size_t num_mat_pt_tensor_vars = 0;
            
        // count the number of material point state vars to write out
        for (auto field : SimulationParameters.output_options.output_mat_pt_state){
            switch(field){
                // scalar vars to write out
                case material_pt_state::density:
                    num_mat_pt_scalar_vars ++;
                    break;
                case material_pt_state::pressure:
                    num_mat_pt_scalar_vars ++;
                    break;
                case material_pt_state::specific_internal_energy:
                    num_mat_pt_scalar_vars ++;
                    break;
                case material_pt_state::sound_speed:
                    num_mat_pt_scalar_vars ++;
                    break;
                case material_pt_state::mass:
                    num_mat_pt_scalar_vars ++;
                    break;
                case material_pt_state::volume_fraction:
                    num_mat_pt_scalar_vars ++; // mat volfrac
                    num_mat_pt_scalar_vars ++; // geometric volfrac
                    break;
                case material_pt_state::eroded_flag:
                    num_mat_pt_scalar_vars ++;
                    break;
                // tensor vars to write out
                case material_pt_state::stress:
                    num_mat_pt_tensor_vars ++;
                    break;
                
                // additional vars for thermal-mechanical solver
                case material_pt_state::thermal_conductivity:
                    num_mat_pt_scalar_vars ++;
                    break;
                
                case material_pt_state::specific_heat:
                    num_mat_pt_scalar_vars ++;
                    break;

                // add other variables here

                // not used
                case material_pt_state::elastic_modulii:
                    break;
                case material_pt_state::shear_modulii:
                    break;
                case material_pt_state::poisson_ratios:
                    break;
                case material_pt_state::heat_flux:
                    break;
                default:
                    std::cout<<"Desired material point state not understood in outputs"<<std::endl;
            } // end switch
        } // end for over mat_pt_states



        size_t num_elem_scalar_vars = 0;
        size_t num_elem_vector_vars = 0;
        size_t num_elem_tensor_vars = 0;

        // count the number of element average fields to write out
        for (auto field : SimulationParameters.output_options.output_elem_state){
            switch(field){
                // scalar vars to write out
                case material_pt_state::density:
                    num_elem_scalar_vars ++;
                    break;
                case material_pt_state::pressure:
                    num_elem_scalar_vars ++;
                    break;
                case material_pt_state::specific_internal_energy:
                    num_elem_scalar_vars ++;
                    break;
                case material_pt_state::sound_speed:
                    num_elem_scalar_vars ++;
                    break;
                case material_pt_state::mass:
                    num_elem_scalar_vars ++;
                    break;
                // tensor vars to write out
                case material_pt_state::stress:
                    num_elem_tensor_vars ++;
                    break;

                // additional vars for thermal-mechanical solver
                case material_pt_state::thermal_conductivity:
                    num_elem_scalar_vars ++;
                    break;
                
                case material_pt_state::specific_heat:
                    num_elem_scalar_vars ++;
                    break;

                // add other variables here

                // not used
                case material_pt_state::volume_fraction:
                    break;
                case material_pt_state::eroded_flag:
                    break;
                case material_pt_state::elastic_modulii:
                    break;
                case material_pt_state::shear_modulii:
                    break;
                case material_pt_state::poisson_ratios:
                    break;
                case material_pt_state::heat_flux:
                    break;
                default:
                    std::cout<<"Desired material point state not understood in outputs"<<std::endl;
            } // end switch
        } // end for over mat_pt_states


        size_t num_gauss_pt_scalar_vars = 0;
        size_t num_gauss_pt_tensor_vars = 0;

        // gauss point values to ouptput
        for (auto field : SimulationParameters.output_options.output_gauss_pt_state){
            switch(field){
                // scalar vars to write out
                case gauss_pt_state::volume:
                    num_gauss_pt_scalar_vars ++;
                    break;
                case gauss_pt_state::level_set:
                    num_gauss_pt_scalar_vars ++;
                    break;
                case gauss_pt_state::divergence_velocity:
                    num_gauss_pt_scalar_vars ++;
                    break;

                // tensor vars to write out
                case gauss_pt_state::gradient_velocity:
                    num_gauss_pt_tensor_vars ++;
                    break;
                default:
                    std::cout<<"Desired Gauss point state not understood in vtk outputs"<<std::endl;

            } // end switch
        } // end loop

        // add the Gauss point state to the element state
        num_elem_scalar_vars += num_gauss_pt_scalar_vars;
        num_elem_tensor_vars += num_gauss_pt_tensor_vars;


        // Scalar, vector, and tensor value names associated with a elem
        std::vector<std::string> elem_scalar_var_names(num_elem_scalar_vars);
        std::vector<std::string> elem_tensor_var_names(num_elem_tensor_vars);

        // Scalar, vector, and tensor values associated with a material in part elems
        std::vector<std::string> mat_elem_scalar_var_names(num_mat_pt_scalar_vars);
        std::vector<std::string> mat_elem_tensor_var_names(num_mat_pt_tensor_vars);


        // the ids to access a variable in the mat_scalar_var_name or tensor list
        int mat_den_id = -1;
        int mat_pres_id = -1;
        int mat_sie_id = -1;
        int mat_sspd_id = -1;
        int mat_mass_id = -1;
        int mat_volfrac_id = -1;  
        int mat_geo_volfrac_id = -1;  // geometric volume fraction of part
        int mat_eroded_id = -1;
        int mat_stress_id = -1;

        int mat_conductivity_id = -1;
        int mat_specific_heat_id = -1;

        // the index for the scalar, vector, and tensor fields
        size_t var = 0;
        size_t vector_var = 0;
        size_t tensor_var = 0;

        // material point state to output
        for (auto field : SimulationParameters.output_options.output_mat_pt_state){
            switch(field){
                // scalar vars
                case material_pt_state::density:
                    mat_elem_scalar_var_names[var] = "mat_den";
                    mat_den_id = var;
                    var++;
                    break;
                case material_pt_state::pressure:
                    mat_elem_scalar_var_names[var] = "mat_pres";
                    mat_pres_id = var;
                    var++;
                    break;
                case material_pt_state::specific_internal_energy:
                    mat_elem_scalar_var_names[var] = "mat_sie";
                    mat_sie_id = var;
                    var++;
                    break;
                case material_pt_state::sound_speed:
                    mat_elem_scalar_var_names[var] = "mat_sspd";
                    mat_sspd_id = var;
                    var++;
                    break;
                case material_pt_state::mass:
                    mat_elem_scalar_var_names[var] = "mat_mass";
                    mat_mass_id = var;
                    var++;
                    break;
                case material_pt_state::volume_fraction:
                    mat_elem_scalar_var_names[var] = "mat_volfrac";
                    mat_volfrac_id = var; 
                    var++;

                    mat_elem_scalar_var_names[var] = "mat_geo_volfrac";
                    mat_geo_volfrac_id = var; 
                    var++;
                    break;
                case material_pt_state::eroded_flag:
                    mat_elem_scalar_var_names[var] = "mat_eroded";
                    mat_eroded_id = var;
                    var++;
                    break;
                // tensor vars
                case material_pt_state::stress:
                    mat_elem_tensor_var_names[tensor_var] = "mat_stress";
                    mat_stress_id = tensor_var;
                    tensor_var++;
                    break;

    
                // additional vars for thermal-mechanical solver
                case material_pt_state::thermal_conductivity:
                    mat_elem_scalar_var_names[var] = "mat_thermal_K";
                    mat_conductivity_id = var;
                    var++;
                    break;
                
                case material_pt_state::specific_heat:
                    mat_elem_scalar_var_names[var] = "mat_Cp";
                    mat_specific_heat_id = var;
                    var++;
                    break;


                // add other variables here

                // not used
                case material_pt_state::elastic_modulii:
                    break;
                case material_pt_state::shear_modulii:
                    break;
                case material_pt_state::poisson_ratios:
                    break;
                case material_pt_state::heat_flux:
                    break;
            } // end switch
        } // end for over mat_pt_states


        // element average fields to output

        // the ids to access a variable in the elem_scalar_var_name or tensor list
        int den_id = -1;
        int pres_id = -1;
        int sie_id = -1;
        int sspd_id = -1;
        int mass_id = -1; 
        int stress_id = -1;

        int conductivity_id = -1;
        int specific_heat_id = -1;

        // reset the counters
        var = 0;
        vector_var = 0;
        tensor_var = 0;

        // element state to output
        for (auto field : SimulationParameters.output_options.output_elem_state){
            switch(field){
                // scalar vars
                case material_pt_state::density:
                    elem_scalar_var_names[var] = "den";
                    den_id = var;
                    var++;
                    break;
                case material_pt_state::pressure:
                    elem_scalar_var_names[var] = "pres";
                    pres_id = var;
                    var++;
                    break;
                case material_pt_state::specific_internal_energy:
                    elem_scalar_var_names[var] = "sie";
                    sie_id = var;
                    var++;
                    break;
                case material_pt_state::sound_speed:
                    elem_scalar_var_names[var] = "sspd";
                    sspd_id = var;
                    var++;
                    break;
                case material_pt_state::mass:
                    elem_scalar_var_names[var] = "mass";
                    mass_id = var;
                    var++;
                    break;
                // tensor vars
                case material_pt_state::stress:
                    elem_tensor_var_names[tensor_var] = "stress";
                    stress_id = tensor_var;
                    tensor_var++;
                    break;

                // heat transfer variables
                case material_pt_state::thermal_conductivity:
                    elem_scalar_var_names[var] = "thermal_K";
                    conductivity_id = var;
                    var++;
                    break;
                
                case material_pt_state::specific_heat:
                    elem_scalar_var_names[var] = "Cp";
                    specific_heat_id = var;
                    var++;
                    break;

                // add other variables here

                // not used
                case material_pt_state::volume_fraction:
                    break;
                case material_pt_state::eroded_flag:
                    break;
                case material_pt_state::elastic_modulii:
                    break;
                case material_pt_state::shear_modulii:
                    break;
                case material_pt_state::poisson_ratios:
                    break;
                case material_pt_state::heat_flux:
                    break;
            } // end switch
        } // end for over mat_pt_states

        // append Gauss point vars to the element arrays
        int vol_id = -1;
        int div_id = -1;
        int level_set_id = -1;
        int vel_grad_id = -1;
        

        for (auto field : SimulationParameters.output_options.output_gauss_pt_state){
            switch(field){
                // scalars
                case gauss_pt_state::volume:
                    elem_scalar_var_names[var] = "vol";
                    vol_id = var;
                    var++;
                    break;
                case gauss_pt_state::divergence_velocity:
                    elem_scalar_var_names[var] = "div";
                    div_id = var;
                    var++;
                    break;

                case gauss_pt_state::level_set:
                    elem_scalar_var_names[var] = "level_set";
                    level_set_id = var;
                    var++;
                    break;

                // tensors
                case gauss_pt_state::gradient_velocity:
                    elem_tensor_var_names[tensor_var] = "vel_grad";
                    vel_grad_id = tensor_var;
                    tensor_var++;
                    break;
            } // end switch
        } // end loop over gauss_pt_states


        // *******************
        //  nodal values
        // *******************

        size_t num_node_scalar_vars = 0;
        size_t num_node_vector_vars = 0;

        for (auto field : SimulationParameters.output_options.output_node_state){
            switch(field){
                // --- scalars
                case node_state::mass:
                    num_node_scalar_vars ++;
                    break;
                case node_state::temp:
                    num_node_scalar_vars ++;
                    break;
                // -- vectors
                case node_state::coords:
                    num_node_vector_vars ++;
                    break;
                case node_state::velocity:
                    num_node_vector_vars ++; // for velocity
                    num_node_vector_vars ++; // for acceleration
                    break;
                case node_state::gradient_level_set:
                    num_node_vector_vars ++;
                    break;                    
                case node_state::force:
                    break;
                
                // heat transer vars
                case node_state::heat_transfer:
                    break;
            } // end switch
        } // end for over 
        Kokkos::fence();


        // Scalar and vector values associated with a node
        std::vector<std::string> node_scalar_var_names(num_node_scalar_vars);
        std::vector<std::string> node_vector_var_names(num_node_vector_vars);

        int node_mass_id = -1;
        int node_vel_id = -1;
        int node_accel_id = -1;
        int node_coord_id = -1;
        int node_temp_id = -1;
        int node_grad_level_set_id = -1;

        // reset counters for node fields
        var = 0;
        vector_var = 0;
        tensor_var = 0;

        for (auto field : SimulationParameters.output_options.output_node_state){
            switch(field){
                // scalars
                case node_state::mass:
                    node_scalar_var_names[var] = "node_mass";
                    node_mass_id = var;
                    var++;
                    break;
                case node_state::temp:
                    node_scalar_var_names[var] = "node_temp";
                    node_temp_id = var;
                    var++;
                    break;

                // vector fields

                case node_state::coords:
                    node_vector_var_names[vector_var] = "node_coords";
                    node_coord_id = vector_var;
                    vector_var++;
                    break;

                case node_state::velocity:
                    node_vector_var_names[vector_var] = "node_vel";
                    node_vel_id = vector_var;
                    vector_var++;

                    node_vector_var_names[vector_var] = "node_accel";
                    node_accel_id = vector_var;
                    vector_var++;
                    break;

                case node_state::gradient_level_set:
                    node_vector_var_names[vector_var] = "node_grad_lvlset";
                    node_grad_level_set_id = vector_var;
                    vector_var++;
                    break;

                // -- not used vars
                case node_state::force:
                    break;

                // heat transer vars
                case node_state::heat_transfer:
                    break;

                // tensors

            } // end switch
        } // end for over 

        // ********************************
        //  Write the collective nodal and elem fields 
        // ********************************

        if (SimulationParameters.output_options.format == output_options::viz ||
            SimulationParameters.output_options.format == output_options::viz_and_state) {

            write_parallel_viz(mesh,
                                State,
                                SimulationParameters,
                                dt,
                                time_value,
                                graphics_times,
                                solver_id,
                                node_states,
                                gauss_pt_states,
                                material_pt_states,
                                num_elem_scalar_vars,
                                num_elem_tensor_vars,
                                num_node_scalar_vars,
                                num_node_vector_vars,
                                num_mat_pt_scalar_vars,
                                num_mat_pt_tensor_vars,
                                elem_scalar_var_names,
                                elem_tensor_var_names,
                                mat_elem_scalar_var_names,
                                mat_elem_tensor_var_names,
                                node_scalar_var_names,
                                node_vector_var_names,
                                den_id,
                                pres_id,
                                sie_id,
                                sspd_id,
                                mass_id,
                                stress_id,
                                vol_id,
                                div_id,
                                level_set_id,
                                vel_grad_id,
                                conductivity_id,
                                specific_heat_id,
                                node_mass_id,
                                node_vel_id,
                                node_accel_id,
                                node_coord_id,
                                node_grad_level_set_id,
                                node_temp_id,
                                mat_den_id,
                                mat_pres_id,
                                mat_sie_id,
                                mat_sspd_id,
                                mat_mass_id,
                                mat_volfrac_id,
                                mat_geo_volfrac_id,  
                                mat_eroded_id,
                                mat_stress_id,
                                mat_conductivity_id,
                                mat_specific_heat_id);
        } // end if viz paraview output is to be written


        // STATE
        if (SimulationParameters.output_options.format == output_options::state ||
            SimulationParameters.output_options.format == output_options::viz_and_state) {

            write_material_point_state(mesh,
                                      State,
                                      SimulationParameters,
                                      time_value,
                                      graphics_times,
                                      node_states,
                                      gauss_pt_states,
                                      material_pt_states);

        } // end if state is to be written


        // will drop ensight outputs in the near future
        if (SimulationParameters.output_options.format == output_options::ensight){
           write_ensight(mesh,
                         State,
                         SimulationParameters,
                         dt,
                         time_value,
                         graphics_times,
                         node_states,
                         gauss_pt_states,
                         material_pt_states);
        }

        return;

    } // end write_mesh

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn write_collective vtm
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
    void write_parallel_viz(Mesh_t& mesh,
        State_t& State,
        SimulationParameters_t& SimulationParameters,
        double dt,
        double time_value,
        CArray<double> graphics_times,
        const size_t solver_id,
        std::vector<node_state> node_states,
        std::vector<gauss_pt_state> gauss_pt_states,
        std::vector<material_pt_state> material_pt_states,
        size_t num_elem_scalar_vars,
        size_t num_elem_tensor_vars,
        size_t num_node_scalar_vars,
        size_t num_node_vector_vars,
        size_t num_mat_pt_scalar_vars,
        size_t num_mat_pt_tensor_vars,
        std::vector<std::string> elem_scalar_var_names,
        std::vector<std::string> elem_tensor_var_names,
        std::vector<std::string> mat_elem_scalar_var_names,
        std::vector<std::string> mat_elem_tensor_var_names,
        std::vector<std::string> node_scalar_var_names,
        std::vector<std::string> node_vector_var_names,
        const int den_id,
        const int pres_id,
        const int sie_id,
        const int sspd_id,
        const int mass_id,
        const int stress_id,
        const int vol_id,
        const int div_id,
        const int level_set_id,
        const int vel_grad_id,
        const int conductivity_id,
        const int specific_heat_id,
        const int node_mass_id,
        const int node_vel_id,
        const int node_accel_id,
        const int node_coord_id,
        const int node_grad_level_set_id,
        const int node_temp_id,
        const int mat_den_id,
        const int mat_pres_id,
        const int mat_sie_id,
        const int mat_sspd_id,
        const int mat_mass_id,
        const int mat_volfrac_id,  
        const int mat_geo_volfrac_id,  
        const int mat_eroded_id,
        const int mat_stress_id,
        const int mat_conductivity_id,
        const int mat_specific_heat_id)
    {
        // **************************************
        //  build and save element average fields
        // **************************************

        // short hand
        const size_t num_nodes = mesh.num_nodes;
        const size_t num_elems = mesh.num_elems;
        const size_t num_dims  = mesh.num_dims;
        const size_t num_nodes_in_elem = mesh.num_nodes_in_elem;
        const size_t num_local_elems = mesh.num_local_elems;
        const size_t num_local_nodes = mesh.num_local_nodes;
        DistributedMap local_element_map = mesh.local_element_map;
        DistributedMap node_map = mesh.node_map;
        DistributedMap all_node_map = mesh.all_node_map;
        DistributedMap nonoverlap_elem_node_map = mesh.nonoverlap_element_node_map;
        const int Pn_order = mesh.Pn;

        const size_t num_mats = State.MaterialPoints.num_material_points.size();
        
        int myrank, nranks;
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        MPI_Comm_size(MPI_COMM_WORLD,&nranks);

        /* save the elem state to an array for exporting to graphics files*/

        //host version of local element map for argument compatibility
        DistributedDFArray<double> elem_scalar_fields(local_element_map, num_elem_scalar_vars, "elem_scalars");
        DistributedDFArray<double> elem_tensor_fields(local_element_map, num_elem_tensor_vars, 3, 3, "elem_tensors");
        DistributedDCArray<size_t> nonoverlap_nodes_in_elem(local_element_map, num_nodes_in_elem, "nonoverlap_nodes_in_elem");
        elem_scalar_fields.set_values(0.0);
        elem_tensor_fields.set_values(0.0);

        // -----------------------------------------------------------------------
        // save the output fields to a single element average array for all state
        // -----------------------------------------------------------------------
        for (int mat_id = 0; mat_id < num_mats; mat_id++) {

            // material point and guass point state are concatenated together
            concatenate_elem_fields(State.MaterialPoints,
                                    State.GaussPoints,
                                    elem_scalar_fields,
                                    elem_tensor_fields,
                                    State.MaterialToMeshMaps.elem_in_mat_elem,
                                    SimulationParameters.output_options.output_elem_state,
                                    SimulationParameters.output_options.output_gauss_pt_state,
                                    State.MaterialToMeshMaps.num_mat_local_elems.host(mat_id),
                                    mat_id,
                                    num_local_elems,
                                    den_id,
                                    pres_id,
                                    sie_id,
                                    sspd_id,
                                    mass_id,
                                    stress_id,
                                    vol_id,
                                    div_id,
                                    level_set_id,
                                    vel_grad_id,
                                    conductivity_id,
                                    specific_heat_id);
        } // end for mats

        // make specific fields for the element average
        if (sie_id>=0){
            FOR_ALL(elem_gid, 0, num_local_elems, {
                // get sie by dividing by the mass
                elem_scalar_fields(elem_gid, sie_id) /= (elem_scalar_fields(elem_gid, mass_id)+1.e-20); 
            });
        } // end if

        Kokkos::fence();
        elem_scalar_fields.update_host();
        elem_tensor_fields.update_host();
        

        // ************************
        //  Build the nodal fields 
        // ************************

        // save the nodal fields to an array for exporting to graphics files
        DistributedDFArray<double> node_scalar_fields(node_map, num_node_scalar_vars, "node_scalars");
        DistributedDFArray<double> node_vector_fields(node_map, num_node_vector_vars,  3, "node_tenors");
        DistributedDFArray<double> nonoverlap_node_scalar_fields(nonoverlap_elem_node_map, num_node_scalar_vars, "node_scalars");
        DistributedDFArray<double> nonoverlap_node_vector_fields(nonoverlap_elem_node_map, num_node_vector_vars,  3, "node_tenors");
    
        concatenate_nodal_fields(State.node,
                                node_scalar_fields,
                                node_vector_fields,
                                SimulationParameters.output_options.output_node_state,
                                dt,
                                num_local_nodes,
                                num_dims,
                                node_mass_id,
                                node_vel_id,
                                node_accel_id,
                                node_coord_id,
                                node_grad_level_set_id,
                                node_temp_id);
                                

        Kokkos::fence();
        node_scalar_fields.update_host();
        node_vector_fields.update_host();

        // ***************************************************************************
        //  Communications for node data from node map to nodes on unique element map
        // ***************************************************************************

        //node data comms
        OutputCommPlan<double> node_scalars_comms(nonoverlap_node_scalar_fields, node_scalar_fields);
        OutputCommPlan<double> node_vectors_comms(nonoverlap_node_vector_fields, node_vector_fields, node_scalars_comms);
        node_scalars_comms.execute_comms();
        node_vectors_comms.execute_comms();

        //nodal coordinates comms
        //convert nodes in elem back to global (convert back to local after we've collected global ids in collective vector)
        DistributedDCArray<double> nonoverlap_node_coords(nonoverlap_elem_node_map, mesh.num_dims);
        CommPlan<double> node_coords_comms(nonoverlap_node_coords, State.node.coords);
        node_coords_comms.execute_comms();

        //convert local node ids in nodes in elem so they correspond to nonoverlap map
        for (size_t elem_id = 0; elem_id < num_local_elems; elem_id++) {
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                nonoverlap_nodes_in_elem.host(elem_id, node_lid) = nonoverlap_elem_node_map.getLocalIndex(all_node_map.getGlobalIndex(mesh.nodes_in_elem.host(elem_id, node_lid)));
            }
        } // end for elem_gid

        //nonoverlap_nodes_in_elem.print();

        if(myrank==0){
            // create the folder structure if it does not exist
            struct stat st;

            if (stat("vtk", &st) != 0) {
                int returnCode = system("mkdir vtk");

                if (returnCode == 1) {
                    std::cout << "Unable to make vtk directory" << std::endl;
                }
            }
            else{
                if(solver_id==0 && graphics_id==0){
                    // delete the existing files inside
                    int returnCode = system("rm vtk/Fierro*");
                    if (returnCode == 1) {
                        std::cout << "Unable to clear vtk/Fierro directory" << std::endl;
                    }
                }
            }

            if (stat("vtk/data", &st) != 0) {
                int returnCode = system("mkdir vtk/data");
                if (returnCode == 1) {
                    std::cout << "Unable to make vtk/data directory" << std::endl;
                }
            }
            else{
                if(solver_id==0 && graphics_id==0){
                    // delete the existing files inside the folder
                    int returnCode = system("rm vtk/data/Fierro*");
                    if (returnCode == 1) {
                        std::cout << "Unable to clear vtk/data directory" << std::endl;
                    }
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        // call the .vtu writer for element fields
        std::string elem_fields_name = "fields";

        write_vtu(nonoverlap_node_coords,
                nonoverlap_nodes_in_elem,
                elem_scalar_fields,
                elem_tensor_fields,
                nonoverlap_node_scalar_fields,
                nonoverlap_node_vector_fields,
                elem_scalar_var_names,
                elem_tensor_var_names,
                node_scalar_var_names,
                node_vector_var_names,
                elem_fields_name,
                graphics_id,
                nonoverlap_elem_node_map.size(),
                mesh.num_local_elems,
                num_nodes_in_elem,
                Pn_order,
                num_dims,
                solver_id);

        // ********************************
        //  Build and write the mat fields 
        // ********************************


        // note: the file path and folder was created in the elem and node outputs
        size_t num_mat_files_written = 0;
        CArray<int> local_mats_in_rank(num_mats);
        int local_num_mats = 0;
        if(num_mat_pt_scalar_vars > 0 || num_mat_pt_tensor_vars >0){
            for (int mat_id = 0; mat_id < num_mats; mat_id++) {

                const size_t num_mat_local_elems = State.MaterialToMeshMaps.num_mat_local_elems.host(mat_id);
                if(num_mat_local_elems){
                    local_mats_in_rank(local_num_mats) = mat_id;
                    local_num_mats++;
                }

                //set global element indices on this rank for this mat
                DistributedMap element_map = mesh.element_map;
                DCArrayKokkos<long long int> global_indices_of_local_mat_elems(num_mat_local_elems, " global_indices_of_local_mat_elems");
                // FOR_ALL(ielem, 0, num_mat_local_elems,{
                //     global_indices_of_local_mat_elems(ielem) = mesh.element_map(State.MaterialToMeshMaps.elem_in_mat_elem(mat_id, ielem));
                // });
                for(int ielem = 0; ielem < num_mat_local_elems; ielem++){
                    global_indices_of_local_mat_elems.host(ielem) = mesh.element_map.getGlobalIndex(State.MaterialToMeshMaps.elem_in_mat_elem.host(mat_id, ielem));
                }
                global_indices_of_local_mat_elems.update_device();
                DistributedMap mat_elem_map = DistributedMap(global_indices_of_local_mat_elems);
                //mat_elem_map.print();

                //allocate arrays for distributed mat elem data
                DistributedDFArray<double> mat_elem_scalar_fields(mat_elem_map, num_mat_pt_scalar_vars, "mat_elem_scalars");
                DistributedDFArray<double> mat_elem_tensor_fields(mat_elem_map, num_mat_pt_tensor_vars, 3, 3, "mat_elem_tensors");

                // set the nodal vars to zero size, we don't write these fields again
                node_scalar_var_names.clear();
                node_vector_var_names.clear();

                // concatenate material fields into a single array
                concatenate_mat_fields(State.MaterialPoints,
                                        mat_elem_scalar_fields,
                                        mat_elem_tensor_fields,
                                        State.MaterialToMeshMaps.elem_in_mat_elem,
                                        SimulationParameters.output_options.output_mat_pt_state,
                                        num_mat_local_elems,
                                        mat_id,
                                        mat_den_id,
                                        mat_pres_id,
                                        mat_sie_id,
                                        mat_sspd_id,
                                        mat_mass_id,
                                        mat_volfrac_id,
                                        mat_geo_volfrac_id,  
                                        mat_eroded_id,
                                        mat_stress_id,
                                        mat_conductivity_id,
                                        mat_specific_heat_id);


                std::string str_mat_val = std::to_string(mat_id);                       
                std::string mat_fields_name = "mat";
                mat_fields_name += str_mat_val;  // add the mat number

                // the number of actual nodes belonging to the part (i.e., the material)
                size_t num_mat_nodes = 0;

                //communicate scalars, tensors, and nodes in elem to collective mat arrays on rank 0
                
                //collect nodes in elem for this material on rank 0
                DistributedDCArray<size_t> mat_nodes_in_mat_elem(mat_elem_map, num_nodes_in_elem, "mat_nodes_in_mat_elem");
                CommPlan<size_t> mat_nodes_in_elem_comms(mat_nodes_in_mat_elem,mesh.nodes_in_elem); //shouldnt do comms since subview of map on this rank
                mat_nodes_in_elem_comms.execute_comms();

                //convert mesh.nodes_in_elem stores local indices and we communicated these in, convert to global
                for (size_t elem_id = 0; elem_id < num_mat_local_elems; elem_id++) {
                    for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                        mat_nodes_in_mat_elem.host(elem_id, node_lid) = mesh.all_node_map.getGlobalIndex(mat_nodes_in_mat_elem.host(elem_id, node_lid));
                    }
                } // end for elem_gid

                //define set of nodes for this mat, collect on rank 0, comms on coords, scalars, and vectors for nodes for this mat
                
                // build a unique mesh (element and nodes) for the material (i.e., the part)
                DCArrayKokkos<long long int> mat_node_indices;
                build_material_node_list(mesh,
                                        mat_node_indices,
                                        mat_nodes_in_mat_elem,
                                        State.MaterialToMeshMaps.elem_in_mat_elem,
                                        mat_id,
                                        num_mat_nodes,
                                        num_mat_local_elems,
                                        num_nodes_in_elem,
                                        num_dims);
                
                //map object for mat node indices
                DistributedMap mat_node_map = DistributedMap(mat_node_indices);

                DistributedDCArray<double> mat_node_coords(mat_node_map, num_dims, "mat_node_coords");
                CommPlan<double> mat_node_coords_comms(mat_node_coords,State.node.coords); //shouldnt do comms since subview of map on this rank
                mat_node_coords_comms.execute_comms();

                DistributedDFArray<double> mat_node_scalar_fields(mat_node_map, num_node_scalar_vars, "mat_node_scalars");
                OutputCommPlan<double> mat_node_scalars_comms(mat_node_scalar_fields,node_scalar_fields); //shouldnt do comms since subview of map on this rank
                mat_node_scalars_comms.execute_comms();

                DistributedDFArray<double> mat_node_vector_fields(mat_node_map, num_node_vector_vars, 3, "mat_node_vectors");
                OutputCommPlan<double> mat_node_vectors_comms(mat_node_vector_fields,node_vector_fields); //shouldnt do comms since subview of map on this rank
                mat_node_vectors_comms.execute_comms();

                //convert mat_nodes_in_mat_elem so it uses contiguous local node ids for this mat portion of the mesh
                for (size_t elem_id = 0; elem_id < num_mat_local_elems; elem_id++) {
                    for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                        mat_nodes_in_mat_elem.host(elem_id, node_lid) = mat_node_map.getLocalIndex(mat_nodes_in_mat_elem.host(elem_id, node_lid));
                    }
                } // end for elem_gid
                
                // only write material data if the mat lives on the mesh, ie. has state allocated
                if (num_mat_local_elems>0){
                    // write out a vtu file this 
                    write_vtu(mat_node_coords,
                                mat_nodes_in_mat_elem,
                                mat_elem_scalar_fields,
                                mat_elem_tensor_fields,
                                mat_node_scalar_fields,
                                mat_node_vector_fields,
                                mat_elem_scalar_var_names,
                                mat_elem_tensor_var_names,
                                node_scalar_var_names,
                                node_vector_var_names,
                                mat_fields_name,
                                graphics_id,
                                num_mat_nodes,
                                num_mat_local_elems,
                                num_nodes_in_elem,
                                Pn_order,
                                num_dims,
                                solver_id);


                    num_mat_files_written++;

                } // end for mat_id

            } // end if material is on the mesh

        } // end if mat variables are to be written


        // *************************************************
        //  write Paraview files to open the graphics files
        // *************************************************

        // save the graphics time
        graphics_times(graphics_id) = time_value;

        // check to see if an mesh state was written 
        bool write_mesh_state = false;
        if( num_elem_scalar_vars > 0 ||
            num_elem_tensor_vars > 0 ||
            num_node_scalar_vars > 0 ||
            num_node_vector_vars > 0)
        {
            write_mesh_state = true;
        }
        // check to see if a mat state was written
        bool write_mat_pt_state = false;
        if( num_mat_pt_scalar_vars > 0 ||
            num_mat_pt_tensor_vars > 0)
        {
                write_mat_pt_state = true;
        }

        // call the vtm file writer
        std::string mat_fields_name = "mat";
        //gather MPI ranks that are writing mat blocks
        write_vtm(graphics_times,
                elem_fields_name,
                mat_fields_name,
                time_value,
                graphics_id,
                num_mats,
                local_num_mats,
                local_mats_in_rank,
                nranks,
                myrank,
                write_mesh_state,
                write_mat_pt_state,
                solver_id);
        
        if(myrank==0){
            // call the pvd file writer
            write_pvd(graphics_times,
                    time_value,
                    graphics_id,
                    solver_id);
        }

        // increment graphics id counter
        graphics_id++; // this is private variable in the class
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn write_collective vtm
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
    void write_collective_viz(Mesh_t& mesh,
        State_t& State,
        SimulationParameters_t& SimulationParameters,
        double dt,
        double time_value,
        CArray<double> graphics_times,
        const size_t solver_id,
        std::vector<node_state> node_states,
        std::vector<gauss_pt_state> gauss_pt_states,
        std::vector<material_pt_state> material_pt_states,
        size_t num_elem_scalar_vars,
        size_t num_elem_tensor_vars,
        size_t num_node_scalar_vars,
        size_t num_node_vector_vars,
        size_t num_mat_pt_scalar_vars,
        size_t num_mat_pt_tensor_vars,
        std::vector<std::string> elem_scalar_var_names,
        std::vector<std::string> elem_tensor_var_names,
        std::vector<std::string> mat_elem_scalar_var_names,
        std::vector<std::string> mat_elem_tensor_var_names,
        std::vector<std::string> node_scalar_var_names,
        std::vector<std::string> node_vector_var_names,
        const int den_id,
        const int pres_id,
        const int sie_id,
        const int sspd_id,
        const int mass_id,
        const int stress_id,
        const int vol_id,
        const int div_id,
        const int level_set_id,
        const int vel_grad_id,
        const int conductivity_id,
        const int specific_heat_id,
        const int node_mass_id,
        const int node_vel_id,
        const int node_accel_id,
        const int node_coord_id,
        const int node_grad_level_set_id,
        const int node_temp_id,
        const int mat_den_id,
        const int mat_pres_id,
        const int mat_sie_id,
        const int mat_sspd_id,
        const int mat_mass_id,
        const int mat_volfrac_id,  
        const int mat_geo_volfrac_id,  
        const int mat_eroded_id,
        const int mat_stress_id,
        const int mat_conductivity_id,
        const int mat_specific_heat_id)
    {
        // // **************************************
        // //  build and save element average fields
        // // **************************************

        // // short hand
        // const size_t num_nodes = mesh.num_nodes;
        // const size_t num_elems = mesh.num_elems;
        // const size_t num_dims  = mesh.num_dims;
        // const size_t num_nodes_in_elem = mesh.num_nodes_in_elem;
        // const size_t num_local_elems = mesh.num_local_elems;
        // const size_t num_local_nodes = mesh.num_local_nodes;
        // DistributedMap local_element_map = mesh.local_element_map;
        // DistributedMap node_map = mesh.node_map;
        // const int Pn_order = mesh.Pn;

        // const size_t num_mats = State.MaterialPoints.num_material_points.size();
        
        // int myrank, nranks;
        // MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        // MPI_Comm_size(MPI_COMM_WORLD,&nranks);

        // /* save the elem state to an array for exporting to graphics files*/

        // //host version of local element map for argument compatibility
        // HostDistributedMap host_local_element_map;
        // DCArrayKokkos<long long int, Kokkos::LayoutLeft , Kokkos::HostSpace> global_indices_of_local_elements(num_local_elems);
        // for(int ielem = 0; ielem < num_local_elems; ielem++){
        //     global_indices_of_local_elements(ielem) = local_element_map.getGlobalIndex(ielem);
        // }
        // host_local_element_map = HostDistributedMap(global_indices_of_local_elements);
        // DistributedDFArray<double> elem_scalar_fields(local_element_map, num_elem_scalar_vars, "elem_scalars");
        // DistributedDFArray<double> elem_tensor_fields(local_element_map, num_elem_tensor_vars, 3, 3, "elem_tensors");
        // elem_scalar_fields.set_values(0.0);
        // elem_tensor_fields.set_values(0.0);
        // //duplicate for now to allow compatibility with comm plan object when using Tpetra (src and dst device type must be equal)
        // //We dont want to make a dual view of the rank 0 collector since that will blow device memory constraints sooner than this duplicate
        // //one other option is to just do the concatenation ops on the host
        // DistributedFArray<double> host_elem_scalar_fields(host_local_element_map, num_elem_scalar_vars, "elem_scalars");
        // DistributedFArray<double> host_elem_tensor_fields(host_local_element_map, num_elem_tensor_vars, 3, 3, "elem_tensors");


        // // -----------------------------------------------------------------------
        // // save the output fields to a single element average array for all state
        // // -----------------------------------------------------------------------
        // for (int mat_id = 0; mat_id < num_mats; mat_id++) {

        //     // material point and guass point state are concatenated together
        //     concatenate_elem_fields(State.MaterialPoints,
        //                             State.GaussPoints,
        //                             elem_scalar_fields,
        //                             elem_tensor_fields,
        //                             State.MaterialToMeshMaps.elem_in_mat_elem,
        //                             SimulationParameters.output_options.output_elem_state,
        //                             SimulationParameters.output_options.output_gauss_pt_state,
        //                             State.MaterialToMeshMaps.num_mat_local_elems.host(mat_id),
        //                             mat_id,
        //                             num_local_elems,
        //                             den_id,
        //                             pres_id,
        //                             sie_id,
        //                             sspd_id,
        //                             mass_id,
        //                             stress_id,
        //                             vol_id,
        //                             div_id,
        //                             level_set_id,
        //                             vel_grad_id,
        //                             conductivity_id,
        //                             specific_heat_id);
        // } // end for mats

        // // make specific fields for the element average
        // if (sie_id>=0){
        //     FOR_ALL(elem_gid, 0, num_local_elems, {
        //         // get sie by dividing by the mass
        //         elem_scalar_fields(elem_gid, sie_id) /= (elem_scalar_fields(elem_gid, mass_id)+1.e-20); 
        //     });
        // } // end if

        // Kokkos::fence();
        // elem_scalar_fields.update_host();
        // elem_tensor_fields.update_host();

        // // -----------------------------------------------------------------------
        // // copy the output fields to host side array compatible with Tpetra Comms
        // // -----------------------------------------------------------------------

        // for (int mat_id = 0; mat_id < num_mats; mat_id++) {
            
        //     // material point and guass point state are concatenated together
        //     copy_elem_fields(elem_scalar_fields,
        //                     elem_tensor_fields,
        //                     host_elem_scalar_fields,
        //                     host_elem_tensor_fields,
        //                     State.MaterialToMeshMaps.elem_in_mat_elem,
        //                     SimulationParameters.output_options.output_elem_state,
        //                     SimulationParameters.output_options.output_gauss_pt_state,
        //                     State.MaterialToMeshMaps.num_mat_local_elems.host(mat_id),
        //                     mat_id,
        //                     num_local_elems,
        //                     den_id,
        //                     pres_id,
        //                     sie_id,
        //                     sspd_id,
        //                     mass_id,
        //                     stress_id,
        //                     vol_id,
        //                     div_id,
        //                     level_set_id,
        //                     vel_grad_id,
        //                     conductivity_id,
        //                     specific_heat_id);
        // } // end for mats
        

        // // ************************
        // //  Build the nodal fields 
        // // ************************

        // //host version of local node map for argument compatibility
        // HostDistributedMap host_node_map;
        // DCArrayKokkos<long long int, Kokkos::LayoutLeft , Kokkos::HostSpace> global_indices_of_local_nodes(mesh.num_local_nodes);
        // for(int inode = 0; inode < mesh.num_local_nodes; inode++){
        //     global_indices_of_local_nodes(inode) = mesh.node_map.getGlobalIndex(inode);
        // }
        // host_node_map = HostDistributedMap(global_indices_of_local_nodes);

        // // save the nodal fields to an array for exporting to graphics files
        // DistributedDFArray<double> node_scalar_fields(node_map, num_node_scalar_vars, "node_scalars");
        // DistributedDFArray<double> node_vector_fields(node_map, num_node_vector_vars,  3, "node_tenors");
        // DistributedFArray<double> host_node_scalar_fields(host_node_map, num_node_scalar_vars, "node_scalars");
        // DistributedFArray<double> host_node_vector_fields(host_node_map, num_node_vector_vars, 3, "node_tenors");
    
        // concatenate_nodal_fields(State.node,
        //                         node_scalar_fields,
        //                         node_vector_fields,
        //                         SimulationParameters.output_options.output_node_state,
        //                         dt,
        //                         num_local_nodes,
        //                         num_dims,
        //                         node_mass_id,
        //                         node_vel_id,
        //                         node_accel_id,
        //                         node_coord_id,
        //                         node_grad_level_set_id,
        //                         node_temp_id);
                                

        // Kokkos::fence();
        // node_scalar_fields.update_host();
        // node_vector_fields.update_host();

        // copy_nodal_fields(host_node_scalar_fields,
        //                 host_node_vector_fields,
        //                 node_scalar_fields,
        //                 node_vector_fields,
        //                 SimulationParameters.output_options.output_node_state,
        //                 dt,
        //                 num_local_nodes,
        //                 num_dims,
        //                 node_mass_id,
        //                 node_vel_id,
        //                 node_accel_id,
        //                 node_coord_id,
        //                 node_grad_level_set_id,
        //                 node_temp_id);

        // // **************************************************
        // //  Collective communications for node and elem data 
        // // **************************************************

        // //elem data collective comms
        // //collective map has all indices on rank 0 and non on other ranks
        // HostDistributedMap collective_elem_map;
        // long long int num_collective_elem_indices = 0;
        // if(myrank==0) num_collective_elem_indices = mesh.global_num_elems;
        // collective_elem_map = HostDistributedMap(mesh.global_num_elems, num_collective_elem_indices);

        // //collective vector and comms to the collective vector for elem fields
        // DistributedFArray<double> collective_elem_scalar_fields(collective_elem_map, num_elem_scalar_vars, "elem_scalars_collective");
        // DistributedFArray<double> collective_elem_tensor_fields(collective_elem_map, num_elem_tensor_vars, 3, 3, "elem_tensors_collective");
        // HostCommPlan<double> collective_elem_scalars_comms(collective_elem_scalar_fields, host_elem_scalar_fields);
        // HostCommPlan<double> collective_elem_tensors_comms(collective_elem_tensor_fields, host_elem_tensor_fields, collective_elem_scalars_comms);
        // collective_elem_scalars_comms.execute_comms();
        // collective_elem_tensors_comms.execute_comms();

        // //host of node in elem for Trilinos template argument compatibility
        // DistributedFArray<size_t> host_local_nodes_in_elem(host_local_element_map, mesh.num_nodes_in_elem);

        // //convert nodes in elem back to global (convert back to local after we've collected global ids in collective vector)
        // for (size_t elem_id = 0; elem_id < num_local_elems; elem_id++) {
        //     for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
        //         host_local_nodes_in_elem(elem_id, node_lid) = mesh.all_node_map.getGlobalIndex(mesh.local_nodes_in_elem(elem_id, node_lid));
        //     }
        // } // end for elem_gid

        // //collect nodes in elem with a conversion back to global node ids
        // DistributedFArray<size_t> collective_nodes_in_elem(collective_elem_map, mesh.num_nodes_in_elem);
        // HostCommPlan<size_t> nodes_in_elem_comms(collective_nodes_in_elem, host_local_nodes_in_elem);

        // nodes_in_elem_comms.execute_comms();

        // //node data collective comms
        // //collective map has all indices on rank 0 and non on other ranks
        // HostDistributedMap collective_node_map;
        // long long int num_collective_node_indices = 0;
        // if(myrank==0) num_collective_node_indices = mesh.global_num_nodes;
        // collective_node_map = HostDistributedMap(mesh.global_num_nodes, num_collective_node_indices);

        // //collective vector and comms to the collective vector for node fields
        // DistributedFArray<double> collective_node_scalar_fields(collective_node_map, num_node_scalar_vars);
        // DistributedFArray<double> collective_node_vector_fields(collective_node_map, num_node_vector_vars, 3);
        // HostCommPlan<double> collective_node_scalars_comms(collective_node_scalar_fields, host_node_scalar_fields);
        // HostCommPlan<double> collective_node_vectors_comms(collective_node_vector_fields, host_node_vector_fields, collective_node_scalars_comms);
        // collective_node_scalars_comms.execute_comms();
        // collective_node_vectors_comms.execute_comms();

        // //collect nodal coordinates
        // //convert nodes in elem back to global (convert back to local after we've collected global ids in collective vector)
        // DistributedFArray<double> host_node_coords(host_node_map, mesh.num_dims);
        // for (size_t node_id = 0; node_id < num_local_nodes; node_id++) {
        //     for (int idim = 0; idim < mesh.num_dims; idim++) {
        //         host_node_coords(node_id, idim) = State.node.coords.host(node_id, idim);
        //     }
        // } // end for elem_gid
        // DistributedFArray<double> collective_node_coords(collective_node_map, mesh.num_dims);
        // HostCommPlan<double> collective_node_coords_comms(collective_node_coords, host_node_coords, collective_node_scalars_comms);
        // collective_node_coords_comms.execute_comms();

        // if(myrank==0){
        //     // create the folder structure if it does not exist
        //     struct stat st;

        //     if (stat("vtk", &st) != 0) {
        //         int returnCode = system("mkdir vtk");

        //         if (returnCode == 1) {
        //             std::cout << "Unable to make vtk directory" << std::endl;
        //         }
        //     }
        //     else{
        //         if(solver_id==0 && graphics_id==0){
        //             // delete the existing files inside
        //             int returnCode = system("rm vtk/Fierro*");
        //             if (returnCode == 1) {
        //                 std::cout << "Unable to clear vtk/Fierro directory" << std::endl;
        //             }
        //         }
        //     }

        //     if (stat("vtk/data", &st) != 0) {
        //         int returnCode = system("mkdir vtk/data");
        //         if (returnCode == 1) {
        //             std::cout << "Unable to make vtk/data directory" << std::endl;
        //         }
        //     }
        //     else{
        //         if(solver_id==0 && graphics_id==0){
        //             // delete the existing files inside the folder
        //             int returnCode = system("rm vtk/data/Fierro*");
        //             if (returnCode == 1) {
        //                 std::cout << "Unable to clear vtk/data directory" << std::endl;
        //             }
        //         }
        //     }
        // }
        // // call the .vtu writer for element fields
        // std::string elem_fields_name = "fields";

        // if(myrank==0){
        //     write_vtu(collective_node_coords,
        //             collective_nodes_in_elem,
        //             collective_elem_scalar_fields,
        //             collective_elem_tensor_fields,
        //             collective_node_scalar_fields,
        //             collective_node_vector_fields,
        //             elem_scalar_var_names,
        //             elem_tensor_var_names,
        //             node_scalar_var_names,
        //             node_vector_var_names,
        //             elem_fields_name,
        //             graphics_id,
        //             mesh.global_num_nodes,
        //             mesh.global_num_elems,
        //             num_nodes_in_elem,
        //             Pn_order,
        //             num_dims,
        //             solver_id);
        // }

        // // ********************************
        // //  Build and write the mat fields 
        // // ********************************


        // // note: the file path and folder was created in the elem and node outputs
        // size_t num_mat_files_written = 0;
        // if(num_mat_pt_scalar_vars > 0 || num_mat_pt_tensor_vars >0){

        //     for (int mat_id = 0; mat_id < num_mats; mat_id++) {

        //         const size_t num_mat_local_elems = State.MaterialToMeshMaps.num_mat_local_elems.host(mat_id);
        //         //array storing number of local elems for this material on each process
        //         CArray<int> processes_num_local_mat_elems, gatherv_displacements;
        //         if(myrank==0){
        //             processes_num_local_mat_elems = CArray<int>(nranks);
        //             gatherv_displacements = CArray<int>(nranks);
        //         }
        //         MPI_Gather(&num_mat_local_elems,1,MPI_INT,processes_num_local_mat_elems.pointer(),1,
        //                     MPI_INT, 0, MPI_COMM_WORLD);

        //         //set global element indices on this rank
        //         HostDistributedMap host_mat_elem_map;
        //         DistributedMap element_map = mesh.element_map;
        //         DCArrayKokkos<long long int, Kokkos::LayoutLeft , Kokkos::HostSpace> global_indices_of_local_mat_elems(num_mat_local_elems);
        //         for(int ielem = 0; ielem < num_mat_local_elems; ielem++){
        //             global_indices_of_local_mat_elems(ielem) = element_map.getGlobalIndex(State.MaterialToMeshMaps.elem_in_mat_elem.host(mat_id, ielem));
        //         }
        //         host_mat_elem_map = HostDistributedMap(global_indices_of_local_mat_elems);

        //         //allocate arrays for distributed mat elem data
        //         DistributedFArray<double> host_mat_elem_scalar_fields(host_mat_elem_map, num_mat_pt_scalar_vars, "mat_elem_scalars");
        //         DistributedFArray<double> host_mat_elem_tensor_fields(host_mat_elem_map, num_mat_pt_tensor_vars, 3, 3, "mat_elem_tensors");

        //         //collect global element indices on rank 0 for this mat
        //         //tally total number of mat elems for rank 0
        //         DCArrayKokkos<long long int, Kokkos::LayoutLeft , Kokkos::HostSpace> global_indices_of_collective_mat_elems;
        //         long long int num_mat_collective_elems = 0;
        //         if(myrank==0){
        //             for(int irank=0; irank < nranks; irank++){
        //                 gatherv_displacements(irank) = num_mat_collective_elems;
        //                 num_mat_collective_elems += processes_num_local_mat_elems(irank);
        //             }
        //             global_indices_of_collective_mat_elems = DCArrayKokkos<long long int, Kokkos::LayoutLeft , Kokkos::HostSpace>(num_mat_collective_elems);
        //         }
        //         // if(myrank==0){
        //         //     for(int irank=0; irank < nranks; irank++){
        //         //         std::cout << "NUM local mat elem on rank " << irank << " is " << processes_num_local_mat_elems(irank) << std::endl;
        //         //         std::cout << "gatherv displacement on rank " << irank << " is " << gatherv_displacements(irank) << std::endl;
        //         //     }
        //         // }
        //         MPI_Gatherv(global_indices_of_local_mat_elems.device_pointer(), num_mat_local_elems, MPI_LONG_LONG_INT,
        //                     global_indices_of_collective_mat_elems.device_pointer(), processes_num_local_mat_elems.pointer(),
        //                     gatherv_displacements.pointer(), MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

        //         //use indices on rank 0 to construct rank 0 collective map for this mat
        //         HostDistributedMap collective_mat_elem_map;
        //         collective_mat_elem_map = HostDistributedMap(global_indices_of_collective_mat_elems);
        //         //collective_mat_elem_map.print();

        //         //collective storage for scalars and tensors using collective elem mat map
        //         DistributedFArray<double> collective_mat_elem_scalar_fields(collective_mat_elem_map, num_mat_pt_scalar_vars, "mat_elem_scalars_collective");
        //         DistributedFArray<double> collective_mat_elem_tensor_fields(collective_mat_elem_map, num_mat_pt_tensor_vars, 3, 3, "mat_elem_tensors_collective");

        //         // set the nodal vars to zero size, we don't write these fields again
        //         node_scalar_var_names.clear();
        //         node_vector_var_names.clear();

        //         // concatenate material fields into a single array
        //         concatenate_mat_fields(State.MaterialPoints,
        //                                 host_mat_elem_scalar_fields,
        //                                 host_mat_elem_tensor_fields,
        //                                 State.MaterialToMeshMaps.elem_in_mat_elem,
        //                                 SimulationParameters.output_options.output_mat_pt_state,
        //                                 num_mat_local_elems,
        //                                 mat_id,
        //                                 mat_den_id,
        //                                 mat_pres_id,
        //                                 mat_sie_id,
        //                                 mat_sspd_id,
        //                                 mat_mass_id,
        //                                 mat_volfrac_id,
        //                                 mat_geo_volfrac_id,  
        //                                 mat_eroded_id,
        //                                 mat_stress_id,
        //                                 mat_conductivity_id,
        //                                 mat_specific_heat_id);


        //         std::string str_mat_val = std::to_string(mat_id);                       
        //         std::string mat_fields_name = "mat";
        //         mat_fields_name += str_mat_val;  // add the mat number

        //         // the number of actual nodes belonging to the part (i.e., the material)
        //         size_t num_mat_nodes = 0;

        //         //communicate scalars, tensors, and nodes in elem to collective mat arrays on rank 0
                
        //         //collect nodes in elem for this material on rank 0
        //         DistributedFArray<size_t> collective_mat_nodes_in_mat_elem(collective_mat_elem_map, num_nodes_in_elem, "collective_mat_nodes_in_mat_elem");
        //         HostCommPlan<size_t> mat_nodes_in_elem_comms(collective_mat_nodes_in_mat_elem,collective_nodes_in_elem); //doesnt really do comms since all on rank 0
        //         mat_nodes_in_elem_comms.execute_comms();

        //         HostCommPlan<double> mat_elem_scalars_comms(collective_mat_elem_scalar_fields,host_mat_elem_scalar_fields);
        //         mat_elem_scalars_comms.execute_comms();

        //         HostCommPlan<double> mat_elem_tensors_comms(collective_mat_elem_tensor_fields,host_mat_elem_tensor_fields);
        //         mat_elem_tensors_comms.execute_comms();

        //         //define set of nodes for this mat, collect on rank 0, comms on coords, scalars, and vectors for nodes for this mat
                
        //         // build a unique mesh (element and nodes) for the material (i.e., the part)
        //         DCArrayKokkos<long long int, Kokkos::LayoutLeft , Kokkos::HostSpace> collective_mat_node_indices;
        //         if(myrank==0){
        //         build_material_node_list(mesh,
        //                                 collective_mat_node_indices,
        //                                 collective_mat_nodes_in_mat_elem,
        //                                 State.MaterialToMeshMaps.elem_in_mat_elem,
        //                                 mat_id,
        //                                 num_mat_nodes,
        //                                 num_mat_collective_elems,
        //                                 num_nodes_in_elem,
        //                                 num_dims);
        //         }
                
        //         //map object for mat node indices
        //         HostDistributedMap collective_mat_node_map = HostDistributedMap(collective_mat_node_indices);

        //         DistributedFArray<double> collective_mat_node_coords(collective_mat_node_map, num_dims, "collective_mat_node_coords");
        //         HostCommPlan<double> mat_node_coords_comms(collective_mat_node_coords,collective_node_coords); //doesnt really do comms since all on rank 0
        //         mat_node_coords_comms.execute_comms();

        //         DistributedFArray<double> collective_mat_node_scalar_fields(collective_mat_node_map, num_node_scalar_vars, "collective_mat_node_scalars");
        //         HostCommPlan<double> mat_node_scalars_comms(collective_mat_node_scalar_fields,collective_node_scalar_fields); //doesnt really do comms since all on rank 0
        //         mat_node_scalars_comms.execute_comms();

        //         DistributedFArray<double> collective_mat_node_vector_fields(collective_mat_node_map, num_node_vector_vars, 3, "collective_mat_node_vectors");
        //         HostCommPlan<double> mat_node_vectors_comms(collective_mat_node_vector_fields,collective_node_vector_fields); //doesnt really do comms since all on rank 0
        //         mat_node_vectors_comms.execute_comms();

        //         //convert collective mat_nodes_in_mat_elem so it uses contiguous node ids for this mat portion of the mesh
        //         for (size_t elem_id = 0; elem_id < num_mat_collective_elems; elem_id++) {
        //             for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
        //                 collective_mat_nodes_in_mat_elem(elem_id, node_lid) = collective_mat_node_map.getLocalIndex(collective_mat_nodes_in_mat_elem(elem_id, node_lid));
        //             }
        //         } // end for elem_gid
                
        //         // only write material data if the mat lives on the mesh, ie. has state allocated
        //         if (num_mat_collective_elems>0&&myrank==0){
        //             // write out a vtu file this 
        //             write_vtu(collective_mat_node_coords,
        //                         collective_mat_nodes_in_mat_elem,
        //                         collective_mat_elem_scalar_fields,
        //                         collective_mat_elem_tensor_fields,
        //                         collective_mat_node_scalar_fields,
        //                         collective_mat_node_vector_fields,
        //                         mat_elem_scalar_var_names,
        //                         mat_elem_tensor_var_names,
        //                         node_scalar_var_names,
        //                         node_vector_var_names,
        //                         mat_fields_name,
        //                         graphics_id,
        //                         num_mat_nodes,
        //                         num_mat_collective_elems,
        //                         num_nodes_in_elem,
        //                         Pn_order,
        //                         num_dims,
        //                         solver_id);


        //             num_mat_files_written++;

        //         } // end for mat_id

        //     } // end if material is on the mesh

        // } // end if mat variables are to be written


        // // *************************************************
        // //  write Paraview files to open the graphics files
        // // *************************************************

        // // save the graphics time
        // graphics_times(graphics_id) = time_value;

        // // check to see if an mesh state was written 
        // bool write_mesh_state = false;
        // if( num_elem_scalar_vars > 0 ||
        //     num_elem_tensor_vars > 0 ||
        //     num_node_scalar_vars > 0 ||
        //     num_node_vector_vars > 0)
        // {
        //     write_mesh_state = true;
        // }
        // // check to see if a mat state was written
        // bool write_mat_pt_state = false;
        // if( num_mat_pt_scalar_vars > 0 ||
        //     num_mat_pt_tensor_vars > 0)
        // {
        //         write_mat_pt_state = true;
        // }

        // // call the vtm file writer
        // std::string mat_fields_name = "mat";
        // if(myrank==0){
        //     write_vtm(graphics_times,
        //             elem_fields_name,
        //             mat_fields_name,
        //             time_value,
        //             graphics_id,
        //             num_mat_files_written,
        //             write_mesh_state,
        //             write_mat_pt_state,
        //             solver_id);

        //     // call the pvd file writer
        //     write_pvd(graphics_times,
        //             time_value,
        //             graphics_id,
        //             solver_id);
        // }

        // // increment graphics id counter
        // graphics_id++; // this is private variable in the class
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
        SimulationParameters_t& SimulationParameters,
        double dt,
        double time_value,
        CArray<double> graphics_times,
        std::vector<node_state> node_states,
        std::vector<gauss_pt_state> gauss_pt_states,
        std::vector<material_pt_state> material_pt_states)
    {
        size_t num_mats = State.MaterialPoints.num_material_points.size();
        
        int myrank, nranks;
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        MPI_Comm_size(MPI_COMM_WORLD,&nranks);
        // ---- Update host data ----

        // material point values
        State.MaterialPoints.den.update_host();
        State.MaterialPoints.pres.update_host();
        State.MaterialPoints.stress.update_host();
        State.MaterialPoints.sspd.update_host();
        State.MaterialPoints.sie.update_host();
        State.MaterialPoints.mass.update_host();
        State.MaterialPoints.eroded.update_host();


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
        
        //host version of local element map for argument compatibility
        HostDistributedMap host_local_element_map;
        DCArrayKokkos<long long int, Kokkos::LayoutLeft , Kokkos::HostSpace> global_indices_of_local_elements(mesh.num_local_elems, "global_indices_of_local_elements");
        for(int ielem = 0; ielem < mesh.num_local_elems; ielem++){
            global_indices_of_local_elements(ielem) = mesh.local_element_map.getGlobalIndex(ielem);
        }
        host_local_element_map = HostDistributedMap(global_indices_of_local_elements);
        DistributedCArray<double> elem_fields = DistributedCArray<double>(host_local_element_map, num_scalar_vars, "elem_fields");
        int  elem_switch = 1;

        DCArrayKokkos<double> speed(num_elems, "speed");
        FOR_ALL(elem_gid, 0, mesh.num_local_elems, {
            double elem_vel[3]; // note:initialization with a list won't work
            elem_vel[0] = 0.0;
            elem_vel[1] = 0.0;
            elem_vel[2] = 0.0;
            // get the coordinates of the element center
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                elem_vel[0] += State.node.vel(mesh.nodes_in_elem(elem_gid, node_lid), 0);
                elem_vel[1] += State.node.vel(mesh.nodes_in_elem(elem_gid, node_lid), 1);
                if (mesh.num_dims == 3) {
                    elem_vel[2] += State.node.vel(mesh.nodes_in_elem(elem_gid, node_lid), 2);
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
            size_t num_mat_elems = State.MaterialToMeshMaps.num_mat_local_elems.host(mat_id);

            for (size_t mat_elem_sid = 0; mat_elem_sid < num_mat_elems; mat_elem_sid++) {
                // 1 material per element

                // get elem gid
                size_t elem_gid = State.MaterialToMeshMaps.elem_in_mat_elem.host(mat_id, mat_elem_sid);

                // save outputs
                elem_fields(elem_gid, 0) = State.MaterialPoints.den.host(mat_id, mat_elem_sid);
                elem_fields(elem_gid, 1) = State.MaterialPoints.pres.host(mat_id, mat_elem_sid);
                elem_fields(elem_gid, 2) = State.MaterialPoints.sie.host(mat_id, mat_elem_sid);
                // 3 is guass point vol
                elem_fields(elem_gid, 4) = State.MaterialPoints.mass.host(mat_id, mat_elem_sid);
                elem_fields(elem_gid, 5) = State.MaterialPoints.sspd.host(mat_id, mat_elem_sid);
                // 6 is elem speed
                elem_fields(elem_gid, 7) = (double)mat_id;
                // 8 is the e_switch
                elem_fields(elem_gid, 9) = (double)State.MaterialPoints.eroded.host(mat_id, mat_elem_sid);
            } // end for mat elems storage
        } // end parallel loop over materials

        // export element centric data
        double e_switch = 1;
        for (size_t elem_gid = 0; elem_gid < mesh.num_local_elems; elem_gid++) {
            elem_fields(elem_gid, 3) = State.GaussPoints.vol.host(elem_gid);
            elem_fields(elem_gid, 6) = speed.host(elem_gid);
            elem_fields(elem_gid, 8) = e_switch;
            elem_switch *= -1;
        } // end for elem_gid

        //collective map has all indices on rank 0 and non on other ranks
        HostDistributedMap collective_elem_map;
        long long int num_collective_elem_indices = 0;
        if(myrank==0) num_collective_elem_indices = mesh.global_num_elems;
        collective_elem_map = HostDistributedMap(mesh.global_num_elems, num_collective_elem_indices);

        //collective vector and comms to the collective vector for elem fields
        DistributedCArray<double> collective_elem_fields(collective_elem_map, num_scalar_vars, "collective_elem_fields");
        HostCommPlanLR<double> collective_elem_comms(collective_elem_fields, elem_fields);
        collective_elem_comms.execute_comms();

        //host of node in elem for Trilinos template argument compatibility
        DistributedCArray<size_t> host_local_nodes_in_elem(host_local_element_map, mesh.num_nodes_in_elem, "host_local_nodes_in_elem");

        //convert nodes in elem back to global (convert back to local after we've collected global ids in collective vector)
        for (size_t elem_id = 0; elem_id < mesh.num_local_elems; elem_id++) {
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                host_local_nodes_in_elem(elem_id, node_lid) = mesh.all_node_map.getGlobalIndex(mesh.local_nodes_in_elem.host(elem_id, node_lid));
            }
        } // end for elem_gid

        //collect nodes in elem with a conversion back to global node ids
        DistributedCArray<size_t> collective_nodes_in_elem(collective_elem_map, mesh.num_nodes_in_elem);
        HostCommPlanLR<size_t> nodes_in_elem_comms(collective_nodes_in_elem, host_local_nodes_in_elem);

        nodes_in_elem_comms.execute_comms();
        
        //NODE DATA COLLECTION
        //host version of local node map for argument compatibility
        HostDistributedMap host_node_map;
        DCArrayKokkos<long long int, Kokkos::LayoutLeft , Kokkos::HostSpace> global_indices_of_local_nodes(mesh.num_local_nodes, "global_indices_of_local_nodes");
        for(int inode = 0; inode < mesh.num_local_nodes; inode++){
            global_indices_of_local_nodes(inode) = mesh.node_map.getGlobalIndex(inode);
        }
        host_node_map = HostDistributedMap(global_indices_of_local_nodes);

        // save the vertex vector fields to an array for exporting to graphics files
        DistributedCArray<double> vec_fields(host_node_map, num_vec_vars, 3, "vec_fields");

        for (size_t node_gid = 0; node_gid < mesh.num_local_nodes; node_gid++) {
            // position, var 0
            vec_fields(node_gid, 0, 0) = State.node.coords.host(node_gid, 0);
            vec_fields(node_gid, 0, 1) = State.node.coords.host(node_gid, 1);
            if (num_dims == 2) {
                vec_fields(node_gid, 0, 2) = 0.0;
            }
            else{
                vec_fields(node_gid, 0, 2) = State.node.coords.host(node_gid, 2);
            }

            // velocity, var 1
            vec_fields(node_gid, 1, 0) = State.node.vel.host(node_gid, 0);
            vec_fields(node_gid, 1, 1) = State.node.vel.host(node_gid, 1);
            if (num_dims == 2) {
                vec_fields(node_gid, 1, 2) = 0.0;
            }
            else{
                vec_fields(node_gid, 1, 2) = State.node.vel.host(node_gid, 2);
            }

            // accelleration, var 2
            vec_fields(node_gid, 2, 0) = (State.node.vel.host(node_gid, 0) - State.node.vel_n0.host(node_gid, 0))/dt;
            vec_fields(node_gid, 2, 1) = (State.node.vel.host(node_gid, 1) - State.node.vel_n0.host(node_gid, 1))/dt;
            if (num_dims == 2) {
                vec_fields(node_gid, 2, 2) = 0.0;
            }
            else{
                vec_fields(node_gid, 2, 2) = (State.node.vel.host(node_gid, 2) - State.node.vel_n0.host(node_gid, 2))/dt;
            }


        } // end for loop over vertices

        //collective map has all indices on rank 0 and non on other ranks
        HostDistributedMap collective_node_map;
        long long int num_collective_node_indices = 0;
        if(myrank==0) num_collective_node_indices = mesh.global_num_nodes;
        collective_node_map = HostDistributedMap(mesh.global_num_nodes, num_collective_node_indices);

        //collective vector and comms to the collective vector for node fields
        DistributedCArray<double> collective_vec_fields(collective_node_map, num_vec_vars, 3, "collective_vec_fields");
        HostCommPlanLR<double> collective_node_comms(collective_vec_fields, vec_fields);
        collective_node_comms.execute_comms();

        if(myrank==0){
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
            fprintf(out[0], "%10lu\n", mesh.global_num_nodes);

            // write all components of the point coordinates
            for (int node_gid = 0; node_gid < mesh.global_num_nodes; node_gid++) {
                fprintf(out[0], "%12.5e\n", collective_vec_fields(node_gid, 0, 0));
            }

            for (int node_gid = 0; node_gid < mesh.global_num_nodes; node_gid++) {
                fprintf(out[0], "%12.5e\n", collective_vec_fields(node_gid, 0, 1));
            }

            for (int node_gid = 0; node_gid < mesh.global_num_nodes; node_gid++) {
                if (num_dims == 3) {
                    fprintf(out[0], "%12.5e\n", collective_vec_fields(node_gid, 0, 2));
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
            fprintf(out[0], "%10lu\n", mesh.global_num_elems);


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
            for (int elem_gid = 0; elem_gid < mesh.global_num_elems; elem_gid++) {
                for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                    fprintf(out[0], "%10lu\t", collective_nodes_in_elem(elem_gid, convert_ijk_to_ensight[node_lid]) + 1); // note: node_gid starts at 1
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

                for (int elem_id = 0; elem_id < mesh.global_num_elems; elem_id++) {
                    fprintf(out[0], "%12.5e\n", collective_elem_fields(elem_id, var));
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

                for (int node_gid = 0; node_gid < mesh.global_num_nodes; node_gid++) {
                    fprintf(out[0], "%12.5e\n", collective_vec_fields(node_gid, var, 0));
                }

                for (int node_gid = 0; node_gid < mesh.global_num_nodes; node_gid++) {
                    fprintf(out[0], "%12.5e\n", collective_vec_fields(node_gid, var, 1));
                }

                for (int node_gid = 0; node_gid < mesh.global_num_nodes; node_gid++) {
                    fprintf(out[0], "%12.5e\n", collective_vec_fields(node_gid, var, 2));
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
        }
        // ---------------------------------------------------------------------------
        // Done writing the graphics dump
        // ---------------------------------------------------------------------------

        // increment graphics id counter
        graphics_id++;

        delete[] name;
        

        return;
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn write_vtk_old
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
    void write_vtk_old(Mesh_t& mesh,
        State_t& State,
        SimulationParameters_t& SimulationParameters,
        double dt,
        double time_value,
        CArray<double> graphics_times,
        std::vector<node_state> node_states,
        std::vector<gauss_pt_state> gauss_pt_states,
        std::vector<material_pt_state> material_pt_states)
    {
        
        int myrank, nranks;
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        MPI_Comm_size(MPI_COMM_WORLD,&nranks);
        size_t num_mats = State.MaterialPoints.num_material_points.size();

        // ---- Update host data ----

        // material point values
        State.MaterialPoints.den.update_host();
        State.MaterialPoints.pres.update_host();
        State.MaterialPoints.stress.update_host();
        State.MaterialPoints.sspd.update_host();
        State.MaterialPoints.sie.update_host();
        State.MaterialPoints.mass.update_host();
        State.MaterialPoints.conductivity.update_host();
        State.MaterialPoints.temp_grad.update_host();
        State.MaterialPoints.eroded.update_host();


        // gauss point values
        State.GaussPoints.vol.update_host();

        // nodal values
        State.node.coords.update_host();
        State.node.vel.update_host();
        State.node.mass.update_host();
        State.node.temp.update_host();

        Kokkos::fence();


        const int num_cell_scalar_vars = 13;
        const int num_cell_vec_vars    = 0;
        const int num_cell_tensor_vars = 0;

        const int num_point_scalar_vars = 1;
        const int num_point_vec_vars = 2;


        // Scalar values associated with a cell
        const char cell_scalar_var_names[num_cell_scalar_vars][15] = {
            "den", "pres", "sie", "vol", "mass", "sspd", "speed", "mat_id", "elem_switch","eroded", "temp_grad_x", "temp_grad_y", "temp_grad_z"
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

        //host version of local element map for argument compatibility
        HostDistributedMap host_local_element_map;
        DCArrayKokkos<long long int, Kokkos::LayoutLeft , Kokkos::HostSpace> global_indices_of_local_elements;
        for(int ielem = 0; ielem < mesh.num_local_elems; ielem++){
            global_indices_of_local_elements(ielem) = mesh.local_element_map.getGlobalIndex(ielem);
        }
        host_local_element_map = HostDistributedMap(global_indices_of_local_elements);

        // save the cell state to an array for exporting to graphics files
        auto elem_fields = DistributedCArray<double>(host_local_element_map, num_cell_scalar_vars);
        int  elem_switch = 1;

        DCArrayKokkos<double> speed(num_elems, "speed");
        FOR_ALL(elem_gid, 0, num_elems, {
            double elem_vel[3]; // note:initialization with a list won't work
            elem_vel[0] = 0.0;
            elem_vel[1] = 0.0;
            elem_vel[2] = 0.0;
            // get the coordinates of the element center
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                elem_vel[0] += State.node.vel(mesh.nodes_in_elem(elem_gid, node_lid), 0);
                elem_vel[1] += State.node.vel(mesh.nodes_in_elem(elem_gid, node_lid), 1);
                if (mesh.num_dims == 3) {
                    elem_vel[2] += State.node.vel(mesh.nodes_in_elem(elem_gid, node_lid), 2);
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
            size_t num_mat_elems = State.MaterialToMeshMaps.num_mat_local_elems.host(mat_id);

            for (size_t mat_elem_sid = 0; mat_elem_sid < num_mat_elems; mat_elem_sid++) {
                // 1 material per element

                // get elem gid
                size_t elem_gid = State.MaterialToMeshMaps.elem_in_mat_elem.host(mat_id, mat_elem_sid);

                // save outputs
                elem_fields(elem_gid, 0) = State.MaterialPoints.den.host(mat_id,mat_elem_sid);
                elem_fields(elem_gid, 1) = State.MaterialPoints.pres.host(mat_id, mat_elem_sid);
                elem_fields(elem_gid, 2) = State.MaterialPoints.sie.host(mat_id, mat_elem_sid);
                // 3 is guass point vol
                elem_fields(elem_gid, 4) = State.MaterialPoints.mass.host(mat_id, mat_elem_sid);
                elem_fields(elem_gid, 5) = State.MaterialPoints.sspd.host(mat_id, mat_elem_sid);
                // 6 is elem speed
                elem_fields(elem_gid, 7) = (double)mat_id;
                // 8 is the e_switch
                elem_fields(elem_gid, 9) = (double)State.MaterialPoints.eroded.host(mat_id, mat_elem_sid);
                elem_fields(elem_gid, 10) = (double)State.MaterialPoints.temp_grad.host(mat_id, elem_gid,0);
                elem_fields(elem_gid, 11) = (double)State.MaterialPoints.temp_grad.host(mat_id, elem_gid,1);
                elem_fields(elem_gid, 12) = (double)State.MaterialPoints.temp_grad.host(mat_id, elem_gid,2);
            } // end for mat elems storage
        } // end parallel loop over materials

        // export element centric data
        double e_switch = 1;
        for (size_t elem_gid = 0; elem_gid < mesh.num_local_elems; elem_gid++) {
            elem_fields(elem_gid, 3) = State.GaussPoints.vol.host(elem_gid);
            elem_fields(elem_gid, 6) = speed.host(elem_gid);
            elem_fields(elem_gid, 8) = State.GaussPoints.div.host(elem_gid);
            elem_switch *= -1;
        } // end for elem_gid

        //collective map has all indices on rank 0 and non on other ranks
        HostDistributedMap collective_elem_map;
        long long int num_collective_elem_indices = 0;
        if(myrank==0) num_collective_elem_indices = mesh.global_num_elems;
        collective_elem_map = HostDistributedMap(mesh.global_num_elems, num_collective_elem_indices);

        //host version of local element map for argument compatibility
        HostDistributedMap host_node_map;
        DCArrayKokkos<long long int, Kokkos::LayoutLeft , Kokkos::HostSpace> global_indices_of_local_nodes;
        for(int inode = 0; inode < mesh.num_local_nodes; inode++){
            global_indices_of_local_nodes(inode) = mesh.node_map.getGlobalIndex(inode);
        }
        host_node_map = HostDistributedMap(global_indices_of_local_nodes);

        // save the vertex vector fields to an array for exporting to graphics files
        DistributedCArray<double> vec_fields(host_node_map, num_point_vec_vars, 3);
        DistributedCArray<double> point_scalar_fields(host_node_map, num_point_scalar_vars);

        for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
            // position, var 0
            vec_fields(node_gid, 0, 0) = State.node.coords.host(node_gid, 0);
            vec_fields(node_gid, 0, 1) = State.node.coords.host(node_gid, 1);
            if (num_dims == 2) {
                vec_fields(node_gid, 0, 2) = 0.0;
            }
            else{
                vec_fields(node_gid, 0, 2) = State.node.coords.host(node_gid, 2);
            }

            // position, var 1
            vec_fields(node_gid, 1, 0) = State.node.vel.host(node_gid, 0);
            vec_fields(node_gid, 1, 1) = State.node.vel.host(node_gid, 1);
            if (num_dims == 2) {
                vec_fields(node_gid, 1, 2) = 0.0;
            }
            else{
                vec_fields(node_gid, 1, 2) = State.node.vel.host(node_gid, 2);
            }

            point_scalar_fields(node_gid, 0) = State.node.temp.host(node_gid);
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

        //sprintf(filename, "vtk/Fierro.%05d.vtk", graphics_id);  // mesh file
        str_output_len = snprintf(filename, max_len, "vtk/Fierro.%05d.vtk", graphics_id);
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
                    State.node.coords.host(node_gid, 0),
                    State.node.coords.host(node_gid, 1),
                    State.node.coords.host(node_gid, 2));
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
        //sprintf(filename, "vtk/Fierro.vtk.series", graphics_id);  // mesh file
        str_output_len = snprintf(filename, max_len, "vtk/Fierro.vtk.series"); 
        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }
        // mesh file

        out[0] = fopen(filename, "w");

        fprintf(out[0], "{\n");
        fprintf(out[0], "  \"file-series-version\" : \"1.0\",\n");
        fprintf(out[0], "  \"files\" : [\n");

        for (int i = 0; i <= graphics_id; i++) {
            fprintf(out[0], "    { \"name\" : \"Fierro.%05d.vtk\", \"time\" : %12.5e },\n", i, graphics_times(i) );
        }

        // fprintf(out[0], "%12.5e\n", graphics_times(i));
        fprintf(out[0], "  ]\n"); // part 4
        fprintf(out[0], "}"); // part 4

        fclose(out[0]);

        // increment graphics id counter
        graphics_id++;


    } // end write vtk old


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn concatenate_elem_fields
    ///
    /// \brief A function to calculate the average of elem fields and concatentate into 1 array
    ///
    ///
    /// \param MaterialPoints a struct containing the material point state arrays
    /// \param elem_scalar_fields the scalar fields
    /// \param elem_tensor_fields the tensor fields
    /// \param elem_in_mat_elem a listing of the element ids the material resides in
    /// \param output_elem_state a std::vector of enums specifying the elem avg outputs
    /// \param num_mat_elems the number of elements the material resides in
    /// \param mat_id the index for the material
    ///
    /////////////////////////////////////////////////////////////////////////////
    void concatenate_elem_fields(const MaterialPoint_t& MaterialPoints,
                                 const GaussPoint_t& GaussPoints,
                                 DistributedDFArray<double>& elem_scalar_fields,
                                 DistributedDFArray<double>& elem_tensor_fields,
                                 const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
                                 const std::vector<material_pt_state>& output_elem_state,
                                 const std::vector<gauss_pt_state>& output_gauss_pt_states,
                                 const size_t num_mat_elems,
                                 const size_t mat_id,
                                 const size_t num_elems,
                                 const int den_id,
                                 const int pres_id,
                                 const int sie_id,
                                 const int sspd_id,
                                 const int mass_id,
                                 const int stress_id,
                                 const int vol_id,
                                 const int div_id,
                                 const int level_set_id,
                                 const int vel_grad_id,
                                 const int conductivity_id,
                                 const int specific_heat_id)
    {

        // --- loop over the material point states

        for (auto field : output_elem_state){
            switch(field){
                // scalar vars
                case material_pt_state::density:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // get elem gid
                        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_lid);

                        // field
                        elem_scalar_fields(elem_gid, den_id) += MaterialPoints.den(mat_id, mat_elem_lid)*
                                                                MaterialPoints.volfrac(mat_id, mat_elem_lid)*
                                                                MaterialPoints.geo_volfrac(mat_id, mat_elem_lid);
                    });
                    break;
                case material_pt_state::pressure:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // get elem gid
                        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_lid);

                        // field
                        elem_scalar_fields(elem_gid, pres_id) += MaterialPoints.pres(mat_id, mat_elem_lid)*
                                                                MaterialPoints.volfrac(mat_id, mat_elem_lid)*
                                                                MaterialPoints.geo_volfrac(mat_id, mat_elem_lid);
                    });
                    break;
                case material_pt_state::specific_internal_energy:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // get elem gid
                        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_lid);

                        // field
                        // extensive ie here, but after this function, it will become specific ie
                        elem_scalar_fields(elem_gid, sie_id) += MaterialPoints.mass(mat_id, mat_elem_lid)*
                                                                MaterialPoints.sie(mat_id, mat_elem_lid);
                    });
                    break;
                case material_pt_state::sound_speed:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // get elem gid
                        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_lid);

                        // field
                        elem_scalar_fields(elem_gid, sspd_id) += MaterialPoints.sspd(mat_id, mat_elem_lid)*
                                                                MaterialPoints.volfrac(mat_id, mat_elem_lid)*
                                                                MaterialPoints.geo_volfrac(mat_id, mat_elem_lid);
                    });
                    break;
                case material_pt_state::mass:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // get elem gid
                        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_lid);

                        // field
                        elem_scalar_fields(elem_gid, mass_id) += MaterialPoints.mass(mat_id, mat_elem_lid);
                    });
                    break;
                // ---------------    
                // tensor vars
                // ---------------
                case material_pt_state::stress:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // get elem gid
                        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_lid);

                        // field
                        // average tensor fields, it is always 3D
                        // note: paraview is row-major, CArray convention
                        for (size_t i=0; i<3; i++){
                            for(size_t j=0; j<3; j++){

                                // stress tensor 
                                elem_tensor_fields(elem_gid, stress_id, i, j) +=
                                                MaterialPoints.stress(mat_id, mat_elem_lid,i,j) *
                                                MaterialPoints.volfrac(mat_id, mat_elem_lid)*
                                                MaterialPoints.geo_volfrac(mat_id, mat_elem_lid);
                            } // end for
                        } // end for
                    });
                    break;

                // thermal solver vars
                case material_pt_state::thermal_conductivity:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // get elem gid
                        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_lid);

                        // field
                        elem_scalar_fields(elem_gid, conductivity_id) += MaterialPoints.conductivity(mat_id, mat_elem_lid)*
                                                                             MaterialPoints.volfrac(mat_id, mat_elem_lid)*
                                                                             MaterialPoints.geo_volfrac(mat_id, mat_elem_lid);
                    });
                    break;

                case material_pt_state::specific_heat:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // get elem gid
                        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_lid);

                        // field
                        elem_scalar_fields(elem_gid, specific_heat_id) += MaterialPoints.specific_heat(mat_id, mat_elem_lid)*
                                                                              MaterialPoints.volfrac(mat_id, mat_elem_lid)*
                                                                              MaterialPoints.geo_volfrac(mat_id, mat_elem_lid);
                    });
                    break;


                // add other variables here

                // not used variables
                case material_pt_state::volume_fraction:
                    break;
                case material_pt_state::eroded_flag:
                    break;
                case material_pt_state::elastic_modulii:
                    break;
                case material_pt_state::shear_modulii:
                    break;
                case material_pt_state::poisson_ratios:
                    break;
                case material_pt_state::heat_flux:
                    break;
            } // end switch
        }// end for over mat point state

        
        // --- add loop over gauss points ---

        // export element centric data
        for (auto field : output_gauss_pt_states){
            switch(field){
                // scalars
                case gauss_pt_state::volume:

                    FOR_ALL(elem_gid, 0, num_elems, {
                        elem_scalar_fields(elem_gid, vol_id) = GaussPoints.vol(elem_gid);
                    });

                    break;
                case gauss_pt_state::divergence_velocity:

                    FOR_ALL(elem_gid, 0, num_elems, {
                        elem_scalar_fields(elem_gid, div_id) = GaussPoints.div(elem_gid);
                    });

                    break;

                case gauss_pt_state::level_set:

                    FOR_ALL(elem_gid, 0, num_elems, {
                        elem_scalar_fields(elem_gid, level_set_id) = GaussPoints.level_set(elem_gid);
                    });

                    break;

                // tensors
                case gauss_pt_state::gradient_velocity:
                    // note: paraview is row-major, CArray convention
                    FOR_ALL(elem_gid, 0, num_elems, {
                        for (size_t i=0; i<3; i++){
                            for(size_t j=0; j<3; j++){
                                elem_tensor_fields(elem_gid, vel_grad_id, i, j) = 
                                                    GaussPoints.vel_grad(elem_gid, i, j);
                            }
                        } // end for
                    });

                    break;

                // add other gauss variables here

            } // end switch
        } // end loop over gauss_pt_states


        // --- add end gauss point loop --

    } // end of function

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn copy_elem_fields
    ///
    /// \brief A function to assign dual values to host only array for tpetra compatibility
    ///
    ///
    /// \param MaterialPoints a struct containing the material point state arrays
    /// \param elem_scalar_fields the scalar fields
    /// \param elem_tensor_fields the tensor fields
    /// \param elem_in_mat_elem a listing of the element ids the material resides in
    /// \param output_elem_state a std::vector of enums specifying the elem avg outputs
    /// \param num_mat_elems the number of elements the material resides in
    /// \param mat_id the index for the material
    ///
    /////////////////////////////////////////////////////////////////////////////
    void copy_elem_fields(DistributedDFArray<double>& elem_scalar_fields,
                          DistributedDFArray<double>& elem_tensor_fields,
                          DistributedFArray<double>& host_elem_scalar_fields,
                          DistributedFArray<double>& host_elem_tensor_fields,
                          const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
                          const std::vector<material_pt_state>& output_elem_state,
                          const std::vector<gauss_pt_state>& output_gauss_pt_states,
                          const size_t num_mat_elems,
                          const size_t mat_id,
                          const size_t num_elems,
                          const int den_id,
                          const int pres_id,
                          const int sie_id,
                          const int sspd_id,
                          const int mass_id,
                          const int stress_id,
                          const int vol_id,
                          const int div_id,
                          const int level_set_id,
                          const int vel_grad_id,
                          const int conductivity_id,
                          const int specific_heat_id)
    {

        // --- loop over the material point states

        for (auto field : output_elem_state){
            switch(field){
                // scalar vars
                case material_pt_state::density:
                    for(int elem_gid = 0; elem_gid < num_elems; elem_gid++){

                        // field
                        host_elem_scalar_fields(elem_gid, den_id) = elem_scalar_fields.host(elem_gid, den_id);
                    }
                    break;
                case material_pt_state::pressure:
                    for(int elem_gid = 0; elem_gid < num_elems; elem_gid++){

                        // field
                        host_elem_scalar_fields(elem_gid, pres_id) = elem_scalar_fields.host(elem_gid, pres_id);
                    }
                    break;
                case material_pt_state::specific_internal_energy:
                    for(int elem_gid = 0; elem_gid < num_elems; elem_gid++){

                        // field
                        // extensive ie here, but after this function, it will become specific ie
                        host_elem_scalar_fields(elem_gid, sie_id) = elem_scalar_fields.host(elem_gid, sie_id);
                    }
                    break;
                case material_pt_state::sound_speed:
                    for(int elem_gid = 0; elem_gid < num_elems; elem_gid++){

                        // field
                        host_elem_scalar_fields(elem_gid, sspd_id) = elem_scalar_fields.host(elem_gid, sspd_id);
                    }
                    break;
                case material_pt_state::mass:
                    for(int elem_gid = 0; elem_gid < num_elems; elem_gid++){

                        // field
                        host_elem_scalar_fields(elem_gid, mass_id) = elem_scalar_fields.host(elem_gid, mass_id);
                    }
                    break;
                // ---------------    
                // tensor vars
                // ---------------
                case material_pt_state::stress:
                    for(int elem_gid = 0; elem_gid < num_elems; elem_gid++){

                        // field
                        // average tensor fields, it is always 3D
                        // note: paraview is row-major, CArray convention
                        for (size_t i=0; i<3; i++){
                            for(size_t j=0; j<3; j++){

                                // stress tensor 
                                host_elem_tensor_fields(elem_gid, stress_id, i, j) = elem_tensor_fields.host(elem_gid, stress_id, i, j);
                            } // end for
                        } // end for
                    }
                    break;

                // thermal solver vars
                case material_pt_state::thermal_conductivity:
                    for(int elem_gid = 0; elem_gid < num_elems; elem_gid++){

                        // field
                        host_elem_scalar_fields(elem_gid, conductivity_id) = elem_scalar_fields.host(elem_gid, conductivity_id);
                    }
                    break;

                case material_pt_state::specific_heat:
                    for(int elem_gid = 0; elem_gid < num_elems; elem_gid++){

                        // field
                        host_elem_scalar_fields(elem_gid, specific_heat_id) = elem_scalar_fields.host(elem_gid, specific_heat_id);
                    }
                    break;


                // add other variables here

                // not used variables
                case material_pt_state::volume_fraction:
                    break;
                case material_pt_state::eroded_flag:
                    break;
                case material_pt_state::elastic_modulii:
                    break;
                case material_pt_state::shear_modulii:
                    break;
                case material_pt_state::poisson_ratios:
                    break;
                case material_pt_state::heat_flux:
                    break;
            } // end switch
        }// end for over mat point state

        
        // --- add loop over gauss points ---

        // export element centric data
        for (auto field : output_gauss_pt_states){
            switch(field){
                // scalars
                case gauss_pt_state::volume:

                    for(int elem_gid = 0; elem_gid < num_elems; elem_gid++){
                        host_elem_scalar_fields(elem_gid, vol_id) = elem_scalar_fields.host(elem_gid, vol_id);
                    }

                    break;
                case gauss_pt_state::divergence_velocity:

                    for(int elem_gid = 0; elem_gid < num_elems; elem_gid++){
                        host_elem_scalar_fields(elem_gid, div_id) = elem_scalar_fields.host(elem_gid, div_id);
                    }

                    break;

                case gauss_pt_state::level_set:

                    for(int elem_gid = 0; elem_gid < num_elems; elem_gid++){
                        host_elem_scalar_fields(elem_gid, level_set_id) = elem_scalar_fields.host(elem_gid, level_set_id);
                    }

                    break;

                // tensors
                case gauss_pt_state::gradient_velocity:
                    // note: paraview is row-major, CArray convention
                    for(int elem_gid = 0; elem_gid < num_elems; elem_gid++){
                        for (size_t i=0; i<3; i++){
                            for(size_t j=0; j<3; j++){
                                host_elem_tensor_fields(elem_gid, vel_grad_id, i, j) = 
                                                    elem_tensor_fields.host(elem_gid, vel_grad_id, i, j);
                            }
                        } // end for
                    }

                    break;

                // add other gauss variables here

            } // end switch
        } // end loop over gauss_pt_states


        // --- add end gauss point loop --

    } // end of function

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn concatenate_mat_fields
    ///
    /// \brief A function to concatentate material fields into 1 array
    ///
    ///
    /// \param MaterialPoints a struct containing the material point state arrays
    /// \param elem_scalar_fields the scalar fields
    /// \param elem_tensor_fields the tensor fields
    /// \param elem_in_mat_elem a listing of the element ids the material resides in
    /// \param output_material_pt_states a std::vector of enums specifying the model
    /// \param num_mat_elems the number of elements the material resides in
    /// \param mat_id the index for the material
    ///
    /////////////////////////////////////////////////////////////////////////////
    void concatenate_mat_fields(const MaterialPoint_t& MaterialPoints,
                                DistributedDFArray<double>& mat_elem_scalar_fields,
                                DistributedDFArray<double>& mat_elem_tensor_fields,
                                const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
                                const std::vector<material_pt_state>& output_material_pt_states,
                                const size_t num_mat_elems,
                                const size_t mat_id,
                                const int mat_den_id,
                                const int mat_pres_id,
                                const int mat_sie_id,
                                const int mat_sspd_id,
                                const int mat_mass_id,
                                const int mat_volfrac_id,  
                                const int mat_geo_volfrac_id,  
                                const int mat_eroded_id,
                                const int mat_stress_id,
                                const int mat_conductivity_id,
                                const int mat_specific_heat_id)
    {
      
        for (auto field : output_material_pt_states){
            switch(field){
                // scalar vars
                case material_pt_state::density:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // field
                        mat_elem_scalar_fields(mat_elem_lid, mat_den_id) = MaterialPoints.den(mat_id, mat_elem_lid);
                    });
                    break;
                case material_pt_state::pressure:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // field
                        mat_elem_scalar_fields(mat_elem_lid, mat_pres_id) = MaterialPoints.pres(mat_id, mat_elem_lid);
                    });
                    break;
                case material_pt_state::specific_internal_energy:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // field
                        // extensive ie here, but after this function, it will become specific ie
                        mat_elem_scalar_fields(mat_elem_lid, mat_sie_id) = MaterialPoints.sie(mat_id, mat_elem_lid);
                    });
                    break;
                case material_pt_state::sound_speed:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // field
                        mat_elem_scalar_fields(mat_elem_lid, mat_sspd_id) = MaterialPoints.sspd(mat_id, mat_elem_lid);
                    });
                    break;
                case material_pt_state::mass:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // field
                        mat_elem_scalar_fields(mat_elem_lid, mat_mass_id) = MaterialPoints.mass(mat_id, mat_elem_lid);
                    });
                    break;
                case material_pt_state::volume_fraction:
                    // material volume fraction
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // field
                        // this is the volume fraction of a material within a part
                        mat_elem_scalar_fields(mat_elem_lid, mat_volfrac_id) = MaterialPoints.volfrac(mat_id, mat_elem_lid);
                    });

                    // geometric volume fraction
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // field
                        // this is the geometric volume fraction (interface reconstruction)
                        mat_elem_scalar_fields(mat_elem_lid, mat_geo_volfrac_id) = MaterialPoints.geo_volfrac(mat_id, mat_elem_lid);
                    });
                    break;
                case material_pt_state::eroded_flag:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // field
                        mat_elem_scalar_fields(mat_elem_lid, mat_eroded_id) = (double)MaterialPoints.eroded(mat_id, mat_elem_lid);
                    });
                    break;
                // ---------------    
                // tensor vars
                // ---------------
                case material_pt_state::stress:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // field
                        // average tensor fields, it is always 3D
                        // note: paraview is row-major, CArray convention
                        for (size_t i=0; i<3; i++){
                            for(size_t j=0; j<3; j++){

                                // stress tensor 
                                mat_elem_tensor_fields(mat_elem_lid, mat_stress_id, i, j) =
                                                MaterialPoints.stress(mat_id, mat_elem_lid,i,j);
                            } // end for
                        } // end for
                    });
                    break;

                // thermal solver vars
                case material_pt_state::thermal_conductivity:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // get elem gid
                        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_lid);

                        // field
                        mat_elem_scalar_fields(elem_gid, mat_conductivity_id) += MaterialPoints.conductivity(mat_id, mat_elem_lid);
                    });
                    break;

                case material_pt_state::specific_heat:
                    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

                        // get elem gid
                        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_lid);

                        // field
                        mat_elem_scalar_fields(elem_gid, mat_specific_heat_id) += MaterialPoints.specific_heat(mat_id, mat_elem_lid);
                    });
                    break;

                // add other variables here

                // not used variables
                case material_pt_state::elastic_modulii:
                    break;
                case material_pt_state::shear_modulii:
                    break;
                case material_pt_state::poisson_ratios:
                    break;
                case material_pt_state::heat_flux:
                    break;
            } // end switch
        }// end for over mat point state
            


    } // end of function
    
    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn concatenate_nodal_fields
    ///
    /// \brief A function to calculate the average of nodal fields
    ///
    ///
    /// \param Node a struct containing the material point state arrays
    /// \param elem_scalar_fields the scalar fields
    /// \param elem_tensor_fields the tensor fields
    /// \param elem_in_mat_elem a listing of the element ids the material resides in
    /// \param output_node_states a std::vector of enums specifying the model
    /// \param num_mat_elems the number of elements the material resides in
    /// \param mat_id the index for the material
    ///
    /////////////////////////////////////////////////////////////////////////////
    void concatenate_nodal_fields(const node_t& Node,
                                  DistributedDFArray<double>& node_scalar_fields,
                                  DistributedDFArray<double>& node_vector_fields,
                                  std::vector<node_state>& output_node_states,
                                  double dt,
                                  const size_t num_nodes,
                                  const size_t num_dims,
                                  const int node_mass_id,
                                  const int node_vel_id,
                                  const int node_accel_id,
                                  const int node_coord_id,
                                  const int node_grad_level_set_id,
                                  const int node_temp_id)
    {
        for (auto field : output_node_states){
            switch(field){
                // scalars
                case node_state::mass:

                    FOR_ALL(node_gid, 0, num_nodes, {
                        node_scalar_fields(node_gid, node_mass_id) = Node.mass(node_gid);
                    });

                    break;
                case node_state::temp:
                    FOR_ALL(node_gid, 0, num_nodes, {
                        node_scalar_fields(node_gid, node_temp_id) = Node.temp(node_gid);
                    });

                    break;

                // vector fields

                case node_state::coords:

                    FOR_ALL(node_gid, 0, num_nodes, {

                        node_vector_fields(node_gid, node_coord_id, 0) = Node.coords(node_gid, 0);
                        node_vector_fields(node_gid, node_coord_id, 1) = Node.coords(node_gid, 1);
                        if (num_dims == 2) {
                            node_vector_fields(node_gid, node_coord_id, 2) = 0.0;
                        }
                        else{
                            node_vector_fields(node_gid, node_coord_id, 2) = Node.coords(node_gid, 2);
                        } // end if

                    }); // end parallel for

                    break;
                case node_state::velocity:

                    FOR_ALL(node_gid, 0, num_nodes, {

                        // velocity, var is node_vel_id 
                        node_vector_fields(node_gid, node_vel_id, 0) = Node.vel(node_gid, 0);
                        node_vector_fields(node_gid, node_vel_id, 1) = Node.vel(node_gid, 1);
                        if (num_dims == 2) {
                            node_vector_fields(node_gid, node_vel_id, 2) = 0.0;
                        }
                        else{
                            node_vector_fields(node_gid, node_vel_id, 2) = Node.vel(node_gid, 2);
                        } // end if

                        // accellerate, var is node_accel_id            
                        node_vector_fields(node_gid, node_accel_id, 0) = (Node.vel(node_gid, 0) - Node.vel_n0(node_gid, 0))/dt;
                        node_vector_fields(node_gid, node_accel_id, 1) = (Node.vel(node_gid, 1) - Node.vel_n0(node_gid, 1))/dt;
                        if (num_dims == 2) {
                            node_vector_fields(node_gid, node_accel_id, 2) = 0.0;
                        }
                        else{
                            node_vector_fields(node_gid, node_accel_id, 2) = (Node.vel(node_gid, 2) - Node.vel_n0(node_gid, 2))/dt;
                        } // end if

                    }); // end parallel for

                    break;
                    
                    
                case node_state::gradient_level_set:

                    FOR_ALL(node_gid, 0, num_nodes, {

                        // velocity, var is node_vel_id 
                        node_vector_fields(node_gid, node_grad_level_set_id, 0) = Node.gradient_level_set(node_gid, 0);
                        node_vector_fields(node_gid, node_grad_level_set_id, 1) = Node.gradient_level_set(node_gid, 1);
                        if (num_dims == 2) {
                            node_vector_fields(node_gid, node_grad_level_set_id, 2) = 0.0;
                        }
                        else{
                            node_vector_fields(node_gid, node_grad_level_set_id, 2) = Node.gradient_level_set(node_gid, 2);
                        } // end if

                    }); // end parallel for

                    break;                
                
                // -- not used vars
                case node_state::force:
                    break;

                // heat transer vars
                case node_state::heat_transfer:
                    break;
                // tensors
            } // end switch
        } // end for over

        

    } // end function

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn copy_nodal_fields
    ///
    /// \brief A function to calculate the average of nodal fields
    ///
    ///
    /// \param Node a struct containing the material point state arrays
    /// \param elem_scalar_fields the scalar fields
    /// \param elem_tensor_fields the tensor fields
    /// \param elem_in_mat_elem a listing of the element ids the material resides in
    /// \param output_node_states a std::vector of enums specifying the model
    /// \param num_mat_elems the number of elements the material resides in
    /// \param mat_id the index for the material
    ///
    /////////////////////////////////////////////////////////////////////////////
    void copy_nodal_fields(DistributedFArray<double>& host_node_scalar_fields,
                            DistributedFArray<double>& host_node_vector_fields,
                            DistributedDFArray<double>& node_scalar_fields,
                            DistributedDFArray<double>& node_vector_fields,
                            std::vector<node_state>& output_node_states,
                            double dt,
                            const size_t num_nodes,
                            const size_t num_dims,
                            const int node_mass_id,
                            const int node_vel_id,
                            const int node_accel_id,
                            const int node_coord_id,
                            const int node_grad_level_set_id,
                            const int node_temp_id)
    {
        for (auto field : output_node_states){
            switch(field){
                // scalars
                case node_state::mass:

                    for(long long int node_gid = 0; node_gid < num_nodes; node_gid++) {
                        host_node_scalar_fields(node_gid, node_mass_id) = node_scalar_fields.host(node_gid, node_mass_id);
                    }

                    break;
                case node_state::temp:
                    for(long long int node_gid = 0; node_gid < num_nodes; node_gid++) {
                        host_node_scalar_fields(node_gid, node_temp_id) = node_scalar_fields.host(node_gid, node_temp_id);
                    }

                    break;

                // vector fields

                case node_state::coords:

                    for(long long int node_gid = 0; node_gid < num_nodes; node_gid++) {

                        host_node_vector_fields(node_gid, node_coord_id, 0) = node_vector_fields.host(node_gid, node_coord_id, 0);
                        host_node_vector_fields(node_gid, node_coord_id, 1) = node_vector_fields.host(node_gid, node_coord_id, 1);
                        if (num_dims == 2) {
                            host_node_vector_fields(node_gid, node_coord_id, 2) = 0.0;
                        }
                        else{
                            host_node_vector_fields(node_gid, node_coord_id, 2) = node_vector_fields.host(node_gid, node_coord_id, 2);
                        } // end if

                    } // end parallel for

                    break;
                case node_state::velocity:

                    for(long long int node_gid = 0; node_gid < num_nodes; node_gid++) {

                        // velocity, var is node_vel_id 
                        host_node_vector_fields(node_gid, node_vel_id, 0) = node_vector_fields.host(node_gid, node_vel_id, 0);
                        host_node_vector_fields(node_gid, node_vel_id, 1) = node_vector_fields.host(node_gid, node_vel_id, 1);
                        if (num_dims == 2) {
                            host_node_vector_fields(node_gid, node_vel_id, 2) = 0.0;
                        }
                        else{
                            host_node_vector_fields(node_gid, node_vel_id, 2) = node_vector_fields.host(node_gid, node_vel_id, 2);
                        } // end if

                        // accellerate, var is node_accel_id            
                        host_node_vector_fields(node_gid, node_accel_id, 0) = node_vector_fields.host(node_gid, node_accel_id, 0);
                        host_node_vector_fields(node_gid, node_accel_id, 1) = node_vector_fields.host(node_gid, node_accel_id, 1);
                        if (num_dims == 2) {
                            host_node_vector_fields(node_gid, node_accel_id, 2) = 0.0;
                        }
                        else{
                            host_node_vector_fields(node_gid, node_accel_id, 2) = node_vector_fields.host(node_gid, node_accel_id, 2);
                        } // end if

                    } // end parallel for

                    break;
                    
                    
                case node_state::gradient_level_set:

                    for(long long int node_gid = 0; node_gid < num_nodes; node_gid++) {

                        // velocity, var is node_vel_id 
                        host_node_vector_fields(node_gid, node_grad_level_set_id, 0) = node_vector_fields.host(node_gid, node_grad_level_set_id, 0);
                        host_node_vector_fields(node_gid, node_grad_level_set_id, 1) = node_vector_fields.host(node_gid, node_grad_level_set_id, 1);
                        if (num_dims == 2) {
                            host_node_vector_fields(node_gid, node_grad_level_set_id, 2) = 0.0;
                        }
                        else{
                            host_node_vector_fields(node_gid, node_grad_level_set_id, 2) = node_vector_fields.host(node_gid, node_grad_level_set_id, 2);
                        } // end if

                    } // end parallel for

                    break;                
                
                // -- not used vars
                case node_state::force:
                    break;

                // heat transer vars
                case node_state::heat_transfer:
                    break;
                // tensors
            } // end switch
        } // end for over

        

    } // end function

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn write_vtu
    ///
    /// \brief Writes a vtu ASCII output file
    ///
    /// \param Simulation mesh
    /// \param State data
    /// \param Simulation parameters
    /// \param current time value
    /// \param Vector of all graphics output times
    ///
    /////////////////////////////////////////////////////////////////////////////
    void write_vtu(
        const DistributedDCArray<double>& node_coords,
        const DistributedDCArray<size_t>& nodes_in_elem,
        const DistributedDFArray<double>& elem_scalar_fields,
        const DistributedDFArray<double>& elem_tensor_fields,
        const DistributedDFArray<double>& node_scalar_fields,
        const DistributedDFArray<double>& node_vector_fields,
        const std::vector<std::string>& elem_scalar_var_names,
        const std::vector<std::string>& elem_tensor_var_names,
        const std::vector<std::string>& node_scalar_var_names,
        const std::vector<std::string>& node_vector_var_names,
        const std::string partname,
        const int graphics_id,
        const size_t num_nodes,
        const size_t num_elems,
        const size_t num_nodes_in_elem,
        const int Pn_order,
        const size_t num_dims,
        const size_t solver_id
        )
    {   
        
        FILE* out[20];   // the output files that are written to
        char  filename[100]; // char string
        int   max_len = sizeof filename;
        int   str_output_len;
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        std::string str_rank_val = std::to_string(myrank);                       
        std::string rank_fields_name = "rank";
        rank_fields_name += str_rank_val;  // add the mat number

        const size_t num_elem_scalar_vars = elem_scalar_var_names.size();
        const size_t num_elem_tensor_vars = elem_tensor_var_names.size();

        const size_t num_node_scalar_vars = node_scalar_var_names.size();
        const size_t num_node_vector_vars = node_vector_var_names.size();


        // create filename
        str_output_len = snprintf(filename, max_len, "vtk/data/Fierro.solver%zu.%s_%s.%05d.vtu", 
                                                                 solver_id, partname.c_str(), rank_fields_name.c_str(), graphics_id);

        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }
        // mesh file
        
        out[0] = fopen(filename, "w");

        fprintf(out[0], "<?xml version=\"1.0\"?>\n");  
        fprintf(out[0], "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"); 
        fprintf(out[0], "  <UnstructuredGrid> \n");
        fprintf(out[0], "    <Piece NumberOfPoints=\"%zu\" NumberOfCells=\"%zu\">\n", num_nodes, num_elems); 

        /*
        ---------------------------------------------------------------------------
        Write the mesh points
        ---------------------------------------------------------------------------
        */
        fprintf(out[0], "\n");
        fprintf(out[0], "      <!-- Define the mesh nodes -->\n");
        fprintf(out[0], "      <Points>\n");
        fprintf(out[0], "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");

        // write all components of the point coordinates
        for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
            double coord_z = 0.0;
            if(num_dims==3){
                coord_z = node_coords.host(node_gid, 2);
            } 
            fprintf(out[0],
                    "          %f %f %f\n",
                    node_coords.host(node_gid, 0),
                    node_coords.host(node_gid, 1),
                    coord_z);
        } // end for
        fprintf(out[0], "        </DataArray>\n");
        fprintf(out[0], "      </Points>\n");

        /*
        ---------------------------------------------------------------------------
        Write the elems
        ---------------------------------------------------------------------------
        */
        fprintf(out[0], "\n");
        fprintf(out[0], "      <!-- Define the elements -->\n");
        fprintf(out[0], "      <Cells>\n");
        fprintf(out[0], "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");  

        // WARNING: look into high-order Pn 2D elements with paraview
        int Pn_order_z = 0;
        if (num_dims == 3){
            Pn_order_z = Pn_order;
        }
        int order[3] = {Pn_order, Pn_order, Pn_order_z};

        // const int num_1D_points = Pn_order+1;

        // write all global point numbers for this elem
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            fprintf(out[0], "          ");  // adding indentation before printing nodes in element
            if (num_dims==3 && Pn_order>1){
                for (int k = 0; k <= Pn_order_z; k++) {
                    for (int j = 0; j <= Pn_order; j++) {
                        for (int i = 0; i <= Pn_order; i++) {
                            size_t node_lid = PointIndexFromIJK(i, j, k, order);
                            fprintf(out[0], "%lu ", nodes_in_elem.host(elem_gid, node_lid));
                        }
                    }
                } // end for
            }
            else if (num_dims == 3 && Pn_order == 1){
               // 3D linear hexahedral elements
                for (int node_lid = 0; node_lid < 8; node_lid++) {
                    fprintf(out[0], "%lu ", nodes_in_elem.host(elem_gid, node_lid));
                } // end for
            }
            else if (num_dims == 2){
                // 2D linear is the only supported option
                for (int node_lid = 0; node_lid < 4; node_lid++) {
                    fprintf(out[0], "%lu ", nodes_in_elem.host(elem_gid, node_lid));
                } // end for
            }
            else {
                std::cout << "ERROR: outputs failed, dimensions and element types are not compatible \n";
            } // end if
            fprintf(out[0], "\n");
        } // end for
        fprintf(out[0], "        </DataArray>\n");

        // Write the element offsets
        fprintf(out[0], "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");  
        size_t count=0;
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            count += num_nodes_in_elem;
            fprintf(out[0], "          %lu\n", count); // num points in this elem + all others before it
        } // end for
        fprintf(out[0], "        </DataArray>\n");


        // Write the element types
        fprintf(out[0], "        <DataArray type=\"Int8\" Name=\"types\" format=\"ascii\">\n"); 
        // ----
        // linear element types
        //   VTK_PIXEL = 8,   linear 2D quad with i,j,k indexing (future format for 2D solver)
        //   VTK_Quad = 9,    linear 2D quad with ensight index ordering (current 2D rz convention)
        //   VTK_VOXEL = 11,  linear 3D hex with i,j,k indexing (current format)
        // arbitrary order types
        //   VTK_LAGRANGE_QUADRILATERAL = 70, use this type when a 2D high-order scheme exists
        //   VTK_LAGRANGE_HEXAHEDRON: 72, this is the current 3D high-order 
        //   VTK_HIGHER_ORDER_HEXAHEDRON: 67
        //   VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33
        // element types: https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
        // element types: https://kitware.github.io/vtk-js/api/Common_DataModel_CellTypes.html
        // vtk format: https://www.kitware.com//modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            if (num_dims==3 && Pn_order>1){
                fprintf(out[0], "          %d \n", 72);
            }
            else if (num_dims == 3 && Pn_order == 1){
                // 3D linear hex
                fprintf(out[0], "          %d \n", 11);
            }
            else {
                // 2D ensight mesh ordering
                fprintf(out[0], "          %d \n", 9);
            }
        }
        fprintf(out[0], "        </DataArray>\n");
        fprintf(out[0], "      </Cells>\n");


        /*
        ---------------------------------------------------------------------------
        Write the nodal variables to file
        ---------------------------------------------------------------------------
        */
        // vtk vector vars = (position, velocity)
        fprintf(out[0], "\n");
        fprintf(out[0], "      <!-- Define the node vector data -->\n");
        if(num_node_vector_vars >0 || num_node_scalar_vars>0){

            fprintf(out[0], "      <PointData>\n");

            // node vectors
            for (int a_var = 0; a_var < num_node_vector_vars; a_var++) {
                fprintf(out[0], "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\">\n", node_vector_var_names[a_var].c_str());
               
                for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
                    fprintf(out[0], "          %f %f %f\n",
                            node_vector_fields.host(node_gid, a_var, 0),
                            node_vector_fields.host(node_gid, a_var, 1),
                            node_vector_fields.host(node_gid, a_var, 2));
                } // end for nodes
                fprintf(out[0], "        </DataArray>\n");

            } // end for vec_vars


            // node scalar vars
            for (int a_var = 0; a_var < num_node_scalar_vars; a_var++) {
                fprintf(out[0], "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", node_scalar_var_names[a_var].c_str());
                for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
                    fprintf(out[0], "          %f\n", node_scalar_fields.host(node_gid, a_var));
                } // end for nodes
                fprintf(out[0], "        </DataArray>\n");
            } // end for vec_vars

            fprintf(out[0], "      </PointData>\n");

        } // end if

        /*
        ---------------------------------------------------------------------------
        Write the elem variables to file
        ---------------------------------------------------------------------------
        */
        fprintf(out[0], "\n");
        fprintf(out[0], "      <!-- Define the cell data -->\n");
        if(num_elem_scalar_vars >0 || num_elem_tensor_vars>0){

            fprintf(out[0], "      <CellData>\n");

            for (int a_var = 0; a_var < num_elem_scalar_vars; a_var++) {

                fprintf(out[0], "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", elem_scalar_var_names[a_var].c_str()); // the 1 is number of scalar components [1:4]

                for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
                    fprintf(out[0], "          %f\n", elem_scalar_fields.host(elem_gid, a_var));
                } // end for elem
                fprintf(out[0], "        </DataArray>\n");
            } // end for elem scalar_vars


            // tensors
            for (int a_var = 0; a_var < num_elem_tensor_vars; a_var++) {
                fprintf(out[0], "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"9\" format=\"ascii\">\n", elem_tensor_var_names[a_var].c_str()); // the 1 is number of scalar components [1:4]
                
                for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
                    // note: paraview is row-major, CArray convention
                    // Txx  Txy  Txz  Tyx  Tyy  Tyz  Tzx  Tzy  Tzz
                    for (size_t i=0; i<3; i++){
                        for(size_t j=0; j<3; j++){
                            fprintf(out[0], "          %f ", elem_tensor_fields.host(elem_gid, a_var, i, j));
                        } // end j
                    } // end i
                } // end for elem
                fprintf(out[0], "\n");
                fprintf(out[0], "        </DataArray>\n");
            } // end for elem scalar_vars

            fprintf(out[0], "      </CellData>\n");
        } // end if

        // end of the vtu file
        fprintf(out[0], "    </Piece>\n");
        fprintf(out[0], "  </UnstructuredGrid>\n");
        fprintf(out[0], "</VTKFile>\n");
        
        //-----------------
        // close the vtu file for element fields
        //-----------------
        fclose(out[0]);

    } // end write vtu


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn write_pvd
    ///
    /// \brief Writes a pvd ASCII output file for the element and nodal fields
    ///
    /// \param Vector of all graphics output times
    /// \param element average field names
    /// \param current time value
    /// \param graphics index
    ///
    /////////////////////////////////////////////////////////////////////////////
    void write_pvd(CArray<double>& graphics_times,
                   double time_value,
                   int graphics_id,
                   const size_t solver_id){

        FILE* out[20];   // the output files that are written to
        char  filename[100]; // char string
        int   max_len = sizeof filename;
        int   str_output_len;

        // Write time series metadata
        str_output_len = snprintf(filename, max_len, "vtk/Fierro.solver%zu.pvd", solver_id); 

        if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }
        // mesh file

        out[0] = fopen(filename, "w");
 
        fprintf(out[0], "<?xml version=\"1.0\"?>\n");
        fprintf(out[0], "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\">\n");
        fprintf(out[0], "  <Collection>\n");

        for (int i = 0; i <= graphics_id; i++) {
            fprintf(out[0], "    <DataSet timestep=\"%12.5e\" file=\"data/Fierro.solver%zu.%05d.vtm\" time= \"%12.5e\" />\n", 
                                                     graphics_times(i), solver_id, i, graphics_times(i) );
            //fprintf(out[0], "    <DataSet timestep=\"%d\" file=\"data/Fierro.solver%zu.%05d.vtm\" time= \"%12.5e\" />\n", 
            //                                         i, solver_id, i, graphics_times(i) );
        }

        fprintf(out[0], "  </Collection>\n");
        fprintf(out[0], "</VTKFile>"); 

        fclose(out[0]);

    } // end pvd


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn write_vtm
    ///
    /// \brief Writes a vtm ASCII output file for all fields -- mesh and material
    ///
    /// \param Vector of all graphics output times
    /// \param element average field names
    /// \param current time value
    /// \param graphics index
    ///
    /////////////////////////////////////////////////////////////////////////////
    void write_vtm(CArray<double>& graphics_times,
                   const  std::string& elem_part_name,
                   const  std::string& mat_part_name,
                   double time_value,
                   int graphics_id,
                   int num_mats,
                   int local_num_mats,
                   CArray<int> local_mats_in_rank,
                   int nranks,
                   int myrank,
                   bool write_mesh_state,
                   bool write_mat_pt_state,
                   const size_t solver_id)
    {   
        //gather number of mats on each rank
        //array storing number of local elems for this material on each process
        CArray<int> num_mats_in_rank, gatherv_displacements;
        CArray<size_t> interface_mats_in_rank;
        if(myrank==0){
            num_mats_in_rank = CArray<int>(nranks);
            interface_mats_in_rank = CArray<size_t>(nranks);
            gatherv_displacements = CArray<int>(nranks);
        }
        MPI_Gather(&local_num_mats,1,MPI_INT,num_mats_in_rank.pointer(),1,
                    MPI_INT, 0, MPI_COMM_WORLD);

        //gather which mats are present on each rank
        long long int length_mats_in_rank = 0;
        RaggedRightArray<int> mats_in_rank;
        if(myrank==0){
            for(int irank=0; irank < nranks; irank++){
                gatherv_displacements(irank) = length_mats_in_rank;
                interface_mats_in_rank(irank) = num_mats_in_rank(irank);
                length_mats_in_rank += num_mats_in_rank(irank);
            }
            mats_in_rank = RaggedRightArray<int>(interface_mats_in_rank);
        }
        // if(myrank==0){
        //     for(int irank=0; irank < nranks; irank++){
        //         std::cout << "NUM local mat elem on rank " << irank << " is " << num_mats_in_rank(irank) << std::endl;
        //         std::cout << "gatherv displacement on rank " << irank << " is " << gatherv_displacements(irank) << std::endl;
        //     }
        // }
        MPI_Gatherv(local_mats_in_rank.pointer(), local_num_mats, MPI_INT,
                    mats_in_rank.pointer(), num_mats_in_rank.pointer(),
                    gatherv_displacements.pointer(), MPI_INT, 0, MPI_COMM_WORLD);


        if(myrank==0){

            //invert map of rank to mat so its mat to rank
            CArray<size_t> num_ranks_in_mat(num_mats);
            num_ranks_in_mat.set_values(0);

            //count how many ranks each material is in
            for(int irank = 0; irank < nranks; irank++){
                for(int imat = 0; imat < num_mats_in_rank(irank); imat++){
                    num_ranks_in_mat(mats_in_rank(irank,imat))++;
                }
            }

            //allocate ragged storage and assign ranks to each mat
            RaggedRightArray<int> ranks_in_mat(num_ranks_in_mat);
            num_ranks_in_mat.set_values(0);

            for(int irank = 0; irank < nranks; irank++){
                for(int imat = 0; imat < num_mats_in_rank(irank); imat++){
                    ranks_in_mat(mats_in_rank(irank,imat),num_ranks_in_mat(mats_in_rank(irank,imat))) = irank;
                    num_ranks_in_mat(mats_in_rank(irank,imat))++;
                }
            }
            
            // loop over all the files that need to be written 
            for(int file_id=0; file_id<=graphics_id; file_id++){

                FILE* out[20];   // the output files that are written to
                char  filename[100]; // char string
                int   max_len = sizeof filename;
                int   str_output_len;


                // Write time series metadata to the data file
                str_output_len = snprintf(filename, max_len, "vtk/data/Fierro.solver%zu.%05d.vtm", solver_id, file_id); 

                if (str_output_len >= max_len) { fputs("Filename length exceeded; string truncated", stderr); }
                // mesh file

                out[0] = fopen(filename, "w");
        
                fprintf(out[0], "<?xml version=\"1.0\"?>\n");
                fprintf(out[0], "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
                fprintf(out[0], "  <vtkMultiBlockDataSet>\n");

                
                // Average mesh fields -- node and elem state written
                size_t block_id = 0;  // this will need to be incremented based on the number of mesh fields written
                if (write_mesh_state){
                    fprintf(out[0], "    <Block index=\"%zu\" name=\"Mesh\">\n", block_id);
                    {
                        block_id++;  // increment block id for material outputs that follow the element avg block

                        // elem and nodal fields are in this file
                        fprintf(out[0], "      <Piece index=\"0\" name=\"Field\">\n");
                        for(int irank = 0; irank < nranks; irank++){
                            fprintf(out[0], "        <DataSet index=\"%d\" file=\"Fierro.solver%zu.%s_rank%d.%05d.vtu\" />\n", 
                                                                    irank, solver_id, elem_part_name.c_str(), irank, file_id );
                        }
                        fprintf(out[0], "      </Piece>\n");

                        // add other Mesh average output Pieces here
                    }
                    fprintf(out[0], "    </Block>\n");
                } // end if write elem and node state is true

                // note: the block_id was incremented if an element average field output was made
                if (write_mat_pt_state){
                    fprintf(out[0], "    <Block index=\"%zu\" name=\"Mat\">\n", block_id);
                    for (size_t mat_id=0; mat_id<num_mats; mat_id++){
                        
                        // output the material specific fields
                        fprintf(out[0], "      <Piece index=\"%zu\" name=\"Mat%zu\">\n", mat_id, mat_id);
                        for(int irank = 0; irank < num_ranks_in_mat(mat_id); irank++){
                            fprintf(out[0], "        <DataSet index=\"%d\" file=\"Fierro.solver%zu.%s%zu_rank%d.%05d.vtu\" />\n", 
                                                                    irank, solver_id, mat_part_name.c_str(), mat_id, ranks_in_mat(mat_id, irank), file_id );
                        }
                        fprintf(out[0], "      </Piece>\n");

                    } // end for loop mat_id
                    fprintf(out[0], "    </Block>\n");
                } // end if write mat satte is true

                // done writing the files to be read by the vtm file
                fprintf(out[0], "  </vtkMultiBlockDataSet>\n");
                fprintf(out[0], "</VTKFile>"); 

                fclose(out[0]);

            } // end for file_id
        }

    } // end vtm


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_material_elem_node_lists
    ///
    /// \brief Creates elems and nodes for a unique mesh of a material (i.e, a part)
    ///
    /// \param Simulation mesh
    /// \param State node data
    /// \param Material node coordinates
    /// \param Material nodes in the material element
    /// \param Material to mesh map for elements
    /// \param number of material nodes
    /// \param number of material elements
    /// \param number of nodes in the element
    /// \param number of dimensions 
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_material_elem_node_lists(
        const Mesh_t& mesh,
        const DistributedDCArray<double>& state_node_coords,
        DCArrayKokkos<double>& mat_node_coords,
        DCArrayKokkos <size_t>& mat_nodes_in_mat_elem,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t mat_id,
        size_t& num_mat_nodes,
        const size_t num_mat_elems,
        const size_t num_nodes_in_elem,
        const size_t num_dims)
    {

        // helper arrays
        DCArrayKokkos <size_t> dummy_counter(mesh.num_nodes, "dummy_counter");
        DCArrayKokkos <size_t> access_mat_node_gids(mesh.num_nodes, "access_mat_node_gids");
        dummy_counter.set_values(0);

        // tag and count the number of nodes in this part
        FOR_ALL (mat_elem_sid, 0, num_mat_elems, {
            // get elem gid
            size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);  // WARNING not GPU compatible
            
            // parallel loop over the nodes in the element
            for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++) {
                size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

                Kokkos::atomic_add(&dummy_counter(node_gid), 1); // values in part will be >0

            } // end for nodes in element
            
        }); // end parallel for
        Kokkos::fence();
        dummy_counter.update_host();

        // loop opperation is not thread safe, must be run serially
        size_t mat_node_gid = 0;
        for(size_t node_gid = 0; node_gid<mesh.num_nodes; node_gid ++) {
            
            // save the nodes on the part (i.e., that belong to the material)
            if (dummy_counter.host(node_gid)>0){
                mat_node_coords.host(mat_node_gid, 0) = state_node_coords.host(node_gid, 0);
                mat_node_coords.host(mat_node_gid, 1) = state_node_coords.host(node_gid, 1);
                if (num_dims == 3){ 
                    mat_node_coords.host(mat_node_gid, 2) = state_node_coords.host(node_gid, 2);
                } // end if on dims

                access_mat_node_gids.host(node_gid) = mat_node_gid; // the part node id

                mat_node_gid ++;

                dummy_counter.host(node_gid) = 0; // set counter to zero, it was accounted for
            } // end if this node is on the part

        } // end loop over all mesh nodes
        mat_node_coords.update_device();
        access_mat_node_gids.update_device();
        dummy_counter.update_device();
        Kokkos::fence();

        // save the number of nodes defining the material region, i.e., the part
        num_mat_nodes = mat_node_gid;
        
        // save the new node id's
        FOR_ALL (mat_elem_sid, 0, num_mat_elems, {
            // get elem gid
            size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);
            
            // parallel loop over the nodes in the element
            for(size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++) {
                size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

                // save the mat_node to the mat elem list
                mat_nodes_in_mat_elem(mat_elem_sid, node_lid) = access_mat_node_gids(node_gid);

            } // end for nodes in element
            
        }); // end parallel for
        Kokkos::fence();
        mat_nodes_in_mat_elem.update_host();

    } // end build part (i.e., material elem and point lists) function

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_material_elem_node_lists
    ///
    /// \brief Creates elems and nodes for a unique mesh of a material (i.e, a part)
    ///
    /// \param Simulation mesh
    /// \param State node data
    /// \param Material node coordinates
    /// \param Material nodes in the material element
    /// \param Material to mesh map for elements
    /// \param number of material nodes
    /// \param number of material elements
    /// \param number of nodes in the element
    /// \param number of dimensions 
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_material_node_list(
        const Mesh_t& mesh,
        DCArrayKokkos<long long int>& mat_node_indices,
        DistributedDCArray<size_t>& mat_nodes_in_mat_elem,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t mat_id,
        size_t& num_mat_nodes,
        const size_t num_mat_elems,
        const size_t num_nodes_in_elem,
        const size_t num_dims)
    {

        
        std::set<long long int> mat_node_set;
        long long int node_gid;
        for (int elem_mat_id = 0; elem_mat_id < num_mat_elems; elem_mat_id++)
        {
            // set nodes per element
            for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++)
            {
                node_gid = mat_nodes_in_mat_elem.host(elem_mat_id, node_lid); //nodes in elem still stores global indices
                mat_node_set.insert(node_gid);
            }
        }

        // save the number of nodes defining the material region, i.e., the part
        num_mat_nodes = mat_node_set.size();

        //copy set to matar view
        int  inode = 0;
        int ighost = 0;
        auto it     = mat_node_set.begin();

        // create a Map for ghost node indices
        mat_node_indices = DCArrayKokkos<long long int>(num_mat_nodes, "mat_nodes");
        while (it != mat_node_set.end()) {
            mat_node_indices.host(ighost++) = *it;
            it++;
        }
        mat_node_indices.update_device();

    } // end build part (i.e., material elem and point lists) function



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
        SimulationParameters_t& SimulationParameters,
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
        size_t num_mats = State.MaterialPoints.num_material_points.size();

        State.MaterialPoints.den.update_host();
        State.MaterialPoints.pres.update_host();
        State.MaterialPoints.stress.update_host();
        State.MaterialPoints.sspd.update_host();
        State.MaterialPoints.sie.update_host();
        State.MaterialPoints.mass.update_host();

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

            size_t num_mat_elems = State.MaterialToMeshMaps.num_mat_elems.host(mat_id);

            for (size_t mat_elem_sid = 0; mat_elem_sid < num_mat_elems; mat_elem_sid++)
            {

                const size_t elem_gid = State.MaterialToMeshMaps.elem_in_mat_elem.host(mat_id, mat_elem_sid);

                double elem_coords[3];
                elem_coords[0] = 0.0;
                elem_coords[1] = 0.0;
                elem_coords[2] = 0.0;

                // get the coordinates of the element center
                for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {

                    elem_coords[0] += State.node.coords.host(mesh.nodes_in_elem.host(elem_gid, node_lid), 0);
                    elem_coords[1] += State.node.coords.host(mesh.nodes_in_elem.host(elem_gid, node_lid), 1);
                    if (num_dims == 3) {
                        elem_coords[2] += State.node.coords.host(mesh.nodes_in_elem.host(elem_gid, node_lid), 2);
                    }
                    else{
                        elem_coords[2] = 0.0;
                    }
                } // end loop over nodes in element

                elem_coords[0] = elem_coords[0] / ((double)mesh.num_nodes_in_elem);
                elem_coords[1] = elem_coords[1] / ((double)mesh.num_nodes_in_elem);
                elem_coords[2] = elem_coords[2] / ((double)mesh.num_nodes_in_elem);

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
                         State.MaterialPoints.den.host(mat_id, mat_elem_sid),
                         State.MaterialPoints.pres.host(mat_id, mat_elem_sid),
                         State.MaterialPoints.sie.host(mat_id, mat_elem_sid),
                         State.MaterialPoints.sspd.host(mat_id, mat_elem_sid),
                         State.GaussPoints.vol.host(elem_gid),
                         State.MaterialPoints.mass.host(mat_id, mat_elem_sid) );

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

            node_coords[0] = State.node.coords.host(node_gid, 0);
            node_coords[1] = State.node.coords.host(node_gid, 1);
            if (num_dims == 3) {
                node_coords[2] = State.node.coords.host(node_gid, 2);
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

           node_vel[0] = State.node.vel.host(node_gid, 0);
           node_vel[1] = State.node.vel.host(node_gid, 1);
            if (num_dims == 3) {
               node_vel[2] = State.node.vel.host(node_gid, 2);
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