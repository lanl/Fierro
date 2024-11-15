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

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // fmin, fmax, abs note: fminl is long
#include <sys/stat.h>
#include <set>
#include <mpi.h>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <Xpetra_MultiVector.hpp>
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "Tpetra_Details_FixedHashTable.hpp"

#include "elements.h"
#include "matar.h"
#include "utilities.h"
#include "node_combination.h"
#include "Solver.h"
#include "FEA_Module.h"
#include "MeshBuilder.h"
#include "MeshIO.h"

// Repartition Package
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_InputTraits.hpp>

// debug and performance includes
#include <sys/time.h>
#include <sys/resource.h>
#include <ctime>

#define BUFFER_LINES 20000
#define MAX_WORD 30
#define MAX_ELEM_NODES 32
#define STRAIN_EPSILON 0.000000001
#define DENSITY_EPSILON 0.0001
#define BC_EPSILON 1.0e-8

Solver::Solver(Simulation_Parameters& _simparam) : simparam(_simparam)
{
    // default flags assume optional routines are off
    setup_flag = finalize_flag = 0;
    communication_time = dev2host_time = host2dev_time = output_time = 0;
    last_print_step    = -1;

    // FEA module data init
    nfea_modules = 0;
    displacement_module = -1;
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn exit_solver
///
/// \brief End the solve
///
/////////////////////////////////////////////////////////////////////////////
void Solver::exit_solver(int status)
{
    Kokkos::finalize();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(status);
}

Solver::~Solver()
{
    // destroy FEA modules
    for (int imodule = 0; imodule < nfea_modules; imodule++)
    {
        delete fea_modules[imodule];
    }
}

namespace elements
{
namespace elem_types
{
/////////////////////////////////////////////////////////////////////////////
///
/// \fn from_vtk
///
/// \brief Check the VTK element type read in
///
/////////////////////////////////////////////////////////////////////////////
elem_type from_vtk(const int& vtk_elem)
{
    switch (vtk_elem)
    {
    case 9:
        return elem_type::Quad4;
    case 12:
        return elem_type::Hex8;
    case 70:
        return elem_type::QuadN;
    case 72:
        return elem_type::HexN;
    default:
        throw std::runtime_error("Unsupported vtk element type: " + std::to_string(vtk_elem));
    }
}
} // namespace elem_types
} // namespace elements

/////////////////////////////////////////////////////////////////////////////
///
/// \fn generate_mesh
///
/// \brief Generate simulation mesh
///
/// \param Options for generating the mesh
///
/////////////////////////////////////////////////////////////////////////////
void Solver::generate_mesh(const std::shared_ptr<MeshBuilderInput>& mesh_generation_options)
{
    auto mesh = MeshBuilder::build_mesh(mesh_generation_options);
    switch (active_node_ordering_convention)
    {
    case ENSIGHT:
        MeshIO::_Impl::reorder_columns(mesh.element_point_index, MeshIO::_Impl::ijk_to_fea().data());
        break;
    case IJK:
        // Already in IJK
        break;
    }

    num_nodes = mesh.points.dims(0);
    num_elem  = mesh.element_point_index.dims(0);
    map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(num_nodes, 0, comm));

    nlocal_nodes = map->getLocalNumElements();

    node_coords_distributed = Teuchos::rcp(new MV(map, mesh.points.dims(1)));
    {
        host_vec_array node_coords = node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
        for (long long int i = 0; i < num_nodes; i++)
        {
            // set global node id (ensight specific order)
            // let map decide if this node id belongs locally; if yes store data
            if (map->isNodeGlobalElement(i))
            {
                // set local node index in this mpi rank
                long long int node_rid = map->getLocalElement(i);
                for (long long int j = 0; j < mesh.points.dims(1); j++)
                {
                    node_coords(node_rid, j) = mesh.points(i, j);
                }
            }
        }
    }
    repartition_nodes();

    max_nodes_per_element = mesh.element_point_index.dims(1);

    // Figure out which elements belong to me.
    std::vector<int> global_indices_temp;
    for (size_t i = 0; i < mesh.element_point_index.dims(0); i++)
    {
        for (size_t j = 0; j < mesh.element_point_index.dims(1); j++)
        {
            if (map->isNodeGlobalElement(mesh.element_point_index(i, j)))
            {
                global_indices_temp.push_back(i);
                break;
            }
        }
    }

    rnum_elem = global_indices_temp.size();

    Element_Types = CArrayKokkos<elements::elem_types::elem_type, array_layout, HostSpace, memory_traits>(global_indices_temp.size());
    auto element_type = elements::elem_types::from_vtk(mesh.element_types(0)); // Not from elements. From up there ^
    switch (mesh.points.dims(1))
    {
    case 2:
        switch (element_type)
        {
        case elements::elem_types::elem_type::Quad4:
            max_nodes_per_patch = 2;
            break;
        default:
            throw std::runtime_error("Higher order meshes are unsupported.");
        }
        break;
    case 3:
        switch (element_type)
        {
        case elements::elem_types::elem_type::Hex8:
            max_nodes_per_patch = 4;
            break;
        default:
            throw std::runtime_error("Higher order meshes are unsupported.");
        }
        break;
    }

    for (size_t i = 0; i < global_indices_temp.size(); i++)
    {
        Element_Types(i) = element_type;
    }

    // copy temporary element storage to multivector storage
    dual_nodes_in_elem = dual_elem_conn_array("dual_nodes_in_elem", global_indices_temp.size(), mesh.element_point_index.dims(1));
    host_elem_conn_array nodes_in_elem = dual_nodes_in_elem.view_host();
    dual_nodes_in_elem.modify_host();

    for (size_t i = 0; i < global_indices_temp.size(); i++)
    {
        for (size_t j = 0; j < mesh.element_point_index.dims(1); j++)
        {
            nodes_in_elem(i, j) = mesh.element_point_index(global_indices_temp[i], j);
        }
    }

    // view storage for all local elements connected to local nodes on this rank
    Kokkos::DualView<GO*, array_layout, device_type, memory_traits> All_Element_Global_Indices("All_Element_Global_Indices", global_indices_temp.size());
    // copy temporary global indices storage to view storage
    for (int i = 0; i < global_indices_temp.size(); i++)
    {
        All_Element_Global_Indices.h_view(i) = global_indices_temp[i];
    }

    // construct overlapping element map (since different ranks can own the same elements due to the local node map)
    All_Element_Global_Indices.modify_host();
    All_Element_Global_Indices.sync_device();

    all_element_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(Teuchos::OrdinalTraits<GO>::invalid(), All_Element_Global_Indices.d_view, 0, comm));
}

/* ----------------------------------------------------------------------
   Read Ensight format mesh file
------------------------------------------------------------------------- */

void Solver::read_mesh_ensight(const char* MESH)
{
    Input_Options input_options = simparam.input_options.value();

    char ch;

    bool zero_index_base = input_options.zero_index_base;

    int negative_index_found = 0;
    int global_negative_index_found = 0;
    int num_dim = simparam.num_dims;
    int p_order = input_options.p_order;
    int local_node_index, current_column_index;
    int buffer_loop, buffer_iteration, buffer_iterations, dof_limit, scan_loop;

    size_t strain_count;
    size_t read_index_start, node_rid, elem_gid;

    std::string skip_line, read_line, substring;
    std::stringstream line_parse;

    CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;

    GO     node_gid;
    real_t dof_value;
    real_t unit_scaling = input_options.unit_scaling;

    // Nodes_Per_Element_Type =  elements::elem_types::Nodes_Per_Element_Type;

    // read the mesh
    // PLACEHOLDER: ensight_format(MESH);
    // abaqus_format(MESH);
    // vtk_format(MESH)

    // task 0 reads file
    if (myrank == 0)
    {
        std::cout << " NUM DIM is " << num_dim << std::endl;
        in = new std::ifstream();
        in->open(MESH);
        if (!(*in))
        {
            throw std::runtime_error(std::string("Can't open ") + MESH);
        }
        // skip 8 lines
        for (int j = 1; j <= 8; j++)
        {
            getline(*in, skip_line);
            std::cout << skip_line << std::endl;
        } // for
    }

    // --- Read the number of nodes in the mesh --- //
    if (myrank == 0)
    {
        getline(*in, read_line);
        line_parse.str(read_line);
        line_parse >> num_nodes;
        std::cout << "declared node count: " << num_nodes << std::endl;
    }

    // broadcast number of nodes
    MPI_Bcast(&num_nodes, 1, MPI_LONG_LONG_INT, 0, world);

    // construct contiguous parallel row map now that we know the number of nodes
    map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(num_nodes, 0, comm));
    // map->describe(*fos,Teuchos::VERB_EXTREME);
    // set the vertices in the mesh read in
    nlocal_nodes = map->getLocalNumElements();
    // populate local row offset data from global data
    global_size_t min_gid    = map->getMinGlobalIndex();
    global_size_t max_gid    = map->getMaxGlobalIndex();
    global_size_t index_base = map->getIndexBase();
    // debug print
    // std::cout << "local node count on task: " << " " << nlocal_nodes << std::endl;

    // allocate node storage with dual view
    // dual_node_coords = dual_vec_array("dual_node_coords", nlocal_nodes,num_dim);

    // local variable for host view in the dual view

    node_coords_distributed = Teuchos::rcp(new MV(map, num_dim));

    // scope ensures view is destroyed for now to avoid calling a device view with an active host view later
    {
        host_vec_array node_coords = node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
        // host_vec_array node_coords = dual_node_coords.view_host();
        // notify that the host view is going to be modified in the file readin
        // dual_node_coords.modify_host();

        // old swage method
        // mesh->init_nodes(local_nrows); // add 1 for index starting at 1

        //std::cout << "Num nodes assigned to task " << myrank << " = " << nlocal_nodes << std::endl;

        // read the initial mesh coordinates
        // x-coords
        /*only task 0 reads in nodes and elements from the input file
        stores node data in a buffer and communicates once the buffer cap is reached
        or the data ends*/

        words_per_line = input_options.words_per_line;
        elem_words_per_line = input_options.elem_words_per_line;

        // allocate read buffer
        read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES, words_per_line, MAX_WORD);

        dof_limit = num_nodes;
        buffer_iterations = dof_limit / BUFFER_LINES;
        if (dof_limit % BUFFER_LINES != 0)
        {
            buffer_iterations++;
        }

        // x-coords
        read_index_start = 0;
        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
        {
            // pack buffer on rank 0
            if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
            {
                for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++)
                {
                    getline(*in, read_line);
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
                while (buffer_iteration * BUFFER_LINES + buffer_loop < num_nodes) {
                    getline(*in, read_line);
                    line_parse.clear();
                    line_parse.str(read_line);
                    for (int iword = 0; iword < words_per_line; iword++)
                    {
                        // read portions of the line into the substring variable
                        line_parse >> substring;
                        // assign the substring variable as a word of the read buffer
                        strcpy(&read_buffer(buffer_loop, iword, 0), substring.c_str());
                    }
                    buffer_loop++;
                }
            }

            // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
            MPI_Bcast(read_buffer.pointer(), BUFFER_LINES * words_per_line * MAX_WORD, MPI_CHAR, 0, world);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, world);

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
                if (map->isNodeGlobalElement(node_gid))
                {
                    // set local node index in this mpi rank
                    node_rid = map->getLocalElement(node_gid);
                    // extract nodal position from the read buffer
                    // for ensight format this is just one coordinate per line
                    dof_value = atof(&read_buffer(scan_loop, 0, 0));
                    node_coords(node_rid, 0) = dof_value * unit_scaling;
                }
            }
            read_index_start += BUFFER_LINES;
        }

        // y-coords
        read_index_start = 0;
        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
        {
            // pack buffer on rank 0
            if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
            {
                for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++)
                {
                    getline(*in, read_line);
                    line_parse.clear();
                    line_parse.str(read_line);
                    for (int iword = 0; iword < words_per_line; iword++)
                    {
                        // read portions of the line into the substring variable
                        line_parse >> substring;
                        // assign the substring variable as a word of the read buffer
                        strcpy(&read_buffer(buffer_loop, iword, 0), substring.c_str());
                    }
                }
            }
            else if (myrank == 0)
            {
                buffer_loop = 0;
                while (buffer_iteration * BUFFER_LINES + buffer_loop < num_nodes) {
                    getline(*in, read_line);
                    line_parse.clear();
                    line_parse.str(read_line);
                    for (int iword = 0; iword < words_per_line; iword++)
                    {
                        // read portions of the line into the substring variable
                        line_parse >> substring;
                        // assign the substring variable as a word of the read buffer
                        strcpy(&read_buffer(buffer_loop, iword, 0), substring.c_str());
                    }
                    buffer_loop++;
                    // std::cout<<" "<< node_coords(node_gid, 0)<<std::endl;
                }
            }

            // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
            MPI_Bcast(read_buffer.pointer(), BUFFER_LINES * words_per_line * MAX_WORD, MPI_CHAR, 0, world);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, world);

            // determine which data to store in the swage mesh members (the local node data)
            // loop through read buffer
            for (scan_loop = 0; scan_loop < buffer_loop; scan_loop++)
            {
                // set global node id (ensight specific order)
                node_gid = read_index_start + scan_loop;
                // let map decide if this node id belongs locally; if yes store data
                if (map->isNodeGlobalElement(node_gid))
                {
                    // set local node index in this mpi rank
                    node_rid = map->getLocalElement(node_gid);
                    // extract nodal position from the read buffer
                    // for ensight format this is just one coordinate per line
                    dof_value = atof(&read_buffer(scan_loop, 0, 0));
                    node_coords(node_rid, 1) = dof_value * unit_scaling;
                }
            }
            read_index_start += BUFFER_LINES;
        }

        // z-coords
        read_index_start = 0;
        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
        {
            // pack buffer on rank 0
            if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
            {
                for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++)
                {
                    getline(*in, read_line);
                    line_parse.clear();
                    line_parse.str(read_line);
                    for (int iword = 0; iword < words_per_line; iword++)
                    {
                        // read portions of the line into the substring variable
                        line_parse >> substring;
                        // assign the substring variable as a word of the read buffer
                        strcpy(&read_buffer(buffer_loop, iword, 0), substring.c_str());
                    }
                }
            }
            else if (myrank == 0)
            {
                buffer_loop = 0;
                while (buffer_iteration * BUFFER_LINES + buffer_loop < num_nodes) {
                    getline(*in, read_line);
                    line_parse.clear();
                    line_parse.str(read_line);
                    for (int iword = 0; iword < words_per_line; iword++)
                    {
                        // read portions of the line into the substring variable
                        line_parse >> substring;
                        // assign the substring variable as a word of the read buffer
                        strcpy(&read_buffer(buffer_loop, iword, 0), substring.c_str());
                    }
                    buffer_loop++;
                    // std::cout<<" "<< node_coords(node_gid, 0)<<std::endl;
                }
            }

            // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
            MPI_Bcast(read_buffer.pointer(), BUFFER_LINES * words_per_line * MAX_WORD, MPI_CHAR, 0, world);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, world);

            // loop through read buffer and store coords in node coords view
            for (scan_loop = 0; scan_loop < buffer_loop; scan_loop++)
            {
                // set global node id (ensight specific order)
                node_gid = read_index_start + scan_loop;
                // let map decide if this node id belongs locally; if yes store data
                if (map->isNodeGlobalElement(node_gid))
                {
                    // set local node index in this mpi rank
                    node_rid = map->getLocalElement(node_gid);
                    // extract nodal position from the read buffer
                    // for ensight format this is just one coordinate per line
                    dof_value = atof(&read_buffer(scan_loop, 0, 0));
                    if (num_dim == 3)
                    {
                        node_coords(node_rid, 2) = dof_value * unit_scaling;
                    }
                }
            }
            read_index_start += BUFFER_LINES;
        }
    } // end active view scope
    // repartition node distribution
    repartition_nodes();

    // synchronize device data
    // dual_node_coords.sync_device();
    // dual_node_coords.modify_device();

    // debug print of nodal data

    // debug print nodal positions and indices
    /*
    std::cout << " ------------NODAL POSITIONS ON TASK " << myrank << " --------------"<<std::endl;
    for (int inode = 0; inode < local_nrows; inode++){
        std::cout << "node: " << map->getGlobalElement(inode) + 1 << " { ";
      for (int istride = 0; istride < num_dim; istride++){
          std::cout << node_coords(inode,istride) << " , ";
      }
      std::cout << " }"<< std::endl;
    }
    */

    // check that local assignments match global total

    // read in element info (ensight file format is organized in element type sections)
    // loop over this later for several element type sections

    num_elem  = 0;
    rnum_elem = 0;
    CArrayKokkos<int, array_layout, HostSpace, memory_traits> node_store(elem_words_per_line);

    if (myrank == 0)
    {
        // skip element type name line
        getline(*in, skip_line);
        std::cout << skip_line << std::endl;
    }

    // --- read the number of cells in the mesh ---
    // --- Read the number of vertices in the mesh --- //
    if (myrank == 0)
    {
        getline(*in, read_line);
        line_parse.clear();
        line_parse.str(read_line);
        line_parse >> num_elem;
        std::cout << "declared element count: " << num_elem << std::endl;
        if (num_elem <= 0)
        {
            std::cout << "ERROR, NO ELEMENTS IN MESH" << std::endl;
        }
    }

    // broadcast number of elements
    MPI_Bcast(&num_elem, 1, MPI_LONG_LONG_INT, 0, world);

    if (myrank == 0)
    {
        std::cout << "before mesh initialization" << std::endl;
    }

    // read in element connectivity
    // we're gonna reallocate for the words per line expected for the element connectivity
    read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES, elem_words_per_line, MAX_WORD);

    // calculate buffer iterations to read number of lines
    buffer_iterations = num_elem / BUFFER_LINES;
    int assign_flag;

    // dynamic buffer used to store elements before we know how many this rank needs
    std::vector<size_t> element_temp(BUFFER_LINES * elem_words_per_line);
    std::vector<size_t> global_indices_temp(BUFFER_LINES);
    size_t buffer_max = BUFFER_LINES * elem_words_per_line;
    size_t indices_buffer_max = BUFFER_LINES;

    if (num_elem % BUFFER_LINES != 0)
    {
        buffer_iterations++;
    }
    read_index_start = 0;
    // std::cout << "ELEMENT BUFFER ITERATIONS: " << buffer_iterations << std::endl;
    rnum_elem = 0;
    for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
    {
        // pack buffer on rank 0
        if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
        {
            for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++)
            {
                getline(*in, read_line);
                line_parse.clear();
                line_parse.str(read_line);
                for (int iword = 0; iword < elem_words_per_line; iword++)
                {
                    // read portions of the line into the substring variable
                    line_parse >> substring;
                    // assign the substring variable as a word of the read buffer
                    strcpy(&read_buffer(buffer_loop, iword, 0), substring.c_str());
                }
            }
        }
        else if (myrank == 0)
        {
            buffer_loop = 0;
            while (buffer_iteration * BUFFER_LINES + buffer_loop < num_elem) {
                getline(*in, read_line);
                line_parse.clear();
                line_parse.str(read_line);
                for (int iword = 0; iword < elem_words_per_line; iword++)
                {
                    // read portions of the line into the substring variable
                    line_parse >> substring;
                    // assign the substring variable as a word of the read buffer
                    strcpy(&read_buffer(buffer_loop, iword, 0), substring.c_str());
                }
                buffer_loop++;
                // std::cout<<" "<< node_coords(node_gid, 0)<<std::endl;
            }
        }

        // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
        MPI_Bcast(read_buffer.pointer(), BUFFER_LINES * elem_words_per_line * MAX_WORD, MPI_CHAR, 0, world);
        // broadcast how many nodes were read into this buffer iteration
        MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, world);

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
                    if (map->isNodeGlobalElement(node_gid) && !assign_flag)
                    {
                        assign_flag = 1;
                        rnum_elem++;
                    }
                }
                else
                {
                    if (map->isNodeGlobalElement(node_gid - 1) && !assign_flag)
                    {
                        assign_flag = 1;
                        rnum_elem++;
                    }
                }
            }

            if (assign_flag)
            {
                for (int inode = 0; inode < elem_words_per_line; inode++)
                {
                    if ((rnum_elem - 1) * elem_words_per_line + inode >= buffer_max)
                    {
                        element_temp.resize((rnum_elem - 1) * elem_words_per_line + inode + BUFFER_LINES * elem_words_per_line);
                        buffer_max = (rnum_elem - 1) * elem_words_per_line + inode + BUFFER_LINES * elem_words_per_line;
                    }
                    element_temp[(rnum_elem - 1) * elem_words_per_line + inode] = node_store(inode);
                    // std::cout << "VECTOR STORAGE FOR ELEM " << rnum_elem << " ON TASK " << myrank << " NODE " << inode+1 << " IS " << node_store(inode) + 1 << std::endl;
                }
                // assign global element id to temporary list
                if (rnum_elem - 1 >= indices_buffer_max)
                {
                    global_indices_temp.resize(rnum_elem - 1 + BUFFER_LINES);
                    indices_buffer_max = rnum_elem - 1 + BUFFER_LINES;
                }
                global_indices_temp[rnum_elem - 1] = elem_gid;
            }
        }
        read_index_start += BUFFER_LINES;
    }

    // Close mesh input file
    if (myrank == 0)
    {
        in->close();
    }

    // std::cout << "RNUM ELEMENTS IS: " << rnum_elem << std::endl;

    Element_Types = CArrayKokkos<elements::elem_types::elem_type, array_layout, HostSpace, memory_traits>(rnum_elem);

    elements::elem_types::elem_type mesh_element_type;

    if (simparam.num_dims == 2)
    {
        if (input_options.element_type == ELEMENT_TYPE::quad4)
        {
            mesh_element_type   = elements::elem_types::Quad4;
            max_nodes_per_patch = 2;
        }
        else if (input_options.element_type == ELEMENT_TYPE::quad8)
        {
            mesh_element_type   = elements::elem_types::Quad8;
            max_nodes_per_patch = 3;
        }
        else if (input_options.element_type == ELEMENT_TYPE::quad12)
        {
            mesh_element_type   = elements::elem_types::Quad12;
            max_nodes_per_patch = 4;
        }
        else
        {
            if (myrank == 0)
            {
                std::cout << "ELEMENT TYPE UNRECOGNIZED" << std::endl;
            }
            exit_solver(0);
        }
        element_select->choose_2Delem_type(mesh_element_type, elem2D);
        max_nodes_per_element = elem2D->num_nodes();
    }

    if (simparam.num_dims == 3)
    {
        if (input_options.element_type == ELEMENT_TYPE::hex8)
        {
            mesh_element_type   = elements::elem_types::Hex8;
            max_nodes_per_patch = 4;
        }
        else if (input_options.element_type == ELEMENT_TYPE::hex20)
        {
            mesh_element_type   = elements::elem_types::Hex20;
            max_nodes_per_patch = 8;
        }
        else if (input_options.element_type == ELEMENT_TYPE::hex32)
        {
            mesh_element_type   = elements::elem_types::Hex32;
            max_nodes_per_patch = 12;
        }
        else
        {
            if (myrank == 0)
            {
                std::cout << "ELEMENT TYPE UNRECOGNIZED" << std::endl;
            }
            exit_solver(0);
        }
        element_select->choose_3Delem_type(mesh_element_type, elem);
        max_nodes_per_element = elem->num_nodes();
    }

    // 1 type per mesh for now
    for (int ielem = 0; ielem < rnum_elem; ielem++)
    {
        Element_Types(ielem) = mesh_element_type;
    }

    // copy temporary element storage to multivector storage
    dual_nodes_in_elem = dual_elem_conn_array("dual_nodes_in_elem", rnum_elem, max_nodes_per_element);
    host_elem_conn_array nodes_in_elem = dual_nodes_in_elem.view_host();
    dual_nodes_in_elem.modify_host();

    for (int ielem = 0; ielem < rnum_elem; ielem++)
    {
        for (int inode = 0; inode < elem_words_per_line; inode++)
        {
            nodes_in_elem(ielem, inode) = element_temp[ielem * elem_words_per_line + inode];
        }
    }

    // view storage for all local elements connected to local nodes on this rank
    // DCArrayKokkos<GO, array_layout, device_type, memory_traits> All_Element_Global_Indices(rnum_elem);
    Kokkos::DualView<GO*, array_layout, device_type, memory_traits> All_Element_Global_Indices("All_Element_Global_Indices", rnum_elem);
    // copy temporary global indices storage to view storage
    for (int ielem = 0; ielem < rnum_elem; ielem++)
    {
        All_Element_Global_Indices.h_view(ielem) = global_indices_temp[ielem];
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
        exit_solver(0);
    }
    // delete temporary element connectivity and index storage
    std::vector<size_t>().swap(element_temp);
    std::vector<size_t>().swap(global_indices_temp);

    All_Element_Global_Indices.modify_host();
    All_Element_Global_Indices.sync_device();

    // debug print
    /*
    Kokkos::View <GO*, array_layout, device_type, memory_traits> All_Element_Global_Indices_pass("All_Element_Global_Indices_pass",rnum_elem);
    deep_copy(All_Element_Global_Indices_pass, All_Element_Global_Indices.h_view);
    std::cout << " ------------ELEMENT GLOBAL INDICES ON TASK " << myrank << " --------------"<<std::endl;
    for (int ielem = 0; ielem < rnum_elem; ielem++){
      std::cout << "elem: " << All_Element_Global_Indices_pass(ielem) + 1;
      std::cout << std::endl;
    }
    */

    // construct overlapping element map (since different ranks can own the same elements due to the local node map)
    all_element_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(Teuchos::OrdinalTraits<GO>::invalid(), All_Element_Global_Indices.d_view, 0, comm));

    // element type selection (subject to change)
    // ---- Set Element Type ---- //
    // allocate element type memory
    // elements::elem_type_t* elem_choice;

    int NE = 1; // number of element types in problem

    // Convert ensight index system to the ijk finite element numbering convention
    // for vertices in cell
    if (active_node_ordering_convention == IJK)
    {
        CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_ensight_to_ijk(max_nodes_per_element);
        CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> tmp_ijk_indx(max_nodes_per_element);
        convert_ensight_to_ijk(0) = 0;
        convert_ensight_to_ijk(1) = 1;
        convert_ensight_to_ijk(2) = 3;
        convert_ensight_to_ijk(3) = 2;
        convert_ensight_to_ijk(4) = 4;
        convert_ensight_to_ijk(5) = 5;
        convert_ensight_to_ijk(6) = 7;
        convert_ensight_to_ijk(7) = 6;

        int nodes_per_element;

        if (num_dim == 2)
        {
            for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++)
            {
                // set nodes per element
                element_select->choose_2Delem_type(Element_Types(cell_rid), elem2D);
                nodes_per_element = elem2D->num_nodes();
                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    tmp_ijk_indx(node_lid) = nodes_in_elem(cell_rid, convert_ensight_to_ijk(node_lid));
                }

                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    nodes_in_elem(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
                }
            }
        }

        if (num_dim == 3)
        {
            for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++)
            {
                // set nodes per element
                element_select->choose_3Delem_type(Element_Types(cell_rid), elem);
                nodes_per_element = elem->num_nodes();
                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    tmp_ijk_indx(node_lid) = nodes_in_elem(cell_rid, convert_ensight_to_ijk(node_lid));
                }

                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    nodes_in_elem(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
                }
            }
        }
    }
    // debug print element edof
    /*
    std::cout << " ------------ELEMENT EDOF ON TASK " << myrank << " --------------"<<std::endl;

    for (int ielem = 0; ielem < rnum_elem; ielem++){
      std::cout << "elem:  " << ielem+1 << std::endl;
      for (int lnode = 0; lnode < 8; lnode++){
          std::cout << "{ ";
            std::cout << lnode+1 << " = " << nodes_in_elem(ielem,lnode) + 1 << " ";

          std::cout << " }"<< std::endl;
      }
      std::cout << std::endl;
    }
    */
} // end read_mesh

/* ----------------------------------------------------------------------
   Read VTK format mesh file
------------------------------------------------------------------------- */

void Solver::read_mesh_vtk(const char* MESH)
{
    Input_Options input_options = simparam.input_options.value();

    char ch;
    std::string skip_line, read_line, substring;
    std::stringstream line_parse;

    int num_dim = simparam.num_dims;
    int p_order = input_options.p_order;
    int local_node_index, current_column_index;
    int buffer_loop, buffer_iteration, buffer_iterations, dof_limit, scan_loop;
    int negative_index_found = 0;
    int global_negative_index_found = 0;

    size_t read_index_start, node_rid, elem_gid;
    size_t strain_count;

    GO     node_gid;
    real_t dof_value;
    real_t unit_scaling                  = input_options.unit_scaling;
    bool   zero_index_base               = input_options.zero_index_base;
    bool   topology_optimization_restart = input_options.topology_optimization_restart;
    bool   shape_optimization_restart    = input_options.shape_optimization_restart;

    CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;

    // Nodes_Per_Element_Type =  elements::elem_types::Nodes_Per_Element_Type;
    simparam.restart_file = topology_optimization_restart;

    // read the mesh
    // --- Read the number of nodes in the mesh --- //
    num_nodes = 0;
    if (myrank == 0)
    {
        std::cout << " NUM DIM is " << num_dim << std::endl;
        in = new std::ifstream();
        in->open(MESH);
        bool found = false;

        while (found == false&&in->good()) {
            std::getline(*in, read_line);
            line_parse.str("");
            line_parse.clear();
            line_parse << read_line;
            line_parse >> substring;

            // looking for the following text:
            //      POINTS %d float
            if (substring == "POINTS")
            {
                line_parse >> num_nodes;
                std::cout << "declared node count: " << num_nodes << std::endl;
                if (num_nodes <= 0)
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
    MPI_Bcast(&num_nodes, 1, MPI_LONG_LONG_INT, 0, world);

    // construct contiguous parallel row map now that we know the number of nodes
    map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(num_nodes, 0, comm));
    // map->describe(*fos,Teuchos::VERB_EXTREME);

    // set the vertices in the mesh read in
    nlocal_nodes = map->getLocalNumElements();
    // populate local row offset data from global data
    global_size_t min_gid    = map->getMinGlobalIndex();
    global_size_t max_gid    = map->getMaxGlobalIndex();
    global_size_t index_base = map->getIndexBase();
    // debug print
    // std::cout << "local node count on task: " << " " << nlocal_nodes << std::endl;

    // allocate node storage with dual view
    // dual_node_coords = dual_vec_array("dual_node_coords", nlocal_nodes,num_dim);

    // local variable for host view in the dual view

    node_coords_distributed = Teuchos::rcp(new MV(map, num_dim));

    // scope ensures view is destroyed for now to avoid calling a device view with an active host view later
    {
        host_vec_array node_coords = node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
        // host_vec_array node_coords = dual_node_coords.view_host();
        // notify that the host view is going to be modified in the file readin
        // dual_node_coords.modify_host();

        // old swage method
        // mesh->init_nodes(local_nrows); // add 1 for index starting at 1

        //std::cout << "Num nodes assigned to task " << myrank << " = " << nlocal_nodes << std::endl;

        // read the initial mesh coordinates
        // x-coords
        /*only task 0 reads in nodes and elements from the input file
        stores node data in a buffer and communicates once the buffer cap is reached
        or the data ends*/

        words_per_line = input_options.words_per_line;
        elem_words_per_line = input_options.elem_words_per_line;

        // allocate read buffer
        read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES, words_per_line, MAX_WORD);

        dof_limit = num_nodes;
        buffer_iterations = dof_limit / BUFFER_LINES;
        if (dof_limit % BUFFER_LINES != 0)
        {
            buffer_iterations++;
        }

        // read coords
        read_index_start = 0;
        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
        {
            // pack buffer on rank 0
            if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
            {
                for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++)
                {
                    getline(*in, read_line);
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
                while (buffer_iteration * BUFFER_LINES + buffer_loop < num_nodes) {
                    getline(*in, read_line);
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
            MPI_Bcast(read_buffer.pointer(), BUFFER_LINES * words_per_line * MAX_WORD, MPI_CHAR, 0, world);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, world);

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
                if (map->isNodeGlobalElement(node_gid))
                {
                    // set local node index in this mpi rank
                    node_rid = map->getLocalElement(node_gid);
                    // extract nodal position from the read buffer
                    // for tecplot format this is the three coords in the same line
                    dof_value = atof(&read_buffer(scan_loop, 0, 0));
                    node_coords(node_rid, 0) = dof_value * unit_scaling;
                    dof_value = atof(&read_buffer(scan_loop, 1, 0));
                    node_coords(node_rid, 1) = dof_value * unit_scaling;
                    if (num_dim == 3)
                    {
                        dof_value = atof(&read_buffer(scan_loop, 2, 0));
                        node_coords(node_rid, 2) = dof_value * unit_scaling;
                    }
                }
            }
            read_index_start += BUFFER_LINES;
        }
    } // end of coordinate readin
    // repartition node distribution
    repartition_nodes(false);

    // synchronize device data
    // dual_node_coords.sync_device();
    // dual_node_coords.modify_device();

    // debug print of nodal data

    // debug print nodal positions and indices
    /*
    std::cout << " ------------NODAL POSITIONS ON TASK " << myrank << " --------------"<<std::endl;
    for (int inode = 0; inode < local_nrows; inode++){
        std::cout << "node: " << map->getGlobalElement(inode) + 1 << " { ";
      for (int istride = 0; istride < num_dim; istride++){
          std::cout << node_coords(inode,istride) << " , ";
      }
      std::cout << " }"<< std::endl;
    }
    */

    // check that local assignments match global total

    // read in element info (ensight file format is organized in element type sections)
    // loop over this later for several element type sections

    num_elem  = 0;
    rnum_elem = 0;
    CArrayKokkos<int, array_layout, HostSpace, memory_traits> node_store(elem_words_per_line);

    // --- read the number of cells in the mesh ---
    // --- Read the number of vertices in the mesh --- //
    if (myrank == 0)
    {
        bool found = false;
        while (found == false&&in->good()) {
            std::getline(*in, read_line);
            line_parse.str("");
            line_parse.clear();
            line_parse << read_line;
            line_parse >> substring;

            // looking for the following text:
            //      CELLS num_cells size
            if (substring == "CELLS")
            {
                line_parse >> num_elem;
                std::cout << "declared element count: " << num_elem << std::endl;
                if (num_elem <= 0)
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
    MPI_Bcast(&num_elem, 1, MPI_LONG_LONG_INT, 0, world);

    if (myrank == 0)
    {
        std::cout << "before mesh initialization" << std::endl;
    }

    // read in element connectivity
    // we're gonna reallocate for the words per line expected for the element connectivity
    read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES, elem_words_per_line, MAX_WORD);

    // calculate buffer iterations to read number of lines
    buffer_iterations = num_elem / BUFFER_LINES;
    int assign_flag;

    // dynamic buffer used to store elements before we know how many this rank needs
    std::vector<size_t> element_temp(BUFFER_LINES * elem_words_per_line);
    std::vector<size_t> global_indices_temp(BUFFER_LINES);
    size_t buffer_max = BUFFER_LINES * elem_words_per_line;
    size_t indices_buffer_max = BUFFER_LINES;

    if (num_elem % BUFFER_LINES != 0)
    {
        buffer_iterations++;
    }
    read_index_start = 0;
    // std::cout << "ELEMENT BUFFER ITERATIONS: " << buffer_iterations << std::endl;
    rnum_elem = 0;
    for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
    {
        // pack buffer on rank 0
        if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
        {
            for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++)
            {
                getline(*in, read_line);
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
            while (buffer_iteration * BUFFER_LINES + buffer_loop < num_elem) {
                getline(*in, read_line);
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
        MPI_Bcast(read_buffer.pointer(), BUFFER_LINES * elem_words_per_line * MAX_WORD, MPI_CHAR, 0, world);
        // broadcast how many nodes were read into this buffer iteration
        MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, world);

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
                    if (map->isNodeGlobalElement(node_gid) && !assign_flag)
                    {
                        assign_flag = 1;
                        rnum_elem++;
                    }
                }
                else
                {
                    if (map->isNodeGlobalElement(node_gid - 1) && !assign_flag)
                    {
                        assign_flag = 1;
                        rnum_elem++;
                    }
                }
            }

            if (assign_flag)
            {
                for (int inode = 0; inode < elem_words_per_line; inode++)
                {
                    if ((rnum_elem - 1) * elem_words_per_line + inode >= buffer_max)
                    {
                        element_temp.resize((rnum_elem - 1) * elem_words_per_line + inode + BUFFER_LINES * elem_words_per_line);
                        buffer_max = (rnum_elem - 1) * elem_words_per_line + inode + BUFFER_LINES * elem_words_per_line;
                    }
                    element_temp[(rnum_elem - 1) * elem_words_per_line + inode] = node_store(inode);
                    // std::cout << "VECTOR STORAGE FOR ELEM " << rnum_elem << " ON TASK " << myrank << " NODE " << inode+1 << " IS " << node_store(inode) + 1 << std::endl;
                }
                // assign global element id to temporary list
                if (rnum_elem - 1 >= indices_buffer_max)
                {
                    global_indices_temp.resize(rnum_elem - 1 + BUFFER_LINES);
                    indices_buffer_max = rnum_elem - 1 + BUFFER_LINES;
                }
                global_indices_temp[rnum_elem - 1] = elem_gid;
            }
        }
        read_index_start += BUFFER_LINES;
    }

    // std::cout << "RNUM ELEMENTS IS: " << rnum_elem << std::endl;

    Element_Types = CArrayKokkos<elements::elem_types::elem_type, array_layout, HostSpace, memory_traits>(rnum_elem);

    elements::elem_types::elem_type mesh_element_type;

    if (simparam.num_dims == 2)
    {
        if (input_options.element_type == ELEMENT_TYPE::quad4)
        {
            mesh_element_type   = elements::elem_types::Quad4;
            max_nodes_per_patch = 2;
        }
        else if (input_options.element_type == ELEMENT_TYPE::quad8)
        {
            mesh_element_type   = elements::elem_types::Quad8;
            max_nodes_per_patch = 3;
        }
        else if (input_options.element_type == ELEMENT_TYPE::quad12)
        {
            mesh_element_type   = elements::elem_types::Quad12;
            max_nodes_per_patch = 4;
        }
        else
        {
            if (myrank == 0)
            {
                std::cout << "ELEMENT TYPE UNRECOGNIZED" << std::endl;
            }
            exit_solver(0);
        }
        element_select->choose_2Delem_type(mesh_element_type, elem2D);
        max_nodes_per_element = elem2D->num_nodes();
    }

    if (simparam.num_dims == 3)
    {
        if (input_options.element_type == ELEMENT_TYPE::hex8)
        {
            mesh_element_type   = elements::elem_types::Hex8;
            max_nodes_per_patch = 4;
        }
        else if (input_options.element_type == ELEMENT_TYPE::hex20)
        {
            mesh_element_type   = elements::elem_types::Hex20;
            max_nodes_per_patch = 8;
        }
        else if (input_options.element_type == ELEMENT_TYPE::hex32)
        {
            mesh_element_type   = elements::elem_types::Hex32;
            max_nodes_per_patch = 12;
        }
        else
        {
            if (myrank == 0)
            {
                std::cout << "ELEMENT TYPE UNRECOGNIZED" << std::endl;
            }
            exit_solver(0);
        }
        element_select->choose_3Delem_type(mesh_element_type, elem);
        max_nodes_per_element = elem->num_nodes();
    }

    // 1 type per mesh for now
    for (int ielem = 0; ielem < rnum_elem; ielem++)
    {
        Element_Types(ielem) = mesh_element_type;
    }

    // copy temporary element storage to multivector storage
    dual_nodes_in_elem = dual_elem_conn_array("dual_nodes_in_elem", rnum_elem, max_nodes_per_element);
    host_elem_conn_array nodes_in_elem = dual_nodes_in_elem.view_host();
    dual_nodes_in_elem.modify_host();

    for (int ielem = 0; ielem < rnum_elem; ielem++)
    {
        for (int inode = 0; inode < elem_words_per_line; inode++)
        {
            nodes_in_elem(ielem, inode) = element_temp[ielem * elem_words_per_line + inode];
        }
    }

    // view storage for all local elements connected to local nodes on this rank
    // DCArrayKokkos<GO, array_layout, device_type, memory_traits> All_Element_Global_Indices(rnum_elem);
    Kokkos::DualView<GO*, array_layout, device_type, memory_traits> All_Element_Global_Indices("All_Element_Global_Indices", rnum_elem);
    // copy temporary global indices storage to view storage
    for (int ielem = 0; ielem < rnum_elem; ielem++)
    {
        All_Element_Global_Indices.h_view(ielem) = global_indices_temp[ielem];
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
        exit_solver(0);
    }

    // delete temporary element connectivity and index storage
    std::vector<size_t>().swap(element_temp);
    std::vector<size_t>().swap(global_indices_temp);

    All_Element_Global_Indices.modify_host();
    All_Element_Global_Indices.sync_device();

    // debug print
    /*
    Kokkos::View <GO*, array_layout, device_type, memory_traits> All_Element_Global_Indices_pass("All_Element_Global_Indices_pass",rnum_elem);
    deep_copy(All_Element_Global_Indices_pass, All_Element_Global_Indices.h_view);
    std::cout << " ------------ELEMENT GLOBAL INDICES ON TASK " << myrank << " --------------"<<std::endl;
    for (int ielem = 0; ielem < rnum_elem; ielem++){
      std::cout << "elem: " << All_Element_Global_Indices_pass(ielem) + 1;
      std::cout << std::endl;
    }
    */

    // construct overlapping element map (since different ranks can own the same elements due to the local node map)
    all_element_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(Teuchos::OrdinalTraits<GO>::invalid(), All_Element_Global_Indices.d_view, 0, comm));

    // element type selection (subject to change)
    // ---- Set Element Type ---- //
    // allocate element type memory
    // elements::elem_type_t* elem_choice;

    int NE = 1; // number of element types in problem

    // Convert ensight index system to the ijk finite element numbering convention
    // for vertices in cell
    if (active_node_ordering_convention == IJK)
    {
        CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_ensight_to_ijk(max_nodes_per_element);
        CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> tmp_ijk_indx(max_nodes_per_element);
        convert_ensight_to_ijk(0) = 0;
        convert_ensight_to_ijk(1) = 1;
        convert_ensight_to_ijk(2) = 3;
        convert_ensight_to_ijk(3) = 2;
        convert_ensight_to_ijk(4) = 4;
        convert_ensight_to_ijk(5) = 5;
        convert_ensight_to_ijk(6) = 7;
        convert_ensight_to_ijk(7) = 6;

        int nodes_per_element;

        if (num_dim == 2)
        {
            for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++)
            {
                // set nodes per element
                element_select->choose_2Delem_type(Element_Types(cell_rid), elem2D);
                nodes_per_element = elem2D->num_nodes();
                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    tmp_ijk_indx(node_lid) = nodes_in_elem(cell_rid, convert_ensight_to_ijk(node_lid));
                }

                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    nodes_in_elem(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
                }
            }
        }

        if (num_dim == 3)
        {
            for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++)
            {
                // set nodes per element
                element_select->choose_3Delem_type(Element_Types(cell_rid), elem);
                nodes_per_element = elem->num_nodes();
                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    tmp_ijk_indx(node_lid) = nodes_in_elem(cell_rid, convert_ensight_to_ijk(node_lid));
                }

                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    nodes_in_elem(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
                }
            }
        }
    }

    //If restarting a topology optimization run, obtain nodal design density data here
    if(topology_optimization_restart){
        design_node_densities_distributed = Teuchos::rcp(new MV(map, 1));
        host_vec_array node_densities = design_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
        if (myrank == 0)
        {
            bool found = false;
            while (found == false&&in->good()) {
                std::getline(*in, read_line);
                //std::cout << read_line << std::endl;
                line_parse.clear();
                line_parse.str(read_line);

                //stop when the design_density string is reached
                while (!line_parse.eof()){
                    line_parse >> substring;
                    //std::cout << substring << std::endl;
                    if(!substring.compare("design_density")){
                        found = true;
                    }
                } //while

            } // end while

            if (!found){
                throw std::runtime_error("ERROR: Failed to find design_density");
            } // end if

            //skip "LOOKUP_TABLE default" line
            std::getline(*in, read_line);
        } // end if(myrank==0)
        
        //read in density of each node
        // allocate read buffer
        words_per_line = 1;
        read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES, words_per_line, MAX_WORD);

        dof_limit = num_nodes;
        buffer_iterations = dof_limit / BUFFER_LINES;
        if (dof_limit % BUFFER_LINES != 0)
        {
            buffer_iterations++;
        }

        // read densities
        read_index_start = 0;
        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
        {
            // pack buffer on rank 0
            if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
            {
                for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++)
                {
                    getline(*in, read_line);
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
                while (buffer_iteration * BUFFER_LINES + buffer_loop < num_nodes) {
                    getline(*in, read_line);
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
            MPI_Bcast(read_buffer.pointer(), BUFFER_LINES * words_per_line * MAX_WORD, MPI_CHAR, 0, world);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, world);

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
                if (map->isNodeGlobalElement(node_gid))
                {
                    // set local node index in this mpi rank
                    node_rid = map->getLocalElement(node_gid);
                    // extract nodal position from the read buffer
                    // for tecplot format this is the three coords in the same line
                    dof_value = atof(&read_buffer(scan_loop, 0, 0));
                    node_densities(node_rid, 0) = dof_value;
                }
            }
            read_index_start += BUFFER_LINES;
        }

        //Find initial objective value to normalize by
        if (myrank == 0 && simparam.optimization_options.normalized_objective)
        {
            bool found = false;
            while (found == false&&in->good()) {
                std::getline(*in, read_line);
                //std::cout << read_line << std::endl;
                line_parse.clear();
                line_parse.str(read_line);

                //stop when the design_density string is reached
                while (!line_parse.eof()){
                    line_parse >> substring;
                    //std::cout << substring << std::endl;
                    if(!substring.compare("Objective_Normalization_Constant")){
                        found = true;
                        line_parse >> substring;
                        simparam.optimization_options.objective_normalization_constant = stod(substring);
                    }
                } //while

            } // end while

            if (!found){
                throw std::runtime_error("ERROR: Failed to find initial objective value for restart");
            } // end if
        } // end if(myrank==0)
    }

    //If restarting a topology optimization run, obtain nodal design density data here
    if(shape_optimization_restart){
        design_node_coords_distributed = Teuchos::rcp(new MV(map, num_dim));
        host_vec_array design_node_coords = design_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
        if (myrank == 0)
        {
            bool found = false;
            while (found == false&&in->good()) {
                std::getline(*in, read_line);
                //std::cout << read_line << std::endl;
                line_parse.clear();
                line_parse.str(read_line);

                //stop when the design_density string is reached
                while (!line_parse.eof()){
                    line_parse >> substring;
                    //std::cout << substring << std::endl;
                    if(!substring.compare("design_coordinates")){
                        found = true;
                    }
                } //while

            } // end while

            if (!found){
                throw std::runtime_error("ERROR: Failed to find design_coordinates");
            } // end if

            //skip "LOOKUP_TABLE default" line
            std::getline(*in, read_line);
        } // end if(myrank==0)
        
        //read in density of each node
        // allocate read buffer
        words_per_line = num_dim;
        read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES, words_per_line, MAX_WORD);

        dof_limit = num_nodes;
        buffer_iterations = dof_limit / BUFFER_LINES;
        if (dof_limit % BUFFER_LINES != 0)
        {
            buffer_iterations++;
        }

        // read densities
        read_index_start = 0;
        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
        {
            // pack buffer on rank 0
            if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
            {
                for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++)
                {
                    getline(*in, read_line);
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
                while (buffer_iteration * BUFFER_LINES + buffer_loop < num_nodes) {
                    getline(*in, read_line);
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
            MPI_Bcast(read_buffer.pointer(), BUFFER_LINES * words_per_line * MAX_WORD, MPI_CHAR, 0, world);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, world);

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
                if (map->isNodeGlobalElement(node_gid))
                {
                    // set local node index in this mpi rank
                    node_rid = map->getLocalElement(node_gid);
                    // extract nodal position from the read buffer
                    // for tecplot format this is the three coords in the same line
                    dof_value = atof(&read_buffer(scan_loop, 0, 0));
                    design_node_coords(node_rid, 0) = dof_value * unit_scaling;
                    dof_value = atof(&read_buffer(scan_loop, 1, 0));
                    design_node_coords(node_rid, 1) = dof_value * unit_scaling;
                    if(num_dim==3){
                        dof_value = atof(&read_buffer(scan_loop, 2, 0));
                        design_node_coords(node_rid, 2) = dof_value * unit_scaling; 
                    }
                }
            }
            read_index_start += BUFFER_LINES;
        }

        //Find initial objective value to normalize by
        if (myrank == 0 && simparam.optimization_options.normalized_objective)
        {
            bool found = false;
            while (found == false&&in->good()) {
                std::getline(*in, read_line);
                //std::cout << read_line << std::endl;
                line_parse.clear();
                line_parse.str(read_line);

                //stop when the design_density string is reached
                while (!line_parse.eof()){
                    line_parse >> substring;
                    //std::cout << substring << std::endl;
                    if(!substring.compare("Objective_Normalization_Constant")){
                        found = true;
                        line_parse >> substring;
                        simparam.optimization_options.objective_normalization_constant = stod(substring);
                    }
                } //while

            } // end while

            if (!found){
                throw std::runtime_error("ERROR: Failed to find initial objective value for restart");
            } // end if
        } // end if(myrank==0)
    }

    // Close mesh input file
    if (myrank == 0)
    {
        in->close();
    }
    // debug print element edof
    /*
    std::cout << " ------------ELEMENT EDOF ON TASK " << myrank << " --------------"<<std::endl;

    for (int ielem = 0; ielem < rnum_elem; ielem++){
      std::cout << "elem:  " << ielem+1 << std::endl;
      for (int lnode = 0; lnode < 8; lnode++){
          std::cout << "{ ";
            std::cout << lnode+1 << " = " << nodes_in_elem(ielem,lnode) + 1 << " ";

          std::cout << " }"<< std::endl;
      }
      std::cout << std::endl;
    }
    */
} // end read_mesh

/* ----------------------------------------------------------------------
   Read Tecplot format mesh file
------------------------------------------------------------------------- */
void Solver::read_mesh_tecplot(const char* MESH)
{
    Input_Options input_options = simparam.input_options.value();

    char ch;
    int  num_dim = simparam.num_dims;
    int  p_order = input_options.p_order;
    int  negative_index_found = 0;
    int  global_negative_index_found = 0;
    int  buffer_loop, buffer_iteration, buffer_iterations, dof_limit, scan_loop;
    int  local_node_index, current_column_index;

    bool zero_index_base = input_options.zero_index_base;
    bool restart_file    = input_options.topology_optimization_restart;
    simparam.restart_file = restart_file;

    size_t strain_count, read_index_start, node_rid, elem_gid;

    std::string skip_line, read_line, substring;
    std::stringstream line_parse;

    CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;

    GO node_gid;

    real_t unit_scaling = input_options.unit_scaling;
    real_t dof_value;

    host_vec_array node_densities;

    // Nodes_Per_Element_Type =  elements::elem_types::Nodes_Per_Element_Type;

    // read the mesh
    // PLACEHOLDER: ensight_format(MESH);
    // abaqus_format(MESH);
    // vtk_format(MESH)

    // task 0 reads file
    if (myrank == 0)
    {
        in = new std::ifstream();
        in->open(MESH);
        // skip 2 lines
        for (int j = 1; j <= 2; j++)
        {
            getline(*in, skip_line);
            std::cout << skip_line << std::endl;
        } // for
    }

    // --- Read the number of nodes in the mesh --- //
    if (myrank == 0)
    {
        getline(*in, read_line);
        line_parse.str(read_line);
        // stop when the NODES= string is reached
        while (!line_parse.eof()) {
            line_parse >> substring;
            if (!substring.compare("NODES="))
            {
                line_parse >> num_nodes;
            }
            if (!substring.compare("ELEMENTS="))
            {
                line_parse >> num_elem;
            }
        } // while
        std::cout << "declared node count: " << num_nodes << std::endl;
        std::cout << "declared element count: " << num_elem << std::endl;
        if (num_elem <= 0)
        {
            std::cout << "ERROR, NO ELEMENTS IN MESH!!!!" << std::endl;
        }
    }

    // broadcast number of nodes
    MPI_Bcast(&num_nodes, 1, MPI_LONG_LONG_INT, 0, world);

    // construct contiguous parallel row map now that we know the number of nodes
    map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(num_nodes, 0, comm));

    // set the vertices in the mesh read in
    nlocal_nodes = map->getLocalNumElements();
    // populate local row offset data from global data
    global_size_t min_gid    = map->getMinGlobalIndex();
    global_size_t max_gid    = map->getMaxGlobalIndex();
    global_size_t index_base = map->getIndexBase();
    // debug print
    // std::cout << "local node count on task: " << " " << nlocal_nodes << std::endl;

    // allocate node storage with dual view
    // dual_node_coords = dual_vec_array("dual_node_coords", nlocal_nodes,num_dim);
    // if(restart_file)
    // dual_node_densities = dual_vec_array("dual_node_densities", nlocal_nodes,1);

    // local variable for host view in the dual view
    node_coords_distributed = Teuchos::rcp(new MV(map, num_dim));
    // active view scrope
    {
        host_vec_array node_coords = node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
        if (restart_file)
        {
            design_node_densities_distributed = Teuchos::rcp(new MV(map, 1));
            node_densities = design_node_densities_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
        }
        // host_vec_array node_coords = dual_node_coords.view_host();
        // if(restart_file)
        // node_densities = dual_node_densities.view_host();
        // notify that the host view is going to be modified in the file readin
        // dual_node_coords.modify_host();
        // if(restart_file)
        // dual_node_densities.modify_host();

        // old swage method
        // mesh->init_nodes(local_nrows); // add 1 for index starting at 1

        //std::cout << "Num nodes assigned to task " << myrank << " = " << nlocal_nodes << std::endl;

        // read the initial mesh coordinates
        // x-coords
        /*only task 0 reads in nodes and elements from the input file
        stores node data in a buffer and communicates once the buffer cap is reached
        or the data ends*/

        words_per_line = input_options.words_per_line;
        if (restart_file)
        {
            words_per_line++;
        }
        elem_words_per_line = input_options.elem_words_per_line;

        // allocate read buffer
        read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES, words_per_line, MAX_WORD);

        dof_limit = num_nodes;
        buffer_iterations = dof_limit / BUFFER_LINES;
        if (dof_limit % BUFFER_LINES != 0)
        {
            buffer_iterations++;
        }

        // read coords, also density if restarting
        read_index_start = 0;
        for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
        {
            // pack buffer on rank 0
            if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
            {
                for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++)
                {
                    getline(*in, read_line);
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
                while (buffer_iteration * BUFFER_LINES + buffer_loop < num_nodes) {
                    getline(*in, read_line);
                    line_parse.clear();
                    line_parse.str(read_line);
                    for (int iword = 0; iword < words_per_line; iword++)
                    {
                        // read portions of the line into the substring variable
                        line_parse >> substring;
                        // assign the substring variable as a word of the read buffer
                        strcpy(&read_buffer(buffer_loop, iword, 0), substring.c_str());
                    }
                    buffer_loop++;
                }
            }

            // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
            MPI_Bcast(read_buffer.pointer(), BUFFER_LINES * words_per_line * MAX_WORD, MPI_CHAR, 0, world);
            // broadcast how many nodes were read into this buffer iteration
            MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, world);

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
                if (map->isNodeGlobalElement(node_gid))
                {
                    // set local node index in this mpi rank
                    node_rid = map->getLocalElement(node_gid);
                    // extract nodal position from the read buffer
                    // for tecplot format this is the three coords in the same line
                    dof_value = atof(&read_buffer(scan_loop, 0, 0));
                    node_coords(node_rid, 0) = dof_value * unit_scaling;
                    dof_value = atof(&read_buffer(scan_loop, 1, 0));
                    node_coords(node_rid, 1) = dof_value * unit_scaling;
                    if (num_dim == 3)
                    {
                        dof_value = atof(&read_buffer(scan_loop, 2, 0));
                        node_coords(node_rid, 2) = dof_value * unit_scaling;
                    }
                    if (restart_file)
                    {
                        dof_value = atof(&read_buffer(scan_loop, num_dim, 0));
                        node_densities(node_rid, 0) = dof_value;
                    }
                    // extract density if restarting
                }
            }
            read_index_start += BUFFER_LINES;
        }
    }
    // repartition node distribution
    repartition_nodes();

    // synchronize device data
    // dual_node_coords.sync_device();
    // dual_node_coords.modify_device();
    // if(restart_file){
    // dual_node_densities.sync_device();
    // dual_node_densities.modify_device();
    // }

    // debug print of nodal data

    // debug print nodal positions and indices

    // std::cout << " ------------NODAL POSITIONS ON TASK " << myrank << " --------------"<<std::endl;
    // for (int inode = 0; inode < local_nrows; inode++){
    // std::cout << "node: " << map->getGlobalElement(inode) + 1 << " { ";
    // for (int istride = 0; istride < num_dim; istride++){
    // std::cout << node_coords(inode,istride) << " , ";
    // }
    // std::cout << node_densities(inode,0);
    // std::cout << " }"<< std::endl;
    // }

    // check that local assignments match global total

    // read in element info (supported tecplot format currently assumes one type)

    CArrayKokkos<int, array_layout, HostSpace, memory_traits> node_store(elem_words_per_line);

    // broadcast number of elements
    MPI_Bcast(&num_elem, 1, MPI_LONG_LONG_INT, 0, world);
    // std::cout<<"before initial mesh initialization"<<std::endl;

    // read in element connectivity
    // we're gonna reallocate for the words per line expected for the element connectivity
    read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES, elem_words_per_line, MAX_WORD);

    // calculate buffer iterations to read number of lines
    buffer_iterations = num_elem / BUFFER_LINES;
    int assign_flag;

    // dynamic buffer used to store elements before we know how many this rank needs
    std::vector<size_t> element_temp(BUFFER_LINES * elem_words_per_line);
    std::vector<size_t> global_indices_temp(BUFFER_LINES);
    size_t buffer_max = BUFFER_LINES * elem_words_per_line;
    size_t indices_buffer_max = BUFFER_LINES;

    if (num_elem % BUFFER_LINES != 0)
    {
        buffer_iterations++;
    }
    read_index_start = 0;
    // std::cout << "ELEMENT BUFFER ITERATIONS: " << buffer_iterations << std::endl;
    rnum_elem = 0;
    // std::cout << "BUFFER ITERATIONS IS: " << buffer_iterations << std::endl;
    for (buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++)
    {
        // pack buffer on rank 0
        if (myrank == 0 && buffer_iteration < buffer_iterations - 1)
        {
            for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++)
            {
                getline(*in, read_line);
                line_parse.clear();
                line_parse.str(read_line);
                for (int iword = 0; iword < elem_words_per_line; iword++)
                {
                    // read portions of the line into the substring variable
                    line_parse >> substring;
                    // assign the substring variable as a word of the read buffer
                    strcpy(&read_buffer(buffer_loop, iword, 0), substring.c_str());
                }
            }
        }
        else if (myrank == 0)
        {
            buffer_loop = 0;
            while (buffer_iteration * BUFFER_LINES + buffer_loop < num_elem) {
                getline(*in, read_line);
                line_parse.clear();
                line_parse.str(read_line);
                for (int iword = 0; iword < elem_words_per_line; iword++)
                {
                    // read portions of the line into the substring variable
                    line_parse >> substring;
                    // assign the substring variable as a word of the read buffer
                    strcpy(&read_buffer(buffer_loop, iword, 0), substring.c_str());
                }
                buffer_loop++;
                // std::cout<<" "<< node_coords(node_gid, 0)<<std::endl;
            }
        }

        // broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
        MPI_Bcast(read_buffer.pointer(), BUFFER_LINES * elem_words_per_line * MAX_WORD, MPI_CHAR, 0, world);
        // broadcast how many nodes were read into this buffer iteration
        MPI_Bcast(&buffer_loop, 1, MPI_INT, 0, world);

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
                    if (map->isNodeGlobalElement(node_gid) && !assign_flag)
                    {
                        assign_flag = 1;
                        rnum_elem++;
                    }
                }
                else
                {
                    if (map->isNodeGlobalElement(node_gid - 1) && !assign_flag)
                    {
                        assign_flag = 1;
                        rnum_elem++;
                    }
                }
            }

            if (assign_flag)
            {
                for (int inode = 0; inode < elem_words_per_line; inode++)
                {
                    if ((rnum_elem - 1) * elem_words_per_line + inode >= buffer_max)
                    {
                        element_temp.resize((rnum_elem - 1) * elem_words_per_line + inode + BUFFER_LINES * elem_words_per_line);
                        buffer_max = (rnum_elem - 1) * elem_words_per_line + inode + BUFFER_LINES * elem_words_per_line;
                    }
                    element_temp[(rnum_elem - 1) * elem_words_per_line + inode] = node_store(inode);
                    // std::cout << "VECTOR STORAGE FOR ELEM " << rnum_elem << " ON TASK " << myrank << " NODE " << inode+1 << " IS " << node_store(inode) + 1 << std::endl;
                }
                // assign global element id to temporary list
                if (rnum_elem - 1 >= indices_buffer_max)
                {
                    global_indices_temp.resize(rnum_elem - 1 + BUFFER_LINES);
                    indices_buffer_max = rnum_elem - 1 + BUFFER_LINES;
                }
                global_indices_temp[rnum_elem - 1] = elem_gid;
            }
        }
        read_index_start += BUFFER_LINES;
    }

    std::cout << "RNUM ELEMENTS IS: " << rnum_elem << std::endl;
    // copy temporary element storage to multivector storage
    Element_Types = CArrayKokkos<elements::elem_types::elem_type, array_layout, HostSpace, memory_traits>(rnum_elem);

    elements::elem_types::elem_type mesh_element_type;

    if (simparam.num_dims == 2)
    {
        if (input_options.element_type == ELEMENT_TYPE::quad4)
        {
            mesh_element_type   = elements::elem_types::Quad4;
            max_nodes_per_patch = 2;
        }
        else if (input_options.element_type == ELEMENT_TYPE::quad8)
        {
            mesh_element_type   = elements::elem_types::Quad8;
            max_nodes_per_patch = 3;
        }
        else if (input_options.element_type == ELEMENT_TYPE::quad12)
        {
            mesh_element_type   = elements::elem_types::Quad12;
            max_nodes_per_patch = 4;
        }
        else
        {
            if (myrank == 0)
            {
                std::cout << "ELEMENT TYPE UNRECOGNIZED" << std::endl;
            }
            exit_solver(0);
        }
        element_select->choose_2Delem_type(mesh_element_type, elem2D);
        max_nodes_per_element = elem2D->num_nodes();
    }

    if (simparam.num_dims == 3)
    {
        if (input_options.element_type == ELEMENT_TYPE::hex8)
        {
            mesh_element_type   = elements::elem_types::Hex8;
            max_nodes_per_patch = 4;
        }
        else if (input_options.element_type == ELEMENT_TYPE::hex20)
        {
            mesh_element_type   = elements::elem_types::Hex20;
            max_nodes_per_patch = 8;
        }
        else if (input_options.element_type == ELEMENT_TYPE::hex32)
        {
            mesh_element_type   = elements::elem_types::Hex32;
            max_nodes_per_patch = 12;
        }
        else
        {
            if (myrank == 0)
            {
                std::cout << "ELEMENT TYPE UNRECOGNIZED" << std::endl;
            }
            exit_solver(0);
        }
        element_select->choose_3Delem_type(mesh_element_type, elem);
        max_nodes_per_element = elem->num_nodes();
    }

    dual_nodes_in_elem = dual_elem_conn_array("dual_nodes_in_elem", rnum_elem, max_nodes_per_element);
    host_elem_conn_array nodes_in_elem = dual_nodes_in_elem.view_host();
    dual_nodes_in_elem.modify_host();

    for (int ielem = 0; ielem < rnum_elem; ielem++)
    {
        for (int inode = 0; inode < elem_words_per_line; inode++)
        {
            nodes_in_elem(ielem, inode) = element_temp[ielem * elem_words_per_line + inode];
        }
    }

    // view storage for all local elements connected to local nodes on this rank
    Kokkos::DualView<GO*, array_layout, device_type, memory_traits> All_Element_Global_Indices("All_Element_Global_Indices", rnum_elem);
    // copy temporary global indices storage to view storage
    for (int ielem = 0; ielem < rnum_elem; ielem++)
    {
        All_Element_Global_Indices.h_view(ielem) = global_indices_temp[ielem];
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
        exit_solver(0);
    }

    // delete temporary element connectivity and index storage
    std::vector<size_t>().swap(element_temp);
    std::vector<size_t>().swap(global_indices_temp);

    All_Element_Global_Indices.modify_host();
    All_Element_Global_Indices.sync_device();

    // construct overlapping element map (since different ranks can own the same elements due to the local node map)
    all_element_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(Teuchos::OrdinalTraits<GO>::invalid(), All_Element_Global_Indices.d_view, 0, comm));

    // element type selection (subject to change)
    // ---- Set Element Type ---- //
    // allocate element type memory
    // elements::elem_type_t* elem_choice;

    int NE = 1; // number of element types in problem

    // Convert ijk index system to the finite element numbering convention
    // for vertices in cell
    if (active_node_ordering_convention == IJK)
    {
        CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_ensight_to_ijk(max_nodes_per_element);
        CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> tmp_ijk_indx(max_nodes_per_element);
        convert_ensight_to_ijk(0) = 0;
        convert_ensight_to_ijk(1) = 1;
        convert_ensight_to_ijk(2) = 3;
        convert_ensight_to_ijk(3) = 2;
        convert_ensight_to_ijk(4) = 4;
        convert_ensight_to_ijk(5) = 5;
        convert_ensight_to_ijk(6) = 7;
        convert_ensight_to_ijk(7) = 6;

        int nodes_per_element;

        if (num_dim == 2)
        {
            for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++)
            {
                // set nodes per element
                element_select->choose_2Delem_type(Element_Types(cell_rid), elem2D);
                nodes_per_element = elem2D->num_nodes();
                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    tmp_ijk_indx(node_lid) = nodes_in_elem(cell_rid, convert_ensight_to_ijk(node_lid));
                }

                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    nodes_in_elem(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
                }
            }
        }

        if (num_dim == 3)
        {
            for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++)
            {
                // set nodes per element
                element_select->choose_3Delem_type(Element_Types(cell_rid), elem);
                nodes_per_element = elem->num_nodes();
                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    tmp_ijk_indx(node_lid) = nodes_in_elem(cell_rid, convert_ensight_to_ijk(node_lid));
                }

                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    nodes_in_elem(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
                }
            }
        }
    }


    //Find initial objective value to normalize by
    if(restart_file){
        if (myrank == 0)
        {
            bool found = false;
            while (found == false&&in->good()) {
                std::getline(*in, read_line);
                //std::cout << read_line << std::endl;
                line_parse.clear();
                line_parse.str(read_line);

                //stop when the design_density string is reached
                while (!line_parse.eof()){
                    line_parse >> substring;
                    //std::cout << substring << std::endl;
                    if(!substring.compare("Objective_Normalization_Constant")){
                        found = true;
                        line_parse >> substring;
                        simparam.optimization_options.objective_normalization_constant = stod(substring);
                        std::cout << "NORMALIZATION CONSTANT FOR OBJECTIVE " << 
                            simparam.optimization_options.objective_normalization_constant << std::endl;
                    }
                } //while

            } // end while

            if (!found){
                throw std::runtime_error("ERROR: Failed to find initial objective value for restart");
            } // end if
        } // end if(myrank==0)
    }

    // Close mesh input file
    if (myrank == 0)
    {
        in->close();
    }

    // debug print element edof

    // std::cout << " ------------ELEMENT EDOF ON TASK " << myrank << " --------------"<<std::endl;

    // for (int ielem = 0; ielem < rnum_elem; ielem++){
    // std::cout << "elem:  " << ielem+1 << std::endl;
    // for (int lnode = 0; lnode < 8; lnode++){
    // std::cout << "{ ";
    // std::cout << lnode+1 << " = " << nodes_in_elem(ielem,lnode) + 1 << " ";

    // std::cout << " }"<< std::endl;
    // }
    // std::cout << std::endl;
    // }
} // end read_mesh

/* ----------------------------------------------------------------------
   Read Abaqus .inp format mesh file
------------------------------------------------------------------------- */
void Solver::read_mesh_abaqus_inp(const char *MESH){
  Input_Options input_options = simparam.input_options.value();
  char ch;
  int num_dim = simparam.num_dims;
  int p_order = input_options.p_order;
  real_t unit_scaling = input_options.unit_scaling;
  bool restart_file = simparam.restart_file;
  int local_node_index, current_column_index;
  size_t strain_count;
  std::string skip_line, read_line, substring, token;
  std::stringstream line_parse;
  CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;
  int buffer_loop, buffer_iteration, buffer_iterations, dof_limit, scan_loop, nodes_per_element;
  size_t read_index_start, node_rid, elem_gid;
  GO node_gid;
  real_t dof_value;
  host_vec_array node_densities;
  bool zero_index_base = input_options.zero_index_base;
  int negative_index_found = 0;
  int global_negative_index_found = 0;

  //Nodes_Per_Element_Type =  elements::elem_types::Nodes_Per_Element_Type;

  //read the mesh
  //PLACEHOLDER: ensight_format(MESH);
  // abaqus_format(MESH);
  // vtk_format(MESH)

  //task 0 reads file
  if(myrank==0){
    in = new std::ifstream();
    in->open(MESH);
  }

  //Abaqus inp file doesn't specify total number of nodes, which is needed for the node map.
  //First pass reads in node section to determine the maximum number of nodes, second pass distributes node data
  //The elements section header does specify element count
  num_nodes = 0;
  if(myrank==0){
    bool searching_for_nodes = true;
    //skip lines at the top with nonessential info; stop skipping when "*Nodes" string is reached
    while (searching_for_nodes&&in->good()) {
      getline(*in, skip_line);
      //std::cout << skip_line << std::endl;
      line_parse.clear();
      line_parse.str(skip_line);
      //stop when the NODES= string is reached
      while (!line_parse.eof()){
        line_parse >> substring;
        //std::cout << substring << std::endl;
        if(!substring.compare("*Node")){
          searching_for_nodes = false;
          break;
        }
      } //while
      
    }
    if(searching_for_nodes){
      std::cout << "FILE FORMAT ERROR" << std::endl;
    }

    //tally node count (bug in dat files seems to print wrong node count so read is done in two passes)
    //stop when apparent "-1" zone delimiter is reacher
    searching_for_nodes = true;
    int node_tally = 0;
    while (searching_for_nodes) {
      getline(*in, read_line);
      //std::cout << read_line << std::endl;
      line_parse.clear();
      line_parse.str(read_line);
      line_parse >> substring;
      
      //std::cout << substring << std::endl;
      if(substring == "*Element,"){
        searching_for_nodes = false;
          break;
      }
      else{
        node_tally++;
      }
      
    }
    num_nodes = node_tally;
    std::cout << "declared node count: " << num_nodes << std::endl;
  }

  //broadcast number of nodes
  MPI_Bcast(&num_nodes,1,MPI_LONG_LONG_INT,0,world);
  
  //construct contiguous parallel row map now that we know the number of nodes
  map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(num_nodes,0,comm));
  
  //close and reopen file for second pass now that global node count is known
  if(myrank==0){
    in->close();
    //in = new std::ifstream();
    in->open(MESH);
  }

  // set the vertices in the mesh read in
  global_size_t local_nrows = map->getLocalNumElements();
  nlocal_nodes = local_nrows;
  //populate local row offset data from global data
  global_size_t min_gid = map->getMinGlobalIndex();
  global_size_t max_gid = map->getMaxGlobalIndex();
  global_size_t index_base = map->getIndexBase();
  //debug print
  //std::cout << "local node count on task: " << " " << nlocal_nodes << std::endl;

  //allocate node storage with dual view
  //dual_node_coords = dual_vec_array("dual_node_coords", nlocal_nodes,num_dim);
  //if(restart_file)
    //dual_node_densities = dual_vec_array("dual_node_densities", nlocal_nodes,1);

  //local variable for host view in the dual view
  node_coords_distributed = Teuchos::rcp(new MV(map, num_dim));
  //view scope
  {
  host_vec_array node_coords = node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  //host_vec_array node_coords = dual_node_coords.view_host();
  if(restart_file){
    design_node_densities_distributed = Teuchos::rcp(new MV(map, 1));
    node_densities = design_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  }
  //notify that the host view is going to be modified in the file readin
  //dual_node_coords.modify_host();
  //if(restart_file)
    //dual_node_densities.modify_host();

  //old swage method
  //mesh->init_nodes(local_nrows); // add 1 for index starting at 1
    
  //std::cout << "Num nodes assigned to task " << myrank << " = " << nlocal_nodes << std::endl;

  // read the initial mesh coordinates
  // x-coords
  /*only task 0 reads in nodes and elements from the input file
  stores node data in a buffer and communicates once the buffer cap is reached
  or the data ends*/

  words_per_line = input_options.words_per_line;
  //if(restart_file) words_per_line++;
  elem_words_per_line = input_options.elem_words_per_line;

  //allocate read buffer
  read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES,words_per_line,MAX_WORD);

  dof_limit = num_nodes;
  buffer_iterations = dof_limit/BUFFER_LINES;
  if(dof_limit%BUFFER_LINES!=0) buffer_iterations++;
  
  //second pass to now read node coords with global node map defines
  if(myrank==0){
    bool searching_for_nodes = true;
    //skip lines at the top with nonessential info; stop skipping when "Nodes for the whole assembly" string is reached
    while (searching_for_nodes&&in->good()) {
      getline(*in, skip_line);
      //std::cout << skip_line << std::endl;
      line_parse.clear();
      line_parse.str(skip_line);
      //stop when the NODES= string is reached
      while (!line_parse.eof()){
        line_parse >> substring;
        //std::cout << substring << std::endl;
        if(!substring.compare("*Node")){
          searching_for_nodes = false;
          break;
        }
      } //while
      
    }
    if(searching_for_nodes){
      std::cout << "FILE FORMAT ERROR" << std::endl;
    }
  }
  
  //read coords, also density if restarting
  read_index_start = 0;
  for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
    //pack buffer on rank 0
    if(myrank==0&&buffer_iteration<buffer_iterations-1){
      for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);

        line_parse >> substring; //skip node index column since coding for sorted inp
        
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //debug print
        //std::cout<<" "<< substring <<std::endl;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
      }
    }
    else if(myrank==0){
      buffer_loop=0;
      while(buffer_iteration*BUFFER_LINES+buffer_loop < num_nodes) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);

        line_parse >> substring; //skip node index column since coding for sorted inp

        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
        buffer_loop++;
      }
      
    }

    //broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
    MPI_Bcast(read_buffer.pointer(),BUFFER_LINES*words_per_line*MAX_WORD,MPI_CHAR,0,world);
    //broadcast how many nodes were read into this buffer iteration
    MPI_Bcast(&buffer_loop,1,MPI_INT,0,world);

    //debug_print
    //std::cout << "NODE BUFFER LOOP IS: " << buffer_loop << std::endl;
    //for(int iprint=0; iprint < buffer_loop; iprint++)
      //std::cout<<"buffer packing: " << std::string(&read_buffer(iprint,0,0)) << std::endl;
    //return;

    //determine which data to store in the swage mesh members (the local node data)
    //loop through read buffer
    for(scan_loop = 0; scan_loop < buffer_loop; scan_loop++){
      //set global node id (ensight specific order)
      node_gid = read_index_start + scan_loop;
      //let map decide if this node id belongs locally; if yes store data
      if(map->isNodeGlobalElement(node_gid)){
        //set local node index in this mpi rank
        node_rid = map->getLocalElement(node_gid);
        //extract nodal position from the read buffer
        //for tecplot format this is the three coords in the same line
        dof_value = atof(&read_buffer(scan_loop,0,0));
        node_coords(node_rid, 0) = dof_value * unit_scaling;
        dof_value = atof(&read_buffer(scan_loop,1,0));
        node_coords(node_rid, 1) = dof_value * unit_scaling;
        if(num_dim==3){
            dof_value = atof(&read_buffer(scan_loop,2,0));
            node_coords(node_rid, 2) = dof_value * unit_scaling;
        }
        //extract density if restarting
      }
    }
    read_index_start+=BUFFER_LINES;
  }
  } //end view scope
  //repartition node distribution
  repartition_nodes();

  //synchronize device data
  //dual_node_coords.sync_device();
  //dual_node_coords.modify_device();
  //if(restart_file){
    //dual_node_densities.sync_device();
    //dual_node_densities.modify_device();
  //}

  //debug print of nodal data
  
  //debug print nodal positions and indices
  /*
  std::cout << " ------------NODAL POSITIONS ON TASK " << myrank << " --------------"<<std::endl;
  for (int inode = 0; inode < local_nrows; inode++){
      std::cout << "node: " << map->getGlobalElement(inode) + 1 << " { ";
    for (int istride = 0; istride < num_dim; istride++){
       std::cout << node_coords(inode,istride) << " , ";
    }
    //std::cout << node_densities(inode,0);
    std::cout << " }"<< std::endl;
  }
  */

  //check that local assignments match global total

  
  //read in element info
  //seek element connectivity zone
  int etype_index = 0;
  if(myrank==0){
    bool searching_for_elements = true;
    //skip lines at the top with nonessential info; stop skipping when "Nodes for the whole assembly" string is reached
    while (searching_for_elements&&in->good()) {
      getline(*in, skip_line);
      //std::cout << skip_line << std::endl;
      line_parse.clear();
      line_parse.str(skip_line);
      //stop when the NODES= string is reached
      while (!line_parse.eof()){
        line_parse >> substring;
        //std::cout << substring << std::endl;
        if(!substring.compare("*Element,")){
          searching_for_elements = false;
          break;
        }
      } //while
      
    }
    if(searching_for_elements){
      std::cout << "FILE FORMAT ERROR" << std::endl;
    }

    if(in->good())
      first_elem_line_streampos = in->tellg();

    //tally node count (bug in dat files seems to print wrong node count so read is done in two passes)
    //stop when apparent "-1" zone delimiter is reacher
    searching_for_elements = true;
    int elem_tally = 0;
    while (searching_for_elements&&in->good()) {
      getline(*in, read_line);
    //   if(elem_tally < 1000)
    //   std::cout << read_line << std::endl;
      line_parse.clear();
      line_parse.str(read_line);
      line_parse >> substring;
      
      //std::cout << substring << std::endl;
      if(substring == "*End"){
        searching_for_elements = false;
          break;
      }
      else{
        elem_tally++;
      }
      
    }
    num_elem = elem_tally;
    
    in->seekg(first_elem_line_streampos);
  }

  //broadcast element type
  etype_index = 1;
  MPI_Bcast(&etype_index,1,MPI_INT,0,world);

  elements::elem_types::elem_type mesh_element_type;
  int elem_words_per_line_no_nodes = 0;
  if(etype_index==1){
    mesh_element_type = elements::elem_types::Hex8;
    nodes_per_element = 8;
    elem_words_per_line = 8;
    max_nodes_per_patch = 4;
  }
  else if(etype_index==2){
    mesh_element_type = elements::elem_types::Hex20;
    nodes_per_element = 20;
    elem_words_per_line = 20;
    max_nodes_per_patch = 8;
  }
  else if(etype_index==3){
    mesh_element_type = elements::elem_types::Hex32;
    nodes_per_element = 32;
    elem_words_per_line = 32;
    max_nodes_per_patch = 12;
  }
  else{
    *fos << "ERROR: ABAQUS ELEMENT TYPE NOT FOUND OR RECOGNIZED" << std::endl;
    exit_solver(0);
  }
  
  //broadcast number of elements
  MPI_Bcast(&num_elem,1,MPI_LONG_LONG_INT,0,world);
  
  *fos << "declared element count: " << num_elem << std::endl;
  //std::cout<<"before initial mesh initialization"<<std::endl;
  
  //read in element connectivity
  //we're gonna reallocate for the words per line expected for the element connectivity
  read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES,elem_words_per_line,MAX_WORD); 
  CArrayKokkos<int, array_layout, HostSpace, memory_traits> node_store(nodes_per_element);

  //calculate buffer iterations to read number of lines
  buffer_iterations = num_elem/BUFFER_LINES;
  int assign_flag;
  //dynamic buffer used to store elements before we know how many this rank needs
  std::vector<size_t> element_temp(BUFFER_LINES*elem_words_per_line);
  std::vector<size_t> global_indices_temp(BUFFER_LINES);
  size_t buffer_max = BUFFER_LINES*elem_words_per_line;
  size_t indices_buffer_max = BUFFER_LINES;

  if(num_elem%BUFFER_LINES!=0) buffer_iterations++;
  read_index_start = 0;
  //std::cout << "ELEMENT BUFFER ITERATIONS: " << buffer_iterations << std::endl;
  rnum_elem = 0;
  //std::cout << "BUFFER ITERATIONS IS: " << buffer_iterations << std::endl;
  for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
    //pack buffer on rank 0
    if(myrank==0&&buffer_iteration<buffer_iterations-1){
      for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
        getline(*in,read_line);
        // if(buffer_iteration==0&&buffer_loop<100){
        //     std::cout << read_line << std::endl;
        // }
        line_parse.clear();
        line_parse.str(read_line);
        std::getline(line_parse, substring, ',');; //skip elem gid since coding for sorted inp
        for(int iword = 0; iword < elem_words_per_line; iword++){
        //read portions of the line into the substring variable
        std::getline(line_parse, substring, ',');
        //line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
      }
    }
    else if(myrank==0){
      buffer_loop=0;
      while(buffer_iteration*BUFFER_LINES+buffer_loop < num_elem) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        std::getline(line_parse, substring, ','); //skip elem gid since coding for sorted inp
        for(int iword = 0; iword < elem_words_per_line; iword++){
        //read portions of the line into the substring variable
        std::getline(line_parse, substring, ',');
        //line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
        buffer_loop++;
        //std::cout<<" "<< node_coords(node_gid, 0)<<std::endl;
      }
    }

    //broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
    MPI_Bcast(read_buffer.pointer(),BUFFER_LINES*elem_words_per_line*MAX_WORD,MPI_CHAR,0,world);
    //broadcast how many nodes were read into this buffer iteration
    MPI_Bcast(&buffer_loop,1,MPI_INT,0,world);
    
    //store element connectivity that belongs to this rank
    //loop through read buffer
    for(scan_loop = 0; scan_loop < buffer_loop; scan_loop++){
      //set global node id (ensight specific order)
      elem_gid = read_index_start + scan_loop;
      //add this element to the local list if any of its nodes belong to this rank according to the map
      //get list of nodes for each element line and check if they belong to the map
      assign_flag = 0;
      for(int inode = elem_words_per_line_no_nodes; inode < elem_words_per_line; inode++){
        //as we loop through the nodes belonging to this element we store them
        //if any of these nodes belongs to this rank this list is used to store the element locally
        node_gid = atoi(&read_buffer(scan_loop,inode,0));
        if(zero_index_base)
          node_store(inode-elem_words_per_line_no_nodes) = node_gid;
        else
          node_store(inode-elem_words_per_line_no_nodes) = node_gid - 1; //subtract 1 since file index start is 1 but code expects 0
        if(node_store(inode-elem_words_per_line_no_nodes) < 0){
          negative_index_found = 1;
        }
        //first we add the elements to a dynamically allocated list
        if(zero_index_base){
          if(map->isNodeGlobalElement(node_gid)&&!assign_flag){
            assign_flag = 1;
            rnum_elem++;
          }
        }
        else{
          if(map->isNodeGlobalElement(node_gid-1)&&!assign_flag){
            assign_flag = 1;
            rnum_elem++;
          }
        }
      }

      if(assign_flag){
        for(int inode = 0; inode < nodes_per_element; inode++){
          if((rnum_elem-1)*nodes_per_element + inode>=buffer_max){ 
            element_temp.resize((rnum_elem-1)*nodes_per_element + inode + BUFFER_LINES*nodes_per_element);
            buffer_max = (rnum_elem-1)*nodes_per_element + inode + BUFFER_LINES*nodes_per_element;
          }
          element_temp[(rnum_elem-1)*nodes_per_element + inode] = node_store(inode); 
          //std::cout << "VECTOR STORAGE FOR ELEM " << rnum_elem << " ON TASK " << myrank << " NODE " << inode+1 << " IS " << node_store(inode) + 1 << std::endl;
        }
        //assign global element id to temporary list
        if(rnum_elem-1>=indices_buffer_max){ 
          global_indices_temp.resize(rnum_elem-1 + BUFFER_LINES);
          indices_buffer_max = rnum_elem-1 + BUFFER_LINES;
        }
        global_indices_temp[rnum_elem-1] = elem_gid;
      }
    }
    read_index_start+=BUFFER_LINES;
  }
  
  //check if ABAQUS file has boundary and loading condition zones
  bool No_Conditions = true;
  if(myrank==0){
    if(in->good())
      before_condition_header = in->tellg();
    bool searching_for_conditions = true;
    //skip lines at the top with nonessential info; stop skipping when "Fixed Supports or Pressure" string is reached
    while (searching_for_conditions&&in->good()) {
      getline(*in, skip_line);
      //std::cout << skip_line << std::endl;
      line_parse.clear();
      line_parse.str(skip_line);
      //stop when the NODES= string is reached
      while (!line_parse.eof()){
        line_parse >> substring;
        //std::cout << substring << std::endl;
        if(!substring.compare("Supports")||!substring.compare("Pressure")){
          No_Conditions = searching_for_conditions = false;
          break;
        }
      } //while
      
    } //while
  }
  
  //broadcast search condition
  MPI_Bcast(&No_Conditions,1,MPI_CXX_BOOL,0,world);

  //flag elasticity fea module for boundary/loading conditions readin that remains
  if(!No_Conditions){
    // check that the input file has configured some kind of acceptable module
    simparam.validate_module_is_specified(FEA_MODULE_TYPE::Elasticity);
    simparam.fea_module_must_read.insert(FEA_MODULE_TYPE::Elasticity);
  }

  // Close mesh input file if no further readin is done by FEA modules for conditions
  if(myrank==0&&No_Conditions){
    in->close();
  }

  std::cout << "RNUM ELEMENTS IS: " << rnum_elem << std::endl;
  //copy temporary element storage to multivector storage
  Element_Types = CArrayKokkos<elements::elem_types::elem_type, array_layout, HostSpace, memory_traits>(rnum_elem);

  //set element object pointer
  if(simparam.num_dims==2){
    element_select->choose_2Delem_type(mesh_element_type, elem2D);
     max_nodes_per_element = elem2D->num_nodes();
  }
  else if(simparam.num_dims==3){
    element_select->choose_3Delem_type(mesh_element_type, elem);
     max_nodes_per_element = elem->num_nodes();
  }

  //1 type per mesh for now
  for(int ielem = 0; ielem < rnum_elem; ielem++)
    Element_Types(ielem) = mesh_element_type;

  dual_nodes_in_elem = dual_elem_conn_array("dual_nodes_in_elem", rnum_elem, max_nodes_per_element);
  host_elem_conn_array nodes_in_elem = dual_nodes_in_elem.view_host();
  dual_nodes_in_elem.modify_host();

  for(int ielem = 0; ielem < rnum_elem; ielem++)
    for(int inode = 0; inode < nodes_per_element; inode++){
      nodes_in_elem(ielem, inode) = element_temp[ielem*nodes_per_element + inode];
    }

  //view storage for all local elements connected to local nodes on this rank
  Kokkos::DualView <GO*, array_layout, device_type, memory_traits> All_Element_Global_Indices("All_Element_Global_Indices",rnum_elem);

  //copy temporary global indices storage to view storage
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    All_Element_Global_Indices.h_view(ielem) = global_indices_temp[ielem];
    if(global_indices_temp[ielem]<0){
      negative_index_found = 1;
    }
  }
  
  MPI_Allreduce(&negative_index_found,&global_negative_index_found,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  if(global_negative_index_found){
    if(myrank==0){
    std::cout << "Node index less than or equal to zero detected; set \"zero_index_base: true\" under \"input_options\" in your yaml file if indices start at 0" << std::endl;
    }
    exit_solver(0);
  }
  
  //debug print element edof
  /*
  std::cout << " ------------ELEMENT EDOF ON TASK " << myrank << " --------------"<<std::endl;
  
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    std::cout << "elem:  " << All_Element_Global_Indices(ielem)+1 << std::endl;
    for (int lnode = 0; lnode < 8; lnode++){
        std::cout << "{ ";
          std::cout << lnode+1 << " = " << nodes_in_elem(ielem,lnode) + 1 << " ";
        
        std::cout << " }"<< std::endl;
    }
    std::cout << std::endl;
  }
  */

  //delete temporary element connectivity and index storage
  std::vector<size_t>().swap(element_temp);
  std::vector<size_t>().swap(global_indices_temp);

  All_Element_Global_Indices.modify_host();
  All_Element_Global_Indices.sync_device();
  
  //construct overlapping element map (since different ranks can own the same elements due to the local node map)
  all_element_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),All_Element_Global_Indices.d_view,0,comm));


  //element type selection (subject to change)
  // ---- Set Element Type ---- //
  // allocate element type memory
  //elements::elem_type_t* elem_choice;

  int NE = 1; // number of element types in problem

  // Convert ijk index system to the finite element numbering convention
  // for vertices in cell
  CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_ensight_to_ijk(max_nodes_per_element);
  CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> tmp_ijk_indx(max_nodes_per_element);
  convert_ensight_to_ijk(0) = 0;
  convert_ensight_to_ijk(1) = 1;
  convert_ensight_to_ijk(2) = 3;
  convert_ensight_to_ijk(3) = 2;
  convert_ensight_to_ijk(4) = 4;
  convert_ensight_to_ijk(5) = 5;
  convert_ensight_to_ijk(6) = 7;
  convert_ensight_to_ijk(7) = 6;
  
  if(num_dim==2)
  for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++) {
    //set nodes per element
    element_select->choose_2Delem_type(Element_Types(cell_rid), elem2D);
    nodes_per_element = elem2D->num_nodes();
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      tmp_ijk_indx(node_lid) = nodes_in_elem(cell_rid, convert_ensight_to_ijk(node_lid));
    }   
        
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      nodes_in_elem(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
    }
  }

  if(num_dim==3)
  for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++) {
    //set nodes per element
    element_select->choose_3Delem_type(Element_Types(cell_rid), elem);
    nodes_per_element = elem->num_nodes();
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      tmp_ijk_indx(node_lid) = nodes_in_elem(cell_rid, convert_ensight_to_ijk(node_lid));
    }   
        
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      nodes_in_elem(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
    }
  }
 
} // end read_mesh

/* ----------------------------------------------------------------------
   Rebalance the initial node decomposition with Zoltan2
------------------------------------------------------------------------- */

void Solver::repartition_nodes(bool repartition_node_densities)
{
    char ch;

    int num_dim = simparam.num_dims;
    int local_node_index, current_column_index, nodes_per_element;

    size_t strain_count;

    std::stringstream line_parse;
    CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;

    GO node_gid;

    // construct input adapted needed by Zoltan2 problem
    typedef Xpetra::MultiVector<real_t, LO, GO, node_type> xvector_t;
    typedef Zoltan2::XpetraMultiVectorAdapter<xvector_t> inputAdapter_t;
    typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;

    Teuchos::RCP<xvector_t> xpetra_node_coords = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t, LO, GO, node_type>(node_coords_distributed));

    Teuchos::RCP<inputAdapter_t> problem_adapter =  Teuchos::rcp(new inputAdapter_t(xpetra_node_coords));

    // Create parameters for an RCB problem

    double tolerance = 1.05;

    Teuchos::ParameterList params("Node Partition Params");
    params.set("debug_level", "basic_status");
    params.set("debug_procs", "0");
    params.set("error_check_level", "debug_mode_assertions");

    // params.set("algorithm", "rcb");
    params.set("algorithm", "multijagged");
    params.set("imbalance_tolerance", tolerance);
    params.set("num_global_parts", nranks);
    params.set("partitioning_objective", "minimize_cut_edge_count");

    Teuchos::RCP<Zoltan2::PartitioningProblem<inputAdapter_t>> problem =
        Teuchos::rcp(new Zoltan2::PartitioningProblem<inputAdapter_t>(&(*problem_adapter), &params));

    // Solve the problem

    problem->solve();

    // create metric object where communicator is Teuchos default

    quality_t* metricObject1 = new quality_t(&(*problem_adapter), &params, // problem1->getComm(),
                                           &problem->getSolution());
    // Check the solution.

    if (myrank == 0)
    {
        metricObject1->printMetrics(std::cout);
    }

    if (myrank == 0)
    {
        real_t imb = metricObject1->getObjectCountImbalance();
        if (imb <= tolerance)
        {
            std::cout << "pass: " << imb << std::endl;
        }
        else
        {
            std::cout << "fail: " << imb << std::endl;
        }
        std::cout << std::endl;
    }
    delete metricObject1;

    // migrate rows of the vector so they correspond to the partition recommended by Zoltan2

    Teuchos::RCP<MV> partitioned_node_coords_distributed = Teuchos::rcp(new MV(map, num_dim));

    Teuchos::RCP<xvector_t> xpartitioned_node_coords_distributed =
        Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t, LO, GO, node_type>(partitioned_node_coords_distributed));

    problem_adapter->applyPartitioningSolution(*xpetra_node_coords, xpartitioned_node_coords_distributed, problem->getSolution());
    *partitioned_node_coords_distributed = Xpetra::toTpetra<real_t, LO, GO, node_type>(*xpartitioned_node_coords_distributed);

    Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> partitioned_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(*(partitioned_node_coords_distributed->getMap())));

    Teuchos::RCP<const Tpetra::Map<LO, GO, node_type>> partitioned_map_one_to_one;
    partitioned_map_one_to_one = Tpetra::createOneToOne<LO, GO, node_type>(partitioned_map);
    Teuchos::RCP<MV> partitioned_node_coords_one_to_one_distributed = Teuchos::rcp(new MV(partitioned_map_one_to_one, num_dim));

    Tpetra::Import<LO, GO> importer_one_to_one(partitioned_map, partitioned_map_one_to_one);
    partitioned_node_coords_one_to_one_distributed->doImport(*partitioned_node_coords_distributed, importer_one_to_one, Tpetra::INSERT);
    node_coords_distributed = partitioned_node_coords_one_to_one_distributed;
    partitioned_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(*partitioned_map_one_to_one));

    // migrate density vector if this is a restart file read
    if (simparam.restart_file&&repartition_node_densities)
    {
        Teuchos::RCP<MV> partitioned_node_densities_distributed = Teuchos::rcp(new MV(partitioned_map, 1));

        // create import object using local node indices map and all indices map
        Tpetra::Import<LO, GO> importer(map, partitioned_map);

        // comms to get ghosts
        partitioned_node_densities_distributed->doImport(*design_node_densities_distributed, importer, Tpetra::INSERT);
        design_node_densities_distributed = partitioned_node_densities_distributed;
    }

    // update nlocal_nodes and node map
    map = partitioned_map;
    nlocal_nodes = map->getLocalNumElements();
}

/* ----------------------------------------------------------------------
   Initialize Ghost and Non-Overlapping Element Maps
------------------------------------------------------------------------- */

void Solver::init_maps()
{
    char ch;
    int  num_dim = simparam.num_dims;
    int  local_node_index, current_column_index;
    int  nodes_per_element;
    GO   node_gid;

    host_elem_conn_array nodes_in_elem = dual_nodes_in_elem.view_host();

    if (rnum_elem >= 1)
    {
        // Construct set of ghost nodes; start with a buffer with upper limit
        size_t buffer_limit = 0;
        if (num_dim == 2)
        {
            for (int ielem = 0; ielem < rnum_elem; ielem++)
            {
                element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
                buffer_limit += elem2D->num_nodes();
            }
        }

        if (num_dim == 3)
        {
            for (int ielem = 0; ielem < rnum_elem; ielem++)
            {
                element_select->choose_3Delem_type(Element_Types(ielem), elem);
                buffer_limit += elem->num_nodes();
            }
        }

        CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> ghost_node_buffer(buffer_limit);

        std::set<GO> ghost_node_set;

        // search through local elements for global node indices not owned by this MPI rank
        if (num_dim == 2)
        {
            for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++)
            {
                // set nodes per element
                element_select->choose_2Delem_type(Element_Types(cell_rid), elem2D);
                nodes_per_element = elem2D->num_nodes();
                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    node_gid = nodes_in_elem(cell_rid, node_lid);
                    if (!map->isNodeGlobalElement(node_gid))
                    {
                        ghost_node_set.insert(node_gid);
                    }
                }
            }
        }

        if (num_dim == 3)
        {
            for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++)
            {
                // set nodes per element
                element_select->choose_3Delem_type(Element_Types(cell_rid), elem);
                nodes_per_element = elem->num_nodes();
                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    node_gid = nodes_in_elem(cell_rid, node_lid);
                    if (!map->isNodeGlobalElement(node_gid))
                    {
                        ghost_node_set.insert(node_gid);
                    }
                }
            }
        }

        // by now the set contains, with no repeats, all the global node indices that are ghosts for this rank
        // now pass the contents of the set over to a CArrayKokkos, then create a map to find local ghost indices from global ghost indices

        nghost_nodes     = ghost_node_set.size();
        ghost_nodes      = Kokkos::DualView<GO*, Kokkos::LayoutLeft, device_type, memory_traits>("ghost_nodes", nghost_nodes);
        ghost_node_ranks = Kokkos::DualView<int*, array_layout, device_type, memory_traits>("ghost_node_ranks", nghost_nodes);
        int  ighost = 0;
        auto it     = ghost_node_set.begin();

        while (it != ghost_node_set.end()) {
            ghost_nodes.h_view(ighost++) = *it;
            it++;
        }

        // debug print of ghost nodes
        // std::cout << " GHOST NODE SET ON TASK " << myrank << std::endl;
        // for(int i = 0; i < nghost_nodes; i++)
        // std::cout << "{" << i + 1 << "," << ghost_nodes(i) + 1 << "}" << std::endl;

        // find which mpi rank each ghost node belongs to and store the information in a CArrayKokkos
        // allocate Teuchos Views since they are the only input available at the moment in the map definitions
        Teuchos::ArrayView<const GO> ghost_nodes_pass(ghost_nodes.h_view.data(), nghost_nodes);

        Teuchos::ArrayView<int> ghost_node_ranks_pass(ghost_node_ranks.h_view.data(), nghost_nodes);

        map->getRemoteIndexList(ghost_nodes_pass, ghost_node_ranks_pass);

        // debug print of ghost nodes
        // std::cout << " GHOST NODE MAP ON TASK " << myrank << std::endl;
        // for(int i = 0; i < nghost_nodes; i++)
        // std::cout << "{" << i + 1 << "," << global2local_map.get(ghost_nodes(i)) + 1 << "}" << std::endl;
    }

    ghost_nodes.modify_host();
    ghost_nodes.sync_device();
    ghost_node_ranks.modify_host();
    ghost_node_ranks.sync_device();
    // create a Map for ghost node indices
    ghost_node_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(Teuchos::OrdinalTraits<GO>::invalid(), ghost_nodes.d_view, 0, comm));

    // Create reference element
    // ref_elem->init(p_order, num_dim, elem->num_basis());
    // std::cout<<"done with ref elem"<<std::endl;

    // communicate ghost node positions; construct multivector distributed object using local node data

    // construct array for all indices (ghost + local)
    nall_nodes = nlocal_nodes + nghost_nodes;
    // CArrayKokkos<GO, array_layout, device_type, memory_traits> all_node_indices(nall_nodes, "all_node_indices");
    Kokkos::DualView<GO*, array_layout, device_type, memory_traits> all_node_indices("all_node_indices", nall_nodes);
    for (int i = 0; i < nall_nodes; i++)
    {
        if (i < nlocal_nodes)
        {
            all_node_indices.h_view(i) = map->getGlobalElement(i);
        }
        else
        {
            all_node_indices.h_view(i) = ghost_nodes.h_view(i - nlocal_nodes);
        }
    }
    all_node_indices.modify_host();
    all_node_indices.sync_device();
    // debug print of node indices
    // for(int inode=0; inode < index_counter; inode++)
    // std::cout << " my_reduced_global_indices " << my_reduced_global_indices(inode) <<std::endl;

    // create a Map for all the node indices (ghost + local)
    all_node_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(Teuchos::OrdinalTraits<GO>::invalid(), all_node_indices.d_view, 0, comm));

    // remove elements from the local set so that each rank has a unique set of global ids

    // local elements belonging to the non-overlapping element distribution to each rank with buffer
    Kokkos::DualView<GO*, array_layout, device_type, memory_traits> Initial_Element_Global_Indices("Initial_Element_Global_Indices", rnum_elem);

    size_t nonoverlapping_count = 0;
    int    my_element_flag;

    // loop through local element set
    if (num_dim == 2)
    {
        for (int ielem = 0; ielem < rnum_elem; ielem++)
        {
            element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
            nodes_per_element = elem2D->num_nodes();

            my_element_flag = 1;
            for (int lnode = 0; lnode < nodes_per_element; lnode++)
            {
                node_gid = nodes_in_elem(ielem, lnode);
                if (ghost_node_map->isNodeGlobalElement(node_gid))
                {
                    local_node_index = ghost_node_map->getLocalElement(node_gid);
                    if (ghost_node_ranks.h_view(local_node_index) < myrank)
                    {
                        my_element_flag = 0;
                    }
                }
            }
            if (my_element_flag)
            {
                Initial_Element_Global_Indices.h_view(nonoverlapping_count++) = all_element_map->getGlobalElement(ielem);
            }
        }
    }

    if (num_dim == 3)
    {
        for (int ielem = 0; ielem < rnum_elem; ielem++)
        {
            element_select->choose_3Delem_type(Element_Types(ielem), elem);
            nodes_per_element = elem->num_nodes();

            my_element_flag = 1;

            for (int lnode = 0; lnode < nodes_per_element; lnode++)
            {
                node_gid = nodes_in_elem(ielem, lnode);
                if (ghost_node_map->isNodeGlobalElement(node_gid))
                {
                    local_node_index = ghost_node_map->getLocalElement(node_gid);
                    if (ghost_node_ranks.h_view(local_node_index) < myrank)
                    {
                        my_element_flag = 0;
                    }
                }
            }
            if (my_element_flag)
            {
                Initial_Element_Global_Indices.h_view(nonoverlapping_count++) = all_element_map->getGlobalElement(ielem);
            }
        }
    }

    // copy over from buffer to compressed storage
    Kokkos::DualView<GO*, array_layout, device_type, memory_traits> Element_Global_Indices("Element_Global_Indices", nonoverlapping_count);
    for (int ibuffer = 0; ibuffer < nonoverlapping_count; ibuffer++)
    {
        Element_Global_Indices.h_view(ibuffer) = Initial_Element_Global_Indices.h_view(ibuffer);
    }
    nlocal_elem_non_overlapping = nonoverlapping_count;
    Element_Global_Indices.modify_host();
    Element_Global_Indices.sync_device();
    // create nonoverlapping element map
    element_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(Teuchos::OrdinalTraits<GO>::invalid(), Element_Global_Indices.d_view, 0, comm));

    // sort element connectivity so nonoverlaps are sequentially found first
    // define initial sorting of global indices

    // element_map->describe(*fos,Teuchos::VERB_EXTREME);

    for (int ielem = 0; ielem < rnum_elem; ielem++)
    {
        Initial_Element_Global_Indices.h_view(ielem) = all_element_map->getGlobalElement(ielem);
    }

    // re-sort so local elements in the nonoverlapping map are first in storage
    CArrayKokkos<GO, array_layout, HostSpace, memory_traits> Temp_Nodes(max_nodes_per_element);

    GO  temp_element_gid, current_element_gid;
    int last_storage_index = rnum_elem - 1;

    for (int ielem = 0; ielem < nlocal_elem_non_overlapping; ielem++)
    {
        current_element_gid = Initial_Element_Global_Indices.h_view(ielem);
        // if this element is not part of the non overlap list then send it to the end of the storage and swap the element at the end
        if (!element_map->isNodeGlobalElement(current_element_gid))
        {
            temp_element_gid = current_element_gid;
            for (int lnode = 0; lnode < max_nodes_per_element; lnode++)
            {
                Temp_Nodes(lnode) = nodes_in_elem(ielem, lnode);
            }
            Initial_Element_Global_Indices.h_view(ielem) = Initial_Element_Global_Indices.h_view(last_storage_index);
            Initial_Element_Global_Indices.h_view(last_storage_index) = temp_element_gid;
            for (int lnode = 0; lnode < max_nodes_per_element; lnode++)
            {
                nodes_in_elem(ielem, lnode) = nodes_in_elem(last_storage_index, lnode);
                nodes_in_elem(last_storage_index, lnode) = Temp_Nodes(lnode);
            }
            last_storage_index--;

            // test if swapped element is also not part of the non overlap map; if so lower loop counter to repeat the above
            temp_element_gid = Initial_Element_Global_Indices.h_view(ielem);
            if (!element_map->isNodeGlobalElement(temp_element_gid))
            {
                ielem--;
            }
        }
    }
    // reset all element map to its re-sorted version
    Initial_Element_Global_Indices.modify_host();
    Initial_Element_Global_Indices.sync_device();

    all_element_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(Teuchos::OrdinalTraits<GO>::invalid(), Initial_Element_Global_Indices.d_view, 0, comm));
    // element_map->describe(*fos,Teuchos::VERB_EXTREME);
    // all_element_map->describe(*fos,Teuchos::VERB_EXTREME);

    // all_element_map->describe(*fos,Teuchos::VERB_EXTREME);
    // construct dof map that follows from the node map (used for distributed matrix and vector objects later)
    Kokkos::DualView<GO*, array_layout, device_type, memory_traits> local_dof_indices("local_dof_indices", nlocal_nodes * num_dim);
    for (int i = 0; i < nlocal_nodes; i++)
    {
        for (int j = 0; j < num_dim; j++)
        {
            local_dof_indices.h_view(i * num_dim + j) = map->getGlobalElement(i) * num_dim + j;
        }
    }

    local_dof_indices.modify_host();
    local_dof_indices.sync_device();
    local_dof_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(num_nodes * num_dim, local_dof_indices.d_view, 0, comm) );

    // construct dof map that follows from the all_node map (used for distributed matrix and vector objects later)
    Kokkos::DualView<GO*, array_layout, device_type, memory_traits> all_dof_indices("all_dof_indices", nall_nodes * num_dim);
    for (int i = 0; i < nall_nodes; i++)
    {
        for (int j = 0; j < num_dim; j++)
        {
            all_dof_indices.h_view(i * num_dim + j) = all_node_map->getGlobalElement(i) * num_dim + j;
        }
    }

    all_dof_indices.modify_host();
    all_dof_indices.sync_device();
    // pass invalid global count so the map reduces the global count automatically
    all_dof_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(Teuchos::OrdinalTraits<GO>::invalid(), all_dof_indices.d_view, 0, comm) );

    // debug print of map
    // debug print

    std::ostream& out = std::cout;

    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    // if(myrank==0)
    // *fos << "Ghost Node Map :" << std::endl;
    // all_node_map->describe(*fos,Teuchos::VERB_EXTREME);
    // *fos << std::endl;
    // std::fflush(stdout);

    // Count how many elements connect to each local node
    node_nconn_distributed = Teuchos::rcp(new MCONN(map, 1));
    // active view scope
    {
        host_elem_conn_array node_nconn = node_nconn_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
        for (int inode = 0; inode < nlocal_nodes; inode++)
        {
            node_nconn(inode, 0) = 0;
        }

        for (int ielem = 0; ielem < rnum_elem; ielem++)
        {
            for (int inode = 0; inode < nodes_per_element; inode++)
            {
                node_gid = nodes_in_elem(ielem, inode);
                if (map->isNodeGlobalElement(node_gid))
                {
                    node_nconn(map->getLocalElement(node_gid), 0)++;
                }
            }
        }
    }

    // create distributed multivector of the local node data and all (local + ghost) node storage

    all_node_coords_distributed   = Teuchos::rcp(new MV(all_node_map, num_dim));
    ghost_node_coords_distributed = Teuchos::rcp(new MV(ghost_node_map, num_dim));

    // create import object using local node indices map and all indices map
    comm_importer_setup();

    // create export objects for reverse comms
    comm_exporter_setup();

    // comms to get ghosts
    all_node_coords_distributed->doImport(*node_coords_distributed, *importer, Tpetra::INSERT);
    // all_node_nconn_distributed->doImport(*node_nconn_distributed, importer, Tpetra::INSERT);

    dual_nodes_in_elem.sync_device();
    dual_nodes_in_elem.modify_device();
    // construct distributed element connectivity multivector
    global_nodes_in_elem_distributed = Teuchos::rcp(new MCONN(all_element_map, dual_nodes_in_elem));

    // construct map of nodes that belong to the non-overlapping element set (contained by ghost + local node set but not all of them)
    std::set<GO> nonoverlap_elem_node_set;
    if (nlocal_elem_non_overlapping)
    {
        // search through local elements for global node indices not owned by this MPI rank
        if (num_dim == 2)
        {
            for (int cell_rid = 0; cell_rid < nlocal_elem_non_overlapping; cell_rid++)
            {
                // set nodes per element
                element_select->choose_2Delem_type(Element_Types(cell_rid), elem2D);
                nodes_per_element = elem2D->num_nodes();
                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    node_gid = nodes_in_elem(cell_rid, node_lid);
                    nonoverlap_elem_node_set.insert(node_gid);
                }
            }
        }

        if (num_dim == 3)
        {
            for (int cell_rid = 0; cell_rid < nlocal_elem_non_overlapping; cell_rid++)
            {
                // set nodes per element
                element_select->choose_3Delem_type(Element_Types(cell_rid), elem);
                nodes_per_element = elem->num_nodes();
                for (int node_lid = 0; node_lid < nodes_per_element; node_lid++)
                {
                    node_gid = nodes_in_elem(cell_rid, node_lid);
                    nonoverlap_elem_node_set.insert(node_gid);
                }
            }
        }
    }

    // by now the set contains, with no repeats, all the global node indices belonging to the non overlapping element list on this MPI rank
    // now pass the contents of the set over to a CArrayKokkos, then create a map to find local ghost indices from global ghost indices
    nnonoverlap_elem_nodes = nonoverlap_elem_node_set.size();
    nonoverlap_elem_nodes  = Kokkos::DualView<GO*, Kokkos::LayoutLeft, device_type, memory_traits>("nonoverlap_elem_nodes", nnonoverlap_elem_nodes);
    if(nnonoverlap_elem_nodes){
        int  inonoverlap_elem_node = 0;
        auto it = nonoverlap_elem_node_set.begin();
        while (it != nonoverlap_elem_node_set.end()) {
            nonoverlap_elem_nodes.h_view(inonoverlap_elem_node++) = *it;
            it++;
        }
        nonoverlap_elem_nodes.modify_host();
        nonoverlap_elem_nodes.sync_device();
    }

    // create a Map for node indices belonging to the non-overlapping set of elements
    nonoverlap_element_node_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(Teuchos::OrdinalTraits<GO>::invalid(), nonoverlap_elem_nodes.d_view, 0, comm));

    // std::cout << "number of patches = " << mesh->num_patches() << std::endl;
    if (myrank == 0)
    {
        std::cout << "End of map setup " << std::endl;
    }
}

/* ----------------------------------------------------------------------
   Find boundary surface segments that belong to this MPI rank
------------------------------------------------------------------------- */

void Solver::Get_Boundary_Patches()
{
    size_t     npatches_repeat, npatches, element_npatches, num_nodes_in_patch, node_gid;
    int        local_node_id;
    int        num_dim = simparam.num_dims;
    CArray<GO> Surface_Nodes;

    const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    // Surface_Nodes = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(4, "Surface_Nodes");

    std::set<Node_Combination> my_patches;
    // inititializes type for the pair variable (finding the iterator type is annoying)
    std::pair<std::set<Node_Combination>::iterator, bool> current_combination;

    std::set<Node_Combination>::iterator it;

    CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_node_order(max_nodes_per_element);
    if ((active_node_ordering_convention == ENSIGHT && num_dim == 3) || (active_node_ordering_convention == IJK && num_dim == 2))
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 3;
        convert_node_order(3) = 2;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 7;
            convert_node_order(7) = 6;
        }
    }
    else if ((active_node_ordering_convention == IJK && num_dim == 3) || (active_node_ordering_convention == ENSIGHT && num_dim == 2))
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 2;
        convert_node_order(3) = 3;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 6;
            convert_node_order(7) = 7;
        }
    }

    // compute the number of patches in this MPI rank with repeats for adjacent cells
    npatches_repeat = 0;

    if (num_dim == 2)
    {
        for (int ielem = 0; ielem < rnum_elem; ielem++)
        {
            element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
            element_npatches = elem2D->nsurfaces;
            npatches_repeat += element_npatches;
        }
    }
    else if (num_dim == 3)
    {
        for (int ielem = 0; ielem < rnum_elem; ielem++)
        {
            element_select->choose_3Delem_type(Element_Types(ielem), elem);
            element_npatches = elem->nsurfaces;
            npatches_repeat += element_npatches;
        }
    }
    // std::cout << "Starting boundary patch allocation of size " << npatches_repeat << std::endl <<std::flush;
    // information for all patches on this rank
    CArrayKokkos<Node_Combination, array_layout, HostSpace, memory_traits> Patch_Nodes(npatches_repeat, "Patch_Nodes");

    CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> Patch_Boundary_Flags(npatches_repeat, "Patch_Boundary_Flags");

    if (myrank == 0)
    {
        std::cout << "Done with boundary patch allocation" << std::endl << std::flush;
    }
    // initialize boundary patch flags
    for (int init = 0; init < npatches_repeat; init++)
    {
        Patch_Boundary_Flags(init) = 1;
    }

    if (myrank == 0)
    {
        std::cout << "Done with boundary patch flags init" << std::endl << std::flush;
    }
    // use set of nodal combinations to find boundary set
    // boundary patches will not try to add nodal combinations twice
    // loop through elements in this rank to find boundary patches
    npatches_repeat = 0;

    if (num_dim == 2)
    {
        for (int ielem = 0; ielem < rnum_elem; ielem++)
        {
            element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
            element_npatches = elem2D->nsurfaces;
            // loop through local surfaces
            for (int isurface = 0; isurface < element_npatches; isurface++)
            {
                num_nodes_in_patch = elem2D->surface_to_dof_lid.stride(isurface);
                Surface_Nodes = CArray<GO>(num_nodes_in_patch);
                for (int inode = 0; inode < num_nodes_in_patch; inode++)
                {
                    local_node_id = elem2D->surface_to_dof_lid(isurface, inode);
                    local_node_id = convert_node_order(local_node_id);
                    Surface_Nodes(inode) = nodes_in_elem(ielem, local_node_id);
                }
                Node_Combination temp(Surface_Nodes);
                // construct Node Combination object for this surface
                Patch_Nodes(npatches_repeat) = temp;
                Patch_Nodes(npatches_repeat).patch_id       = npatches_repeat;
                Patch_Nodes(npatches_repeat).element_id     = ielem;
                Patch_Nodes(npatches_repeat).local_patch_id = isurface;
                // test if this patch has already been added; if yes set boundary flags to 0
                current_combination = my_patches.insert(Patch_Nodes(npatches_repeat));
                // if the set determines this is a duplicate, access the original element's patch id and set flag to 0
                if (current_combination.second == false)
                {
                    // set original element flag to 0
                    Patch_Boundary_Flags((*current_combination.first).patch_id) = 0;
                    // set this current flag to 0 for the duplicate as well
                    Patch_Boundary_Flags(npatches_repeat) = 0;
                }
                npatches_repeat++;
            }
        }
    }

    if (num_dim == 3)
    {
        for (int ielem = 0; ielem < rnum_elem; ielem++)
        {
            element_select->choose_3Delem_type(Element_Types(ielem), elem);
            element_npatches = elem->nsurfaces;
            // loop through local surfaces
            for (int isurface = 0; isurface < element_npatches; isurface++)
            {
                num_nodes_in_patch = elem->surface_to_dof_lid.stride(isurface);
                // debug print
                // std::cout << "NUMBER OF PATCH NODES FOR ELEMENT " << ielem+1 << " ON LOCAL SURFACE " << isurface+1 << " IS " << num_nodes_in_patch << std::endl;
                Surface_Nodes = CArray<GO>(num_nodes_in_patch);
                for (int inode = 0; inode < num_nodes_in_patch; inode++)
                {
                    local_node_id = elem->surface_to_dof_lid(isurface, inode);
                    local_node_id = convert_node_order(local_node_id);
                    Surface_Nodes(inode) = nodes_in_elem(ielem, local_node_id);
                }
                Node_Combination temp(Surface_Nodes);
                // construct Node Combination object for this surface
                Patch_Nodes(npatches_repeat) = temp;
                Patch_Nodes(npatches_repeat).patch_id       = npatches_repeat;
                Patch_Nodes(npatches_repeat).element_id     = ielem;
                Patch_Nodes(npatches_repeat).local_patch_id = isurface;
                // test if this patch has already been added; if yes set boundary flags to 0
                current_combination = my_patches.insert(Patch_Nodes(npatches_repeat));
                // if the set determines this is a duplicate access the original element's patch id and set flag to 0

                if (current_combination.second == false)
                {
                    // set original element flag to 0
                    Patch_Boundary_Flags((*current_combination.first).patch_id) = 0;
                    // set this current flag to 0 for the duplicate as well
                    Patch_Boundary_Flags(npatches_repeat) = 0;
                }
                npatches_repeat++;
            }
        }
    }
    // debug print of all patches
    /*
    std::cout << " ALL PATCHES " << npatches_repeat <<std::endl;
    for(int iprint = 0; iprint < npatches_repeat; iprint++){
      std::cout << "Patch " << iprint + 1 << " ";
      for(int j = 0; j < Patch_Nodes(iprint).node_set.size(); j++)
        std::cout << Patch_Nodes(iprint).node_set(j) << " ";
      std::cout << std::endl;
    }
    */
    if (myrank == 0)
    {
        std::cout << "Done with boundary patch loop" << std::endl << std::flush;
    }
    // loop through patch boundary flags to isolate boundary patches
    nboundary_patches = 0;
    for (int iflags = 0 ; iflags < npatches_repeat; iflags++)
    {
        if (Patch_Boundary_Flags(iflags))
        {
            nboundary_patches++;
        }
    }

    // std::cout << " BOUNDARY PATCHES PRE COUNT ON TASK " << myrank << " = " << nboundary_patches <<std::endl;
    // upper bound that is not much larger
    Boundary_Patches  = CArrayKokkos<Node_Combination, array_layout, HostSpace, memory_traits>(nboundary_patches, "Boundary_Patches");
    nboundary_patches = 0;

    bool my_rank_flag;

    size_t remote_count;
    for (int ipatch = 0 ; ipatch < npatches_repeat; ipatch++)
    {
        if (Patch_Boundary_Flags(ipatch))
        {
            /*check if Nodes in the combination for this patch belong to this MPI rank.
              If all are local then this is a boundary patch belonging to this rank.
              If all nodes are remote then another rank must decide if that patch is a boundary.
              If only a subset of the nodes are local it must be a boundary patch; this
              case assigns the patch to the lowest mpi rank index the nodes in this patch belong to */
            num_nodes_in_patch = Patch_Nodes(ipatch).node_set.size();
            my_rank_flag = true;
            remote_count = 0;

            // assign as a local boundary patch if any of the nodes on the patch are local
            // only the local nodes on the patch will contribute to the equation assembly on this rank
            for (int inode = 0; inode < num_nodes_in_patch; inode++)
            {
                node_gid = Patch_Nodes(ipatch).node_set(inode);
                if (!map->isNodeGlobalElement(node_gid))
                {
                    // if(ghost_node_ranks.h_view(global2local_map.get(node_gid))<myrank)
                    // my_rank_flag = false;
                    remote_count++;
                    // test
                }
            }

            if (remote_count == num_nodes_in_patch)
            {
                my_rank_flag = false;
            }
            // all nodes were remote
            // if(remote_count == num_nodes_in_patch) my_rank_flag = false;

            // if at least one node was local
            if (my_rank_flag)
            {
                Boundary_Patches(nboundary_patches++) = Patch_Nodes(ipatch);
                boundary_patch_to_index[Patch_Nodes(ipatch)] = nboundary_patches - 1;
            }
        }
    }

    // debug print of boundary patches
    /*
    std::cout << " BOUNDARY PATCHES ON TASK " << myrank << " = " << nboundary_patches <<std::endl;
    for(int iprint = 0; iprint < nboundary_patches; iprint++){
      std::cout << "Patch " << iprint + 1 << " ";
      for(int j = 0; j < Boundary_Patches(iprint).node_set.size(); j++)
        std::cout << Boundary_Patches(iprint).node_set(j) << " ";
      std::cout << std::endl;
    }
    */
    // std::fflush(stdout);
}

/* ----------------------------------------------------------------------
  Parameters to feed ROL Trilinos package
------------------------------------------------------------------------- */

void Solver::set_rol_params(Teuchos::RCP<Teuchos::ParameterList> parlist)
{
    //set defaults here
    parlist->sublist("General").set("Variable Objective Function", false);
    parlist->sublist("General").set("Scale for Epsilon Active Sets", (double) 1.0);
    parlist->sublist("General").set("Output Level", (int) 1);
    parlist->sublist("General").set("Inexact Objective Function", false);
    parlist->sublist("General").set("Inexact Gradient", false);
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector", false);
    parlist->sublist("General").set("Projected Gradient Criticality Measure", false);

    parlist->sublist("General").sublist("Secant").set("Type", "Limited-Memory BFGS");
    parlist->sublist("General").sublist("Secant").set("Use as Preconditioner", false);
    parlist->sublist("General").sublist("Secant").set("Use as Hessian", false);
    parlist->sublist("General").sublist("Secant").set("Maximum Storage", (int) 5);
    parlist->sublist("General").sublist("Secant").set("Use Default Scaling", false);
    parlist->sublist("General").sublist("Secant").set("Initial Hessian Scale", (double) 1e-16);
    parlist->sublist("General").sublist("Secant").set("Barzilai-Borwein Type", (int) 1);

    parlist->sublist("General").sublist("Krylov").set("Type", "Conjugate Gradients");
    parlist->sublist("General").sublist("Krylov").set("Absolute Tolerance", (double) 1e-4);
    parlist->sublist("General").sublist("Krylov").set("Relative Tolerance", (double) 1e-2);
    parlist->sublist("General").sublist("Krylov").set("Iteration Limit", (int) 50);

    parlist->sublist("General").sublist("Polyhedral Projection").set("Type", "Dai-Fletcher");
    parlist->sublist("General").sublist("Polyhedral Projection").set("Iteration Limit", (int) 1000);
    parlist->sublist("General").sublist("Polyhedral Projection").set("Absolute Tolerance", (double) 1e-4);
    parlist->sublist("General").sublist("Polyhedral Projection").set("Relative Tolerance", (double) 1e-2);

    //Line search settings
    parlist->sublist("Step").sublist("Line Search").set("Function Evaluation Limit", (int) 20);
    parlist->sublist("Step").sublist("Line Search").set("Sufficient Decrease Tolerance", (double) 1.e-2);
    parlist->sublist("Step").sublist("Line Search").set("Initial Step Size", (double) 5e0);
    parlist->sublist("Step").sublist("Line Search").set("User Defined Initial Step Size", true);
    parlist->sublist("Step").sublist("Line Search").set("Normalize Initial Step Size", false);
    parlist->sublist("Step").sublist("Line Search").set("Accept Last Alpha", false);
    parlist->sublist("Step").sublist("Line Search").set("Use Previous Step Length as Initial Guess", false);
    parlist->sublist("Step").sublist("Line Search").set("Maximum Step Size", (double) 5e3);
    parlist->sublist("Step").sublist("Line Search").set("Use Adaptive Step Size Selection", true);

    parlist->sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", "Quasi-Newton Method");
    parlist->sublist("Step").sublist("Line Search").sublist("Descent Method").set("Nonlinear CG Type", "Hestenes-Stiefel");

    parlist->sublist("Step").sublist("Line Search").sublist("Curvature Condition").set("Type", "Strong Wolfe Conditions");
    parlist->sublist("Step").sublist("Line Search").sublist("Curvature Condition").set("General Parameter", (double) 0.9);
    parlist->sublist("Step").sublist("Line Search").sublist("Curvature Condition").set("Generalized Wolfe Parameter", (double) 0.6);

    parlist->sublist("Step").sublist("Line Search").sublist("Line-Search Method").set("Type", "Cubic Interpolation");
    parlist->sublist("Step").sublist("Line Search").sublist("Line-Search Method").set("Increase Rate", (double) 5e0);
    parlist->sublist("Step").sublist("Line Search").sublist("Line-Search Method").set("Backtracking Rate" , (double) 0.5);
    parlist->sublist("Step").sublist("Line Search").sublist("Line-Search Method").set("Bracketing Tolerance" , (double) 1e-8);

    parlist->sublist("Step").sublist("Line Search").sublist("Line-Search Method").sublist("Path-Based Target Level").set("Target Relaxation Parameter" , (double) 1.0);
    parlist->sublist("Step").sublist("Line Search").sublist("Line-Search Method").sublist("Path-Based Target Level").set("Upper Bound on Path Length" , (double) 1.0);

    //Trust region settings
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Solver", "Truncated CG");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model", "SPG");
    parlist->sublist("Step").sublist("Trust Region").set("Initial Radius", (double) 2e1);
    parlist->sublist("Step").sublist("Trust Region").set("Maximum Radius", (double) 5e8);
    parlist->sublist("Step").sublist("Trust Region").set("Step Acceptance Threshold", (double) 0.05);
    parlist->sublist("Step").sublist("Trust Region").set("Radius Shrinking Threshold", (double) 0.05);
    parlist->sublist("Step").sublist("Trust Region").set("Radius Growing Threshold", (double) 0.9);
    parlist->sublist("Step").sublist("Trust Region").set("Radius Shrinking Rate (Negative rho)", (double) 0.0625);
    parlist->sublist("Step").sublist("Trust Region").set("Radius Shrinking Rate (Positive rho)", (double) 0.25);
    parlist->sublist("Step").sublist("Trust Region").set("Radius Growing Rate", (double) 2.5);
    parlist->sublist("Step").sublist("Trust Region").set("Safeguard Size", (double) 1e1);

    //Trust region Lin-More
    parlist->sublist("Step").sublist("Trust Region").sublist("Lin-More").set("Maximum Number of Minor Iterations", (int) 10);
    parlist->sublist("Step").sublist("Trust Region").sublist("Lin-More").set("Sufficient Decrease Parameter", (double) 1e-2);
    parlist->sublist("Step").sublist("Trust Region").sublist("Lin-More").set("Relative Tolerance Exponent", (double) 1.1);

    parlist->sublist("Step").sublist("Trust Region").sublist("Lin-More").sublist("Cauchy Point").set("Maximum Number of Reduction Steps", (int) 10);
    parlist->sublist("Step").sublist("Trust Region").sublist("Lin-More").sublist("Cauchy Point").set("Maximum Number of Expansion Steps", (int) 10);
    parlist->sublist("Step").sublist("Trust Region").sublist("Lin-More").sublist("Cauchy Point").set("Initial Step Size", (double) 1.0);
    parlist->sublist("Step").sublist("Trust Region").sublist("Lin-More").sublist("Cauchy Point").set("Normalize Initial Step Size", true);
    parlist->sublist("Step").sublist("Trust Region").sublist("Lin-More").sublist("Cauchy Point").set("Reduction Rate", (double) 0.1);
    parlist->sublist("Step").sublist("Trust Region").sublist("Lin-More").sublist("Cauchy Point").set("Expansion Rate", (double) 5.0);
    parlist->sublist("Step").sublist("Trust Region").sublist("Lin-More").sublist("Cauchy Point").set("Decrease Tolerance", (double) 1e-8);

    parlist->sublist("Step").sublist("Trust Region").sublist("Lin-More").sublist("Projected Search").set("Backtracking Rate", (double) 0.5);
    parlist->sublist("Step").sublist("Trust Region").sublist("Lin-More").sublist("Projected Search").set("Maximum Number of Steps", (int) 20);

    //Trust region SPG (Spectral Projected Gradient)
    parlist->sublist("Step").sublist("Trust Region").sublist("SPG").set("Use Nonmonotone Trust Region", false);
    parlist->sublist("Step").sublist("Trust Region").sublist("SPG").set("Maximum Storage Size", (int) 10);
    parlist->sublist("Step").sublist("Trust Region").sublist("SPG").sublist("Solver").set("Iteration Limit", (int) 25);
    parlist->sublist("Step").sublist("Trust Region").sublist("SPG").sublist("Solver").set("Minimum Spectral Step Size", (double) 1e-12);
    parlist->sublist("Step").sublist("Trust Region").sublist("SPG").sublist("Solver").set("Maximum Spectral Step Size", (double) 1e12);
    parlist->sublist("Step").sublist("Trust Region").sublist("SPG").sublist("Solver").set("Use Smallest Model Iterate", false);

    //Controls for Inexactness
    parlist->sublist("Step").sublist("Trust Region").sublist("Inexact").sublist("Value").set("Tolerance Scaling", (double) 1e-1);
    parlist->sublist("Step").sublist("Trust Region").sublist("Inexact").sublist("Value").set("Exponent", (double) 0.9);
    parlist->sublist("Step").sublist("Trust Region").sublist("Inexact").sublist("Value").set("Forcing Sequence Initial Value", (double) 1.0);
    parlist->sublist("Step").sublist("Trust Region").sublist("Inexact").sublist("Value").set("Forcing Sequence Update Frequency", (int) 10);
    parlist->sublist("Step").sublist("Trust Region").sublist("Inexact").sublist("Value").set("Forcing Sequence Reduction Factor", (double) 0.1);

    parlist->sublist("Step").sublist("Trust Region").sublist("Inexact").sublist("Gradient").set("Tolerance Scaling", (double) 1e-1);
    parlist->sublist("Step").sublist("Trust Region").sublist("Inexact").sublist("Gradient").set("Relative Tolerance", (double) 2.0);

    //Spectral gradient options
    parlist->sublist("Step").sublist("Spectral Gradient").set("Minimum Spectral Step Size", (double) 1e-12);
    parlist->sublist("Step").sublist("Spectral Gradient").set("Maximum Spectral Step Size", (double) 1e12);

    //Primal dual active set options
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Dual Scaling", (double) 1.0);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Iteration Limit", (int) 10);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Step Tolerance", (double) 1e-8);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Gradient Tolerance", (double) 1e-6);
    
    parlist->sublist("Step").sublist("Composite Step").set("Output Level", (int) 0);
    parlist->sublist("Step").sublist("Composite Step").sublist("Optimality System Solver").set("Nominal Relative Tolerance", (double) 1e-10);
    parlist->sublist("Step").sublist("Composite Step").sublist("Optimality System Solver").set("Fix Tolerance", true);

    parlist->sublist("Step").sublist("Composite Step").sublist("Tangential Subproblem Solver").set("Iteration Limit", (int) 20);
    parlist->sublist("Step").sublist("Composite Step").sublist("Tangential Subproblem Solver").set("Relative Tolerance", (double) 1e-2);

    //Augmented Lagrangian Options
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Level of Hessian Approximation", (int) 0);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Use Default Problem Scaling", true);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Objective Scaling", (double) 1e0);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Constraint Scaling", (double) 1e0);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Use Default Initial Penalty Parameter", false);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Initial Penalty Parameter", simparam.optimization_options.rol_params.initial_constraint_penalty);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Penalty Parameter Growth Factor", (double) 1e1);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Minimum Penalty Parameter Reciprocal", (double) 0.1);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Initial Optimality Tolerance", (double) 1.0);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Optimality Tolerance Update Exponent", (double) 1.0);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Optimality Tolerance Decrease Exponent", (double) 1.0);
    
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Initial Feasibility Tolerance", (double) 1000.0);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Feasibility Tolerance Update Exponent", (double) 0.1);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Feasibility Tolerance Decrease Exponent", (double) 0.9);
    
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Print Intermediate Optimization History", false);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Step Type", simparam.optimization_options.rol_params.subproblem_algorithm_string);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Iteration Limit", simparam.optimization_options.rol_params.subproblem_iteration_limit);
    
    parlist->sublist("Step").sublist("Moreau-Yosida Penalty").set("Initial Penalty Parameter", (double) 1e-9);
    parlist->sublist("Step").sublist("Moreau-Yosida Penalty").set("Penalty Parameter Growth Factor", (double) 1.5);
    parlist->sublist("Step").sublist("Moreau-Yosida Penalty").sublist("Subproblem").set("Optimality Tolerance", (double) 1e-12);
    parlist->sublist("Step").sublist("Moreau-Yosida Penalty").sublist("Subproblem").set("Feasibility Tolerance", (double) 1e-12);
    parlist->sublist("Step").sublist("Moreau-Yosida Penalty").sublist("Subproblem").set("Print History", true);
    parlist->sublist("Step").sublist("Moreau-Yosida Penalty").sublist("Subproblem").set("Iteration Limit", (int) 200);
    
    parlist->sublist("Step").sublist("Bundle").set("Initial Trust-Region Parameter", (double) 1e1);
    parlist->sublist("Step").sublist("Bundle").set("Maximum Trust-Region Parameter", (double) 1e8);
    parlist->sublist("Step").sublist("Bundle").set("Tolerance for Trust-Region Parameter", (double) 1e-4);
    parlist->sublist("Step").sublist("Bundle").set("Epsilon Solution Tolerance", (double) 1e-8);
    parlist->sublist("Step").sublist("Bundle").set("Upper Threshold for Serious Step", (double) 1e-1);
    parlist->sublist("Step").sublist("Bundle").set("Lower Threshold for Serious Step", (double) 2e-1);
    parlist->sublist("Step").sublist("Bundle").set("Upper Threshold for Null Step", (double) 9e-1);
    parlist->sublist("Step").sublist("Bundle").set("Distance Measure Coefficient", (double) 1e-6);
    parlist->sublist("Step").sublist("Bundle").set("Maximum Bundle Size", (int) 50);
    parlist->sublist("Step").sublist("Bundle").set("Removal Size for Bundle Update", (int) 2);
    parlist->sublist("Step").sublist("Bundle").set("Cutting Plane Tolerance", (double) 1e-8);
    parlist->sublist("Step").sublist("Bundle").set("Cutting Plane Iteration Limit", (int) 1000);

    
    parlist->sublist("Status Test").set("Gradient Tolerance", simparam.optimization_options.rol_params.gradient_tolerance);
    parlist->sublist("Status Test").set("Constraint Tolerance", simparam.optimization_options.rol_params.constraint_tolerance);
    parlist->sublist("Status Test").set("Step Tolerance", simparam.optimization_options.rol_params.step_tolerance);
    parlist->sublist("Status Test").set("Iteration Limit", simparam.optimization_options.rol_params.iteration_limit);
    parlist->sublist("Status Test").set("Use Relative Tolerances", true);
}

/* ----------------------------------------------------------------------
  Setup Tpetra importers for comms
------------------------------------------------------------------------- */

void Solver::comm_importer_setup()
{
    // create import object using local node indices map and ghost indices map
    importer = Teuchos::rcp(new Tpetra::Import<LO, GO>(map, all_node_map));
    ghost_importer = Teuchos::rcp(new Tpetra::Import<LO, GO>(map, ghost_node_map));
    dof_importer   = Teuchos::rcp(new Tpetra::Import<LO, GO>(local_dof_map, all_dof_map));

    // output map and importers
    sorted_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(num_nodes, 0, comm));
    node_sorting_importer = Teuchos::rcp(new Tpetra::Import<LO, GO>(map, sorted_map));
    // sorted element mapping
    sorted_element_map = Teuchos::rcp(new Tpetra::Map<LO, GO, node_type>(num_elem, 0, comm));
    element_sorting_importer = Teuchos::rcp(new Tpetra::Import<LO, GO>(all_element_map, sorted_element_map));;
}

/* ----------------------------------------------------------------------
  Setup Tpetra exporters for reverse comms
------------------------------------------------------------------------- */

void Solver::comm_exporter_setup()
{
    // create import object using local node indices map and ghost indices map
    exporter = Teuchos::rcp(new Tpetra::Export<LO, GO>(all_node_map, map));
}

/* ----------------------------------------------------------------------
  Communicate updated nodal coordinates to ghost nodes
------------------------------------------------------------------------- */

void Solver::comm_coordinates()
{
    // create import object using local node indices map and ghost indices map
    // Tpetra::Import<LO, GO> importer(map, ghost_node_map);

    // comms to get ghosts
    ghost_node_coords_distributed->doImport(*node_coords_distributed, *ghost_importer, Tpetra::INSERT);
    // all_node_map->describe(*fos,Teuchos::VERB_EXTREME);
    // all_node_coords_distributed->describe(*fos,Teuchos::VERB_EXTREME);
}

/* ----------------------------------------------------------------------
   Return the CPU time for the current process in seconds very
   much in the same way as MPI_Wtime() returns the wall time.
------------------------------------------------------------------------- */

double Solver::CPU_Time()
{
    std::chrono::system_clock::time_point zero_time;

    auto zero_time_duration = zero_time.time_since_epoch();
    auto time = std::chrono::system_clock::now();
    auto time_duration = time.time_since_epoch();

    // double calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count();
    double calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(time_duration - zero_time_duration).count();
    calc_time *= 1e-09;

    return calc_time;
}

/* ----------------------------------------------------------------------
   Clock variable initialization
------------------------------------------------------------------------- */

void Solver::init_clock()
{
    double current_cpu = 0;
    initial_CPU_time = CPU_Time();
}
