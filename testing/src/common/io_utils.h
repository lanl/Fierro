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

#include <map>
#include <memory>

// Reads in mesh connectivity data
class MeshReader
{
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
    /// \param Element state struct
    /// \param Node state struct
    /// \param Corner state struct
    /// \param Number of dimensions
    /// \param Number of RK bins
    ///
    ///
    /////////////////////////////////////////////////////////////////////////////
    void read_mesh(mesh_t& mesh, elem_t& elem, node_t& node, corner_t& corner, int num_dims, int rk_num_bins)
    {
        if (mesh_file_ == NULL) {
            printf("No mesh given\n");
            exit(0);
        }

        // Check mesh file extension
        // and read based on extension
        read_ensight_mesh(mesh, elem, node, corner, num_dims, rk_num_bins);
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
    void read_ensight_mesh(mesh_t& mesh, elem_t& elem, node_t& node, corner_t& corner, int num_dims, int rk_num_bins)
    {
        const size_t rk_level = 0;

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
                // printf("%c",ch);
            }
            // printf("\n");
        }

        // --- Read in the nodes in the mesh ---

        size_t num_nodes = 0;

        fscanf(in, "%lu", &num_nodes);
        printf("Num nodes read in %lu\n", num_nodes);

        // intialize node variables
        mesh.initialize_nodes(num_nodes);
        node.initialize(rk_num_bins, num_nodes, num_dims);

        // read the initial mesh coordinates
        // x-coords
        for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
            fscanf(in, "%le", &node.coords(rk_level, node_id, 0));
        }

        // y-coords
        for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
            fscanf(in, "%le", &node.coords(rk_level, node_id, 1));
        }

        // z-coords
        for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
            if (num_dims == 3) {
                fscanf(in, "%le", &node.coords(rk_level, node_id, 2));
            }
            else{
                double dummy;
                fscanf(in, "%le", &dummy);

                // printf("dummy = %le\n", dummy);
            }
        } // end for

        ch = (char)fgetc(in);
        // printf("%c",ch);

        // skip 1 line
        for (int j = 1; j <= 1; j++) {
            int i = 0;
            while ((ch = (char)fgetc(in)) != '\n') {
                i++;
                // printf("%c",ch);
            }
            // printf("\n");
        }

        // --- read in the elements in the mesh ---
        size_t num_elem = 0;

        fscanf(in, "%lu", &num_elem);
        printf("Num elements read in %lu\n", num_elem);

        // intialize elem variables
        mesh.initialize_elems(num_elem, num_dims);
        elem.initialize(rk_num_bins, num_nodes, 3); // always 3D here, even for 2D

        // for each cell read the list of associated nodes
        for (int elem_gid = 0; elem_gid < num_elem; elem_gid++) {
            for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                fscanf(in, "%lu", &mesh.nodes_in_elem.host(elem_gid, node_lid));  // %d vs zu

                // shift to start node index space at 0
                mesh.nodes_in_elem.host(elem_gid, node_lid) -= 1;
            }
        }
        // update device side
        mesh.nodes_in_elem.update_device();

        // intialize corner variables
        int num_corners = num_elem * mesh.num_nodes_in_elem;
        mesh.initialize_corners(num_corners);
        corner.initialize(num_corners, num_dims);

        // save the node coords to the current RK value
        for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
            for (int rk = 1; rk < rk_num_bins; rk++) {
                for (int dim = 0; dim < num_dims; dim++) {
                    node.coords(rk, node_gid, dim) = node.coords(0, node_gid, dim);
                } // end for dim
            } // end for rk
        } // end parallel for

        // Close mesh input file
        fclose(in);

        return;
    }
};

#endif // end Header Guard