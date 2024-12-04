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

#include <math.h>  // fmin, fmax, abs note: fminl is long
#include "FEA_Module.h"
#include "Solver.h"
#include "Simulation_Parameters/Simulation_Parameters.h"

#define BC_EPSILON 1.0e-8
using namespace utils;

FEA_Module::FEA_Module(Solver* Solver_Pointer)
{
    Solver_Pointer_ = Solver_Pointer;
    simparam = &Solver_Pointer->simparam;

    num_dim = simparam->num_dims;
    num_gauss_points = simparam->num_gauss_points;

    // obtain global and local node and element counts
    num_nodes    = Solver_Pointer->num_nodes;
    num_elem     = Solver_Pointer->num_elem;
    nall_nodes   = Solver_Pointer->nall_nodes;
    rnum_elem    = Solver_Pointer->rnum_elem;
    nlocal_nodes = Solver_Pointer->nlocal_nodes;
    nghost_nodes = Solver_Pointer->nghost_nodes;
    max_nodes_per_element = Solver_Pointer->max_nodes_per_element;

    hessvec_count     = update_count = 0;
    linear_solve_time = hessvec_time = hessvec_linear_time = 0;
    last_compute_step = -1;

    Matrix_alloc = 0;
    gradient_print_sync = 0;
    // RCP pointer to *this (Parallel Nonlinear Solver Object)
    // FEM_pass = Teuchos::rcp(this);

    // Trilinos output stream
    std::ostream& out = std::cout;
    fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    (*fos).setOutputToRootOnly(0);

    // MPI Data copy
    myrank   = Solver_Pointer->myrank;
    nranks   = Solver_Pointer->nranks;
    world    = Solver_Pointer->world;
    importer = Solver_Pointer->importer;
    ghost_importer = Solver_Pointer->ghost_importer;
    node_sorting_importer    = Solver_Pointer->node_sorting_importer;
    element_sorting_importer = Solver_Pointer->element_sorting_importer;
    dof_importer = Solver_Pointer->dof_importer;

    // obtain node and element maps
    comm = Solver_Pointer->comm;
    map  = Solver_Pointer->map; // map of node indices
    ghost_node_map  = Solver_Pointer->ghost_node_map; // map of node indices with ghosts on each rank
    all_node_map    = Solver_Pointer->all_node_map; // map of node indices with ghosts on each rank
    element_map     = Solver_Pointer->element_map; // non overlapping map of elements owned by each rank used in reduction ops
    all_element_map = Solver_Pointer->all_element_map; // overlapping map of elements connected to the local nodes in each rank
    local_dof_map   = Solver_Pointer->local_dof_map; // map of local dofs (typically num_node_local*num_dim)
    all_dof_map     = Solver_Pointer->all_dof_map; // map of local and ghost dofs (typically num_node_all*num_dim)

    // obtain mesh coordinates, densities, and element connectivity
    global_nodes_in_elem_distributed = Solver_Pointer->global_nodes_in_elem_distributed; // element to node connectivity table
    node_nconn_distributed                  = Solver_Pointer->node_nconn_distributed; // how many elements a node is connected to
    node_coords_distributed                 = Solver_Pointer->node_coords_distributed;
    all_node_coords_distributed             = Solver_Pointer->all_node_coords_distributed;
    initial_node_coords_distributed         = Solver_Pointer->initial_node_coords_distributed;
    all_initial_node_coords_distributed     = Solver_Pointer->all_initial_node_coords_distributed;
    design_node_densities_distributed       = Solver_Pointer->design_node_densities_distributed;
    design_node_coords_distributed          = Solver_Pointer->design_node_coords_distributed;
    filtered_node_densities_distributed     = Solver_Pointer->filtered_node_densities_distributed;
    test_node_densities_distributed         = Solver_Pointer->test_node_densities_distributed;
    all_node_densities_distributed          = Solver_Pointer->all_node_densities_distributed;
    all_design_node_coords_distributed      = Solver_Pointer->all_design_node_coords_distributed;
    all_filtered_node_densities_distributed = Solver_Pointer->all_filtered_node_densities_distributed;
    Global_Element_Densities                = Solver_Pointer->Global_Element_Densities;
    Element_Types = Solver_Pointer->Element_Types;

    // element select data
    element_select = Solver_Pointer->element_select;
    element_select->choose_3Delem_type(Element_Types(0), elem);

    // obtain boundary condition and loading data
    nboundary_patches  = Solver_Pointer->nboundary_patches;
    Boundary_Patches   = Solver_Pointer->Boundary_Patches;
    node_specified_bcs = false;
    // initialize for default
    num_boundary_conditions = 0;
    bcs_initialized = false;

    // flag init
    body_term_flag = nonzero_bc_flag = false;

    // output data
    noutput = 0;
}

FEA_Module::~FEA_Module()
{
}

/* ----------------------------------------------------------------------
   find which boundary patches correspond to the given BC.
   bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
   val = plane value, cylinder radius, shell radius
------------------------------------------------------------------------- */

void FEA_Module::tag_boundaries(int bc_tag, real_t val, int bdy_set, real_t* patch_limits)
{
    int is_on_set;
    /*
    if (bdy_set == num_bdy_sets_){
      std::cout << " ERROR: number of boundary sets must be increased by "
        << bdy_set-num_bdy_sets_+1 << std::endl;
      exit(0);
    }
    */

    // test patch limits for feasibility
    if (patch_limits != NULL)
    {
        // test for upper bounds being greater than lower bounds
        if (patch_limits[1] <= patch_limits[0])
        {
            std::cout << " Warning: patch limits for boundary condition are infeasible " << patch_limits[0] << " and " << patch_limits[1] << std::endl;
        }
        if (patch_limits[3] <= patch_limits[2])
        {
            std::cout << " Warning: patch limits for boundary condition are infeasible " << patch_limits[2] << " and " << patch_limits[3] << std::endl;
        }
    }

    // save the boundary vertices to this set that are on the plane
    int counter = 0;
    for (int iboundary_patch = 0; iboundary_patch < nboundary_patches; iboundary_patch++)
    {
        // check to see if this patch is on the specified plane
        is_on_set = check_boundary(Boundary_Patches(iboundary_patch), bc_tag, val, patch_limits); // no=0, yes=1

        if (is_on_set == 1)
        {
            Boundary_Condition_Patches(bdy_set, counter) = iboundary_patch;
            counter++;
        }
    } // end for bdy_patch

    // save the number of bdy patches in the set
    NBoundary_Condition_Patches(bdy_set) = counter;

    *fos << " tagged boundary patches " << std::endl;
}

/* -------------------------------------------------------------------------------------------
   Communicate ghosts using the current optimization design data
---------------------------------------------------------------------------------------------- */

void FEA_Module::comm_densities(Teuchos::RCP<const MV> zp)
{
    // set density vector to the current value chosen by the optimizer
    test_node_densities_distributed = zp;

    // debug print of design vector
    // std::ostream &out = std::cout;
    // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    // if(myrank==0)
    // *fos << "Density data :" << std::endl;
    // node_densities_distributed->describe(*fos,Teuchos::VERB_EXTREME);
    // *fos << std::endl;
    // std::fflush(stdout);

    // communicate design densities
    // create import object using local node indices map and all indices map
    Tpetra::Import<LO, GO> importer(map, all_node_map);

    // comms to get ghosts
    all_node_densities_distributed->doImport(*test_node_densities_distributed, importer, Tpetra::INSERT);

    // update_count++;
    // if(update_count==1){
    // MPI_Barrier(world);
    // MPI_Abort(world,4);
    // }
}

void FEA_Module::comm_filtered_densities()
{
    // debug print of design vector
    // std::ostream &out = std::cout;
    // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    // if(myrank==0)
    // *fos << "Density data :" << std::endl;
    // node_densities_distributed->describe(*fos,Teuchos::VERB_EXTREME);
    // *fos << std::endl;
    // std::fflush(stdout);

    // communicate design densities
    // create import object using local node indices map and all indices map
    Tpetra::Import<LO, GO> importer(map, all_node_map);

    // comms to get ghosts
    all_filtered_node_densities_distributed->doImport(*filtered_node_densities_distributed, importer, Tpetra::INSERT);

    // update_count++;
    // if(update_count==1){
    // MPI_Barrier(world);
    // MPI_Abort(world,4);
    // }
}

/* ----------------------------------------------------------------------
   routine for checking to see if a patch is on a boundary set
   bc_tag = 0 xplane, 1 yplane, 3 zplane, 4 cylinder, 5 is shell
   val = plane value, radius, radius
------------------------------------------------------------------------- */

int FEA_Module::check_boundary(Node_Combination& Patch_Nodes, int bc_tag, real_t val, real_t* patch_limits)
{
    int is_on_set = 1;
    const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);

    // Nodes on the Patch
    auto   node_list = Patch_Nodes.node_set;
    int    num_dim   = simparam->num_dims;
    size_t nnodes    = node_list.size();
    size_t node_rid;
    real_t node_coord[num_dim];
    int    dim_other1, dim_other2;
    CArrayKokkos<size_t, array_layout, device_type, memory_traits> node_on_flags(nnodes, "node_on_flags");

    // initialize
    for (int inode = 0; inode < nnodes; inode++)
    {
        node_on_flags(inode) = 0;
    }

    if (bc_tag == 0)
    {
        dim_other1 = 1;
        dim_other2 = 2;
    }
    else if (bc_tag == 1)
    {
        dim_other1 = 0;
        dim_other2 = 2;
    }
    else if (bc_tag == 2)
    {
        dim_other1 = 0;
        dim_other2 = 1;
    }

    // test for planes
    if (bc_tag < 3)
    {
        for (int inode = 0; inode < nnodes; inode++)
        {
            node_rid = all_node_map->getLocalElement(node_list(inode));
            for (int init = 0; init < num_dim; init++)
            {
                node_coord[init] = all_node_coords(node_rid, init);
            }
            if (fabs(node_coord[bc_tag] - val) <= BC_EPSILON)
            {
                node_on_flags(inode) = 1;

                // test if within patch segment if user specified
                if (patch_limits != NULL)
                {
                    if (node_coord[dim_other1] - patch_limits[0] <= -BC_EPSILON)
                    {
                        node_on_flags(inode) = 0;
                    }
                    if (node_coord[dim_other1] - patch_limits[1] >= BC_EPSILON)
                    {
                        node_on_flags(inode) = 0;
                    }
                    if (node_coord[dim_other2] - patch_limits[2] <= -BC_EPSILON)
                    {
                        node_on_flags(inode) = 0;
                    }
                    if (node_coord[dim_other2] - patch_limits[3] >= BC_EPSILON)
                    {
                        node_on_flags(inode) = 0;
                    }
                }
            }
            // debug print of node id and node coord
            // std::cout << "node coords on task " << myrank << " for node " << node_rid << std::endl;
            // std::cout << "coord " <<node_coord << " flag " << node_on_flags(inode) << " bc_tag " << bc_tag << std::endl;
        }
    }

    /*
    // cylinderical shell where radius = sqrt(x^2 + y^2)
    else if (this_bc_tag == 3){

        real_t R = sqrt(these_patch_coords[0]*these_patch_coords[0] +
                        these_patch_coords[1]*these_patch_coords[1]);

        if ( fabs(R - val) <= 1.0e-8 ) is_on_bdy = 1;


    }// end if on type

    // spherical shell where radius = sqrt(x^2 + y^2 + z^2)
    else if (this_bc_tag == 4){

        real_t R = sqrt(these_patch_coords[0]*these_patch_coords[0] +
                        these_patch_coords[1]*these_patch_coords[1] +
                        these_patch_coords[2]*these_patch_coords[2]);

        if ( fabs(R - val) <= 1.0e-8 ) is_on_bdy = 1;

    } // end if on type
    */
    // check if all nodes lie on the boundary set
    for (int inode = 0; inode < nnodes; inode++)
    {
        if (!node_on_flags(inode))
        {
            is_on_set = 0;
        }
    }

    // debug print of return flag
    // std::cout << "patch flag on task " << myrank << " is " << is_on_set << std::endl;
    return is_on_set;
} // end method to check bdy