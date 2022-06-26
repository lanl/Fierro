// -----------------------------------------------------------------------------
// This code reads the mesh in different formats
//------------------------------------------------------------------------------
#include "mesh.h"
#include "state.h"
#include "Implicit_Solver.h"

//MPI data
extern int myrank; //index of this mpi rank in the world communicator
extern int nranks; //number of mpi ranks in the world communicator
extern MPI_Comm world; //stores the default communicator object (MPI_COMM_WORLD)
extern Implicit_Solver implicit_solver_object; //current hack to get parallel file readin and spatial decomposition

// -----------------------------------------------------------------------------
// Reads an ensight .geo mesh file
//------------------------------------------------------------------------------
void read_mesh_ensight(char* MESH,
                       mesh_t &mesh,
                       node_t &node,
                       elem_t &elem,
                       corner_t &corner,
                       const size_t num_dims,
                       const size_t rk_num_bins){

    const size_t rk_level = 0;

	//FILE *in;
    //char ch;
    

    implicit_solver_object.nranks = nranks;
    implicit_solver_object.myrank = myrank;
    implicit_solver_object.world = world;

    implicit_solver_object.read_mesh_ensight(MESH);
    implicit_solver_object.init_maps();
    
    size_t num_nodes_in_elem = 1;
    for (int dim=0; dim<num_dims; dim++){
        num_nodes_in_elem *= 2;
    }


    // --- Read in the nodes in the mesh ---

    size_t num_nodes = implicit_solver_object.nlocal_nodes;
    
    printf("Num nodes assigned to MPI rank %lu is %lu\n" , myrank, num_nodes);
    
    // intialize node variables
    mesh.initialize_nodes(num_nodes);
    node.initialize(rk_num_bins, num_nodes, num_dims);
    std::cout << "Bin counts " << rk_num_bins << " Node counts " << num_nodes << " Num dim " << num_dims << std::endl;
    
    CArrayKokkos<double, DefaultLayout, HostSpace> host_node_coords_state(rk_num_bins, num_nodes, num_dims);
    Implicit_Solver::host_vec_array interface_node_coords = implicit_solver_object.node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    host_node_coords_state.get_kokkos_view() = node.coords.get_kokkos_dual_view().view_host();
    //host_node_coords_state = CArrayKokkos<double, DefaultLayout, HostSpace>(rk_num_bins, num_nodes, num_dims);
    host_node_coords_state.get_kokkos_view() = Kokkos::View<double*,DefaultLayout, HostSpace>("debug", rk_num_bins*num_nodes*num_dims);
    //save node data to node.coords
    for(int inode = 0; inode < num_nodes; inode++){
      host_node_coords_state(0,inode,0) = interface_node_coords(inode,0);
      host_node_coords_state(0,inode,1) = interface_node_coords(inode,1);
      host_node_coords_state(0,inode,2) = interface_node_coords(inode,2);
    }
    node.coords.update_device();


    // --- read in the elements in the mesh ---
    size_t num_elem = 0;
    
    num_elem = implicit_solver_object.rnum_elem;
    printf("Num elems assigned to MPI rank %lu is %lu\n" , myrank, num_elem);

    // intialize elem variables
    mesh.initialize_elems(num_elem, num_dims);
    elem.initialize(rk_num_bins, num_nodes, 3); // always 3D here, even for 2D

    //save data to mesh.nodes_in_elem.host
    CArrayKokkos<size_t, DefaultLayout, HostSpace> host_mesh_nodes_in_elem(num_elem, num_nodes_in_elem);
    Implicit_Solver::host_elem_conn_array interface_nodes_in_elem = implicit_solver_object.nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    host_mesh_nodes_in_elem.get_kokkos_view() = mesh.nodes_in_elem.get_kokkos_dual_view().view_host();
    //save node data to node.coords
    for(int ielem = 0; ielem < num_elem; ielem++){
        for(int inode = 0; inode < num_nodes_in_elem; inode++){
            host_mesh_nodes_in_elem(ielem,inode) = interface_nodes_in_elem(ielem,inode);
        }
    }

    // update device side
    mesh.nodes_in_elem.update_device();
    
    
    // intialize corner variables
    int num_corners = num_elem*mesh.num_nodes_in_elem;
    mesh.initialize_corners(num_corners);
    corner.initialize(num_corners, num_dims);

    
    // save the node coords to the current RK value
    for (size_t node_gid=0; node_gid<num_nodes; node_gid++){
        
        for(int rk=1; rk<rk_num_bins; rk++){
            for (int dim = 0; dim < num_dims; dim++){
                node.coords(rk, node_gid, dim) = node.coords(0, node_gid, dim);
            } // end for dim
        } // end for rk
        
    } // end parallel for
    
    size_t nall_nodes = implicit_solver_object.nall_nodes;
    node.all_coords = DCArrayKokkos <double> (rk_num_bins, nall_nodes, num_dims);
    node.all_vel    = DCArrayKokkos <double> (rk_num_bins, nall_nodes, num_dims);
    node.all_mass   = DCArrayKokkos <double> (nall_nodes);

    // Close mesh input file
    //fclose(in);

    return;
    
}
