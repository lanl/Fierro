// -----------------------------------------------------------------------------
// This code reads the mesh in different formats
//------------------------------------------------------------------------------
#include "mesh.h"
#include "state.h"
#include "Explicit_Solver_SGH.h"

// -----------------------------------------------------------------------------
// Reads an ensight .geo mesh file
//------------------------------------------------------------------------------
void sgh_interface_setup(Explicit_Solver_SGH *explicit_solver_pointer,
                       mesh_t &mesh,
                       node_t &node,
                       elem_t &elem,
                       corner_t &corner,
                       const size_t num_dims,
                       const size_t rk_num_bins){

    const size_t rk_level = 0;
    
    size_t num_nodes_in_elem = 1;
    for (int dim=0; dim<num_dims; dim++){
        num_nodes_in_elem *= 2;
    }

    // --- Read in the nodes in the mesh ---

    size_t num_nodes = explicit_solver_pointer->nall_nodes;
    int myrank = explicit_solver_pointer->myrank;
    int nranks = explicit_solver_pointer->nranks;
    //printf("Num nodes assigned to MPI rank %lu is %lu\n" , myrank, num_nodes);

    // intialize node variables
    mesh.initialize_nodes(num_nodes);
    mesh.initialize_local_nodes(explicit_solver_pointer->nlocal_nodes);
    node.initialize(rk_num_bins, num_nodes, num_dims);
    //std::cout << "Bin counts " << rk_num_bins << " Node counts " << num_nodes << " Num dim " << num_dims << std::endl;

    //view scope
    {
      Explicit_Solver_SGH::host_vec_array interface_node_coords = explicit_solver_pointer->all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
      //save node data to node.coords
      //std::cout << "NODE DATA ON RANK " << myrank << std::endl;
      for(int inode = 0; inode < num_nodes; inode++){
        //std::cout << "Node index " << inode+1 << " ";
        node.coords.host(0,inode,0) = interface_node_coords(inode,0);
        //std::cout << host_node_coords_state(0,inode,0)+1<< " ";
        node.coords.host(0,inode,1) = interface_node_coords(inode,1);
        //std::cout << host_node_coords_state(0,inode,1)+1<< " ";
        node.coords.host(0,inode,2) = interface_node_coords(inode,2);
        //std::cout << host_node_coords_state(0,inode,2)+1<< std::endl;
      }
    } //end view scope
    // --- read in the elements in the mesh ---
    size_t num_elem = 0;
    
    num_elem = explicit_solver_pointer->rnum_elem;
    //printf("Num elems assigned to MPI rank %lu is %lu\n" , myrank, num_elem);

    // intialize elem variables
    mesh.initialize_elems(num_elem, num_dims);
    elem.initialize(rk_num_bins, num_nodes, 3); // always 3D here, even for 2D

    //save data to mesh.nodes_in_elem.host
    //CArrayKokkos<size_t, DefaultLayout, HostSpace> host_mesh_nodes_in_elem(num_elem, num_nodes_in_elem);
    //view scope
    {
      Explicit_Solver_SGH::host_elem_conn_array interface_nodes_in_elem = explicit_solver_pointer->nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
      //save node data to node.coords
      //std::cout << "ELEMENT CONNECTIVITY ON RANK " << myrank << std::endl;
      for(int ielem = 0; ielem < num_elem; ielem++){
        //std::cout << "Element index " << ielem+1 << " ";
        for(int inode = 0; inode < num_nodes_in_elem; inode++){
            mesh.nodes_in_elem.host(ielem,inode) = explicit_solver_pointer->all_node_map->getLocalElement(interface_nodes_in_elem(ielem,inode));
            //debug print
            //std::cout << mesh.nodes_in_elem.get_kokkos_dual_view().h_view(ielem*num_nodes_in_elem + inode)+1<< " ";
        }
        //std::cout << std::endl;
      }
    }
    // update device side
    mesh.nodes_in_elem.update_device();

    //debug print
    
    //CArrayKokkos<size_t> device_mesh_nodes_in_elem(num_elem, num_nodes_in_elem);
    //device_mesh_nodes_in_elem.get_kokkos_view() = mesh.nodes_in_elem.get_kokkos_dual_view().d_view;
    //host_mesh_nodes_in_elem.get_kokkos_view() = mesh.nodes_in_elem.get_kokkos_dual_view().view_host();
    /*
    if(myrank==1){
    std::cout << "ELEMENT CONNECTIVITY ON RANK 1 in LOCAL INDICES" << myrank << std::endl;
    for(int ielem = 0; ielem < num_elem; ielem++){
        std::cout << "Element index " << ielem+1 << " ";
        for(int inode = 0; inode < num_nodes_in_elem; inode++){
            //debug print
            //device_mesh_nodes_in_elem(ielem,inode) = explicit_solver_pointer->all_node_map->getLocalElement(interface_nodes_in_elem(ielem,inode));
            std::cout << mesh.nodes_in_elem(ielem, inode)+1<< " ";
        }
        std::cout << std::endl;
    }
    }
    */
    /*
    std::cout.flush();
    if(myrank==1){
    std::cout << "ELEMENT CONNECTIVITY ON RANK 1 in GLOBAL INDICES" << myrank << std::endl;
    std::cout << "local node index of global index 275 on rank 1 " << explicit_solver_pointer->all_node_map->getLocalElement(275) << std::endl;
    for(int ielem = 0; ielem < num_elem; ielem++){
        std::cout << ielem << " ";
        for(int inode = 0; inode < num_nodes_in_elem; inode++){
            //debug print
            //device_mesh_nodes_in_elem(ielem,inode) = explicit_solver_pointer->all_node_map->getLocalElement(interface_nodes_in_elem(ielem,inode));
            std::cout << explicit_solver_pointer->all_node_map->getGlobalElement(mesh.nodes_in_elem(ielem, inode))<< " ";
        }
        std::cout << std::endl;
    }
    }
    std::cout.flush();
    */
    /*
    size_t nall_nodes = explicit_solver_pointer->nall_nodes;
    node.all_coords = DCArrayKokkos <double> (rk_num_bins, nall_nodes, num_dims);
    node.all_vel    = DCArrayKokkos <double> (rk_num_bins, nall_nodes, num_dims);
    node.all_mass   = DCArrayKokkos <double> (nall_nodes);

    //save all data (nlocal +nghost)
    CArrayKokkos<double, DefaultLayout, HostSpace> host_all_node_coords_state(rk_num_bins, nall_nodes, num_dims);
    Explicit_Solver_SGH::host_vec_array interface_all_node_coords = explicit_solver_pointer->all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    host_all_node_coords_state.get_kokkos_view() = node.all_coords.get_kokkos_dual_view().view_host();
    //host_node_coords_state = CArrayKokkos<double, DefaultLayout, HostSpace>(rk_num_bins, num_nodes, num_dims);
    //host_all_node_coords_state.get_kokkos_view() = Kokkos::View<double*,DefaultLayout, HostSpace>("debug", rk_num_bins*nall_nodes*num_dims);
    //save node data to node.coords
    
    //std::cout << "ALL NODE DATA ON RANK " << myrank << std::endl;
    for(int inode = 0; inode < nall_nodes; inode++){
        //std::cout << "Node index " << inode+1 << " ";
        node.all_coords.host(0,inode,0) = interface_all_node_coords(inode,0);
        //std::cout << host_all_node_coords_state(0,inode,0)+1<< " ";
        node.all_coords.host(0,inode,1) = interface_all_node_coords(inode,1);
        //std::cout << host_all_node_coords_state(0,inode,1)+1<< " ";
        node.all_coords.host(0,inode,2) = interface_all_node_coords(inode,2);
        //std::cout << host_all_node_coords_state(0,inode,2)+1<< std::endl;
    }
    */

    // save the node coords to the current RK value
    for (size_t node_gid=0; node_gid<num_nodes; node_gid++){
        
        for(int rk=1; rk<rk_num_bins; rk++){
            for (int dim = 0; dim < num_dims; dim++){
                node.coords.host(rk, node_gid, dim) = node.coords.host(0, node_gid, dim);
            } // end for dim
        } // end for rk
        
    } // end parallel for
    
    /*
    // save the node coords to the current RK value
    for (size_t node_gid=0; node_gid<nall_nodes; node_gid++){
        
        for(int rk=1; rk<rk_num_bins; rk++){
            for (int dim = 0; dim < num_dims; dim++){
                node.all_coords.host(rk, node_gid, dim) = node.all_coords.host(0, node_gid, dim);
            } // end for dim
        } // end for rk
        
    } // end parallel for
    */
    
    node.coords.update_device();
    //node.all_coords.update_device();

    
    // intialize corner variables
    int num_corners = num_elem*mesh.num_nodes_in_elem;
    mesh.initialize_corners(num_corners);
    corner.initialize(num_corners, num_dims);
    
    /*
    for(int inode = 0; inode < num_nodes; inode++){
        std::cout << "Node index " << inode+1 << " ";
        for(int rk=0; rk<rk_num_bins; rk++){
          std::cout << "rk index " << rk+1 << " ";
          std::cout << node.coords(rk,inode,0)+1<< " ";
          std::cout << node.coords(rk,inode,1)+1<< " ";
          std::cout << node.coords(rk,inode,2)+1<< std::endl;
        }
    }
    */
    // Close mesh input file
    //fclose(in);

    return;
    
}
