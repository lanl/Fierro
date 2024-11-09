// -----------------------------------------------------------------------------
// This code reads the mesh in different formats
//------------------------------------------------------------------------------
#include "mesh.h"
#include "state.h"
#include "rdh.h"


#include <fstream>

int PointIndexFromIJK(int i, int j, int k, const int* order);



// -----------------------------------------------------------------------------
// Reads an ensight .geo mesh file
//------------------------------------------------------------------------------
#if 0
void read_mesh_ensight(char* MESH,
                       mesh_t &mesh,
                       node_t &node,
                       elem_t &elem,
                       corner_t &corner,
                       const size_t num_dims,
                       const size_t rk_num_bins){

    const size_t rk_level = 0;

	FILE *in;
    char ch;
    
    
    size_t num_nodes_in_elem = 1;
    for (int dim=0; dim<num_dims; dim++){
        num_nodes_in_elem *= 2;
    }


    //read the mesh    WARNING: assumes a .geo file
    in = fopen(MESH,"r");  
    
    //skip 8 lines
    for (int j=1; j<=8;j++) {
        int i=0;
        while ((ch=(char)fgetc(in))!='\n') {
            i++;
            //printf("%c",ch);
        }
        //printf("\n");
    }  

    

    // --- Read in the nodes in the mesh ---
    
    size_t num_nodes = 0;
    
    fscanf(in,"%lu",&num_nodes);
    printf("Num nodes read in %lu\n" , num_nodes);
    
    // intialize node variables
    mesh.initialize_nodes(num_nodes);
    node.initialize(rk_num_bins, num_nodes, num_dims);
    

    // read the initial mesh coordinates
    // x-coords
    for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
        fscanf(in,"%le",&node.coords(rk_level,node_id, 0));
    }

    // y-coords
    for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
        fscanf(in,"%le",&node.coords(rk_level,node_id, 1));
    }  

    // z-coords
    for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
        if(num_dims==3){
            fscanf(in,"%le",&node.coords(rk_level,node_id, 2));
        } else
        {
            double dummy;
            fscanf(in,"%le",&dummy);
            
            //printf("dummy = %le\n", dummy);
        }
    } // end for

    
    //ch = (char)fgetc(in);
    ////printf("%c",ch);
//
    ////skip 1 line
    //for (int j=1; j<=1; j++) {
    //    int i=0;
    //    while ((ch=(char)fgetc(in))!='\n') {
    //        i++;
    //        //printf("%c",ch);
    //    }
    //    //printf("\n");
    //}

    std::string elem_str_name;
    char temp_name[100];
    fscanf(in,"%101s",&temp_name);
    
    elem_str_name = temp_name;
    printf("element kind = %s \n", elem_str_name.c_str());
    
    if (elem_str_name == "hexa8"){
        mesh.elem_kind = mesh_init::linear_tensor_element;
    }
    else if (elem_str_name == "quad4"){
        mesh.elem_kind = mesh_init::linear_tensor_element;
    }
    else {
        printf("\n ERROR: element kind is unknown \n");
    }
    //...
    

    // --- read in the elements in the mesh ---
    size_t num_elem = 0;
    
    fscanf(in,"%lu",&num_elem);
    printf("Num elements read in %lu\n" , num_elem);

    // intialize elem variables
    mesh.initialize_elems(num_elem, num_dims);
    elem.initialize(rk_num_bins, num_elem, 3); // always 3D here, even for 2D

    // for each cell read the list of associated nodes
    for (int elem_gid = 0; elem_gid < num_elem; elem_gid++) {
        for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
            
            fscanf(in,"%lu",&mesh.nodes_in_elem.host(elem_gid, node_lid));  // %d vs zu

            // shift to start node index space at 0
            mesh.nodes_in_elem.host(elem_gid, node_lid) -= 1;
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




    // Close mesh input file
    fclose(in);

    return;
    
}

#endif

// This function reads a VTKPn file
//--------------------------------------------------------
//
void readVTKPn(char* MESH,
                 mesh_t &mesh,
                 node_t &node,
                 elem_t &elem,
                 zone_t &zone,
                 mat_pt_t &mat_pt,
                 corner_t &corner,
                 fe_ref_elem_t &ref_elem,
                 const size_t num_dims,
                 const size_t rk_num_bins)
{

    
    printf("READING arbitrary-order element mesh \n");
    
    const size_t rk_level = 0;
    
    size_t i;  // a counter on searching a mesh file
    
    size_t num_nodes;
    size_t num_elems;
    size_t num_nodes_in_elem;
    size_t num_zones_in_elem;
    size_t num_surfs_in_elem;

    std::string token;
    
    bool found = false;
    
    std::ifstream in;  // FILE *in;
    in.open(MESH);
    

    // look for nodes
    i = 0;
    while (found==false) {
        std::string str;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split (str, delimiter);
        
        // looking for the following text:
        //      POINTS %d float
        if(v[0] == "POINTS"){
            num_nodes = std::stoi(v[1]);
            printf("Num nodes read in %zu\n", num_nodes);
            
            found=true;
        } // end if
        
        
        if (i>1000){
            printf("ERROR: Failed to find POINTS \n");
            break;
        } // end if
        
        i++;
    } // end while
    
    // intialize node variables
    mesh.initialize_nodes(num_nodes);
    // printf("here-node initialize\n");
    node.initialize(rk_num_bins, num_nodes, num_dims);
    // printf("here-node finish initialize\n");
    
    
    // read the point coordinates
    for (size_t node_gid=0; node_gid<num_nodes; node_gid++){
        
        
        std::string str;
        std::getline(in, str);
        
        std::string delimiter = " ";
        std::vector<std::string> v = split (str, delimiter);
        
        for (size_t dim=0; dim<3; dim++){
            node.coords(rk_level, node_gid, dim) = std::stod(v[dim]);
        }
        
    } // end for points
    found=false;
    
    
    // look for CELLS
    i = 0;
    while (found==false) {
        std::string str;
        std::getline(in, str);
        
        std::string delimiter = " ";
        std::vector<std::string> v = split (str, delimiter);
        
        // looking for the following text:
        //      CELLS num_elems size
        if(v[0] == "CELLS"){
            num_elems = std::stoi(v[1]);
            printf("Num elements read in %zu\n", num_elems);
            
            found=true;
        } // end if
        
        
        if (i>1000){
            printf("ERROR: Failed to find CELLS \n");
            break;
        } // end if
        
        i++;
    } // end while
    
    
    // next line has the number of nodes in an element
    {
        std::streampos oldpos = in.tellg();  // stores the position
        
        std::string str;
        std::getline(in, str);
        
        std::string delimiter = " ";
        std::vector<std::string> v = split (str, delimiter);
        num_nodes_in_elem = std::stoi(v[0]); // setting the num_nodes_in_elem
        
        // go back to prior line
        in.seekg(oldpos);
        
    }
    
    size_t p = int(std::cbrt(num_nodes_in_elem)-1);
    //num_zones_in_elem = int(std::pow(p, 3)); // one order lower than nodal index space
    num_surfs_in_elem = 2*num_dims; // 4 (2D) or 6 (3D)
    
    // initialize reference element //
    //printf("polynomial order = %d \n", p);
    //printf("num_zones_in_elem = %d \n", num_zones_in_elem);
    //printf("num_surfs_in_elem = %d \n", num_surfs_in_elem);
    
    ref_elem.init(p, num_dims);
    
    size_t num_gauss_leg_in_elem = ref_elem.num_gauss_leg_in_elem;
    num_zones_in_elem = ref_elem.num_elem_basis;
    //printf("num thermo basis functions = %zu \n", ref_elem.num_elem_basis);
    size_t num_nodes_in_zone = 8;
    // intialize elem mesh
    mesh.initialize_elems_Pn(num_elems, num_nodes_in_elem, num_gauss_leg_in_elem, num_zones_in_elem, num_nodes_in_zone, num_surfs_in_elem, num_dims);
    elem.initialize_Pn(rk_num_bins, num_elems, num_nodes_in_elem, num_zones_in_elem, num_dims); 
    zone.initialize_Pn(rk_num_bins, num_elems, num_zones_in_elem);
    mat_pt.initialize_Pn(rk_num_bins, num_elems, num_nodes_in_elem, num_zones_in_elem, num_dims, p);
    
    // intialize corner variables
    size_t num_corners = num_elems*num_nodes_in_elem;
    mesh.initialize_corners(num_corners);
    corner.initialize(num_corners, num_dims);
    
    
    // read the point ids in the element
    
    // first step is to build a node ordering map
    const int num_1D_points = std::cbrt(num_nodes_in_elem);  // cube root
    const int Pn_order = num_1D_points - 1;
    
    mesh.Pn = Pn_order;
    
    printf("Pn_order = %d \n", Pn_order);
    
    CArray <size_t> get_ijk_from_vtk(num_nodes_in_elem, 3);
    mesh.convert_vtk_to_fierro = CArray <size_t> (num_nodes_in_elem);
    mesh.convert_fierro_to_vtk = CArray <size_t> (num_nodes_in_elem);
    
    // p_order   = 1, 2, 3, 4, 5
    // num_nodes = 2, 3, 4, 5, 6
    const int order[3] = {Pn_order, Pn_order, Pn_order};
    
    // re-order the nodes to be in i,j,k format of Fierro
    size_t this_node_lid = 0;
    for (size_t k=0; k<num_1D_points; k++){
        for (size_t j=0; j<num_1D_points; j++){
            for (size_t i=0; i<num_1D_points; i++){
                
                // convert this_point index to the FE index convention
                size_t vtk_index = PointIndexFromIJK(i, j, k, order);
                
                // store the points in this elem according the the finite
                // element numbering convention
                mesh.convert_vtk_to_fierro(vtk_index) = this_node_lid;
                mesh.convert_fierro_to_vtk(this_node_lid) = vtk_index;
                
                get_ijk_from_vtk(vtk_index, 0) = i;
                get_ijk_from_vtk(vtk_index, 1) = j;
                get_ijk_from_vtk(vtk_index, 2) = k;
                
                // increment the point counting index
                this_node_lid ++;
                
            } // end for icount
        } // end for jcount
    }  // end for kcount
    

    for (size_t elem_gid=0; elem_gid<num_elems; elem_gid++) {
        
        std::string str;
        std::getline(in, str);
        
        
        std::string delimiter = " ";
        std::vector<std::string> v = split (str, delimiter);
        for (size_t node_lid=0; node_lid<num_nodes_in_elem; node_lid++){
            
            // covert the local indexing to the ijk convention used by Fierro
            int vtk_index = mesh.convert_fierro_to_vtk(node_lid);
            
            // must add 1 because 1st value is num_nodes_in_elem
            mesh.nodes_in_elem.host(elem_gid, node_lid) = std::stoi(v[vtk_index+1]);
            
            mesh.convert_vtk_to_fierro(vtk_index) = node_lid;
        }
        
        
    } // end for
    

    
    
    // update device side
    mesh.nodes_in_elem.update_device();
    
    found=false;

    
    // look for CELL_TYPE
    i = 0;
    size_t elem_type = 0;
    while (found==false) {
        std::string str;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split (str, delimiter);
        
        // looking for the following text:
        //      CELLS num_elems size
        if(v[0] == "CELL_TYPES"){

            std::getline(in, str);
            elem_type = std::stoi(str);
            printf("elem type found is %d \n", elem_type);
            
            found=true;
        } // end if
        
        
        if (i>1000){
            printf("ERROR: Failed to find CELL_TYPE \n");
            break;
        } // end if
        
        i++;
    } // end while
    //printf("elem type = %d \n", elem_type);
    // elem type:
    // VTK_LAGRANGE_HEXAHEDRON: 72,
    if (elem_type == 72){
        printf("abitrary-order tensor product element \n");
        mesh.elem_kind = mesh_init::arbitrary_tensor_element;
    }
    else {
        printf("\n ERROR: element kind is unknown \n");
    }
    
    printf("done reading vtk mesh \n");
    
    found=false;
    
    in.close();
    
    
    // save the node coords to the current RK value
    for (size_t node_gid=0; node_gid<num_nodes; node_gid++){
        
        for(int rk=1; rk<rk_num_bins; rk++){
            for (int dim = 0; dim < num_dims; dim++){
                node.coords(rk, node_gid, dim) = node.coords(0, node_gid, dim);
            } // end for dim
        } // end for rk
        
    } // end parallel for
    
}



/**\brief Given (i,j,k) coordinates within the Lagrange hex, return an offset into the local connectivity (PointIds) array.
  *
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

  if (nbdy == 3) // Vertex DOF
    { // ijk is a corner node. Return the proper index (somewhere in [0,7]):
    return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
    }

  int offset = 8;
  if (nbdy == 2) // Edge DOF
    {
    if (!ibdy)
      { // On i axis
      return (i - 1) +
        (j ? order[0] - 1 + order[1] - 1 : 0) +
        (k ? 2 * (order[0] - 1 + order[1] - 1) : 0) +
        offset;
      }
    if (!jbdy)
      { // On j axis
      return (j - 1) +
        (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) +
        (k ? 2 * (order[0] - 1 + order[1] - 1) : 0) +
        offset;
      }
    // !kbdy, On k axis
    offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
    return (k - 1) + (order[2] - 1) * (i ? (j ? 3 : 1) : (j ? 2 : 0)) + offset;
    }

  offset += 4 * (order[0] - 1 + order[1] - 1 + order[2] - 1);
  if (nbdy == 1) // Face DOF
    {
    if (ibdy) // On i-normal face
      {
      return (j - 1) + ((order[1] - 1) * (k - 1)) + (i ? (order[1] - 1) * (order[2] - 1) : 0) + offset;
      }
    offset += 2 * (order[1] - 1) * (order[2] - 1);
    if (jbdy) // On j-normal face
      {
      return (i - 1) + ((order[0] - 1) * (k - 1)) + (j ? (order[2] - 1) * (order[0] - 1) : 0) + offset;
      }
    offset += 2 * (order[2] - 1) * (order[0] - 1);
    // kbdy, On k-normal face
    return (i - 1) + ((order[0] - 1) * (j - 1)) + (k ? (order[0] - 1) * (order[1] - 1) : 0) + offset;
    }

  // nbdy == 0: Body DOF
  offset += 2 * (
    (order[1] - 1) * (order[2] - 1) +
    (order[2] - 1) * (order[0] - 1) +
    (order[0] - 1) * (order[1] - 1));
  return offset +
    (i - 1) + (order[0] - 1) * (
      (j - 1) + (order[1] - 1) * (
        (k - 1)));
}





