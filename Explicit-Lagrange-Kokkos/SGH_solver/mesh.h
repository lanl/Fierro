#ifndef MESH_H
#define MESH_H


#include "matar.h"
#include "state.h"

#define PI 3.141592653589793

// mesh sizes and connectivity data structures
struct mesh_t {
    
    size_t num_dims;
    size_t num_nodes;
    size_t num_nodes_in_elem;
    size_t num_elems;
    size_t num_corners;

    
    // node ids in elem
    DCArrayKokkos <size_t> nodes_in_elem;
    
    // corner ids in elem
    CArrayKokkos <size_t> corners_in_elem;
    
    // corner ids in node
    RaggedRightArrayKokkos <size_t> corners_in_node;
    CArrayKokkos <size_t> num_corners_in_node;
    
    // initialization methods
    void initialize_nodes(const size_t num_nodes_inp)
    {
        this->num_nodes = num_nodes_inp;
    }; // end method
    
    
    // initialization methods
    void initialize_elems(const size_t num_elems_inp, const size_t num_dims_inp)
    {
        this->num_dims = num_dims_inp;
        this->num_nodes_in_elem = 1;
        for (int dim=0; dim<num_dims; dim++){
            num_nodes_in_elem *= 2;
        }
        this->num_elems = num_elems_inp;
        this->nodes_in_elem = DCArrayKokkos <size_t> (num_elems, num_nodes_in_elem);
        this->corners_in_elem = CArrayKokkos <size_t> (num_elems, num_nodes_in_elem);
    }; // end method
    
    
    // initialization methods
    void initialize_corners(const size_t num_corners_inp)
    {
        this->num_corners = num_corners_inp;
    }; // end method
    
    
    // build the corner mesh connectivity arrays
    void build_corner_connectivity(){
        
        this -> num_corners_in_node = CArrayKokkos <size_t> (num_nodes); // stride sizes
        
        // initializing the number of corners (node-cell pair) to be zero
        FOR_ALL(node_gid, 0, num_nodes, {
            num_corners_in_node(node_gid) = 0;
        });
        
        
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++){
            FOR_ALL(node_lid, 0, num_nodes_in_elem,{
                
                // get the global_id of the node
                size_t node_gid = nodes_in_elem(elem_gid, node_lid);

                // increment the number of corners attached to this point
                num_corners_in_node(node_gid) = num_corners_in_node(node_gid) + 1;
                
            });  // end FOR_ALL over nodes in element
        } // end for elem_gid
        
        
        // the stride sizes are the num_corners_in_node at the node
        this -> corners_in_node = RaggedRightArrayKokkos <size_t> (num_corners_in_node);

        CArrayKokkos <size_t> count_saved_corners_in_node(num_nodes);

        // reset num_corners to zero
        FOR_ALL(node_gid, 0, num_nodes, {
            count_saved_corners_in_node(node_gid) = 0;
        });
        
        
        // populate the cells connected to a node list and corners in a node
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++){
            FOR_ALL(node_lid, 0, num_nodes_in_elem, {
                
                // get the global_id of the node
                size_t node_gid = nodes_in_elem(elem_gid, node_lid);
                
                // the column index is the num corners saved
                size_t j = count_saved_corners_in_node(node_gid);

                // Save corner index to this node_gid
                size_t corner_gid = node_lid + elem_gid*num_nodes_in_elem;
                corners_in_node(node_gid, j) = corner_gid;

                // Save corner index to element
                size_t corner_lid = node_lid;
                corners_in_elem(elem_gid, corner_lid) = corner_gid;

                // increment the number of corners saved to this node_gid
                count_saved_corners_in_node(node_gid) = count_saved_corners_in_node(node_gid) + 1;

            });  // end FOR_ALL over nodes in element
        } // end for elem_gid


    } // end of build_corner_connectivity
  

}; // end mesh_t


namespace region
{

    // for tagging boundary faces
    enum vol_tag
    {
        global = 0,     // tag every cell in the mesh
        box = 1,        // tag all cells inside a box
        cylinder = 2,   // tag all cells inside a cylinder
        sphere = 3      // tag all cells inside a sphere
    };

} // end of namespace


namespace init_conds
{
    
    // applying initial conditions
    enum init_velocity_conds
    {
        // uniform
        cartesian = 0,   // cart velocity
        radial = 1,      // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
        spherical = 2,   // spherical
    
        // linear variation
        radial_linear = 3,     // linear variation from 0,0,0
        spherical_linear = 4,   // linear variation from 0,0,0
    
        // vortical initial conditions
        tg_vortex = 5
    };
    
} // end of initial conditions namespace


// fill instructions
struct mat_fill_t {
    
    // type
    region::vol_tag volume; // 1 is global, 2 are planes, 3 is a sphere
    
    // material id
    size_t mat_id;
    
    // planes
    double x1;
    double x2;
    double y1;
    double y2;
    double z1;
    double z2;
    
    // radius
    double radius1;
    double radius2;

    
    // initial conditions
    init_conds::init_velocity_conds velocity;
    
    // velocity coefficients by component
    double u,v,w;
    
    // velocity magnitude for radial velocity initialization
    double speed;
    
    double sie;  // specific internal energy
    double den;  // density
};


namespace bdy
{
    
    // for tagging boundary faces
    enum bdy_tag
    {
        x_plane  = 0,   // tag an x-plane
        y_plane  = 1,   // tag an y-plane
        z_plane  = 2,   // tag an z-plane
        cylinder = 3,   // tag an cylindrical surface
        sphere   = 4,   // tag a spherical surface
        readFile = 5    // read from a file
    };
    
    
    
    // for enforcing boundary conditions
    enum bdy_hydro_conds
    {
        fixed = 0,          // zero velocity
        reflected = 1,      // reflected or wall condition
        velocity = 2,       // constant velocity
        pressure = 3,       // constant pressure
        acceleration = 4,   // constant acceleration
        contact = 5         // contact surface
    };

} // end of bdy namespace


// tag mesh points on bdy's and set the BC type
struct boundary_t {

    // tag surface type
    bdy::bdy_tag surface;    // 0=xplane, 1=yplane, 2=zplane, 3=cylinder, 4=sphere, 5=read file
    
    // tag surface value or radius
    real_t value;
    
    // BC type
    bdy::bdy_hydro_conds hydro_bc;
    
};


void read_mesh_ensight(char* MESH,
                       mesh_t &mesh,
                       node_t &node,
                       elem_t &elem,
                       corner_t &corner,
                       size_t num_dims);


void input(CArrayKokkos <material_t> &material,
           CArrayKokkos <mat_fill_t> &mat_fill,
           CArrayKokkos <boundary_t> &boundary,
           CArrayKokkos <double> &state_vars);


KOKKOS_FUNCTION
void get_vol_hex(const DViewCArrayKokkos <double> &elem_vol,
                 const size_t elem_gid,
                 const DViewCArrayKokkos <double> &node_coords,
                 const mesh_t &mesh);


KOKKOS_FUNCTION
void get_bmatrix(const ViewCArrayKokkos <double> &B_matrix,
                 const size_t elem_gid,
                 const DViewCArrayKokkos <double> &node_coords,
                 const mesh_t &mesh);


void setup( const CArrayKokkos <material_t> &material,
            const CArrayKokkos <mat_fill_t> &mat_fill,
            const CArrayKokkos <boundary_t> &boundary,
            const mesh_t &mesh,
            const DViewCArrayKokkos <double> &node_coords,
            const DViewCArrayKokkos <double> &node_vel,
            const DViewCArrayKokkos <double> &node_mass,      
            const DViewCArrayKokkos <double> &elem_den,
            const DViewCArrayKokkos <double> &elem_pres,
            const DViewCArrayKokkos <double> &elem_stress,
            const DViewCArrayKokkos <double> &elem_sspd,       
            const DViewCArrayKokkos <double> &elem_sie,
            const DViewCArrayKokkos <double> &elem_vol,
            const DViewCArrayKokkos <double> &elem_mass,
            const DViewCArrayKokkos <size_t> &elem_mat_id,
            const DViewCArrayKokkos <double> &elem_statev,
            const CArrayKokkos <double> &state_vars);

void ensight( mesh_t &mesh,
              DViewCArrayKokkos <double> &node_coords,
              DViewCArrayKokkos <double> &node_vel,
              DViewCArrayKokkos <double> &node_mass,
              DViewCArrayKokkos <double> &elem_den,
              DViewCArrayKokkos <double> &elem_pres,
              DViewCArrayKokkos <double> &elem_stress,
              DViewCArrayKokkos <double> &elem_sspd, 
              DViewCArrayKokkos <double> &elem_sie,
              DViewCArrayKokkos <double> &elem_vol,
              DViewCArrayKokkos <double> &elem_mass,
              DViewCArrayKokkos <size_t> &elem_mat_id);
#endif 
