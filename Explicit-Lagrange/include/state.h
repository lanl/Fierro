#ifndef STATE_H
#define STATE_H  

#include <iostream>
#include <string>

#include "utilities.h"
#include "geometry.h"

using namespace utils;


// should be node_state
class node_t {

private:
    int num_dim_;
    int num_rk_;

    // **** Node State **** //
    int num_nodes_;

    // Position
    real_t *node_coords_ = NULL; 

    // velocity 
    real_t *node_vel_ = NULL;     

    // velocity normal
    real_t *node_vel_norm_ = NULL; 

    // stress
    real_t *node_stress_ = NULL;
    
    // forces on the nodes
    real_t *node_force_ = NULL;

    // mass at nodes
    real_t *node_mass_ = NULL;

    real_t *corner_norm_sum_ = NULL;

public:

    void init_node_state (int num_dim, swage::mesh_t& mesh, int num_rk)
    {

        num_dim_   = num_dim;
        num_rk_    = num_rk;

        // **** Node State **** //
        num_nodes_ = mesh.num_nodes();

        node_coords_ = new real_t[num_rk_*num_nodes_*num_dim_]();
        node_vel_ = new real_t[num_rk_*num_nodes_*num_dim_]();
        node_vel_norm_ = new real_t[num_nodes_*num_dim_]();
        // node_stress_ = new real_t[num_rk_*num_nodes_*num_dim_*num_dim_]();
        node_force_ = new real_t[num_nodes_*num_dim_]();
        node_mass_ = new real_t[num_nodes_]();
        corner_norm_sum_  = new real_t[num_nodes_*num_dim_]();

    }


    // **** Node State **** //
    inline real_t& coords(int rk_stage, int node_gid, int this_dim) const
    {
        return node_coords_[rk_stage*num_nodes_*num_dim_ + node_gid*num_dim_ + this_dim];
    }

    inline real_t& norm_sum(int node_gid, int this_dim) const
    {
        return corner_norm_sum_[node_gid*num_dim_ + this_dim];
    }

    inline real_t& vel(int rk_stage, int node_gid, int this_dim) const
    {
        return node_vel_[rk_stage*num_nodes_*num_dim_ + node_gid*num_dim_ + this_dim];
    }

    inline real_t& vel_norm(int node_gid, int this_dim) const
    {
        return node_vel_norm_[node_gid*num_dim_ + this_dim];
    }

    // inline real_t& stress(int node_gid, int dim_i, int dim_j) const
    // {
    // return node_stress_[node_gid*num_dim_*num_dim_ + dim_i*num_dim_ + dim_j];
    // }
    
    inline real_t& force(int node_gid, int this_dim) const
    {
        return node_force_[node_gid*num_dim_ + this_dim];
    }

    inline real_t& mass(int node_gid) const
    {
        return node_mass_[node_gid];
    }



    // deconstructor
    ~node_t ( ) {

        delete[] node_coords_;
        delete[] node_vel_;
        delete[] node_vel_norm_;
        delete[] corner_norm_sum_;
        delete[] node_force_;
        delete[] node_mass_;
    }

};



class corner_t {

private:
    int num_dim_;
    int num_rk_;

    // **** Node State **** //
    int num_corners_;

    // Impeadance
    real_t *corn_impedance_ = NULL;

    // Density
    real_t *corn_density_ = NULL;

    // Sound speed
    real_t *corn_sspd_ = NULL;

    // velocity 
    real_t *corn_vel_ = NULL; 

    // velocity difference 
    real_t *corn_vel_diff_ = NULL; 
    
    // forces on the corner
    real_t *corn_force_ = NULL; 

    // power on the corner
    real_t *corn_power_ = NULL; 

    // Corner normal
    real_t *corn_normal_ = NULL; 

    // Corner normals
    real_t *corn_normals_ = NULL; 

    // stress at corner
    real_t *corn_stress_ = NULL;  

    // Pressure at the corner
    real_t *corn_press_ = NULL;

    // Shock direction
    real_t *corn_shock_dir_ = NULL;


    

public:

    void init_corner_state (int num_dim, swage::mesh_t& mesh, int num_rk)
    {

        num_dim_   = num_dim;
        num_rk_    = num_rk;

        // **** Node State **** //
        num_corners_ = mesh.num_corners();

        corn_impedance_ = new real_t[num_corners_]();
    
        corn_density_ = new real_t[num_corners_]();

        corn_press_ = new real_t[num_corners_]();
    
        corn_sspd_ = new real_t[num_corners_]();

        corn_vel_ = new real_t[num_corners_*num_dim_]();
    
        corn_vel_diff_ = new real_t[num_corners_]();
    
        corn_shock_dir_ = new real_t[num_corners_*num_dim_]();

        corn_force_ = new real_t[num_corners_*num_dim_]();
    
        corn_power_ = new real_t[num_corners_]();

        corn_normal_ = new real_t[num_corners_*num_dim_]();

        corn_normals_ = new real_t[num_corners_*num_dim_*num_dim_]();

        corn_stress_ = new real_t[num_corners_*num_dim_*num_dim_]();

    }


    // **** Corner State **** //

    inline real_t& impedance(int corner_gid) const
    {
        return corn_impedance_[corner_gid];
    }

    inline real_t& density(int corner_gid) const
    {
        return corn_density_[corner_gid];
    }

    inline real_t& sspd(int corner_gid) const
    {
        return corn_sspd_[corner_gid];
    }


    inline real_t& normal(int corner_gid, int this_dim) const
    {
        return corn_normal_[corner_gid*num_dim_ + this_dim];
    }

    inline real_t& normals(int corner_gid, int this_corn_patch, int this_dim) const
    {
        return corn_normals_[corner_gid*num_dim_*num_dim_ + this_corn_patch*num_dim_ + this_dim];
    }

    inline real_t& vel(int corner_gid, int this_dim) const
    {
        return corn_vel_[corner_gid*num_dim_ + this_dim];
    }

    inline real_t& vel_diff(int corner_gid) const
    {
        return corn_vel_diff_[corner_gid];
    }

    inline real_t& shock_dir(int corner_gid, int this_dim) const
    {
        return corn_shock_dir_[corner_gid*num_dim_ + this_dim];
    }
    
    inline real_t& force(int corner_gid, int this_dim) const
    {
        return corn_force_[corner_gid*num_dim_ + this_dim];
    }

    inline real_t& power(int corner_gid) const
    {
        return corn_power_[corner_gid];
    }

    inline real_t& stress(int corner_gid, int dim_i, int dim_j) const
    {
        return corn_stress_[corner_gid*num_dim_*num_dim_ + dim_i*num_dim_ + dim_j];
    }

    inline real_t& pressure(int corner_gid) const
    {
        return corn_press_[corner_gid];
    }

    // deconstructor
    ~corner_t ( ) {
        
        delete[] corn_impedance_;
        delete[] corn_density_;
        delete[] corn_press_;
        delete[] corn_sspd_;
        delete[] corn_vel_;
        delete[] corn_vel_diff_;
        delete[] corn_shock_dir_;
        delete[] corn_force_;
        delete[] corn_power_;
        delete[] corn_normal_;
        delete[] corn_normals_;
        delete[] corn_stress_;
    }

};





class cell_state_t {

    private:

        int num_dim_;
        int num_rk_;
        int num_nodes_;

        int *mat_id_ = NULL;        // the material region id  //****NOT IN MESH CLASSS******
        
        // **** Cell State **** //
        int num_cells_;

        real_t *mass_ = NULL;

        real_t *density_ = NULL;

        real_t *pressure_ = NULL;

        real_t *cs_ = NULL;

        real_t *power_ = NULL;

        real_t *divergence_ = NULL;

        real_t *ie_ = NULL;

        real_t *ke_ = NULL;

        real_t *total_energy_ = NULL;

        real_t *den_phi_ = NULL;
        real_t *vel_phi_ = NULL;
        real_t *te_phi_ = NULL;

        real_t *stress_ = NULL;

        real_t *force_ = NULL;

        // force related variables
        real_t *f_matrix_ = NULL; // corner forces in the cell
    
        // ---- the modal fields ---- //

        real_t *velocity_ = NULL;
        // real_t *total_energy_ = NULL;


        // -- Variables for SGH only -- //

        real_t *vel_grad_ = NULL;
        real_t *b_mat_ = NULL;



        

    public:
    
        void init_cell_state (int num_dim, swage::mesh_t& mesh, int num_rk)
        {

            num_dim_   = num_dim;
            num_rk_    = num_rk;

            // **** Element State **** //
            num_cells_ = mesh.num_cells();
            num_nodes_ = mesh.num_nodes_in_cell();

            mat_id_ = new int[num_cells_]();

            // vol and mass for the FV method, DG requires a vol_matrix, and mass_matrix
            
            mass_ = new real_t[num_cells_]();

            density_ = new real_t[num_cells_]();

            pressure_ = new real_t[num_cells_]();

            cs_ = new real_t[num_cells_]();

            power_ = new real_t [num_cells_]();

            divergence_ = new real_t [num_cells_]();
            
            ie_ = new real_t [num_rk_*num_cells_]();
            ke_ = new real_t [num_rk_*num_cells_]();
            total_energy_ = new real_t [num_rk_*num_cells_]();

            den_phi_ = new real_t [num_cells_]();
            vel_phi_ = new real_t [num_cells_]();
            te_phi_  = new real_t [num_cells_]();

            stress_ = new real_t [num_rk_*num_cells_*num_dim_*num_dim]();

            force_ = new real_t [num_cells_*num_dim_]();
        
            // corner forces
            f_matrix_ = new real_t[num_cells_*8*3]();
        
            // modal Taylor-series fields
            velocity_ = new real_t[num_rk_*num_cells_*3]();

            vel_grad_ = new real_t[num_cells_ *num_dim_*num_dim_]();
            b_mat_ = new real_t[num_cells_ * num_nodes_ * num_dim_]();

        }

        inline int& mat_id(int cell_gid) const
        {
            return mat_id_[cell_gid];
        }

        // **** Cell State **** //
    
        inline real_t& mass(int cell_gid) const
        {
            return mass_[cell_gid];
        }

        inline real_t& density(int cell_gid) const
        {
            return density_[cell_gid];
        }

        inline real_t& pressure(int cell_gid) const
        {
            return pressure_[cell_gid];
        }

        inline real_t& cs(int cell_gid) const
        {
            return cs_[cell_gid];
        }

        inline real_t& power(int cell_gid) const
        {
            return power_[cell_gid];
        }

        inline real_t& divergence(int cell_gid) const
        {
            return divergence_[cell_gid];
        }

        inline real_t& ie(int rk_stage, int cell_gid) const
        {
            return ie_[rk_stage*num_cells_ + cell_gid];
        }

        inline real_t& ke(int rk_stage, int cell_gid) const
        {
            return ke_[rk_stage*num_cells_ + cell_gid];
        }

        inline real_t& total_energy(int rk_stage, int cell_gid) const
        {
            return total_energy_[rk_stage*num_cells_ + cell_gid];
        }

        inline real_t& den_phi(int cell_gid) const
        {
            return den_phi_[cell_gid];
        }
    
        inline real_t& vel_phi(int cell_gid) const
        {
            return vel_phi_[cell_gid];
        }
    
        inline real_t& te_phi(int cell_gid) const
        {
            return te_phi_[cell_gid];
        }

        inline real_t& stress(int rk_stage, int cell_gid, int dim_i, int dim_j) const
        {
            return stress_[rk_stage*num_cells_*num_dim_*num_dim_ + cell_gid*num_dim_*num_dim_ + dim_i*num_dim_ + dim_j];
        }


        inline real_t& force(int cell_gid, int this_dim) const
        {
            // i*num_j*num_k + j*num_k + k
            return force_[cell_gid*num_dim_ + this_dim];
        }


    
        inline real_t& f_mat(int cell_gid, int node_lid, int this_dim) const
        {
            int idx = cell_gid*24 + node_lid*num_dim_ + this_dim;
            return f_matrix_[idx];
        } 

        inline real_t& velocity(int rk_stage, int cell_gid, int this_dim) const
        {
            // i*num_j*num_k + j*num_k + k
            return velocity_[rk_stage*num_cells_*num_dim_ + cell_gid*num_dim_ + this_dim];
        }


        inline real_t& vel_grad(int cell_gid, int dim_i, int dim_j) const
        {
            return vel_grad_[cell_gid*num_dim_*num_dim_ + dim_i*num_dim_ + dim_j];
        } 

        inline real_t& b_mat(int cell_gid, int node_lid, int dim) const
        {
            return b_mat_[cell_gid*num_nodes_*num_dim_ + node_lid*num_dim_ + dim];
        } 


        // deconstructor
        ~cell_state_t ( ) {

            delete[] mat_id_;
            delete[] mass_;
            delete[] density_;
            delete[] pressure_;

            delete[] cs_;
            delete[] power_;
            delete[] divergence_;
            delete[] ie_;
            delete[] ke_;
            delete[] total_energy_;
            delete[] force_;

            delete[] den_phi_;
            delete[] vel_phi_;
            delete[] te_phi_;

            delete[] stress_;

            delete[] f_matrix_;

            delete[] velocity_;

            delete[] vel_grad_;
            delete[] b_mat_;
        }
};



class elem_state_t {

    private:

        int num_dim_;
        int num_rk_;

        int num_basis_;

        // **** Element State **** //
        int num_elem_;
        int *bad_ = NULL;

        int *mat_id_;

        // Nodal fields
        real_t *inverse_mass_matrix_ = NULL;

        // modal polynomial fields
        real_t *velocity_ = NULL;
        real_t *specific_total_energy_ = NULL;
        real_t *specific_volume_ = NULL;
    
        // element averages
        real_t *avg_velocity_ = NULL;
        real_t *avg_specific_total_energy_ = NULL;
        real_t *avg_density_ = NULL;
        real_t *avg_specific_volume_ = NULL;


    public:
    
        void init_elem_state (int num_dim, swage::mesh_t& mesh, int num_rk, elements::ref_element& ref_elem)
        {

            num_dim_   = num_dim;
            num_rk_    = num_rk;

            // element state
            num_elem_ = mesh.num_elems();
            num_basis_ = ref_elem.num_basis();  // this num_basis is only needed for modal DG.

            mat_id_ = new int[num_elem_]();
            bad_ = new int[num_elem_]();

            inverse_mass_matrix_ = new real_t[num_elem_*num_basis_*num_basis_]();
        
            // modal polynomial (Taylor-series, Legendre, etc.) fields
            velocity_ = new real_t[num_rk_*num_elem_*num_basis_*num_dim_]();
            specific_total_energy_ = new real_t[num_rk_*num_elem_*num_basis_]();
            specific_volume_ = new real_t[num_elem_]();
        
            avg_velocity_ = new real_t[num_elem_*num_dim_]();
            avg_specific_total_energy_ = new real_t[num_elem_]();
            avg_density_ = new real_t[num_elem_]();
            avg_specific_volume_ = new real_t[num_elem_]();
        
        }

 
        // **** Element State **** //

        inline int& mat_id(int elem_gid) const
        {
            return mat_id_[elem_gid];
        }

        inline int& bad(int elem_gid) const
        {
            return bad_[elem_gid];
        }
    
        inline real_t& mass_mat_inv(int elem_gid, int basis_m, int basis_n) const
        {
            return inverse_mass_matrix_[elem_gid*num_basis_*num_basis_ + basis_m * num_basis_ + basis_n];
        }
    
        // were are the dims?
        inline real_t& velocity(int rk_stage, int elem_gid, int this_basis, int dim) const
        {
            return velocity_[rk_stage*num_elem_*num_basis_*num_dim_ + elem_gid*num_basis_*num_dim_ + this_basis*num_dim_ + dim];
        }

        inline real_t& specific_total_energy(int rk_stage, int elem_gid, int this_basis) const
        {
            return specific_total_energy_[rk_stage*num_elem_*num_basis_ + elem_gid*num_basis_ + this_basis];
        }

        inline real_t& specific_volume(int elem_gid, int this_basis) const
        {
            return specific_volume_[elem_gid*num_basis_ + this_basis];
        }
    
        inline real_t& avg_density(int elem_gid) const
        {
        return avg_density_[elem_gid];
        }

        inline real_t& avg_velocity(int elem_gid, int dim) const
        {
        return avg_velocity_[elem_gid*num_dim_ + dim];
        }
    
        inline real_t& avg_specific_total_energy(int elem_gid) const
        {
        return avg_specific_total_energy_[elem_gid];
        }
    
        inline real_t& avg_specific_volume(int elem_gid) const
        {
        return avg_specific_volume_[elem_gid];
        }

        // deconstructor
        ~elem_state_t ( ) {

            delete[] mat_id_;
            delete[] inverse_mass_matrix_;
            delete[] velocity_;
            delete[] specific_total_energy_;
            delete[] specific_volume_;
            delete[] avg_density_;
            delete[] avg_velocity_;
            delete[] avg_specific_total_energy_;
            delete[] avg_specific_volume_;
            delete[] bad_;

        }
};





class mat_pt_t {

private:

    int num_dim_;
    int num_rk_;

    // **** Material Point State **** //
    int num_matpt_;
    int *mat_id_ = NULL;        // the material region id  //****NOT IN MESH CLASSS******

    // physics state

    // Polynomial fields
    real_t *velocity_ = NULL; 
    real_t *specific_total_energy_ = NULL; 
    real_t *density_ = NULL; 
    real_t *specific_volume_ = NULL; 
    
	//limiter values
	real_t *energy_limiter_ = NULL;
	real_t *velocity_limiter_ = NULL;

    // Gauss point fields
    real_t *stress_=NULL;
    real_t *pressure_ = NULL; 
    real_t *sspd_ = NULL; 
    real_t *mass_ = NULL; 
    real_t *ie_ = NULL;
    real_t *ke_ = NULL;

    real_t *div_vel_ = NULL;
    real_t *grad_vel_ = NULL;  // symmetric part of velocity gradient

    real_t *q_visc_ = NULL;
    
    real_t *den_phi_ = NULL;
    real_t *vel_phi_ = NULL;
    real_t *te_phi_ = NULL;

    // Quadrature weights
    real_t *weights_ = NULL; 


public:

    void init_mat_pt_state (int num_dim, swage::mesh_t& mesh, int num_rk)
    {

        num_dim_   = num_dim;
        num_rk_    = num_rk;
        
        // **** Material Point State **** //
        num_matpt_ = mesh.num_gauss_pts();
        
        mat_id_ = new int[num_matpt_]();

        velocity_ = new real_t[num_rk_*num_matpt_*num_dim_]();
     
		//limiting values
		//for velocity_limiter, copy velocity_
		//for energy_limiter, copy specific_total_energy
   
		energy_limiter_ = new real_t[num_rk_*num_matpt_]();
		velocity_limiter_ = new real_t[num_rk_*num_matpt_]();

        specific_total_energy_ = new real_t[num_rk_*num_matpt_]();
        specific_volume_ = new real_t[num_rk_*num_matpt_]();
        density_ = new real_t[num_matpt_]();
        pressure_ = new real_t[num_matpt_]();
        stress_ = new real_t[num_matpt_*num_dim_*num_dim_]();
        sspd_ = new real_t[num_matpt_]();
        mass_ = new real_t[num_matpt_]();
        ie_ = new real_t[num_matpt_]();
        ke_ = new real_t[num_matpt_]();
    
        den_phi_ = new real_t[num_matpt_]();
        vel_phi_ = new real_t[num_matpt_]();
        te_phi_  = new real_t[num_matpt_]();

        div_vel_ = new real_t[num_matpt_]();
        grad_vel_ = new real_t[num_matpt_*num_dim_*num_dim_]();
    
        q_visc_ = new real_t[num_matpt_]();

        weights_ = new real_t[num_matpt_]();

    }

    inline int& mat_id(int mat_pt_gid) const
    {
        return mat_id_[mat_pt_gid];
    }


	//function for returning the velocity and specific total energy limiting values
	//velocity_limiter, _
	//energy_limiter, copy specific_total_energy
	
	inline real_t& velocity_limiter(int rk_stage, int mat_pt_gid) const
	{
		return velocity_limiter_[rk_stage*num_matpt_*num_dim_ + mat_pt_gid];
	}

	inline real_t& energy_limiter(int rk_stage, int mat_pt_gid) const
	{
		return energy_limiter_[rk_stage*num_matpt_ + mat_pt_gid];
	}

	//-------end of adding new functions; rest is og---------	

    inline real_t& velocity(int rk_stage, int mat_pt_gid, int this_dim) const
    {
        return velocity_[rk_stage*num_matpt_*num_dim_ + mat_pt_gid*num_dim_ + this_dim];
    }

    inline real_t& specific_total_energy(int rk_stage, int mat_pt_gid) const
    {
        return specific_total_energy_[rk_stage*num_matpt_+ mat_pt_gid];
    }

    inline real_t& specific_volume(int rk_stage, int mat_pt_gid) const
    {
        return specific_volume_[rk_stage*num_matpt_+ mat_pt_gid];
    }

    inline real_t& density(int mat_pt_gid) 
    {
        return density_[mat_pt_gid];
    }

    inline real_t& pressure(int mat_pt_gid) const
    {
        return pressure_[mat_pt_gid];
    }

    inline real_t& stress(int mat_pt_gid, int dim_i, int dim_j) const
    {
        return stress_[mat_pt_gid*num_dim_*num_dim_ + dim_i*num_dim_ + dim_j];
    }

    inline real_t& sspd(int mat_pt_gid) const
    {
        return sspd_[mat_pt_gid];
    }

    inline real_t& mass(int mat_pt_gid) const
    {
        return mass_[mat_pt_gid];
    }

    inline real_t& ie(int mat_pt_gid) const
    {
        return ie_[mat_pt_gid];
    }

    inline real_t& ke(int mat_pt_gid) const
    {
        return ke_[mat_pt_gid];
    }

    inline real_t& div_vel(int mat_pt_gid) const
    {
        return div_vel_[mat_pt_gid];
    }
    
    inline real_t& grad_vel(int mat_pt_gid, int dim_i, int dim_j) const
    {
        return grad_vel_[mat_pt_gid*num_dim_*num_dim_ + dim_i*num_dim_ + dim_j];
    }

    inline real_t& q_visc(int mat_pt_gid) const
    {
        return q_visc_[mat_pt_gid];
    }

    inline real_t& den_phi(int mat_pt_gid) const
    {
        return den_phi_[mat_pt_gid];
    }
    inline real_t& vel_phi(int mat_pt_gid) const
    {
        return vel_phi_[mat_pt_gid];
    }
    inline real_t& te_phi(int mat_pt_gid) const
    {
        return te_phi_[mat_pt_gid];
    }
    
    

    inline real_t& weight(int mat_pt_gid) const
    {
        return weights_[mat_pt_gid];
    }

    // deconstructor
    ~mat_pt_t ( ) {

        delete[] mat_id_;

        delete[] velocity_;
        delete[] specific_total_energy_;
        delete[] specific_volume_;
        delete[] density_;
        delete[] pressure_;
        delete[] stress_;
        delete[] sspd_;
        delete[] mass_;
        delete[] weights_;

        delete[] ie_;
        delete[] ke_;
        delete[] div_vel_;
        delete[] grad_vel_;
        
        delete[] den_phi_;
        delete[] vel_phi_;
        delete[] te_phi_;
        
        
        delete[] q_visc_;

    }
};













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
    int mat_id;
    
    // planes
    real_t x1;
    real_t x2;
    real_t y1;
    real_t y2;
    real_t z1;
    real_t z2;
    
    // radius
    real_t radius1;
    real_t radius2;

    
    // initial conditions
    init_conds::init_velocity_conds velocity;
    
    // velocity coefficients by component
    real_t u,v,w;
    
    // velocity magnitude for radial velocity initialization
    real_t speed;
    
    real_t ie; // internal energy
    real_t r;  // density
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


// eos variable to return from eos fcn
enum eos_return_t{
    p_of_de = 1,   // return pressure as a fcn of den and energy
    e_of_dp = 2,   // return energy as a fcn of den and pres
    cv = 3,        // return the specific heat
    sspd = 4,      // return sound speed
    t_of_cve = 5   // return temperature as a fcn of cv and energy
};

real_t ideal_gas(int kind, int k, real_t d, real_t e);



// CCH Code
void get_force_cch();
void update_position_cch(real_t rk_alpha);
void update_mom_energy_cch(real_t rk_alpha);
void setup_cch(char *MESH);
void update_state_cch();
void track_cch(real_t &x, real_t &y);
void cch_hydro();

// General code
void input();
void read_mesh_ensight(char* MESH);
void ensight();
void run_info(int cycle);

void boundary_force();
void boundary_velocity();
void setup_material();

void get_timestep();
void rk_init();



void test_corner_normals();


// SGH code
void update_energy_sgh(real_t rk_alpha);
void get_force_sgh();
void get_bmatrix();
void update_position_sgh(real_t rk_alpha);
void update_velocity_sgh(real_t rk_alpha);
void update_state_sgh();
void setup_sgh(char *MESH);
void track_sgh(real_t &x, real_t &y);
void get_velgrad();
void sgh_hydro();

// DG code
void dg_hydro();
void build_corner_normals();
void build_RHS();
void mass_mat_inverse();
void riemann();
void setup_dgh(char *MESH);
void strong_mass_dg();
void specific_vol_dg(real_t rk_alpha);
void gradvel_dg();
void momentum_dg(real_t rk_alpha,  int cycle);
void energy_dg(real_t rk_alpha, int cycle);
void limit_density(swage::mesh_t& mesh, elements::ref_element& ref_elem, std::string type, int elem_id );
void limit_vel(swage::mesh_t& mesh, elements::ref_element& ref_elem, std::string type, int elem_id );
void limit_energy(swage::mesh_t& mesh, elements::ref_element& ref_elem, std::string type, int elem_id);
void limit_specific_volume(swage::mesh_t& mesh, elements::ref_element& ref_elem, std::string type, int elem_id);
void track_dgh(real_t &x, real_t &y);
void calc_average_density();
void calc_average_velocity();
void calc_average_specific_total_energy();
void calc_average_specific_vol();
#endif 
