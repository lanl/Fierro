#ifndef ELEMENTS_H
#define ELEMENTS_H 
/*****************************************************************************
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utilities.h"
#include "matar.h"

using namespace utils;




namespace elements{

//==============================================================================
//   Function Declaration
//==============================================================================

    // Used by Lobatto 1D/2D to set Lobatto quadrature points
    void lobatto_nodes_1D(
        CArray <real_t> &lob_nodes_1D,
        const int &num);

    void lobatto_weights_1D(
        CArray <real_t> &lob_weights_1D,  // Labbatto weights
        const int &num);

    void length_weights(
        CArray <real_t> &len_weights_1D,  // Labbatto weights
        CArray <real_t> &lab_weights_1D,  // Labbatto weights
        CArray <real_t> &lab_nodes_1D,
        const int &order);

    void sub_weights(
        CArray <real_t> &sub_weights_1D,  // Labbatto weights
        CArray <real_t> &lab_weights_1D,  // Labbatto weights
        CArray <real_t> &lab_nodes_1D,
        const int &order);

    void mat_inverse(
        ViewCArray <real_t> &mat_inv,
        ViewCArray <real_t> &matrix);

    // void mat_mult(
    //     CArray <real_t> &result,
    //     CArray <real_t> &A,
    //     CArray <real_t> &B);

    // void mat_trans(
    //     CArray <real_t> &trans,
    //     CArray <real_t> &mat);

    void set_nodes_wgts(
        CArray <real_t> &lab_nodes_1D,
        CArray <real_t> &lab_weights_1D,
        CArray <real_t> &len_weights_1D,
        CArray <real_t> &sub_weights_1D, 
        int p_order);

    // void sub_cells(
    //     CArray <real_t> &lab_nodes_1D,
    //     int &p_order, 
    //     int &dim);

    // Create the unit normals for the corner patches in reference space
    void set_unit_normals(CArray <real_t> &unit_normals);


    // Used by Gauss2/3D to set quadrature points
    void line_gauss_info(
        real_t &x, 
        real_t &w, 
        int &m,  
        int &p);

    // Used by Lovatto 1D/2D to set Lobatto quadrature points
    void line_lobatto_info(
        real_t &x, 
        real_t &w, 
        int &m, 
        int &p);

    // setting gauss quadrature points for 2D elements
    void gauss_2d(
        ViewCArray <real_t> &these_g_pts,     // gauss points
        ViewCArray <real_t> &these_weights,   // gauss weights
        ViewCArray <real_t> &tot_g_weight,    // 2D product of gauss weights
        int &quad_order);                     // quadrature order (n)

    // setting gauss quadrature points for 2D elements
    void gauss_3d(
        ViewCArray <real_t> &these_g_pts,   // gauss points
        ViewCArray <real_t> &these_weights, // gauss weights
        ViewCArray <real_t> &tot_g_weight,  // 3D product of gauss weights
        int &quad_order);                     // quadrature order (n)

    // setting gauss quadrature points for 4D elements
    void gauss_4d(
        ViewCArray <real_t> &these_g_pts,     // gauss points
        ViewCArray <real_t> &these_weights,   // gauss weights
        int &quad_order,                        // quadrature order (n)
        const int &dim);

    // setting Gauss-Lobatto quadrature points for 2D elements
    void lobatto_2d(
        ViewCArray <real_t> &these_L_pts,     // gauss points
        ViewCArray <real_t> &these_weights,   // gauss weights
        int &quad_order);                       // quadrature order (n)

    // setting Gauss-Lobatto quadrature points for 3D elements
    void lobatto_3d(
        ViewCArray <real_t> &these_L_pts,     // gauss points
        ViewCArray <real_t> &these_weights,   // gauss weights
        int &quad_order); 

    // setting gauss quadrature points for 4D elements
    void lobatto_4d(
        ViewCArray <real_t> &these_L_pts,     // gauss points
        ViewCArray <real_t> &these_weights,   // gauss weights
        int &quad_order,                        // quadrature order (n)
        const int &dim);

    //defining the jacobian for 2D elements
    void jacobian_2d(
        ViewCArray <real_t> &J_matrix, 
        real_t &det_J,
        const ViewCArray <real_t> &vertices, 
        const ViewCArray <real_t> &this_partial,
        const int &num_nodes);

    //defining the jacobian for 3D elements
    void jacobian_3d(
        ViewCArray <real_t> &J_matrix, 
        real_t &det_J,
        const ViewCArray <real_t> &vertices, 
        const ViewCArray <real_t> &this_partial,
        const int &num_nodes);

    //defining the jacobian for 4D elements
    void jacobian_4d(
        ViewCArray <real_t> &J_matrix, 
        real_t &det_J,
        const ViewCArray <real_t> &vertices, 
        const ViewCArray <real_t> &this_partial,
        const int &num_nodes,
        const int &dim);

    //defining the inverse jacobian for 2D element    
    void jacobian_inverse_2d(
        ViewCArray <real_t> &J_inverse, 
        const ViewCArray <real_t> &jacobian);

    //defining the inverse jacobian for 2D element    
    void jacobian_inverse_3d(
        ViewCArray <real_t> &J_inverse_matrix,
        const ViewCArray <real_t> &jacobian);

    // defining  the inverse of the Jacobian for 4D elements
    void jacobian_inverse_4d(
        ViewCArray <real_t> &J_inverse_matrix,
        const ViewCArray <real_t> &jacobian,
        const real_t &det_J);

    // creates nodal positions with Chebyshev spacing
    void chebyshev_nodes_1D(
        ViewCArray <real_t> &cheb_nodes_1D,   // Chebyshev nodes
        const int &order);                      // Interpolation order


// Reference Element Informations


class ref_element{
    private:
        
        int num_dim_;
        
        int num_ref_nodes_1D_;
        int num_ref_cells_1D_;
        int num_ref_corners_1D_;
    
        int num_ref_surface_nodes_in_elem_;
        int num_ref_inside_nodes_in_elem_;

        // Zones
        int num_zones_1d_;
        int num_zones_in_elem_;
        
        // cells
        int num_ref_cells_in_elem_;

        CArray <int>cells_in_zone_list_;
        
        // nodes
        int num_ref_nodes_in_elem_;
        int num_ref_nodes_in_cell_;

        CArray <int> ref_nodes_in_cell_;
        CArray <int> ref_surface_nodes_in_elem_;
        CArray <int> ref_inside_nodes_in_elem_;

        CArray <real_t> ref_node_positions_;
        CArray <real_t> ref_node_g_weights_;

        CArray <int> cell_nodes_in_elem_list_;

        // Vertices
        int num_ref_verts_1d_;
        int num_ref_verts_in_elem_;

        // corners
        int num_ref_corners_in_cell_;
        int num_ref_corners_in_elem_;

        CArray <int> ref_corners_in_cell_;

        CArray <real_t> ref_corner_surf_normals_;
        CArray <real_t> ref_corner_g_weights_;
        CArray <real_t> ref_corner_surf_g_weights_;

        // Num basis functions
        int num_basis_;

        // Basis evaluation at nodes
        CArray <real_t> ref_nodal_basis_;
    
    public:
        // DANIELLOOK
        // Gradient of basis
        CArray <real_t> ref_nodal_gradient_;
        //real_t * ref_nodal_gradient_;
    
        // Function Declarations

        // Initialize reference element information
        void init(int poly_order, int num_dim);

        int num_dim() const;

        int num_basis() const;

        int num_ref_nodes() const;

        int num_ref_cells_in_elem() const;
        int num_ref_corners_in_cell() const;
        
        int node_rid(int i, int j, int k) const;
        int cell_rid(int i, int j, int k) const;
        int corner_rid(int i, int j, int k) const;
        
        int ref_corners_in_cell(int cell_rid, int corner_rlid) const;
        int ref_nodes_in_cell(int cell_rid, int node_rlid) const;
    
        int ref_surface_nodes_in_elem(int node_rlid) const;
        int ref_inside_nodes_in_elem(int node_rlid) const;
        int num_ref_surface_nodes_in_elem() const;
        int num_ref_inside_nodes_in_elem() const;
    
        real_t ref_node_positions(int node_rid, int dim) const;

        real_t ref_corner_surface_normals(int corner_rid, int surf_rlid, int dim) const;
        
        real_t ref_corner_g_surface_weights(int corner_rid, int surf_rlid) const;
        
        real_t ref_node_g_weights(int node_rid) const;

        real_t ref_corner_g_weights(int corner_rid) const;

        real_t &ref_nodal_gradient(int node_rid, int basis_id, int dim) const;

        real_t &ref_nodal_basis(int node_rid, int basis_id) const;
        // Nodal jacobian and determinant should be in mesh class, I think....
        //real_t ref_nodal_jacobian(int node_rid, int basis_id, int dim) const;
        
        int& cell_lid_in_zone(int zone_lid, int cell_lid) const;

        int& cell_nodes_in_elem(int cell_lid, int node_lid) const;

        int vert_node_map(int vert_lid);

        // Deconstructor
        ~ref_element();

    };


    class Element2D {
        
        protected:
            const static int num_dim_ = 2;

        public:
        
        virtual int num_verts() = 0;
        virtual int num_nodes() = 0;
        virtual int num_basis() = 0;

        // calculate a physical position in an element for a given xi,eta
        virtual void physical_position(
            ViewCArray <real_t>  &x_point,
            const ViewCArray <real_t>  &xi_point,
            const ViewCArray <real_t> &vertices) = 0;

        // calculate the value for the basis at each node for a given xi,eta
        virtual void basis(
            ViewCArray <real_t>  &basis,
            const ViewCArray <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Xi
        virtual void  partial_xi_basis(
            ViewCArray <real_t>  &partial_xi, 
            const ViewCArray <real_t>  &xi_point) = 0;


        // Partial derivative of shape functions with respect to Xi
        virtual void  partial_eta_basis(
            ViewCArray <real_t> &partial_eta, 
            const ViewCArray <real_t>  &xi_point) = 0;

        // Map from vertex to node
        virtual int vert_node_map( const int vert_lid) = 0;


    }; // end of 2D element class

    class Element3D {
        
        protected:
            const static int num_dim_ = 3;

        public:
        
        virtual int num_verts() = 0;
        virtual int num_nodes() = 0;
        virtual int num_basis() = 0;

        // calculate a physical position in an element for a given xi,eta,mu
        virtual void physical_position(
            ViewCArray <real_t>  &x_point,
            const ViewCArray <real_t>  &xi_point,
            const ViewCArray <real_t> &vertices) = 0;

        // calculate the value for the basis at each node for a given xi,eta, mu
        virtual void basis(
            ViewCArray <real_t>  &basis,
            const ViewCArray <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Xi at Xi_point
        virtual void partial_xi_basis(
            ViewCArray <real_t>  &partial_xi, 
            const ViewCArray <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Eta
        virtual void partial_eta_basis(
            ViewCArray <real_t> &partial_eta, 
            const ViewCArray <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Mu
        virtual void partial_mu_basis(
            ViewCArray <real_t> &partial_mu, 
            const ViewCArray <real_t>  &xi_point) = 0;

        // Map from vertex to node
        virtual inline int vert_node_map( const int vert_lid) = 0;

        // Reference vertices location
        virtual real_t& ref_locs(const int vert_lid, const int dim) = 0;

    }; // end of 3D parent class

    class Element4D {

        public:

        // calculate a physical position in an element for a given xi,eta,mu,tau
        virtual void physical_position(
            ViewCArray <real_t>  &x_point,
            const ViewCArray <real_t>  &xi_point,
            const ViewCArray <real_t> &vertices) = 0;

        // calculate the value for the basis at each node for a given xi,eta,mu,tau
        virtual void basis(
            ViewCArray <real_t>  &basis,
            const ViewCArray <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Xi at Xi_point
        virtual void partial_xi_basis(
            ViewCArray <real_t>  &partial_xi, 
            const ViewCArray <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Eta
        virtual void partial_eta_basis(
            ViewCArray <real_t> &partial_eta, 
            const ViewCArray <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Mu
        virtual void partial_mu_basis(
            ViewCArray <real_t> &partial_mu, 
            const ViewCArray <real_t>  &xi_point) = 0;

        // Partial derivative of shape functions with respect to Tau
        virtual void partial_tau_basis(
            ViewCArray <real_t> &partial_tau, 
            const ViewCArray <real_t>  &xi_point) = 0;


    }; // end of 3D parent class




    /*
     .-------------------------------. 
    | .----------------------------. |
    | |    _____       ________    | |
    | |   / ___ `.    |_   ___ `.  | |
    | |  |_/___) |      | |   `. \ | |
    | |   .'____.'      | |    | | | |
    | |  / /____       _| |___.' / | |
    | |  |_______|    |________.'  | |
    | |                            | |
    | '----------------------------' |
     '-------------------------------' 
    */
    /*
    ===========================
    2D Quad 4 Elements
    ===========================


    The finite element local point numbering for a 4 node quadralateral is
    as follows

          Eta
           ^
           |
    3------+-----2
    |      |     |
    |      |     |
    |      |     |
    |      ------+------> Xi
    |            |
    |            |
    0------------1
    */




    class Quad4: public Element2D {
        
        protected:

            static const int num_verts_ = 4;
            static const int num_nodes_ = 4;
            static const int num_basis_ = 4;

            static real_t ref_vert[num_verts_*num_dim_];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts_];

        public:


            int num_verts();
            int num_nodes();
            int num_basis();
            
            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                ViewCArray <real_t>  &x_point, 
                const ViewCArray <real_t>  &xi_point, 
                const ViewCArray <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                ViewCArray <real_t>  &basis,
                const ViewCArray <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void  partial_xi_basis(
                ViewCArray <real_t>  &partial_xi, 
                const ViewCArray <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void  partial_eta_basis(
                ViewCArray <real_t> &partial_eta, 
                const ViewCArray <real_t>  &xi_point); 

            // Map from vertex to node
            int vert_node_map( const int vert_lid);

    }; // end of quad_4_2D class


    /*
    ===========================
    2D Quad 8 Elements
    ===========================


    The finite element local point numbering for a 8 node Hexahedral is
    as follows

           Eta
            ^
            |
    3-------6------2
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    7       +------5-----> Xi   
    |              |
    |              |
    |              |
    0------4-------1
    */

    class Quad8: public Element2D {
        
        protected:

            static const int num_verts_ = 8;
            static const int num_nodes_ = 8;
            static const int num_basis_ = 8;

            static real_t ref_vert[num_verts_*num_dim_];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts_];

        public:


            int num_verts();
            int num_nodes();
            int num_basis();

            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                ViewCArray <real_t>  &x_point, 
                const ViewCArray <real_t>  &xi_point, 
                const ViewCArray <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                ViewCArray <real_t>  &basis,
                const ViewCArray <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void  partial_xi_basis(
                ViewCArray <real_t>  &partial_xi, 
                const ViewCArray <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void  partial_eta_basis(
                ViewCArray <real_t> &partial_eta, 
                const ViewCArray <real_t>  &xi_point);

            // Map from vertex to node
            int vert_node_map( const int vert_lid);
    
    }; // end of quad_8_2D class

    /*
    ===========================
    2D Quad 12 Elements
    ===========================


    The finite element local point numbering for a 8 node Hexahedral is
    as follows (NEED TO DEFINE)

             Eta
              ^
              |
      3---7------6---2
      |       |      |
      |       |      |
     11       |      10
      |       |      |
      |       +------|-----> Xi   
      |              |
      8              9
      |              |
      0----4-----5---1

    */

    class Quad12: public Element2D {
        
        protected:

            static const int num_verts_ = 12;
            static const int num_nodes_ = 12;
            static const int num_basis_ = 12;

            static real_t ref_vert[num_verts_*num_dim_];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts_];

        public:


            int num_verts();
            int num_nodes();
            int num_basis();

            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                ViewCArray <real_t>  &x_point, 
                const ViewCArray <real_t>  &xi_point, 
                const ViewCArray <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                ViewCArray <real_t>  &basis,
                const ViewCArray <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void  partial_xi_basis(
                ViewCArray <real_t>  &partial_xi, 
                const ViewCArray <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void  partial_eta_basis(
                ViewCArray <real_t> &partial_eta, 
                const ViewCArray <real_t>  &xi_point);

            // Map from vertex to node
            int vert_node_map( const int vert_lid);

    }; // end of quad_8_2D class

    /*
    ==========================
     Arbitrary Order Elements
    ==========================

       __                   _ _   _
     / __ \                | | \ | |
    | |  | |_   _  __ _  __| |  \| |
    | |  | | | | |/ _` |/ _` | . ` |
    | |__| | |_| | (_| | (_| | |\  |
     \___\_\\__,_|\__,_|\__,_|_| \_| 

    Representative linear element for visualization
     
           Eta (j)
            ^
            |
    3--------------2
    |       |      |
    |       |      |
    |       |      |
    |       |      |
    |       +------|-----> Xi (i) 
    |              |
    |              |
    |              |
    0--------------1
    */


    class QuadN{
        public:

            const static int num_dim = 2;

            int num_basis;
            int num_verts;

            // calculates the basis values and derivatives in 1D
            // used in the basis_partials functiosn to build the 3D element
            void lagrange_1D(
                ViewCArray <real_t> &interp,          // interpolant
                ViewCArray <real_t> &Dinterp,         // derivative of function
                const real_t &x_point,                  // point of interest in element
                const ViewCArray <real_t> &xi_point,  // nodal positions in 1D, normally chebyshev
                const int &orderN);                     // order of element

            void corners (
                ViewCArray <real_t> &lag_nodes,   // Nodes of Lagrange elements 
                ViewCArray <real_t> &lag_corner,  // corner nodes of QuadN element
                const int &orderN);                 // Element order

            void physical_position (
                ViewCArray <real_t> &x_point,             // location in real space
                const ViewCArray <real_t> &lag_nodes,     // Nodes of Lagrange elements 
                const ViewCArray <real_t> &lag_basis_2d,  // 2D basis values 
                const int &orderN);                         // order of the element

            void basis_partials (
                ViewCArray <real_t> &lag_nodes,       // Nodes of Lagrange elements (to be filled in)
                ViewCArray <real_t> &nodes_1d,        // Nodal spacing in 1D, any spacing is accepted
                ViewCArray <real_t> &val_1d,          // Interpolant Value in 1D
                ViewCArray <real_t> &DVal_1d,         // Derivateive of basis in 1D
                ViewCArray <real_t> &val_2d,          // for holding the interpolant in each direction
                ViewCArray <real_t> &DVal_2d,         // for holding the derivatives in each direction
                ViewCArray <real_t> &lag_basis_2d,    // 2D basis values 
                ViewCArray <real_t> &lag_partial,     // Partial of basis 
                const ViewCArray <real_t> &xi_point,  // point of interest
                const int &orderN);                     // Element order
    };



    /* 
     .-------------------------------. 
    | .----------------------------. |
    | |    ______      ________    | |
    | |   / ____ `.   |_   ___ `.  | |
    | |   `'  __) |     | |   `. \ | |
    | |   _  |__ '.     | |    | | | |
    | |  | \____) |    _| |___.' / | |
    | |   \______.'   |________.'  | |
    | |                            | |
    | '----------------------------' |
     '------------------------------' 
    */

/*
==========================
  Hex 8
==========================

 The finite element local vertex numbering for a 8 node Hexahedral is
 as follows

         Mu (k)
         |     Eta (j)    
         |    /
         |   /
     6---+----7
    /|   |   /|
   / |   |  / |
  4--------5  |
  |  |    -|--+---> Xi (i)
  |  |     |  |
  |  2-----|--3
  | /      | /       
  |/       |/
  0----*----1
 
*/

class Hex8: public Element3D {
        
        protected:

            static const int num_verts_ = 8;
            static const int num_nodes_ = 8;
            static const int num_basis_ = 8;

            static real_t ref_vert[num_verts_*num_dim_];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts_];

        public:


            int num_verts();
            int num_nodes();
            int num_basis();

            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                ViewCArray <real_t>  &x_point, 
                const ViewCArray <real_t>  &xi_point, 
                const ViewCArray <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                ViewCArray <real_t>  &basis,
                const ViewCArray <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void partial_xi_basis(
                ViewCArray <real_t>  &partial_xi, 
                const ViewCArray <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void partial_eta_basis(
                ViewCArray <real_t> &partial_eta, 
                const ViewCArray <real_t>  &xi_point);

            // with repsect to Mu
            void partial_mu_basis(
                ViewCArray <real_t> &partial_mu, 
                const ViewCArray <real_t>  &xi_point);

            // Map from vertex to node
            inline int vert_node_map( const int vert_lid);

            real_t& ref_locs(const int vert_lid, const int dim);

    }; // end of hex 8 class


    /*
    ==========================
      Hex 20
    ==========================

    The finite element local point numbering for a 20 node Hexahedral is 
    as follows

         Mu (k)
             |     Eta (j)
             |    /
             |   /
        7----14----6
       /|         /|
     15 |       13 |
     / 19       /  18
    4----12----5   |
    |   |      |   |  --> Xi (i)
    |   |      |   |
    |   3---10-|---2
    16 /      17  /
    | 11       | 9         
    |/         |/
    0-----8----1
    */

    class Hex20: public Element3D {
        protected:

            static const int num_verts_ = 20;
            static const int num_nodes_ = 20;
            static const int num_basis_ = 20;

            static real_t ref_vert[num_verts_*num_dim_];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts_];

        public:


            int num_verts();
            int num_nodes();
            int num_basis();

            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                ViewCArray <real_t>  &x_point, 
                const ViewCArray <real_t>  &xi_point, 
                const ViewCArray <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                ViewCArray <real_t>  &basis,
                const ViewCArray <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void  partial_xi_basis(
                ViewCArray <real_t>  &partial_xi, 
                const ViewCArray <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void  partial_eta_basis(
                ViewCArray <real_t> &partial_eta, 
                const ViewCArray <real_t>  &xi_point);

            // with repsect to Mu
            void partial_mu_basis(
                ViewCArray <real_t> &partial_mu, 
                const ViewCArray <real_t>  &xi_point);

            // Map from vertex to node
            inline int vert_node_map( const int vert_lid);

            real_t& ref_locs(const int vert_lid, const int dim);

    }; // end of hex20 class

    /* 
    ==========================
      Hex 32
    ==========================

    The finite element local point numbering for a 32 node Hexahedral is 
    shown below


                 Mu (k)
                  ^         Eta (j)
                  |        /
                  |       /
                         /
          7----23------22----6
         /|                 /|
       15 |               14 |
       /  |               /  |
     12  31             13   30 
     /    |             /    |
    4-----20-----21----5     |
    |     |            |     |   ----> Xi (i)
    |    27            |     26  
    |     |            |     |
    28    |           29     |
    |     3----19------|18---2
    |    /             |    /
    |  11              |   10
    24 /              25  /
    | 8                | 9         
    |/                 |/
    0----16------17----1
    */

    class Hex32: public Element3D {
        
        protected:

            static const int num_verts_ = 32;
            static const int num_nodes_ = 32;
            static const int num_basis_ = 32;

            static real_t ref_vert[num_verts_*num_dim_];  // listed as {Xi, Eta, Mu}
            static const int vert_to_node[num_verts_];

        public:


            int num_verts();
            int num_nodes();
            int num_basis();

            // calculate a physical position in an element for a given xi,eta
            void physical_position(
                ViewCArray <real_t>  &x_point, 
                const ViewCArray <real_t>  &xi_point, 
                const ViewCArray <real_t>  &vertices);


            // calculate the value for the basis at each node for a given xi,eta
            void basis(
                ViewCArray <real_t>  &basis,
                const ViewCArray <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi
            void  partial_xi_basis(
                ViewCArray <real_t>  &partial_xi, 
                const ViewCArray <real_t>  &xi_point);


            // Partial derivative of shape functions with respect to Eta
            void  partial_eta_basis(
                ViewCArray <real_t> &partial_eta, 
                const ViewCArray <real_t>  &xi_point);

            // with repsect to Mu
            void partial_mu_basis(
                ViewCArray <real_t> &partial_mu, 
                const ViewCArray <real_t>  &xi_point);

            // Map from vertex to node
            int vert_node_map( const int vert_lid);

            real_t& ref_locs(const int vert_lid, const int dim);

    }; // end of hex32 class


    /*
    ==========================
     Arbitrary Order Elements
    ==========================
     _   _           _   _ 
    | | | | _____  _| \ | |
    | |_| |/ _ \ \/ /  \| |
    |  _  |  __/>  <| |\  |
    |_| |_|\___/_/\_\_| \_|
                            
    Representative linear element for visualization
       
           k
           |     j    
           |    /
           |   /
       6---+----7
      /|   |   /|
     / |   |  / |
    4--------5  |
    |  |    -|--+---> i
    |  |     |  |
    |  2-----|--3
    | /      | /       
    |/       |/
    0--------1
    */

    class HexN{
            
        protected:
            
            const static int num_dim_ = 3;

            // Nodes
            int num_nodes_1d_;
            int num_nodes_;

            CArray <real_t> HexN_Nodes_1d_;
            CArray <real_t> HexN_Nodes_;

            // Vertices
            int num_verts_1d_;
            int num_verts_;
            int num_basis_;

            CArray <real_t> HexN_Verts_1d_;
            CArray <real_t> HexN_Verts_;

            CArray <int> Vert_Node_map_;
            
            int order_;


        public:

            void setup_HexN(int elem_order);

            int num_verts();
            int num_nodes();
            int num_basis();
            int node_rid(int i, int j, int k) const;
            int vert_rid(int i, int j, int k) const;

            // Return the noda coordinates in reference space
            real_t &node_coords(int node_rlid, int dim);


            int vert_node_map(int vert_rid) const;
            
            // Evaluate the basis at a given point
            void basis(
                CArray <real_t> &basis,
                CArray <real_t> &point);

            void build_nodal_gradient(
                CArray <real_t> &gradient);

            // calculate the partial of the basis w.r.t xi at a given point
            void partial_xi_basis(
                CArray <real_t> &partial_xi, 
                CArray <real_t> &point);

            // calculate the partial of the basis w.r.t eta at a given point
            void partial_eta_basis(
                CArray <real_t> &partial_eta, 
                CArray <real_t> &point);

            // calculate the partial of the basis w.r.t mu at a given point
            void partial_mu_basis(
                CArray <real_t> &partial_mu, 
                CArray <real_t> &point);

            void lagrange_basis_1D(
                CArray <real_t> &interp,    // interpolant
                const real_t &x_point);     // point of interest in element
            
            void lagrange_derivative_1D(
                CArray <real_t> &partials,  //derivative
                const real_t &x_point);     // point of interest in element

            void create_lobatto_nodes(int element_order);

    };



    /*
    ==========================
     4D Tesseract element
    ==========================
     
    The finite element local point numbering for a 16 node Tesseract is
    based on the 3D Hex8 Ensight element
     

                     _.15--------------------------------------14
                _.+<    |\                              . >-"/ |
          _ .+>         | \                         .>"" ./    |
      .>""              |  \                     .""    /      |
    12----------------------+------------------13    ./        |
    | )<=               |    \               / | _/""          |
    |     )\+           |     \            / /"|               |
    |         (\=       |   _. 7---------+--6  |               |
    |             \>   .|+<    |       / . "|  |               |
    |               '4--+------+------5'    |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |      |      |     |  |               |
    |                |  |   _ .3------+-----2_ |               |
    |                |  | "   /       |   /'  '| \= _          |
    |                0--+---+---------1*"      |     '\        |
    |             ./'   |  "           \       |       ""\     |
    |            /      |/              \      |           '". |
    |         /        11----------------+-----+---------------10
    |      ./    .+<""                    )    |           .</
    |    /(   /(                            \  |     _.+</
    | ./  /"                                 \ |  >(
    8------------------------------------------9'

                      j
                      ^        k
                      |      /
                      |    / 
                      |  /
                      |/
                      +---------->i

       i = Xi
       j = Eta
       k = Mu
       t = Tau 

    */

    class Tess16: public Element4D {
        public:
            const static int num_verts = 16;
            const static int num_dim = 4;
            const static int num_basis = 16;

        protected:
            static real_t ref_vert[num_verts*num_dim];  // listed as {Xi, Eta, Mu, Tau}


        public:

            // calculate a physical position in an element for a given xi,eta,mu
            void physical_position(
                ViewCArray <real_t> &x_point,
                const ViewCArray <real_t> &xi_point,
                const ViewCArray <real_t> &vertices);
            
            // calculate the value for the basis at each node for a given xi,eta,mu,tau
            void basis(
                ViewCArray <real_t>  &basis,
                const ViewCArray <real_t>  &xi_point);

            // Partial derivative of shape functions with respect to Xi at Xi_point
            void partial_xi_basis(
                ViewCArray <real_t> &partial_xi, 
                const ViewCArray <real_t> &xi_point);

            // Partial derivative of shape functions with respect to Eta
            void partial_eta_basis(
                ViewCArray <real_t> &partial_eta, 
                const ViewCArray <real_t> &xi_point);

            // Partial derivative of shape functions with respect to Mu
            void partial_mu_basis(
                ViewCArray <real_t> &partial_mu, 
                const ViewCArray <real_t> &xi_point);

            // Partial derivative of shape functions with respect to Tau
            void partial_tau_basis(
                ViewCArray <real_t> &partial_tau, 
                const ViewCArray <real_t> &xi_point);                                          
    }; // End of Tess16 Element Class



/* Add enumerated list of element types to choose from */

// namespace elem_types
// {

//     enum elem_type
//     {
//         // 2D element types
//         Quad4  = 0,   // 4 node serendipity quadrilateral element 
//         Quad8  = 1,   // 8 node serendipity quadrilateral element 
//         Quad12 = 2,   // 12 node serendipity quadrilateral element 
//         QuadN  = 3,   // N node lagrangian quadrilateral element 
        
//         // 3D element types
//         Hex8   = 4,   // 8 node serendipity hexahedral element 
//         Hex20  = 5,   // 20 node serendipity hexahedral element 
//         Hex32  = 6,   // 32 node serendipity hexahedral element 
//         HexN   = 7    // N node lagrangian hexahedral element 
    
//         // add tesseract element

//     };

// } // end of elem_type namespace


// // Choose the type of element to use for the calculation
// struct elem_type_t {
//     // element type 
//     elem_types::elem_type type;   
// };


// void choose_elem_type(elem_type_t elem_type);


} // end namespace elements


// --- Element type choice ---

// 2D Element Types

// extern elements::Element2D* elem2D;

// extern elements::Quad4      Quad4_elem;
// extern elements::Quad8      Quad8_elem;
// extern elements::Quad12     Quad12_elem;
// extern elements::QuadN      QuadN_elem;

// // 3D element types
// extern elements::Element3D* elem;

// extern elements::Hex8      Hex8_elem;
// extern elements::Hex20     Hex20_elem;
// extern elements::Hex32     Hex32_elem;
// extern elements::HexN      HexN_elem;



// extern elements::elem_type_t*   elem_choice;

// Everything is HexN
extern elements::HexN      elem;



// Reference element
extern elements::ref_element  ref_elem;




#endif //ELEMENTS_H
