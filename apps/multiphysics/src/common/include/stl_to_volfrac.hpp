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
#ifndef STL_TO_VOLFRAC_H
#define STL_TO_VOLFRAC_H

#include "matar.h"

// -----------------------------------------------
// types
// -----------------------------------------------

struct bin_keys_t{
    size_t i,j,k;
};

// a vector type with 3 components
struct vec_t{
    double x;
    double y;
    double z;
    
    // default constructor
    KOKKOS_INLINE_FUNCTION
    vec_t (){};
    
    // overloaded constructor
    KOKKOS_INLINE_FUNCTION
    vec_t(const double x_in, const double y_in, const double z_in){
        x = x_in;
        y = y_in;
        z = z_in;
    };
    
}; // end vec_t


// a triangle data type
struct triangle_t {
    
    vec_t normal; // surface normal
    
    vec_t p[3];   // three nodes with x,y,z coords
    
    // default constructor
    KOKKOS_INLINE_FUNCTION
    triangle_t(){};
    
    // overloaded constructor to accept 3 vectors
    KOKKOS_INLINE_FUNCTION
    triangle_t (const vec_t p_in[3])
    {
        // store the coords
        p[0]=p_in[0]; 
        p[1]=p_in[1]; 
        p[2]=p_in[2];

        // calculate the normal to this surface

        //A = p1 - p0;
        //B = p2 - p0;
        vec_t A;
        A.x = p[1].x - p[0].x;
        A.y = p[1].y - p[0].y;
        A.z = p[1].z - p[0].z;
        
        vec_t B;
        B.x = p[2].x - p[0].x;
        B.y = p[2].y - p[0].y;
        B.z = p[2].z - p[0].z;
        
        normal.x = A.y * B.z - A.z * B.y;
        normal.y = A.z * B.x - A.x * B.z;
        normal.z = A.x * B.y - A.y * B.x;
        
        const double mag = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
        
        // save the unit normal
        normal.x /= mag;
        normal.y /= mag;
        normal.z /= mag;
    };
    
}; // end triangle_t

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_bernstein_basis_fcns
///
/// \brief a function to calculate the bernstein basis at (xi, eta, mu)
///
/// \param bern_basis the bernstein basis
/// \param xi first reference coordinate to evaluate at
/// \param eta second reference coordinate to evaluate at
/// \param mu third reference coordinate to evaluate at
/// \param eval_pnt_rid the local storage id for basis value at (xi, eta, mu)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_INLINE_FUNCTION
void get_bernstein_basis_fcns(CArrayKokkos <double> &bern_basis,
                              const double xi, 
                              const double eta, 
                              const double mu,
                              const size_t eval_pnt_rid);


/////////////////////////////////////////////////////////////////////////////
///
/// \fn calc_scalar_in_elem
///
/// \brief a function that returns a scalar at (xi, eta, mu) in the element
///
/// \param node_scalar the scalar field at the element nodes
/// \param node_basis the basis functions from each node
/// \param nodes_in_elems a list of nodes in the element
/// \param elem_gid the element grid index
/// \param eval_pnt_rid the local storage id for basis value at (xi, eta, mu)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_INLINE_FUNCTION
double calc_scalar_in_elem(const CArrayKokkos <double> &node_scalar,
                           const CArrayKokkos <double> &node_basis, 
                           const CArrayKokkos <size_t> &nodes_in_elems,
                           const size_t elem_gid,
                           const size_t eval_pnt_rid);


/////////////////////////////////////////////////////////////////////////////
///
/// \fn calc_vector_in_elem
///
/// \brief a function to calculate a vector at (xi, eta, mu) in the element
///
/// \param vec_pnt the vector field at (xi, eta, mu) is returned by reference
/// \param node_vec the vector field at the element nodes
/// \param node_basis the basis functions from each node
/// \param nodes_in_elems a list of nodes in the element
/// \param elem_gid the element grid index
/// \param eval_pnt_rid the local storage id for basis value at (xi, eta, mu)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_INLINE_FUNCTION
void calc_vector_in_elem(CArrayKokkos <double> &vec_pnt,
                         const CArrayKokkos <double> &node_vec,
                         const CArrayKokkos <double> &node_basis, 
                         const CArrayKokkos <size_t> &nodes_in_elems,
                         const size_t elem_gid,
                         const size_t eval_pnt_rid);


/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_sdf_to_tri
///
/// \brief a function to calculate a signed distance function (sdf) to a
///        triangular facet of a part surface
///
/// \param node_sdf the sdf at the nodes of the mesh
/// \param num_tris_in_bin the number of triangles in a bin mesh
/// \param i_bin the i index of the bin mesh
/// \param j_bin the j index of the bin mesh
/// \param k_bin the k index of the bin mesh
/// \param num_bins_x the number of bins in the x direction (i.e. i direction)
/// \param num_bins_y the number of bins in the x direction (i.e. j direction)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_INLINE_FUNCTION
void get_sdf_to_tri(CArrayKokkos<double> &node_sdf,
                    const CArrayKokkos <size_t> &num_tris_in_bin,
                    const size_t i_bin, 
                    const size_t j_bin, 
                    const size_t k_bin, 
                    const size_t num_bins_x, 
                    const size_t num_bins_y);


/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_bin_keys
///
/// \brief a function to that returns the bin keys, which are (i_bin, j_bin, k_bin),
///        at the (x_pt, y_pt, z_pt)
///
/// \param x_pt a point in physical space on the mesh
/// \param y_pt a point in physical space on the mesh
/// \param z_pt a point in physical space on the mesh
///
/////////////////////////////////////////////////////////////////////////////                    
KOKKOS_INLINE_FUNCTION
bin_keys_t get_bin_keys(const double x_pt, 
                        const double y_pt, 
                        const double z_pt);


/////////////////////////////////////////////////////////////////////////////
///
/// \fn build_3D_zone_nodes
///
/// \brief a function to calculate the ref local indices (rid) to access the 
///        nodes of the zones in an element
///
/// \param node_rids_in_zone_lids an array storing the re
/// \param num_zones_1d number of zones in each direction of the element
/// \param num_nodes_1d number of nodes in each direction of the element
///
/////////////////////////////////////////////////////////////////////////////  
void build_3D_zone_nodes(DCArrayKokkos <size_t> &node_rids_in_zone_lids, 
                         const size_t num_zones_1d,
                         const size_t num_nodes_1d);


/////////////////////////////////////////////////////////////////////////////
///
/// \fn cross product
///
/// \brief returns the cross product of two vectors
///
/// \param vec_t the first vector having coordinates in x, y, z space
/// \param vec_t the second vector having coordinates in x, y, z space
///
///////////////////////////////////////////////////////////////////////////// 
KOKKOS_INLINE_FUNCTION
vec_t cross(const vec_t &a, const vec_t &b);


/////////////////////////////////////////////////////////////////////////////
///
/// \fn dot product
///
/// \brief the dot product of two vectors
///
/// \param vec_t the first vector having coordinates in x, y, z space
/// \param vec_t the second vector having coordinates in x, y, z space
///
///////////////////////////////////////////////////////////////////////////// 
KOKKOS_INLINE_FUNCTION
double dot(const vec_t &a, const vec_t &b);


/////////////////////////////////////////////////////////////////////////////
///
/// \fn magnitude
///
/// \brief the magnitude of a vector
///
/// \param vec_t the vector having coordinates in x, y, z space
///
///////////////////////////////////////////////////////////////////////////// 
KOKKOS_INLINE_FUNCTION
double magnitude(const vec_t &a);


/////////////////////////////////////////////////////////////////////////////
///
/// \fn distance
///
/// \brief returns the magnitude of a the difference tween two vectors
///
/// \param vec_t the first vector having coordinates in x, y, z space
/// \param vec_t the second vector having coordinates in x, y, z space
///
///////////////////////////////////////////////////////////////////////////// 
KOKKOS_INLINE_FUNCTION
double distance(const vec_t &a, const vec_t &b);


/////////////////////////////////////////////////////////////////////////////
///
/// \fn paint_stl_on_mesh
///
/// \brief paints volume fraction on a mesh based on an interface defined 
//         by a STL CAD file
///
/// \param elem_vol_frac the element volume fraction
/// \param node_coords the (x,y,z) coordinate values at the mesh nodes
/// \param nodes_in_elems the nodes in the element
/// \param Pn_order the element Pn order based on nodes
/// \param filename the file name of the STL file
///
///////////////////////////////////////////////////////////////////////////// 
int paint_stl_on_mesh(DCArrayKokkos <double> &elem_vol_frac, 
                      const DCArrayKokkos <double> &node_coords,
                      const DCArrayKokkos <size_t> &nodes_in_elems,
                      const size_t Pn_order,
                      const std::string &filename);

#endif