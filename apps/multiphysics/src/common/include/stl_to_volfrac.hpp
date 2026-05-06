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
                           const DCArrayKokkos <size_t> &nodes_in_elems,
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
                         const DCArrayKokkos <size_t> &nodes_in_elems,
                         const size_t elem_gid,
                         const size_t eval_pnt_rid);


/////////////////////////////////////////////////////////////////////////////
///
/// \fn paint_stl_on_mesh
///
/// \brief paints volume fraction on a mesh based on an interface defined 
//         by a STL CAD file
///
/// \param elem_geo_volfrac_fill is the element volume fraction for a fill 
/// \param node_coords the (x,y,z) coordinate values at the mesh nodes
/// \param nodes_in_elems the nodes in the element
/// \param filename the file name of the STL file
///
///////////////////////////////////////////////////////////////////////////// 
int paint_stl_on_mesh(DCArrayKokkos <double> &elem_geo_volfrac_fill, 
                      const DCArrayKokkos <double> &node_coords,
                      const DCArrayKokkos <size_t> &nodes_in_elem,
                      const size_t num_nodes,
                      const std::string &file_path);


/////////////////////////////////////////////////////////////////////////////
///
/// \fn binary_stl_reader
///
/// \brief a function to read a binary STL file, exporting triangular 
//         facet coordinates of the surface and the number of facets.
///
/// \param filepath to the STL file
///
/////////////////////////////////////////////////////////////////////////////
inline
std::tuple<
    CArray<double>,   // normal
    CArray<double>, CArray<double>, CArray<double>,   // v0X, v0Y, v0Z 
    CArray<double>, CArray<double>, CArray<double>,   // v1X, v1Y, v1Z 
    CArray<double>, CArray<double>, CArray<double>,   // v2X, v2Y, v2Z
    size_t // n_facets
>
binary_stl_reader(const std::string& path);




// end of hpp file
#endif