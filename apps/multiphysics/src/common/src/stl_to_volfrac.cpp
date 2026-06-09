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


#include <chrono>   // for timing

#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>

#include "stl_utils.hpp"
#include "stl_to_volfrac.hpp"


using namespace mtr;

const double PI = 3.14159265358979323846;


//
// -----------------------------------------------
//


// This function is for seriallizing the integation locations in the reference element
KOKKOS_INLINE_FUNCTION
size_t get_id_of_ijk(size_t i, size_t j, size_t k, size_t num_x, size_t num_y){
    return i + (j + k*num_y)*num_x;
}


KOKKOS_INLINE_FUNCTION
void get_bernstein_basis_fcns(const CArrayKokkos <double> &bern_basis,
                              const double xi, 
                              const double eta, 
                              const double mu,
                              const size_t eval_pnt_rid){

    // Output: bern_basis[8] = shape functions

    double Bx[2];
    Bx[0] = 1.0 - xi;
    Bx[1] = xi;

    double By[2];
    By[0] = 1.0 - eta;
    By[1] = eta;

    double Bz[2];
    Bz[0] = 1.0 - mu;
    Bz[1] = mu;

    size_t node_rid = 0;
    for(size_t k=0; k<=1; ++k) {
        for(size_t j=0; j<=1; ++j) {
            for(size_t i=0; i<=1; ++i) {

                // the order of nodes in the elem follows i,j,k ordering
                bern_basis(eval_pnt_rid, node_rid) = Bx[i] * By[j] * Bz[k];
                node_rid++;

            } // end for i
        } // end for j
    } // end for k

    return;

} // end function


KOKKOS_INLINE_FUNCTION
double calc_scalar_in_elem(const CArrayKokkos <double> &node_scalar,
                           const CArrayKokkos <double> &node_basis, 
                           const DCArrayKokkos <size_t> &nodes_in_elem,
                           const size_t elem_gid,
                           const size_t eval_pnt_rid){

    // bern_basis(eval_pnt_rid, num_basis)
    const size_t num_basis = node_basis.dims(1);

    // the physical location (x,y,z) in the element is vec_pt
    double scalar_pnt = 0;

    // calculate x,y,z location in the element using basis functions evaluated at (xi,eta,mu)
    for(size_t node_rid=0; node_rid<num_basis; node_rid++) {

        // get the node index for this node_rid
        size_t node_gid = nodes_in_elem(elem_gid, node_rid);

        scalar_pnt += node_basis(eval_pnt_rid,node_rid)*node_scalar(node_gid);

    } // end for nodes

    return scalar_pnt;

} // end function


KOKKOS_INLINE_FUNCTION
void calc_vector_in_elem(const CArrayKokkos <double> &vec_pnt,
                         const CArrayKokkos <double> &node_vec,
                         const CArrayKokkos <double> &node_basis, 
                         const CArrayKokkos <size_t> &nodes_in_elem,
                         const size_t elem_gid,
                         const size_t eval_pnt_rid){

    // bern_basis(eval_pnt_rid, num_basis)
    const size_t num_basis = node_basis.dims(1);

    // the vector value in the element is vec_pnt at this eval pnt
    vec_pnt(0) = 0;
    vec_pnt(1) = 0;
    vec_pnt(2) = 0;

    // calculate x,y,z location in the element using basis functions evaluated at (xi,eta,mu)
    for(size_t node_rid=0; node_rid<num_basis; node_rid++) {

        // get the node index for this node_rid
        size_t node_gid = nodes_in_elem(elem_gid, node_rid);

        for(size_t dim=0; dim<3; dim++){
            vec_pnt(dim) += node_basis(eval_pnt_rid,node_rid)*node_vec(node_gid,dim);
        } // end for dim

    } // end for nodes

    return;

} // end function




//------------------------------------------------------------------------
//
// Function that takes a stl file and paints it on a mesh 
//
// This function works as follows here:
//  1) Read in an STL file
//  2) loop over the nodes in the FE mesh then calculate the signed distance
//     function to the surface
//
//------------------------------------------------------------------------
int paint_stl_on_mesh(DCArrayKokkos <double> &elem_geo_volfrac_fill, 
                      const MPICArrayKokkos <double> &node_coords,
                      const DCArrayKokkos <size_t> &nodes_in_elem,
                      const size_t num_nodes,
                      const double scale_x,
                      const double scale_y,
                      const double scale_z,
                      const double origin_x,
                      const double origin_y,
                      const double origin_z,
                      const std::string &file_path)
{

    // -----------------
    // read .STL file


    printf("Reading STL file \n");

    auto [normal_host, 
            v0X_host, v0Y_host, v0Z_host, 
            v1X_host, v1Y_host, v1Z_host, 
            v2X_host, v2Y_host, v2Z_host, 
            num_inp_triangles_host] = binary_stl_reader(file_path);
    
    // Warning on C++ support:
    // At this time with C++, the contents from a tuple cannot 
    // be used inside a lambda function.  The parallel loops use 
    // lambda functions. To overcome this C++ limitation, all 
    // contents in the tuple will be copied or pointed to (Using 
    // a MATAR dual view) allowing the data to be used in parallel.
    const size_t num_inp_triangles = num_inp_triangles_host;
    DViewCArrayKokkos <double> normal(&normal_host(0,0), num_inp_triangles, 3);
    DViewCArrayKokkos <double> v0X(&v0X_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v0Y(&v0Y_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v0Z(&v0Z_host(0),num_inp_triangles);
    DViewCArrayKokkos <double> v1X(&v1X_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v1Y(&v1Y_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v1Z(&v1Z_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v2X(&v2X_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v2Y(&v2Y_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v2Z(&v2Z_host(0),num_inp_triangles); 

    normal.update_device(); 
    v0X.update_device(); 
    v0Y.update_device(); 
    v0Z.update_device();
    v1X.update_device(); 
    v1Y.update_device(); 
    v1Z.update_device(); 
    v2X.update_device(); 
    v2Y.update_device(); 
    v2Z.update_device(); 
    Kokkos::fence();



    // check for a closed surface and calculate volume
    double sum0 = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;

    // no MATAR reductions for multiple values yet, thus harding coding the kokkos loop
    Kokkos::parallel_reduce(
        "stl_mesh_sum_checks",
        num_inp_triangles,
    
        // this is the for loop coding
        KOKKOS_LAMBDA(const int tri,         
                      double& sum_lcl_0,
                      double& sum_lcl_1,
                      double& sum_lcl_2,
                      double& sum_lcl_3) {

            //A = p1 - p0;
            //B = p2 - p0;
            vec_t A;
            A.x =  v1X(tri) - v0X(tri);
            A.y =  v1Y(tri) - v0Y(tri);
            A.z =  v1Z(tri) - v0Z(tri);
            
            vec_t B;
            B.x =  v2X(tri) - v0X(tri);
            B.y =  v2Y(tri) - v0Y(tri);
            B.z =  v2Z(tri) - v0Z(tri);
            
            // area normal of triangle
            vec_t N;
            N.x = 0.5*(A.y * B.z - A.z * B.y); 
            N.y = 0.5*(A.z * B.x - A.x * B.z); 
            N.z = 0.5*(A.x * B.y - A.y * B.x); 

            sum_lcl_0 += N.x; // checking area normal in x sums to zero
            sum_lcl_1 += N.y; // checking area normal in y sums to zero
            sum_lcl_2 += N.z; // checking area normal in z sums to zero

            const double area = magnitude(N); // surface area
            if(fabs(area)<1.e-15) Kokkos::abort("ERROR: STL triangle has zero surface area \n");

            // get the volume

            double mid_point[3];
            mid_point[0] = v0X(tri);
            mid_point[1] = v0Y(tri);
            mid_point[2] = v0Z(tri);
            
            mid_point[0] += v1X(tri);
            mid_point[1] += v1Y(tri);
            mid_point[2] += v1Z(tri);

            mid_point[0] += v2X(tri);
            mid_point[1] += v2Y(tri);
            mid_point[2] += v2Z(tri);
   
            // need to divide by 3 to get mid point location
            mid_point[0] /= 3.0;
            mid_point[1] /= 3.0;
            mid_point[2] /= 3.0;

            // n dot x_mid_point, this is for getting the volume
            sum_lcl_3 += N.x*mid_point[0] + N.y*mid_point[1] + N.z*mid_point[2]; 

        },
     sum0,
     sum1,
     sum2,
     sum3); // end parallel sum over all STL triangles

    
    // DEBUG:
    //printf("closed surface checks, nx=%f, ny=%f, nz=%f \n", sum0, sum1, sum2);

    if(fabs(sum0)>1.e-8 || 
       fabs(sum1)>1.e-8 ||
       fabs(sum2)>1.e-8 )
    {
       printf("*** Error in STL geometry, it is not a closed part surface \n");
       printf("error = (%f, %f, %f) \n", sum0, sum1, sum2);
       std::runtime_error("**** Failed to Paint STL file on mesh ****");
    } // end if on normal check
    printf("STL volume = %f \n", sum3/3.0); // Vol = (1/3* sum N dot x_mid_point)


    // apply scaling and shifts if they are set
    if (fabs(origin_x) > 1.0e-13 ||
        fabs(origin_y) > 1.0e-13 ||
        fabs(origin_z) > 1.0e-13 || 
        fabs(scale_x - 1.0) > 1.0e-13 ||
        fabs(scale_y - 1.0) > 1.0e-13 ||
        fabs(scale_z - 1.0) > 1.0e-13)
    {
        // loop over all facets
        FOR_ALL(tri,0, num_inp_triangles, {

            // update x coordinates of triangle vertices
            v0X(tri) = scale_x*(v0X(tri) - origin_x);
            v1X(tri) = scale_x*(v1X(tri) - origin_x);
            v2X(tri) = scale_x*(v2X(tri) - origin_x);

            // update y coordinates of triangle vertices
            v0Y(tri) = scale_y*(v0Y(tri) - origin_y);
            v1Y(tri) = scale_y*(v1Y(tri) - origin_y);
            v2Y(tri) = scale_y*(v2Y(tri) - origin_y);

            // update z coordinates of triangle vertices
            v0Z(tri) = scale_z*(v0Z(tri) - origin_z);
            v1Z(tri) = scale_z*(v1Z(tri) - origin_z);
            v2Z(tri) = scale_z*(v2Z(tri) - origin_z);

        }); // end for all over triangles
        Kokkos::fence();

    } // end if

    // -----------------
    // Getting SDF at the mesh nodes
    printf("Getting SDF at the mesh nodes \n");

    CArrayKokkos <double> node_sdf(num_nodes,  "node_sdf");

    // loop over the nodes and calculate signed distance function (sdf)
    FOR_ALL(node_gid, 0, num_nodes, {

        // initialize to a big value
        node_sdf(node_gid) = 1.e32;
        
        // loop over all triangles in the file
        // Note: future work should make this an oct-tree search
        for(size_t tri=0; tri<num_inp_triangles; tri++){

            // node vec
            vec_t node_p(node_coords(node_gid,0), node_coords(node_gid,1), node_coords(node_gid,2));

            // points defining the triangular facet
            const vec_t tri_v0(v0X(tri), v0Y(tri), v0Z(tri));
            const vec_t tri_v1(v1X(tri), v1Y(tri), v1Z(tri));
            const vec_t tri_v2(v2X(tri), v2Y(tri), v2Z(tri));

            double sdf_lcl = signed_distance_to_triangle(node_p, 
                                                         tri_v0, 
                                                         tri_v1, 
                                                         tri_v2);

            // check to see if this triangle is closer to the point
            if(fabs(sdf_lcl) < fabs(node_sdf(node_gid))){
                node_sdf(node_gid) = sdf_lcl;
            } // end if

        } // end for tri in STL file

    }); // end parallel for  



    // -----------------
    // reference element integration
    // WARNING: move to Element routines, this is a place-holder
    printf("building reference element \n");

    const size_t num_eval_pnts_1D = 9;  
    const double ref_pnts_1D[9]={0.0, 0.125, 0.24, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};

    const size_t num_eval_pnts = num_eval_pnts_1D*num_eval_pnts_1D*num_eval_pnts_1D;
    const size_t num_nodes_in_elem = 8; // WARNING WARNING WARNING

    // build the reference element
    CArrayKokkos <double> bern_basis(num_eval_pnts,num_nodes_in_elem);

    // loop over the eval reference element
    FOR_ALL(k, 0, num_eval_pnts_1D,
            j, 0, num_eval_pnts_1D,
            i, 0, num_eval_pnts_1D, {

        // coords
        const double xi  = ref_pnts_1D[i];
        const double eta = ref_pnts_1D[j];
        const double mu  = ref_pnts_1D[k];

        // get the eval ref elem local id
        const size_t eval_pnt_rid = get_id_of_ijk(i, j, k, num_eval_pnts_1D, num_eval_pnts_1D);

        // get bern basis (WARNING: update to use ELEMENTS library arbitrary order library)
        get_bernstein_basis_fcns(bern_basis, xi, eta, mu, eval_pnt_rid);

    }); // end for all


    const size_t num_elems = nodes_in_elem.dims(0); // nodes_in_elem(num_elems, num_nodes_in_elem) 

    printf("calculating volfrac using SDF vals at nodes \n");

    // evaluate SDF at this eval point
    FOR_ALL(elem_gid, 0, num_elems, {

        double num_inside_part = 0.0;

        // loop over the eval points in this element
        for(size_t k=0; k<num_eval_pnts_1D; k++){
            for(size_t j=0; j<num_eval_pnts_1D; j++){
                for(size_t i=0; i<num_eval_pnts_1D; i++){

                    // get the eval ref elem local id
                    const size_t eval_pnt_rid = get_id_of_ijk(i, j, k, num_eval_pnts_1D, num_eval_pnts_1D);

                    // evaluate SDF at this point
                    const double sdf_val = calc_scalar_in_elem(node_sdf, bern_basis, nodes_in_elem, elem_gid, eval_pnt_rid);

                    if(sdf_val<=0){
                        // we are inside the part
                        num_inside_part++;

                        //printf("elem= %d, sdf = %f \n", elem_gid, sdf_val);
                    } // end if

                } // end for i
            } // end for j
        } // end for k


        //  The ratio of hits to number of points is vol frac
        elem_geo_volfrac_fill(elem_gid) = num_inside_part/((double)num_eval_pnts); // coded for 3D

    }); // end parallel for


    
    printf("Finished STL to Volfrac Calculation \n\n");


    return 1;

} // end function


