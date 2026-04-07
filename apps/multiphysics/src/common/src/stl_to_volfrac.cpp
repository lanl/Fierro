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

// -----------------------------------------------
// pointcloud reproducing kernels in C++
//  Nathaniel Morgan
// -----------------------------------------------

#include <chrono>   // for timing

#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cmath>

#include <cstdlib> // For rand() and srand()


#include "matar.h"

#include "lu_solver.hpp"


#include <set> // for unorded map testing 



#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

using namespace mtr;

const double PI = 3.14159265358979323846;

// -----------------------------------------------
// inputs:


// the bin sizes for finding neighboring points
const double bin_dx = 0.05; // bins in x
const double bin_dy = 0.05; // bins in y
const double bin_dz = 0.05; // bins in z


//
// -----------------------------------------------


KOKKOS_INLINE_FUNCTION
void get_bernstein_basis_fcns(CArrayKokkos <double> &bern_basis,
                              const double xi, 
                              const double eta, 
                              const double mu,
                              const size_t eval_pnt_lid){

    // Output: bern_basis[8] = shape functions

    const double Bx[2];
    Bx[0] = 1.0 - xi;
    Bx[1] = xi;

    const double By[2];
    By[0] = 1.0 - eta;
    By[1] = eta;

    const double Bz[2];
    Bz[0] = 1.0 - mu;
    Bz[1] = mu;

    size_t node_rid = 0;
    for(size_t k=0; k<=1; ++k) {
        for(size_t j=0; j<=1; ++j) {
            for(size_t i=0; i<=1; ++i) {

                // the order of nodes in the elem follows i,j,k ordering
                bern_basis(eval_pnt_lid, node_rid) = Bx[i] * By[j] * Bz[k];
                node_rid++;

            } // end for i
        } // end for j
    } // end for k

    return;

} // end function


KOKKOS_INLINE_FUNCTION
double calc_scalar_in_elem(const CArrayKokkos <double> &node_scalar,
                           const CArrayKokkos <double> &node_basis, 
                           const CArrayKokkos <size_t> &nodes_in_elems,
                           const size_t elem_gid,
                           const size_t eval_pnt_lid){

    // bern_basis(eval_pnt_lid, num_basis)
    const size_t num_basis = basis.dims(1);

    // the physical location (x,y,z) in the element is vec_pt
    double scalar_pnt = 0;

    // calculate x,y,z location in the element using basis functions evaluated at (xi,eta,mu)
    for(size_t node_rid=0; node_rid<num_basis; node_rid++) {

        // get the node index for this node_rid
        size_t node_gid = nodes_in_elems(elem_gid, node_rid);

        scalar_pnt += node_basis(eval_pnt_lid,node_rid)*node_scalar(node_gid);

    } // end for nodes

    return scalar_pnt;

} // end function



KOKKOS_INLINE_FUNCTION
void calc_vector_in_elem(CArrayKokkos <double> &vec_pnt,
                         const CArrayKokkos <double> &node_vec,
                         const CArrayKokkos <double> &node_basis, 
                         const CArrayKokkos <size_t> &nodes_in_elems,
                         const size_t elem_gid,
                         const size_t eval_pnt_lid){

    // bern_basis(eval_pnt_lid, num_basis)
    const size_t num_basis = basis.dims(1);

    // the vector value in the element is vec_pnt at this eval pnt
    vec_pnt(0) = 0;
    vec_pnt(1) = 0;
    vec_pnt(2) = 0;

    // calculate x,y,z location in the element using basis functions evaluated at (xi,eta,mu)
    for(size_t node_rid=0; node_rid<num_basis; node_rid++) {

        // get the node index for this node_rid
        size_t node_gid = nodes_in_elems(elem_gid, node_rid);

        for(size_t dim=0; dim<3; dim++){
            vec_pnt(dim) += node_basis(eval_pnt_lid,node_rid)*node_vec(node_gid,dim);
        } // end for dim

    } // end for nodes

    return;

} // end function


KOKKOS_INLINE_FUNCTION
void get_sdf_to_tri(CArrayKokkos<double> &node_sdf,
                    const CArrayKokkos <size_t> &num_tris_in_bin,
                    const size_t icount, 
                    const size_t jcount, 
                    const size_t kcount, 
                    const size_t num_bins_x, 
                    const size_t num_bins_y){
    
    // get bin neighbor gid on this search boundary
    const size_t neighbor_bin_gid = get_gid(icount, jcount, kcount, num_bins_x, num_bins_y);

    // loop over all the triangles in this bin
    for(size_t tri_lid=0; tri_lid<num_tris_in_bin(neighbor_bin_gid); tri_lid++){

        // get the triangle index
        const size_t tri_gid = tris_in_bin(neighbor_bin_gid,tri_lid);

        // get the x,y,z position of the triangle 
        const double x_tri = tri_coords(tri_gid, 0);
        const double y_tri = tri_coords(tri_gid, 1);
        const double z_tri = tri_coords(tri_gid, 2);

        // store tri coords in a vector
        vec_t vec_of_tri_coords(x_tri, y_tri, z_tri);

        // create the triangle object
        triangle_t triangle(vec_of_tri_coords);

        // distance from node position to the triangular facet
        const double dx_tri_pt = x_tri-x_pt;
        const double dy_tri_pt = y_tri-y_pt;
        const double dz_tri_pt = z_tri-z_pt;
        vec_t distance_vec_to_tri(dx_tri_pt, dy_tri_pt, dz_tri_pt);

        // get the signed distance for this triangle
        double tri_distance = magnitude(distance_vec_to_tri);
        double sign = dot(triangle.normal, distance_vec_to_tri);

        if(fabs(tri_distance) < fabs(node_sdf(node_gid))){
            node_sdf(node_gid) = fabs(tri_distance)*sign;
        }

    } // end for
} // end function 

///////////


struct bin_keys_t{
    size_t i,j,k;
};

KOKKOS_INLINE_FUNCTION
size_t get_gid(size_t i, size_t j, size_t k, size_t num_x, size_t num_y){
    return i + (j + k*num_y)*num_x;
}

KOKKOS_INLINE_FUNCTION
bin_keys_t get_bin_keys(const double x_pt, 
                        const double y_pt, 
                        const double z_pt){
            

    double i_dbl = fmax(0, round((x_pt - X0 - bin_dx*0.5)/bin_dx - 1.0e-10)); // x = ih + X0 + dx_bin*0.5
    double j_dbl = fmax(0, round((y_pt - Y0 - bin_dy*0.5)/bin_dy - 1.0e-10));
    double k_dbl = fmax(0, round((z_pt - Z0 - bin_dz*0.5)/bin_dz - 1.0e-10));

    bin_keys_t bin_keys; // save i,j,k to the bin keys

    // get the integer for the bins
    bin_keys.i = (size_t)i_dbl;
    bin_keys.j = (size_t)j_dbl;
    bin_keys.k = (size_t)k_dbl;

    return bin_keys;

} // end function

KOKKOS_INLINE_FUNCTION
size_t get_bin_gid(const double x_pt, 
                   const double y_pt, 
                   const double z_pt, 
                   const size_t num_bins_x,
                   const size_t num_bins_y,
                   const size_t num_bins_z){
            

    double i_dbl = fmin(num_bins_x-1, fmax(0.0, round((x_pt - X0)/bin_dx - 1.0e-8))); // x = ih + X0
    double j_dbl = fmin(num_bins_y-1, fmax(0.0, round((y_pt - Y0)/bin_dy - 1.0e-8)));
    double k_dbl = fmin(num_bins_z-1, fmax(0.0, round((z_pt - Z0)/bin_dz - 1.0e-8)));

    // get the integers for the bins
    size_t i = (size_t)i_dbl;
    size_t j = (size_t)j_dbl;
    size_t k = (size_t)k_dbl;
    
    // get the 1D index for this bin                               
    return get_gid(i, j, k, num_bins_x, num_bins_y);

} // end function



void build_3D_zone_nodes(DCArrayKokkos <size_t> &node_lids_in_zone_lids, 
                         const size_t num_sub_zones_1d,
                         const size_t num_points_1d){  

    // this is 3D
    const size_t num_nodes_in_zone = 8;

    // running on CPU as there is little parallelism

    // loop over i,j,k of the sub-zones
    FOR_LOOP(k, 0, num_sub_zones_1d,
             j, 0, num_sub_zones_1d,
             i, 0, num_sub_zones_1d, {
                
                // get the sub_zone local index
                const size_t sub_zone_lid = get_gid(i, j, k, num_sub_zones_1d, num_sub_zones_1d);

                // get node lids for this sub-zone
                node_lids_in_zone_lids.host(sub_zone_lid,0) = get_gid(i,   j,   k, num_points_1d, num_points_1d); 
                node_lids_in_zone_lids.host(sub_zone_lid,1) = get_gid(i+1, j,   k, num_points_1d, num_points_1d); 
                node_lids_in_zone_lids.host(sub_zone_lid,2) = get_gid(i,   j+1, k, num_points_1d, num_points_1d); 
                node_lids_in_zone_lids.host(sub_zone_lid,3) = get_gid(i+1, j+1, k, num_points_1d, num_points_1d); 

                node_lids_in_zone_lids.host(sub_zone_lid,4) = get_gid(i,   j,   k+1, num_points_1d, num_points_1d); 
                node_lids_in_zone_lids.host(sub_zone_lid,5) = get_gid(i+1, j,   k+1, num_points_1d, num_points_1d); 
                node_lids_in_zone_lids.host(sub_zone_lid,6) = get_gid(i,   j+1, k+1, num_points_1d, num_points_1d); 
                node_lids_in_zone_lids.host(sub_zone_lid,7) = get_gid(i+1, j+1, k+1, num_points_1d, num_points_1d); 

    }); // end parallel k,j,i
    node_lids_in_zone_lids.update_device();

} // end function


//
// -----------------------------------------------


std::tuple<
    CArray<double>,   // normal
    CArray<double>, CArray<double>, CArray<double>,   // v1X, v1Y, v1Z
    CArray<double>, CArray<double>, CArray<double>,   // v2X, v2Y, v2Z
    CArray<double>, CArray<double>, CArray<double>,   // v3X, v3Y, v3Z
    size_t // n_facets
>
binary_stl_reader(const std::string& path)
{
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in) { std::perror("open"); std::exit(EXIT_FAILURE); }

    const std::streamoff filesize = in.tellg();
    if (filesize < 100) {
        std::cerr << "ERROR: File too small to be a valid STL\n";
        std::exit(EXIT_FAILURE);
    }
    in.seekg(0);

    // ---- check if ASCII -------------------------------------------------
    char magic[6] = { 0 };
    in.read(magic, 5);          // read first 5 chars
    in.seekg(0);               // rewind
    if (std::strncmp(magic, "solid", 5) == 0) {
        std::cerr
            << "ERROR: \"" << path
            << "\" looks like an **ASCII** STL (starts with \"solid\").\n"
            << "Re‑export it as *binary* or implement an ASCII parser.\n";
        std::exit(EXIT_FAILURE);        // or call ascii_stl_reader();
    }

    // ---- read 80‑byte header + nominal facet count ----------------------
    char header[80];                in.read(header, 80);
    size_t n_facets_nominal;  in.read(reinterpret_cast<char*>(&n_facets_nominal), 4);

    // ---- compute expected count from file size to sanity‑check ----------
    // binary facet record = 50 bytes (12×4 + 12×4 + 12×4 + 2)
    const size_t n_facets_from_size =
        static_cast<size_t>((filesize - 84) / 50);

    size_t n_facets = n_facets_nominal;
    if (n_facets_nominal != n_facets_from_size) {
        std::cout << "WARNING: facet count in header (" << n_facets_nominal
            << ") disagrees with file size (" << n_facets_from_size
            << ").  Using size‑derived value.\n";
        n_facets = n_facets_from_size;
    }
    std::cout << "STL facets: " << n_facets << '\n';

    // ---- allocate MATAR arrays -----------------------------------------
    CArray<double> normal(n_facets, 3);
    CArray<double> v1X(n_facets), v1Y(n_facets), v1Z(n_facets);
    CArray<double> v2X(n_facets), v2Y(n_facets), v2Z(n_facets);
    CArray<double> v3X(n_facets), v3Y(n_facets), v3Z(n_facets);

    // ---- read facet records --------------------------------------------
    double nrm[3], v1[3], v2[3], v3[3];
    for (unsigned int i = 0; i < n_facets; ++i) {
        in.read(reinterpret_cast<char*>(nrm), 12);
        in.read(reinterpret_cast<char*>(v1), 12);
        in.read(reinterpret_cast<char*>(v2), 12);
        in.read(reinterpret_cast<char*>(v3), 12);
        in.ignore(2);                        // attribute byte count

        for (int d = 0; d < 3; ++d) normal(i, d) = nrm[d];
        v1X(i) = v1[0]; v1Y(i) = v1[1]; v1Z(i) = v1[2];
        v2X(i) = v2[0]; v2Y(i) = v2[1]; v2Z(i) = v2[2];
        v3X(i) = v3[0]; v3Y(i) = v3[1]; v3Z(i) = v3[2];
    }
    return { normal,v1X,v1Y,v1Z,v2X,v2Y,v2Z,v3X,v3Y,v3Z,n_facets };

} // end of function




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
    
    // overloaded constructor
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



// cross prodcut
KOKKOS_INLINE_FUNCTION
vec_t cross(const vec_t &a, const vec_t &b) {
    return {a.y*b.z - a.z*b.y,
            a.z*b.x - a.x*b.z,
            a.x*b.y - a.y*b.x};
}

KOKKOS_INLINE_FUNCTION
double dot(const vec_t &a, const vec_t &b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

KOKKOS_INLINE_FUNCTION
double magnitude(const vect_t &a){
    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

KOKKOS_INLINE_FUNCTION
double distance(const vect_t &a, , const vec_t &b){
    return sqrt((a.x-b.x)*(a.x-b.x) + 
                (a.y-b.y)*(a.y-b.y) + 
                (a.z-b.z)*(a.z-b.z));
}

//---------------------------------------------------------
//
// Function that takes a stl file and paints it on a mesh 
//
//---------------------------------------------------------
int paint_stl_on_mesh(DCArrayKokkos <double> &elem_vol_frac, 
                      const DCArrayKokkos <double> &node_coords,
                      const DCArrayKokkos <size_t> &nodes_in_elems,
                      const size_t Pn_order,
                      const std::string &filename)
{


    // -----------------
    // read .STL file
    // -----------------

    printf("Reading STL file \n\n");

    auto [normal_host, 
            v1X_host, v1Y_host, v1Z_host, 
            v2X_host, v2Y_host, v2Z_host, 
            v3X_host, v3Y_host, v3Z_host, 
            num_inp_triangles_host] = binary_stl_reader(filename);
    
    // Warning on C++ support:
    // At this time with C++, the contents from a tuple cannot 
    // be used inside a lambda function.  The parallel loops use 
    // lambda functions. To overcome this C++ limitation, all 
    // contents in the tuple will be copied or pointed to (Using 
    // a MATAR dual view) allowing the data to be used in parallel.
    const size_t num_inp_triangles = num_inp_triangles_host;
    DViewCArrayKokkos <double> normal(&normal_host(0,0), num_inp_triangles, 3);
    DViewCArrayKokkos <double> v1X(&v1X_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v1Y(&v1Y_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v1Z(&v1Z_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v2X(&v2X_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v2Y(&v2Y_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v2Z(&v2Z_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v3X(&v3X_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v3Y(&v3Y_host(0),num_inp_triangles); 
    DViewCArrayKokkos <double> v3Z(&v3Z_host(0),num_inp_triangles);

    normal.update_device(); 
    v1X.update_device(); 
    v1Y.update_device(); 
    v1Z.update_device(); 
    v2X.update_device(); 
    v2Y.update_device(); 
    v2Z.update_device(); 
    v3X.update_device(); 
    v3Y.update_device(); 
    v3Z.update_device();

    
    // array storing the geo-center of triangles
    DCArrayKokkos <double> tri_coords(num_inp_triangles, 3, "tri_coords");

    // point values are the geo-center of triangles
    FOR_ALL(tri, 0, num_inp_triangles, {

        // point on surface
        tri_coords(tri, 0) =  1.0/3.0*(v1X(tri) + v2X(tri) + v3X(tri));
        tri_coords(tri, 1) =  1.0/3.0*(v1Y(tri) + v2Y(tri) + v3Y(tri));
        tri_coords(tri, 2) =  1.0/3.0*(v1Z(tri) + v2Z(tri) + v3Z(tri));

    }); // end parallel for tri's in the file


    // -----------------
    // create bin mesh
    // -----------------

    printf("making bins \n");

    // find (xmin, ymin, zmin) and (xmax, ymax, zmax) for building bin mesh
    double xmin_lcl, ymin_lcl, zmin_lcl;
    double xmin, ymin, zmin;

    Kokkos::parallel_reduce("MultiMin", num_inp_triangles, KOKKOS_LAMBDA(int tri, double& xmin_lcl, double& ymin_lcl, double& zmin_lcl) { 
        if (tri_coords(tri, 0) > xmin_lcl) xmin_lcl = tri_coords(tri, 0); 
        if (tri_coords(tri, 1) > ymin_lcl) ymin_lcl = tri_coords(tri, 1); 
        if (tri_coords(tri, 2) > zmin_lcl) zmin_lcl = tri_coords(tri, 2); 
    }, Kokkos::Min<double>(xmin), Kokkos::Min<double>(ymin), Kokkos::Min<double>(zmin));
    
    Kokkos::parallel_reduce("MultiMax", num_inp_triangles, KOKKOS_LAMBDA(int tri, double& xmax_lcl, double& ymax_lcl, double& zmax_lcl) { 
        if (tri_coords(tri, 0) > xmax_lcl) xmax_lcl = tri_coords(tri, 0); 
        if (tri_coords(tri, 1) > ymax_lcl) ymax_lcl = tri_coords(tri, 1); 
        if (tri_coords(tri, 2) > zmax_lcl) zmax_lcl = tri_coords(tri, 2); 
    }, Kokkos::Max<double>(xmax), Kokkos::Max<double>(ymax), Kokkos::Max<double>(zmax));

    
    // the number of nodes in the bin mesh
    size_t num_bins_x = 20;//(size_t)( round( (xmax - xmin)/bin_dx) + 1 );  
    size_t num_bins_y = 20;//(size_t)( round( (ymax - ymin)/bin_dy) + 1 );  
    size_t num_bins_z = 20;//(size_t)( round( (zmax - zmin)/bin_dz) + 1 );  

    size_t num_bins = num_bins_x*num_bins_y*num_bins_z;
    printf("num bins x=%zu, y=%zu, z=%zu \n", num_bins_x, num_bins_y, num_bins_z);

    // bins and their connectivity to each other and points
    DCArrayKokkos <bin_keys_t> keys_in_bin(num_bins, "keys_in_bin"); // mapping from gid to (i,j,k)
    DCArrayKokkos <size_t> num_tris_in_bin(num_bins, "num_bins");
    num_tris_in_bin.set_values(0);
    DRaggedRightArrayKokkos <size_t> tris_in_bin; // allocated later

    // connectivity from points to bins
    DCArrayKokkos <size_t> tri_bin_gid(num_inp_triangles, "points_in_gid");
    CArrayKokkos <size_t>  tri_bin_lid_storage(num_inp_triangles, "bin_lid_storage");  // only used to create storage
    DCArrayKokkos <int> bin_tri_stencil(num_bins, 6, "bin_tri_stencil");   // how imin,imax,jmin,jmax,kmin,kmax range for bins in stencil

        
    printf("Starting timers \n\n");


    // -------------
    // start timer
    auto time_1 = std::chrono::high_resolution_clock::now();

    // build reverse mapping between gid and i,j,k
    FOR_ALL(i, 0, num_bins_x,
            j, 0, num_bins_y,
            k, 0, num_bins_z, {
        

        // get bin gid for this i,j,k
        size_t bin_gid = get_gid(i, j, k, num_bins_x, num_bins_y);

        // the i,j,k for this bin
        bin_keys_t bin_keys;
        bin_keys.i = i;
        bin_keys.j = j;
        bin_keys.k = k;

        // save mapping from bin_gid to bin_keys i,j,k
        keys_in_bin(bin_gid) = bin_keys;

    });
    Kokkos::fence();
    keys_in_bin.update_host();

    // end timer
    auto time_2 = std::chrono::high_resolution_clock::now();


    // -------------------------------------------------------------------
    // below here, making dual maps between bins and triangles
    // -------------------------------------------------------------------

    // start timer
    auto time_3 = std::chrono::high_resolution_clock::now();

    printf("building maps between triangles and bins \n");

    // save bin id to triangles
    FOR_ALL(tri_gid, 0, num_inp_triangles, {

        // get the 1D index for this bin
        size_t bin_gid = get_bin_gid(tri_coords(tri_gid,0), 
                                     tri_coords(tri_gid,1), 
                                     tri_coords(tri_gid,2),
                                     num_bins_x, 
                                     num_bins_y,
                                     num_bins_z);

        size_t storage_lid = Kokkos::atomic_fetch_add(&num_tris_in_bin(bin_gid), 1);
        tri_bin_gid(tri_gid) = bin_gid; // the id of the bin
        tri_bin_lid_storage(tri_gid) = storage_lid; // the storage place in the bin

    }); // end for all
    Kokkos::fence();
    tri_bin_gid.update_host();
    num_tris_in_bin.update_host();


    // allocate tris in bin connectivity
    tris_in_bin = DRaggedRightArrayKokkos <size_t> (num_tris_in_bin, "num_tris_in_bin");

    // save tris in bin
    FOR_ALL(tri_gid, 0, num_inp_triangles, {

        // get bin gid
        size_t bin_gid = tri_bin_gid(tri_gid);

        // get it's storage location in the ragged right compressed storage
        size_t storage_lid = tri_bin_lid_storage(tri_gid);

        // save the point to this bin
        tris_in_bin(bin_gid, storage_lid) = tri_gid;

    }); // end for all



    // ------------------------------------------------
    // Find the closest tri neighbors on the bin mesh
    // ------------------------------------------------
    
    FOR_ALL(i, 0, num_bins_x,
            j, 0, num_bins_y,
            k, 0, num_bins_z, {

        
        // get bin gid for this i,j,k
        size_t bin_gid = get_gid(i, j, k, num_bins_x, num_bins_y);  
        

        // initialize imin stencil to -1, used for knowing if CAD was in range of bin mesh
        for(int id=0; id<6; id++){
            bin_tri_stencil(bin_gid,id) = -1-id;
        }


        // get bin gid at i,j,k
        size_t num_tris_found = num_tris_in_bin(bin_gid);

        // commenting this out to remove thread divergence
        // if (num_tris_found > 0){
        //     // Once I find at least 1 triangle, I can exit as this is the closest one found

        //     bin_tri_stencil(bin_gid,0) = i; //imin
        //     bin_tri_stencil(bin_gid,1) = i; //imax
        //     bin_tri_stencil(bin_gid,2) = j; //jmin
        //     bin_tri_stencil(bin_gid,3) = j; //jmax
        //     bin_tri_stencil(bin_gid,4) = k; //kmin
        //     bin_tri_stencil(bin_gid,5) = k; //kmax

        // } else 

        // establish the stencil size to find at least 1 triangle
        for(size_t stencil=1; stencil<1000000; stencil++){

            // i-1:i+1
            const size_t imin = MAX(0, i-stencil);
            const size_t imax = MIN(num_bins_x-1, i+stencil);

            // j-1:j+1
            const size_t jmin = MAX(0, j-stencil);
            const size_t jmax = MIN(num_bins_y-1, j+stencil);

            // k-1:k+1
            const size_t kmin = MAX(0, k-stencil);
            const size_t kmax = MIN(num_bins_z-1, k+stencil);

            // i-search boundaries
            for (size_t kcount = kmin; kcount <= kmax; kcount++) {
                for (size_t jcount = jmin; jcount <= jmax; jcount++) {

                    size_t icount = imin;     
                    // get bin neighbor gid on this search boundary
                    const size_t neighbor_bin_gid = get_gid(icount, jcount, kcount, num_bins_x, num_bins_y);
                    num_tris_found += num_tris_in_bin(neighbor_bin_gid);

                    icount = imax;
                    // get bin neighbor gid on this search boundary
                    const size_t neighbor_bin_gid = get_gid(icount, jcount, kcount, num_bins_x, num_bins_y);
                    num_tris_found += num_tris_in_bin(neighbor_bin_gid);

                } // end for jcount
            } // end for kcount

            // j-search boundaries
            for (size_t kcount = kmin; kcount <= kmax; kcount++) {
                for (size_t icount = imin; icount <= imax; icount++) {

                    size_t jcount = jmin;
                    // get bin neighbor gid on this search boundary
                    const size_t neighbor_bin_gid = get_gid(icount, jcount, kcount, num_bins_x, num_bins_y);
                    num_tris_found += num_tris_in_bin(neighbor_bin_gid);

                    jcount = jmax;
                    // get bin neighbor gid on this search boundary
                    const size_t neighbor_bin_gid = get_gid(icount, jcount, kcount, num_bins_x, num_bins_y);
                    num_tris_found += num_tris_in_bin(neighbor_bin_gid);

                } // end for icount
            } // end for kcount
            
            // k-search boundaries
            for (size_t jcount = jmin; jcount <= jmax; jcount++) {
                for (size_t icount = imin; icount <= imax; icount++) {

                    size_t kcount = kmin;
                    // get bin neighbor gid on this search boundary
                    const size_t neighbor_bin_gid = get_gid(icount, jcount, kcount, num_bins_x, num_bins_y);
                    num_tris_found += num_tris_in_bin(neighbor_bin_gid);

                    kcount = kmax;
                    // get bin neighbor gid on this search boundary
                    const size_t neighbor_bin_gid = get_gid(icount, jcount, kcount, num_bins_x, num_bins_y);
                    num_tris_found += num_tris_in_bin(neighbor_bin_gid);

                } // end icount
            } // end jcount
                

            // Once I find at least 1 triangle, I can exit as this is the closest one found
            if (num_tris_found > 0){

                bin_tri_stencil(bin_gid,0) = imin;
                bin_tri_stencil(bin_gid,1) = imax;
                bin_tri_stencil(bin_gid,2) = jmin;
                bin_tri_stencil(bin_gid,3) = jmax;
                bin_tri_stencil(bin_gid,4) = kmin;
                bin_tri_stencil(bin_gid,5) = kmax;

                break;
            }
            // else increase stencil size, the for loop interator increases stencil range

            
        } // end for stencil 

    }); // end for all
    Kokkos::fence();
    bin_tri_stencil.update_host();

     
    printf("building signed distance function using maps and the bins \n");


    // ----------------------------------------------------------------
    // calculating the signed distance function (SDF) at element nodes
    // ----------------------------------------------------------------

    // memory to store SDF at the nodes
    CArrayKokkos <double> node_sdf(num_nodes);
    node_sdf.set_values(1.0e32); // initialize distance to a large value

    // find closest triangles (using bins) to each node and save SDF
    FOR_ALL(node_gid, 0, num_nodes, {

        // get physical position for this node
        const double x_pt = node_coords(node_gid,0);
        const double y_pt = node_coords(node_gid,1);
        const double z_pt = node_coords(node_gid,2);

        // make vector of node coords
        vec_t vec_node_position(x_pt, y_pt, z_pt);

        // get i,j,k for this bin
        bin_keys_t bin_keys = get_bin_keys(x_pt, y_pt, z_pt);

        const size_t i = bin_keys.i;
        const size_t j = bin_keys.j;
        const size_t k = bin_keys.k;

        // get bin gid for this i,j,k
        size_t bin_gid = get_gid(i, j, k, num_bins_x, num_bins_y);  

        // get imin, imax, ..., kmin, kmax for closest triangles in the bin
        imin = bin_tri_stencil(bin_gid,0);
        imax = bin_tri_stencil(bin_gid,1);
        jmin = bin_tri_stencil(bin_gid,2);
        jmax = bin_tri_stencil(bin_gid,3);
        kmin = bin_tri_stencil(bin_gid,4);
        kmax = bin_tri_stencil(bin_gid,5);

        // try if check on imin=imax=jmin=jmax=kmin=kmax

        // i-search boundaries
        for (size_t kcount = kmin; kcount <= kmax; kcount++) {
            for (size_t jcount = jmin; jcount <= jmax; jcount++) {

                size_t icount = imin;   
                
                // get the signed distance function
                get_sdf_to_tri(node_sdf,
                               num_tris_in_bin,
                               icount, 
                               jcount, 
                               kcount, 
                               num_bins_x, 
                               num_bins_y);

                icount = imax;

                // get the signed distance function
                get_sdf_to_tri(node_sdf,
                               num_tris_in_bin,
                               icount, 
                               jcount, 
                               kcount, 
                               num_bins_x, 
                               num_bins_y);

            } // end for jcount
        } // end for kcount

        // j-search boundaries
        for (size_t kcount = kmin; kcount <= kmax; kcount++) {
            for (size_t icount = imin; icount <= imax; icount++) {

                size_t jcount = jmin;

                // get the signed distance function
                get_sdf_to_tri(node_sdf,
                               num_tris_in_bin,
                               icount, 
                               jcount, 
                               kcount, 
                               num_bins_x, 
                               num_bins_y);

                jcount = jmax;

                // get the signed distance function
                get_sdf_to_tri(node_sdf,
                               num_tris_in_bin,
                               icount, 
                               jcount, 
                               kcount, 
                               num_bins_x, 
                               num_bins_y);

            } // end for icount
        } // end for kcount
        
        // k-search boundaries
        for (size_t jcount = jmin; jcount <= jmax; jcount++) {
            for (size_t icount = imin; icount <= imax; icount++) {

                size_t kcount = kmin;

                // get the signed distance function
                get_sdf_to_tri(node_sdf,
                               num_tris_in_bin,
                               icount, 
                               jcount, 
                               kcount, 
                               num_bins_x, 
                               num_bins_y);

                kcount = kmax;

                // get the signed distance function
                get_sdf_to_tri(node_sdf,
                               num_tris_in_bin,
                               icount, 
                               jcount, 
                               kcount, 
                               num_bins_x, 
                               num_bins_y);

            } // end icount
        } // end jcount

    }); // end parallel for nodes in mesh


    // --------------------------------------------
    // calculating SDF in subzones using nodal SDF
    // --------------------------------------------

    // build details on the element
    const size_t num_sub_zones_1d = Pn_order;
    const size_t num_points_1d = Pn_order+1;
    const size_t num_nodes_in_elem = num_points_1d*num_points_1d*num_points_1d; 

    // PLIC (using SDF) can be applied on the zones in the element
    //const size_t num_sub_zones_in_elem = Pn_order*Pn_order*Pn_order; // if using PLIC

    // storage for the local ids of the nodes in the sub zone, needed for PLIC
    //DCArrayKokkos <size_t> node_lids_in_zone_lids(num_sub_zones_in_elem,8); // nodes of the sub zones in an element
    //build_3D_zone_nodes(node_lids_in_zone_lids, num_sub_zones_1d, num_points_1d);
    
    const size_t num_eval_pnts = 9;  
    const double ref_pnts_1D[9]={0.0, 0.125, 0.24, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};

    // build the reference element
    CArrayKokkos <size_t> bern_basis(num_eval_pnts,num_nodes_in_elem);

    // loop over the eval reference element
    FOR_ALL(k, 0, num_eval_pnts,
            j, 0, num_eval_pnts,
            i, 0, num_eval_pnts, {

        // coords
        const double xi  = ref_pnts_1D[i];
        const double eta = ref_pnts_1D[j];
        const double mu  = ref_pnts_1D[k];

        // get the eval_lid
        const size_t eval_pnt_lid = get_gid(i, j, k, num_eval_pnts, num_eval_pnts);

        // get bern basis (WARNING: update to use ELEMENTS library arbitrary order library)
        get_bernstein_basis_fcns(bern_basis, xi, eta, mu, eval_pnt_lid);

    }); // end for all


    // evaluate SDF at this eval point
    FOR_ALL(elem_gid, 0, num_elems, {

        double num_inside_part = 0;

        // loop over the eval points in this element
        for(size_t k=0; k<num_eval_pnts; k++){
            for(size_t j=0; j<num_eval_pnts; j++){
                for(size_t i=0; i<num_eval_pnts; i++){

                    // evaluate SDF at this point
                    const double sdf_val = calc_scalar_in_elem(node_sdf, bern_basis, nodes_in_elems, elem_gid, eval_pnt_lid);

                    if(sdf_val<0){
                        // we are inside the part
                        num_inside_part++
                    } // end if

                } // end for i
            } // end for j
        } // end for k


        //  The ratio of hits to number of points is vol frac
        elem_vol_frac(elem_gid) = num_inside_part/(num_eval_pnts*num_eval_pnts*num_eval_pnts); // coded for 3D

    }); // end parallel for
        

        


    // end timer
    auto time_4 = std::chrono::high_resolution_clock::now();

    printf("done calculating volume fraction in elements \n");




    // -----------------
    //  Timers
    // -----------------
    printf("\n");
    std::chrono::duration <double, std::milli> ms = time_2 - time_1;
    std::cout << "runtime to create bins = " << ms.count() << "ms\n\n";

    ms = time_4 - time_3;
    std::cout << "runtime to find and save neighbors = " << ms.count() << "ms\n\n";





    printf("Writing VTK Graphics File \n\n");

    std::ofstream out("cloud.vtk");

    out << "# vtk DataFile Version 3.0\n";
    out << "3D point cloud\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << num_points << " float\n";
    for (size_t tri_gid = 0; tri_gid < num_points; ++tri_gid) {
        out << tri_coords.host(tri_gid,0) << " " 
            << tri_coords.host(tri_gid,1) << " " 
            << tri_coords.host(tri_gid,2) << "\n";
    }

    out << "\nPOINT_DATA " << num_points << "\n";
    out << "SCALARS lvlset float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (size_t tri_gid = 0; tri_gid < num_points; ++tri_gid) {
        out << 0 << "\n";
    }

    
    printf("Finished \n\n");



    return 0;
    
} // end main
