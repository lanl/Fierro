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


//
// -----------------------------------------------
//


std::tuple<
    CArray<double>,   // normal
    CArray<double>, CArray<double>, CArray<double>,   // v0X, v0Y, v0Z 
    CArray<double>, CArray<double>, CArray<double>,   // v1X, v1Y, v1Z
    CArray<double>, CArray<double>, CArray<double>,   // v2X, v2Y, v2Z
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

    uint32_t n_facets = n_facets_nominal;
    if (n_facets_nominal != n_facets_from_size) {
        std::cout << "WARNING: facet count in header (" << n_facets_nominal
            << ") disagrees with file size (" << n_facets_from_size
            << ").  Using size‑derived value.\n";
        n_facets = n_facets_from_size;
    }
    std::cout << "STL facets: " << n_facets << '\n';

    // ---- allocate MATAR arrays -----------------------------------------
    CArray<double> normal(n_facets, 3);
    CArray<double> v0X(n_facets), v0Y(n_facets), v0Z(n_facets);
    CArray<double> v1X(n_facets), v1Y(n_facets), v1Z(n_facets);
    CArray<double> v2X(n_facets), v2Y(n_facets), v2Z(n_facets);

    // ---- read facet records --------------------------------------------
    float nrm[3], v0[3], v1[3], v2[3];
    for (size_t i = 0; i < n_facets; ++i) {
        in.read(reinterpret_cast<char*>(nrm), 12);
        in.read(reinterpret_cast<char*>(v0), 12);
        in.read(reinterpret_cast<char*>(v1), 12);
        in.read(reinterpret_cast<char*>(v2), 12);
        in.ignore(2);                        // attribute byte count

        for (int dim = 0; dim < 3; ++dim){ 
             normal(i, dim) = (double)nrm[dim]; 
        } // end for


        v0X(i) = (double)v0[0]; 
        v0Y(i) = (double)v0[1]; 
        v0Z(i) = (double)v0[2];

        v1X(i) = (double)v1[0]; 
        v1Y(i) = (double)v1[1]; 
        v1Z(i) = (double)v1[2];

        v2X(i) = (double)v2[0]; 
        v2Y(i) = (double)v2[1]; 
        v2Z(i) = (double)v2[2];
        
    } // end for facet

    return { normal,v0X,v0Y,v0Z,v1X,v1Y,v1Z,v2X,v2Y,v2Z,n_facets };

} // end of function to read STL file


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
                      const DCArrayKokkos <double> &node_coords,
                      const DCArrayKokkos <size_t> &nodes_in_elem,
                      const size_t num_nodes,
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

        printf("node_sdf %d = %f \n", node_gid, node_sdf(node_gid));

    }); // end parallel for  


    // ----- testing coding above ------

    // -----------------
    //  Viz for checking STL reader
    // -----------------
    printf("Writing VTK STL Point File \n\n");

    std::ofstream out("tri_points.vtk");

    out << "# vtk DataFile Version 3.0\n";
    out << "3D point cloud\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << num_inp_triangles << " float\n";
    for (size_t tri_gid = 0; tri_gid < num_inp_triangles; tri_gid++) {

        const double x_tri = (v0X_host(tri_gid) + v1X_host(tri_gid) + v2X_host(tri_gid))/3.0;
        const double y_tri = (v0Y_host(tri_gid) + v1Y_host(tri_gid) + v2Y_host(tri_gid))/3.0;
        const double z_tri = (v0Z_host(tri_gid) + v1Z_host(tri_gid) + v2Z_host(tri_gid))/3.0;

        out << x_tri << " " 
            << y_tri << " " 
            << z_tri << "\n";
    }

    out << "\nPOINT_DATA " << num_inp_triangles << "\n";
    out << "SCALARS lvlset float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (size_t tri_gid = 0; tri_gid < num_inp_triangles; ++tri_gid) {
        out << 0 << "\n";
    }

    // -----

    printf("Writing VTK SDF File \n\n");

    std::ofstream out_sdf("sdf_nodes.vtk");

    out_sdf << "# vtk DataFile Version 3.0\n";
    out_sdf << "3D point cloud\n";
    out_sdf << "ASCII\n";
    out_sdf << "DATASET POLYDATA\n";
    out_sdf << "POINTS " << num_nodes << " float\n";
    for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {

        out_sdf << node_coords.host(node_gid,0) << " " 
                << node_coords.host(node_gid,1) << " " 
                << node_coords.host(node_gid,2) << "\n";

    }

    out_sdf << "\nPOINT_DATA " << num_nodes << "\n";
    out_sdf << "SCALARS SDF float 1\n";
    out_sdf << "LOOKUP_TABLE default\n";
    for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
        out_sdf << node_sdf(node_gid) << "\n";
    }


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

        printf("vol frac in elem %d = %f \n", elem_gid, elem_geo_volfrac_fill(elem_gid));
    }); // end parallel for


    
    printf("Finished STL to Volfrac Calculation \n\n");


    return 1;

} // end function


