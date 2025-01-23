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

#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // fmin, fmax, abs note: fminl is long
#include <sys/stat.h>
#include <mpi.h>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "Tpetra_Details_FixedHashTable.hpp"
#include "Tpetra_Import.hpp"

#include "elements.h"
#include "matar.h"
#include "utilities.h"
#include "FEA_Module_Inertial.h"
#include "Simulation_Parameters/FEA_Module/Inertial_Parameters.h"
#include "Simulation_Parameters/Simulation_Parameters.h"
#include "Solver.h"

#define MAX_ELEM_NODES 8
#define DENSITY_EPSILON 0.0001

using namespace utils;

FEA_Module_Inertial::FEA_Module_Inertial(
    Inertial_Parameters& params, Solver* Solver_Pointer,
    const int my_fea_module_index)
    : FEA_Module(Solver_Pointer)
{
    // assign interfacing index
    my_fea_module_index_ = my_fea_module_index;
    Module_Type = FEA_MODULE_TYPE::Inertial;

    // acquire base class data from existing simparam in solver (gets yaml options etc.)
    module_params = &params;
    simparam      = &(Solver_Pointer->simparam);

    // TO parameters
    nodal_density_flag = simparam->nodal_density_flag;

    // property initialization flags
    mass_init   = false;
    com_init[0] = com_init[1] = com_init[2] = false;

    // property update counters
    mass_update = com_update[0] = com_update[1] = com_update[2] = -1;

    // RCP initialization
    mass_gradients_distributed = Teuchos::null;
    center_of_mass_gradients_distributed = Teuchos::null;

    // construct per element inertial property vectors
    Global_Element_Masses    = Teuchos::rcp(new MV(element_map, 1));
    Global_Element_Volumes   = Teuchos::rcp(new MV(element_map, 1));
    Global_Element_Moments_x = Teuchos::rcp(new MV(element_map, 1));
    Global_Element_Moments_y = Teuchos::rcp(new MV(element_map, 1));
    Global_Element_Moments_z = Teuchos::rcp(new MV(element_map, 1));
    Global_Element_Moments_of_Inertia_xx = Teuchos::rcp(new MV(element_map, 1));
    Global_Element_Moments_of_Inertia_yy = Teuchos::rcp(new MV(element_map, 1));
    Global_Element_Moments_of_Inertia_zz = Teuchos::rcp(new MV(element_map, 1));
    Global_Element_Moments_of_Inertia_xy = Teuchos::rcp(new MV(element_map, 1));
    Global_Element_Moments_of_Inertia_xz = Teuchos::rcp(new MV(element_map, 1));
    Global_Element_Moments_of_Inertia_yz = Teuchos::rcp(new MV(element_map, 1));
}

FEA_Module_Inertial::~FEA_Module_Inertial()
{
}

/* ----------------------------------------------------------------------
   Compute the mass of each element; estimated with quadrature
------------------------------------------------------------------------- */

void FEA_Module_Inertial::compute_element_masses(const_host_vec_array design_variables, bool max_flag, bool use_initial_coords)
{
    if(simparam->topology_optimization_on){
        compute_element_masses_TO(design_variables, max_flag, use_initial_coords);
    }
    else if(simparam->shape_optimization_on){
        compute_element_masses_SO(design_variables, max_flag, use_initial_coords);
    }

}

/* ------------------------------------------------------------------------------------------
   Compute the mass of each element when using shape optimization; estimated with quadrature
--------------------------------------------------------------------------------------------- */

void FEA_Module_Inertial::compute_element_masses_SO(const_host_vec_array design_coords, bool max_flag, bool use_initial_coords)
{
    // local number of uniquely assigned elements
    size_t nonoverlap_nelements = element_map->getLocalNumElements();
    // initialize memory for volume storage
    host_vec_array Element_Masses = Global_Element_Masses->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    // local variable for host view in the dual view
    const_host_vec_array all_node_coords;
    if (use_initial_coords)
    {
        all_node_coords = all_initial_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    else
    {
        all_node_coords = all_design_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    const_host_vec_array all_design_densities;
    if (nodal_density_flag)
    {
        all_design_densities = all_node_densities_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    const_host_elem_conn_array nodes_in_elem  = global_nodes_in_elem_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    int                        nodes_per_elem = elem->num_basis();
    int                        z_quad, y_quad, x_quad, direct_product_count;
    size_t                     local_node_id;
    LO                         ielem;
    GO                         global_element_index;

    real_t Jacobian, current_density, weight_multiply;
    // CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
    // CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
    CArray<real_t>     legendre_nodes_1D(num_gauss_points);
    CArray<real_t>     legendre_weights_1D(num_gauss_points);
    
    real_t pointer_quad_coordinate[num_dim];
    real_t pointer_quad_coordinate_weight[num_dim];
    real_t pointer_interpolated_point[num_dim];
    real_t pointer_JT_row1[num_dim];
    real_t pointer_JT_row2[num_dim];
    real_t pointer_JT_row3[num_dim];
    
    ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate, num_dim);
    ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight, num_dim);
    ViewCArray<real_t> interpolated_point(pointer_interpolated_point, num_dim);
    ViewCArray<real_t> JT_row1(pointer_JT_row1, num_dim);
    ViewCArray<real_t> JT_row2(pointer_JT_row2, num_dim);
    ViewCArray<real_t> JT_row3(pointer_JT_row3, num_dim);

    real_t pointer_basis_values[elem->num_basis()];
    real_t pointer_basis_derivative_s1[elem->num_basis()];
    real_t pointer_basis_derivative_s2[elem->num_basis()];
    real_t pointer_basis_derivative_s3[elem->num_basis()];
    
    ViewCArray<real_t> basis_values(pointer_basis_values, elem->num_basis());
    ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1, elem->num_basis());
    ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2, elem->num_basis());
    ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3, elem->num_basis());
    
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(), num_dim);
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());

    // initialize weights
    elements::legendre_nodes_1D(legendre_nodes_1D, num_gauss_points);
    elements::legendre_weights_1D(legendre_weights_1D, num_gauss_points);

    Solver::node_ordering_convention active_node_ordering_convention = Solver_Pointer_->active_node_ordering_convention;
    
    CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_node_order(max_nodes_per_element);
    
    if ((active_node_ordering_convention == Solver::ENSIGHT && num_dim == 3) || (active_node_ordering_convention == Solver::IJK && num_dim == 2))
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 3;
        convert_node_order(3) = 2;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 7;
            convert_node_order(7) = 6;
        }
    }
    else
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 2;
        convert_node_order(3) = 3;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 6;
            convert_node_order(7) = 7;
        }
    }

    // loop over elements and use quadrature rule to compute volume from Jacobian determinant
    for (int nonoverlapping_ielem = 0; nonoverlapping_ielem < nonoverlap_nelements; nonoverlapping_ielem++)
    {
        global_element_index = element_map->getGlobalElement(nonoverlapping_ielem);
        ielem = all_element_map->getLocalElement(global_element_index);
        // acquire set of nodes for this local element
        for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
        {
            local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, convert_node_order(node_loop)));
            nodal_positions(node_loop, 0) = all_node_coords(local_node_id, 0);
            nodal_positions(node_loop, 1) = all_node_coords(local_node_id, 1);
            nodal_positions(node_loop, 2) = all_node_coords(local_node_id, 2);
            if (nodal_density_flag)
            {
                nodal_density(node_loop) = all_design_densities(local_node_id, 0);
            }
            /*
            if(myrank==1&&nodal_positions(node_loop,2)>10000000){
                std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
                std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
                std::fflush(stdout);
            }
            */
            // std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " "
            // << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) << " " << nodal_density(node_loop) <<std::endl;
        }

        // debug print of index
        // std::cout << "nonoverlap element id on TASK " << myrank << " is " << nonoverlapping_ielem << std::endl;
        // std::fflush(stdout);

        // initialize element mass
        Element_Masses(nonoverlapping_ielem, 0) = 0;

        direct_product_count = std::pow(num_gauss_points, num_dim);

        // loop over quadrature points
        for (int iquad = 0; iquad < direct_product_count; iquad++)
        {
            // set current quadrature point
            if (num_dim == 3)
            {
                z_quad = iquad / (num_gauss_points * num_gauss_points);
            }
            y_quad = (iquad % (num_gauss_points * num_gauss_points)) / num_gauss_points;
            x_quad = iquad % num_gauss_points;
            quad_coordinate(0) = legendre_nodes_1D(x_quad);
            quad_coordinate(1) = legendre_nodes_1D(y_quad);
            if (num_dim == 3)
            {
                quad_coordinate(2) = legendre_nodes_1D(z_quad);
            }

            // set current quadrature weight
            quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
            quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
            if (num_dim == 3)
            {
                quad_coordinate_weight(2) = legendre_weights_1D(z_quad);
            }
            else
            {
                quad_coordinate_weight(2) = 1;
            }
            weight_multiply = quad_coordinate_weight(0) * quad_coordinate_weight(1) * quad_coordinate_weight(2);

            // compute shape functions at this point for the element type
            elem->basis(basis_values, quad_coordinate);

            // compute all the necessary coordinates and derivatives at this point

            // compute shape function derivatives
            elem->partial_xi_basis(basis_derivative_s1, quad_coordinate);
            elem->partial_eta_basis(basis_derivative_s2, quad_coordinate);
            elem->partial_mu_basis(basis_derivative_s3, quad_coordinate);

            // compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
            // derivative of x,y,z w.r.t s
            JT_row1(0) = 0;
            JT_row1(1) = 0;
            JT_row1(2) = 0;
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                JT_row1(0) += nodal_positions(node_loop, 0) * basis_derivative_s1(node_loop);
                JT_row1(1) += nodal_positions(node_loop, 1) * basis_derivative_s1(node_loop);
                JT_row1(2) += nodal_positions(node_loop, 2) * basis_derivative_s1(node_loop);
            }

            // derivative of x,y,z w.r.t t
            JT_row2(0) = 0;
            JT_row2(1) = 0;
            JT_row2(2) = 0;
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                JT_row2(0) += nodal_positions(node_loop, 0) * basis_derivative_s2(node_loop);
                JT_row2(1) += nodal_positions(node_loop, 1) * basis_derivative_s2(node_loop);
                JT_row2(2) += nodal_positions(node_loop, 2) * basis_derivative_s2(node_loop);
            }

            // derivative of x,y,z w.r.t w
            JT_row3(0) = 0;
            JT_row3(1) = 0;
            JT_row3(2) = 0;
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                JT_row3(0) += nodal_positions(node_loop, 0) * basis_derivative_s3(node_loop);
                JT_row3(1) += nodal_positions(node_loop, 1) * basis_derivative_s3(node_loop);
                JT_row3(2) += nodal_positions(node_loop, 2) * basis_derivative_s3(node_loop);
                // debug print
                /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
                std::cout << " ELEMENT VOLUME JACOBIAN DEBUG ON TASK " << myrank << std::endl;
                std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
                std::fflush(stdout);
                }*/
            }

            // compute the determinant of the Jacobian
            Jacobian = JT_row1(0) * (JT_row2(1) * JT_row3(2) - JT_row3(1) * JT_row2(2)) -
                        JT_row1(1) * (JT_row2(0) * JT_row3(2) - JT_row3(0) * JT_row2(2)) +
                        JT_row1(2) * (JT_row2(0) * JT_row3(1) - JT_row3(0) * JT_row2(1));
            if (Jacobian < 0)
            {
                Jacobian = -Jacobian;
            }

            // compute density
            current_density = 0;
            if (max_flag)
            {
                current_density = 1;
            }
            else
            {
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    current_density += nodal_density(node_loop) * basis_values(node_loop);
                }
            }

            Element_Masses(nonoverlapping_ielem, 0) += current_density * weight_multiply * Jacobian;
        }
    
    }
}

/* ----------------------------------------------------------------------
   Compute the mass of each element; estimated with quadrature
------------------------------------------------------------------------- */

void FEA_Module_Inertial::compute_element_masses_TO(const_host_vec_array design_densities, bool max_flag, bool use_initial_coords)
{
    // local number of uniquely assigned elements
    size_t nonoverlap_nelements = element_map->getLocalNumElements();
    // initialize memory for volume storage
    host_vec_array Element_Masses = Global_Element_Masses->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    if (!nodal_density_flag)
    {
        compute_element_volumes();
    }
    const_host_vec_array Element_Volumes = Global_Element_Volumes->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    // local variable for host view in the dual view
    const_host_vec_array all_node_coords;
    if (use_initial_coords)
    {
        all_node_coords = all_initial_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    else
    {
        all_node_coords = all_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    const_host_vec_array all_design_densities;
    if (nodal_density_flag)
    {
        all_design_densities = all_node_densities_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    const_host_elem_conn_array nodes_in_elem  = global_nodes_in_elem_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    int                        nodes_per_elem = elem->num_basis();
    int                        z_quad, y_quad, x_quad, direct_product_count;
    size_t                     local_node_id;
    LO                         ielem;
    GO                         global_element_index;

    real_t Jacobian, current_density, weight_multiply;
    // CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
    // CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
    CArray<real_t>     legendre_nodes_1D(num_gauss_points);
    CArray<real_t>     legendre_weights_1D(num_gauss_points);
    
    real_t pointer_quad_coordinate[num_dim];
    real_t pointer_quad_coordinate_weight[num_dim];
    real_t pointer_interpolated_point[num_dim];
    real_t pointer_JT_row1[num_dim];
    real_t pointer_JT_row2[num_dim];
    real_t pointer_JT_row3[num_dim];
    
    ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate, num_dim);
    ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight, num_dim);
    ViewCArray<real_t> interpolated_point(pointer_interpolated_point, num_dim);
    ViewCArray<real_t> JT_row1(pointer_JT_row1, num_dim);
    ViewCArray<real_t> JT_row2(pointer_JT_row2, num_dim);
    ViewCArray<real_t> JT_row3(pointer_JT_row3, num_dim);

    real_t pointer_basis_values[elem->num_basis()];
    real_t pointer_basis_derivative_s1[elem->num_basis()];
    real_t pointer_basis_derivative_s2[elem->num_basis()];
    real_t pointer_basis_derivative_s3[elem->num_basis()];
    
    ViewCArray<real_t> basis_values(pointer_basis_values, elem->num_basis());
    ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1, elem->num_basis());
    ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2, elem->num_basis());
    ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3, elem->num_basis());
    
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(), num_dim);
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());

    // initialize weights
    elements::legendre_nodes_1D(legendre_nodes_1D, num_gauss_points);
    elements::legendre_weights_1D(legendre_weights_1D, num_gauss_points);

    Solver::node_ordering_convention active_node_ordering_convention = Solver_Pointer_->active_node_ordering_convention;
    
    CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_node_order(max_nodes_per_element);
    
    if ((active_node_ordering_convention == Solver::ENSIGHT && num_dim == 3) || (active_node_ordering_convention == Solver::IJK && num_dim == 2))
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 3;
        convert_node_order(3) = 2;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 7;
            convert_node_order(7) = 6;
        }
    }
    else
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 2;
        convert_node_order(3) = 3;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 6;
            convert_node_order(7) = 7;
        }
    }

    // loop over elements and use quadrature rule to compute volume from Jacobian determinant
    for (int nonoverlapping_ielem = 0; nonoverlapping_ielem < nonoverlap_nelements; nonoverlapping_ielem++)
    {
        global_element_index = element_map->getGlobalElement(nonoverlapping_ielem);
        ielem = all_element_map->getLocalElement(global_element_index);
        if (nodal_density_flag)
        {
            // acquire set of nodes for this local element
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, convert_node_order(node_loop)));
                nodal_positions(node_loop, 0) = all_node_coords(local_node_id, 0);
                nodal_positions(node_loop, 1) = all_node_coords(local_node_id, 1);
                nodal_positions(node_loop, 2) = all_node_coords(local_node_id, 2);
                if (nodal_density_flag)
                {
                    nodal_density(node_loop) = all_design_densities(local_node_id, 0);
                }
                /*
                if(myrank==1&&nodal_positions(node_loop,2)>10000000){
                  std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
                  std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
                  std::fflush(stdout);
                }
                */
                // std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " "
                // << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) << " " << nodal_density(node_loop) <<std::endl;
            }

            // debug print of index
            // std::cout << "nonoverlap element id on TASK " << myrank << " is " << nonoverlapping_ielem << std::endl;
            // std::fflush(stdout);

            // initialize element mass
            Element_Masses(nonoverlapping_ielem, 0) = 0;

            direct_product_count = std::pow(num_gauss_points, num_dim);

            // loop over quadrature points
            for (int iquad = 0; iquad < direct_product_count; iquad++)
            {
                // set current quadrature point
                if (num_dim == 3)
                {
                    z_quad = iquad / (num_gauss_points * num_gauss_points);
                }
                y_quad = (iquad % (num_gauss_points * num_gauss_points)) / num_gauss_points;
                x_quad = iquad % num_gauss_points;
                quad_coordinate(0) = legendre_nodes_1D(x_quad);
                quad_coordinate(1) = legendre_nodes_1D(y_quad);
                if (num_dim == 3)
                {
                    quad_coordinate(2) = legendre_nodes_1D(z_quad);
                }

                // set current quadrature weight
                quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
                quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
                if (num_dim == 3)
                {
                    quad_coordinate_weight(2) = legendre_weights_1D(z_quad);
                }
                else
                {
                    quad_coordinate_weight(2) = 1;
                }
                weight_multiply = quad_coordinate_weight(0) * quad_coordinate_weight(1) * quad_coordinate_weight(2);

                // compute shape functions at this point for the element type
                elem->basis(basis_values, quad_coordinate);

                // compute all the necessary coordinates and derivatives at this point

                // compute shape function derivatives
                elem->partial_xi_basis(basis_derivative_s1, quad_coordinate);
                elem->partial_eta_basis(basis_derivative_s2, quad_coordinate);
                elem->partial_mu_basis(basis_derivative_s3, quad_coordinate);

                // compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
                // derivative of x,y,z w.r.t s
                JT_row1(0) = 0;
                JT_row1(1) = 0;
                JT_row1(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row1(0) += nodal_positions(node_loop, 0) * basis_derivative_s1(node_loop);
                    JT_row1(1) += nodal_positions(node_loop, 1) * basis_derivative_s1(node_loop);
                    JT_row1(2) += nodal_positions(node_loop, 2) * basis_derivative_s1(node_loop);
                }

                // derivative of x,y,z w.r.t t
                JT_row2(0) = 0;
                JT_row2(1) = 0;
                JT_row2(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row2(0) += nodal_positions(node_loop, 0) * basis_derivative_s2(node_loop);
                    JT_row2(1) += nodal_positions(node_loop, 1) * basis_derivative_s2(node_loop);
                    JT_row2(2) += nodal_positions(node_loop, 2) * basis_derivative_s2(node_loop);
                }

                // derivative of x,y,z w.r.t w
                JT_row3(0) = 0;
                JT_row3(1) = 0;
                JT_row3(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row3(0) += nodal_positions(node_loop, 0) * basis_derivative_s3(node_loop);
                    JT_row3(1) += nodal_positions(node_loop, 1) * basis_derivative_s3(node_loop);
                    JT_row3(2) += nodal_positions(node_loop, 2) * basis_derivative_s3(node_loop);
                    // debug print
                    /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
                    std::cout << " ELEMENT VOLUME JACOBIAN DEBUG ON TASK " << myrank << std::endl;
                    std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
                    std::fflush(stdout);
                    }*/
                }

                // compute the determinant of the Jacobian
                Jacobian = JT_row1(0) * (JT_row2(1) * JT_row3(2) - JT_row3(1) * JT_row2(2)) -
                           JT_row1(1) * (JT_row2(0) * JT_row3(2) - JT_row3(0) * JT_row2(2)) +
                           JT_row1(2) * (JT_row2(0) * JT_row3(1) - JT_row3(0) * JT_row2(1));
                if (Jacobian < 0)
                {
                    Jacobian = -Jacobian;
                }

                // compute density
                current_density = 0;
                if (max_flag)
                {
                    current_density = 1;
                }
                else
                {
                    for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                    {
                        current_density += nodal_density(node_loop) * basis_values(node_loop);
                    }
                }

                Element_Masses(nonoverlapping_ielem, 0) += current_density * weight_multiply * Jacobian;
            }
        }
        else
        {
            Element_Masses(nonoverlapping_ielem, 0) = Element_Volumes(nonoverlapping_ielem, 0) * design_densities(nonoverlapping_ielem, 0);
        }
    }

    // std::ostream &out = std::cout;
    // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    // if(myrank==0)
    // *fos << "Global Element Masses:" << std::endl;
    // Global_Element_Masses->describe(*fos,Teuchos::VERB_EXTREME);
    // *fos << std::endl;
}

/* ----------------------------------------------------------------------
   Compute the gradients of mass function with respect to nodal densities
------------------------------------------------------------------------- */

void FEA_Module_Inertial::compute_nodal_gradients(const_host_vec_array design_variables, host_vec_array design_gradients, bool use_initial_coords){
    if(simparam->topology_optimization_on){
        compute_TO_gradients(design_variables, design_gradients, use_initial_coords);
    }
    else if(simparam->shape_optimization_on){
        compute_shape_gradients(design_variables, design_gradients, use_initial_coords);
    }

}

/* --------------------------------------------------------------------------------
   Compute the gradients of mass function with respect to nodal design coordinates
----------------------------------------------------------------------------------- */

void FEA_Module_Inertial::compute_shape_gradients(const_host_vec_array design_coords, host_vec_array design_gradients, bool use_initial_coords)
{

}

/* ----------------------------------------------------------------------
   Compute the gradients of mass function with respect to nodal densities
------------------------------------------------------------------------- */

void FEA_Module_Inertial::compute_TO_gradients(const_host_vec_array design_variables, host_vec_array design_gradients, bool use_initial_coords)
{
    // local number of uniquely assigned elements
    size_t nonoverlap_nelements = element_map->getLocalNumElements();
    // local variable for host view in the dual view
    const_host_vec_array all_node_coords;
    if (use_initial_coords)
    {
        all_node_coords = all_initial_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    else
    {
        all_node_coords = all_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    const_host_vec_array       all_node_densities;
    if (nodal_density_flag)
    {
        all_node_densities = all_node_densities_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    int    nodes_per_elem = elem->num_basis();
    int    z_quad, y_quad, x_quad, direct_product_count;
    size_t local_node_id;
    LO     ielem;
    GO     global_element_index;

    real_t Jacobian, weight_multiply;
    // CArrayKokkos<real_t> legendre_nodes_1D(num_gauss_points);
    // CArrayKokkos<real_t> legendre_weights_1D(num_gauss_points);
    CArray<real_t> legendre_nodes_1D(num_gauss_points);
    CArray<real_t> legendre_weights_1D(num_gauss_points);
    real_t pointer_quad_coordinate[num_dim];
    real_t pointer_quad_coordinate_weight[num_dim];
    real_t pointer_interpolated_point[num_dim];
    real_t pointer_JT_row1[num_dim];
    real_t pointer_JT_row2[num_dim];
    real_t pointer_JT_row3[num_dim];
    ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate, num_dim);
    ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight, num_dim);
    ViewCArray<real_t> interpolated_point(pointer_interpolated_point, num_dim);
    ViewCArray<real_t> JT_row1(pointer_JT_row1, num_dim);
    ViewCArray<real_t> JT_row2(pointer_JT_row2, num_dim);
    ViewCArray<real_t> JT_row3(pointer_JT_row3, num_dim);

    real_t pointer_basis_values[elem->num_basis()];
    real_t pointer_basis_derivative_s1[elem->num_basis()];
    real_t pointer_basis_derivative_s2[elem->num_basis()];
    real_t pointer_basis_derivative_s3[elem->num_basis()];
    ViewCArray<real_t> basis_values(pointer_basis_values, elem->num_basis());
    ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1, elem->num_basis());
    ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2, elem->num_basis());
    ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3, elem->num_basis());
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(), num_dim);
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());

    // initialize weights
    elements::legendre_nodes_1D(legendre_nodes_1D, num_gauss_points);
    elements::legendre_weights_1D(legendre_weights_1D, num_gauss_points);

    Solver::node_ordering_convention                             active_node_ordering_convention = Solver_Pointer_->active_node_ordering_convention;
    CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_node_order(max_nodes_per_element);
    if ((active_node_ordering_convention == Solver::ENSIGHT && num_dim == 3) || (active_node_ordering_convention == Solver::IJK && num_dim == 2))
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 3;
        convert_node_order(3) = 2;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 7;
            convert_node_order(7) = 6;
        }
    }
    else
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 2;
        convert_node_order(3) = 3;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 6;
            convert_node_order(7) = 7;
        }
    }

    // initialize design gradients to 0
    for (int init = 0; init < nlocal_nodes; init++)
    {
        design_gradients(init, 0) = 0;
    }

    // loop over elements and use quadrature rule to compute volume from Jacobian determinant
    for (int ielem = 0; ielem < rnum_elem; ielem++)
    {
        // acquire set of nodes for this local element
        for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
        {
            local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, convert_node_order(node_loop)));
            nodal_positions(node_loop, 0) = all_node_coords(local_node_id, 0);
            nodal_positions(node_loop, 1) = all_node_coords(local_node_id, 1);
            nodal_positions(node_loop, 2) = all_node_coords(local_node_id, 2);
            if (nodal_density_flag)
            {
                nodal_density(node_loop) = all_node_densities(local_node_id, 0);
            }
            /*
            if(myrank==1&&nodal_positions(node_loop,2)>10000000){
              std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
              std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
              std::fflush(stdout);
            }
            */
            // std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " " << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) <<std::endl;
        }

        direct_product_count = std::pow(num_gauss_points, num_dim);

        // loop over quadrature points
        for (int iquad = 0; iquad < direct_product_count; iquad++)
        {
            // set current quadrature point
            if (num_dim == 3)
            {
                z_quad = iquad / (num_gauss_points * num_gauss_points);
            }
            y_quad = (iquad % (num_gauss_points * num_gauss_points)) / num_gauss_points;
            x_quad = iquad % num_gauss_points;
            quad_coordinate(0) = legendre_nodes_1D(x_quad);
            quad_coordinate(1) = legendre_nodes_1D(y_quad);
            if (num_dim == 3)
            {
                quad_coordinate(2) = legendre_nodes_1D(z_quad);
            }

            // set current quadrature weight
            quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
            quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
            if (num_dim == 3)
            {
                quad_coordinate_weight(2) = legendre_weights_1D(z_quad);
            }
            else
            {
                quad_coordinate_weight(2) = 1;
            }
            weight_multiply = quad_coordinate_weight(0) * quad_coordinate_weight(1) * quad_coordinate_weight(2);

            // compute shape functions at this point for the element type
            elem->basis(basis_values, quad_coordinate);

            // compute all the necessary coordinates and derivatives at this point

            // compute shape function derivatives
            elem->partial_xi_basis(basis_derivative_s1, quad_coordinate);
            elem->partial_eta_basis(basis_derivative_s2, quad_coordinate);
            elem->partial_mu_basis(basis_derivative_s3, quad_coordinate);

            // compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
            // derivative of x,y,z w.r.t s
            JT_row1(0) = 0;
            JT_row1(1) = 0;
            JT_row1(2) = 0;
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                JT_row1(0) += nodal_positions(node_loop, 0) * basis_derivative_s1(node_loop);
                JT_row1(1) += nodal_positions(node_loop, 1) * basis_derivative_s1(node_loop);
                JT_row1(2) += nodal_positions(node_loop, 2) * basis_derivative_s1(node_loop);
            }

            // derivative of x,y,z w.r.t t
            JT_row2(0) = 0;
            JT_row2(1) = 0;
            JT_row2(2) = 0;
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                JT_row2(0) += nodal_positions(node_loop, 0) * basis_derivative_s2(node_loop);
                JT_row2(1) += nodal_positions(node_loop, 1) * basis_derivative_s2(node_loop);
                JT_row2(2) += nodal_positions(node_loop, 2) * basis_derivative_s2(node_loop);
            }

            // derivative of x,y,z w.r.t w
            JT_row3(0) = 0;
            JT_row3(1) = 0;
            JT_row3(2) = 0;
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                JT_row3(0) += nodal_positions(node_loop, 0) * basis_derivative_s3(node_loop);
                JT_row3(1) += nodal_positions(node_loop, 1) * basis_derivative_s3(node_loop);
                JT_row3(2) += nodal_positions(node_loop, 2) * basis_derivative_s3(node_loop);
                // debug print
                /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
                std::cout << " ELEMENT VOLUME JACOBIAN DEBUG ON TASK " << myrank << std::endl;
                std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
                std::fflush(stdout);
                }*/
            }

            // compute the determinant of the Jacobian
            Jacobian = JT_row1(0) * (JT_row2(1) * JT_row3(2) - JT_row3(1) * JT_row2(2)) -
                       JT_row1(1) * (JT_row2(0) * JT_row3(2) - JT_row3(0) * JT_row2(2)) +
                       JT_row1(2) * (JT_row2(0) * JT_row3(1) - JT_row3(0) * JT_row2(1));
            if (Jacobian < 0)
            {
                Jacobian = -Jacobian;
            }

            // assign contribution to every local node this element has
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                if (map->isNodeGlobalElement(nodes_in_elem(ielem, node_loop)))
                {
                    local_node_id = map->getLocalElement(nodes_in_elem(ielem, node_loop));
                    design_gradients(local_node_id, 0) += weight_multiply * basis_values(node_loop) * Jacobian;
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------------
   Compute the moment of each element for a specified component; estimated with quadrature
--------------------------------------------------------------------------------------------------------------------------- */

void FEA_Module_Inertial::compute_element_moments(const_host_vec_array design_densities, bool max_flag, int moment_component, bool use_initial_coords)
{
    // local number of uniquely assigned elements
    size_t nonoverlap_nelements = element_map->getLocalNumElements();
    // initialize memory for volume storage
    host_vec_array Element_Masses = Global_Element_Masses->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    host_vec_array Element_Moments;

    if (moment_component == 0)
    {
        Element_Moments = Global_Element_Moments_x->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    }
    if (moment_component == 1)
    {
        Element_Moments = Global_Element_Moments_y->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    }
    if (moment_component == 2)
    {
        Element_Moments = Global_Element_Moments_z->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    }

    const_host_vec_array Element_Volumes = Global_Element_Volumes->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    // local variable for host view in the dual view
    const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    const_host_vec_array all_initial_node_coords;
    if (use_initial_coords)
    {
        all_initial_node_coords = all_initial_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    const_host_vec_array all_design_densities;
    if (nodal_density_flag)
    {
        all_design_densities = all_node_densities_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    const_host_elem_conn_array nodes_in_elem  = global_nodes_in_elem_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    int                        nodes_per_elem = elem->num_basis();
    int                        z_quad, y_quad, x_quad, direct_product_count;
    size_t                     local_node_id;
    LO                         ielem;
    GO                         global_element_index;

    Solver::node_ordering_convention                             active_node_ordering_convention = Solver_Pointer_->active_node_ordering_convention;
    CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_node_order(max_nodes_per_element);
    if ((active_node_ordering_convention == Solver::ENSIGHT && num_dim == 3) || (active_node_ordering_convention == Solver::IJK && num_dim == 2))
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 3;
        convert_node_order(3) = 2;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 7;
            convert_node_order(7) = 6;
        }
    }
    else
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 2;
        convert_node_order(3) = 3;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 6;
            convert_node_order(7) = 7;
        }
    }

    real_t Jacobian, current_density, weight_multiply;
    // CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
    // CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
    CArray<real_t>     legendre_nodes_1D(num_gauss_points);
    CArray<real_t>     legendre_weights_1D(num_gauss_points);
    real_t             pointer_quad_coordinate[num_dim];
    real_t             pointer_quad_coordinate_weight[num_dim];
    real_t             pointer_interpolated_point[num_dim];
    real_t             pointer_JT_row1[num_dim];
    real_t             pointer_JT_row2[num_dim];
    real_t             pointer_JT_row3[num_dim];
    ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate, num_dim);
    ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight, num_dim);
    ViewCArray<real_t> interpolated_point(pointer_interpolated_point, num_dim);
    ViewCArray<real_t> JT_row1(pointer_JT_row1, num_dim);
    ViewCArray<real_t> JT_row2(pointer_JT_row2, num_dim);
    ViewCArray<real_t> JT_row3(pointer_JT_row3, num_dim);

    real_t                                                         pointer_basis_values[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s1[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s2[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s3[elem->num_basis()];
    ViewCArray<real_t>                                             basis_values(pointer_basis_values, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s1(pointer_basis_derivative_s1, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s2(pointer_basis_derivative_s2, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s3(pointer_basis_derivative_s3, elem->num_basis());
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(), num_dim);
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> initial_nodal_positions(elem->num_basis(), num_dim);
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_position(num_dim);

    // initialize weights
    elements::legendre_nodes_1D(legendre_nodes_1D, num_gauss_points);
    elements::legendre_weights_1D(legendre_weights_1D, num_gauss_points);

    // loop over elements and use quadrature rule to compute volume from Jacobian determinant
    for (int nonoverlapping_ielem = 0; nonoverlapping_ielem < nonoverlap_nelements; nonoverlapping_ielem++)
    {
        global_element_index = element_map->getGlobalElement(nonoverlapping_ielem);
        ielem = all_element_map->getLocalElement(global_element_index);
        // acquire set of nodes for this local element
        for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
        {
            local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, convert_node_order(node_loop)));
            nodal_positions(node_loop, 0) = all_node_coords(local_node_id, 0);
            nodal_positions(node_loop, 1) = all_node_coords(local_node_id, 1);
            nodal_positions(node_loop, 2) = all_node_coords(local_node_id, 2);
            if (use_initial_coords)
            {
                initial_nodal_positions(node_loop, 0) = all_initial_node_coords(local_node_id, 0);
                initial_nodal_positions(node_loop, 1) = all_initial_node_coords(local_node_id, 1);
                initial_nodal_positions(node_loop, 2) = all_initial_node_coords(local_node_id, 2);
            }

            if (nodal_density_flag)
            {
                nodal_density(node_loop) = all_design_densities(local_node_id, 0);
            }
            /*
            if(myrank==1&&nodal_positions(node_loop,2)>10000000){
              std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
              std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
              std::fflush(stdout);
            }
            */
            // std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " "
            // << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) << " " << nodal_density(node_loop) <<std::endl;
        }

        // debug print of index
        // std::cout << "nonoverlap element id on TASK " << myrank << " is " << nonoverlapping_ielem << std::endl;
        // std::fflush(stdout);

        // initialize element mass
        Element_Moments(nonoverlapping_ielem, 0) = 0;

        direct_product_count = std::pow(num_gauss_points, num_dim);

        // loop over quadrature points
        for (int iquad = 0; iquad < direct_product_count; iquad++)
        {
            // set current quadrature point
            if (num_dim == 3)
            {
                z_quad = iquad / (num_gauss_points * num_gauss_points);
            }
            y_quad = (iquad % (num_gauss_points * num_gauss_points)) / num_gauss_points;
            x_quad = iquad % num_gauss_points;
            quad_coordinate(0) = legendre_nodes_1D(x_quad);
            quad_coordinate(1) = legendre_nodes_1D(y_quad);
            if (num_dim == 3)
            {
                quad_coordinate(2) = legendre_nodes_1D(z_quad);
            }

            // set current quadrature weight
            quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
            quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
            if (num_dim == 3)
            {
                quad_coordinate_weight(2) = legendre_weights_1D(z_quad);
            }
            else
            {
                quad_coordinate_weight(2) = 1;
            }
            weight_multiply = quad_coordinate_weight(0) * quad_coordinate_weight(1) * quad_coordinate_weight(2);

            // compute shape functions at this point for the element type
            elem->basis(basis_values, quad_coordinate);

            // compute all the necessary coordinates and derivatives at this point

            // compute shape function derivatives
            elem->partial_xi_basis(basis_derivative_s1, quad_coordinate);
            elem->partial_eta_basis(basis_derivative_s2, quad_coordinate);
            elem->partial_mu_basis(basis_derivative_s3, quad_coordinate);

            // compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
            // derivative of x,y,z w.r.t s
            JT_row1(0) = 0;
            JT_row1(1) = 0;
            JT_row1(2) = 0;
            if (use_initial_coords)
            {
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row1(0) += initial_nodal_positions(node_loop, 0) * basis_derivative_s1(node_loop);
                    JT_row1(1) += initial_nodal_positions(node_loop, 1) * basis_derivative_s1(node_loop);
                    JT_row1(2) += initial_nodal_positions(node_loop, 2) * basis_derivative_s1(node_loop);
                }

                // derivative of x,y,z w.r.t t
                JT_row2(0) = 0;
                JT_row2(1) = 0;
                JT_row2(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row2(0) += initial_nodal_positions(node_loop, 0) * basis_derivative_s2(node_loop);
                    JT_row2(1) += initial_nodal_positions(node_loop, 1) * basis_derivative_s2(node_loop);
                    JT_row2(2) += initial_nodal_positions(node_loop, 2) * basis_derivative_s2(node_loop);
                }

                // derivative of x,y,z w.r.t w
                JT_row3(0) = 0;
                JT_row3(1) = 0;
                JT_row3(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row3(0) += initial_nodal_positions(node_loop, 0) * basis_derivative_s3(node_loop);
                    JT_row3(1) += initial_nodal_positions(node_loop, 1) * basis_derivative_s3(node_loop);
                    JT_row3(2) += initial_nodal_positions(node_loop, 2) * basis_derivative_s3(node_loop);
                    // debug print
                    /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
                    std::cout << " ELEMENT VOLUME JACOBIAN DEBUG ON TASK " << myrank << std::endl;
                    std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
                    std::fflush(stdout);
                    }*/
                }
            }
            else
            {
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row1(0) += nodal_positions(node_loop, 0) * basis_derivative_s1(node_loop);
                    JT_row1(1) += nodal_positions(node_loop, 1) * basis_derivative_s1(node_loop);
                    JT_row1(2) += nodal_positions(node_loop, 2) * basis_derivative_s1(node_loop);
                }

                // derivative of x,y,z w.r.t t
                JT_row2(0) = 0;
                JT_row2(1) = 0;
                JT_row2(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row2(0) += nodal_positions(node_loop, 0) * basis_derivative_s2(node_loop);
                    JT_row2(1) += nodal_positions(node_loop, 1) * basis_derivative_s2(node_loop);
                    JT_row2(2) += nodal_positions(node_loop, 2) * basis_derivative_s2(node_loop);
                }

                // derivative of x,y,z w.r.t w
                JT_row3(0) = 0;
                JT_row3(1) = 0;
                JT_row3(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row3(0) += nodal_positions(node_loop, 0) * basis_derivative_s3(node_loop);
                    JT_row3(1) += nodal_positions(node_loop, 1) * basis_derivative_s3(node_loop);
                    JT_row3(2) += nodal_positions(node_loop, 2) * basis_derivative_s3(node_loop);
                    // debug print
                    /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
                    std::cout << " ELEMENT VOLUME JACOBIAN DEBUG ON TASK " << myrank << std::endl;
                    std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
                    std::fflush(stdout);
                    }*/
                }
            }

            // compute the determinant of the Jacobian
            Jacobian = JT_row1(0) * (JT_row2(1) * JT_row3(2) - JT_row3(1) * JT_row2(2)) -
                       JT_row1(1) * (JT_row2(0) * JT_row3(2) - JT_row3(0) * JT_row2(2)) +
                       JT_row1(2) * (JT_row2(0) * JT_row3(1) - JT_row3(0) * JT_row2(1));
            if (Jacobian < 0)
            {
                Jacobian = -Jacobian;
            }

            // compute density
            if (max_flag)
            {
                current_density = 1;
            }
            else
            {
                if (nodal_density_flag)
                {
                    current_density = 0;
                    for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                    {
                        current_density += nodal_density(node_loop) * basis_values(node_loop);
                    }
                } // if
                else
                {
                    current_density = design_densities(nonoverlapping_ielem, 0);
                }
            }

            // compute current position
            current_position(0) = current_position(1) = current_position(2) = 0;
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                current_position(moment_component) += nodal_positions(node_loop, moment_component) * basis_values(node_loop);
            }

            Element_Moments(nonoverlapping_ielem, 0) += current_density * current_position(moment_component) * weight_multiply * Jacobian;
        }
    }

    // std::ostream &out = std::cout;
    // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    // if(myrank==0)
    // *fos << "Global Element Masses:" << std::endl;
    // Global_Element_Masses->describe(*fos,Teuchos::VERB_EXTREME);
    // *fos << std::endl;
}

/* ---------------------------------------------------------------------------------------------------
   Compute the gradients of the specified moment component with respect to design densities
------------------------------------------------------------------------------------------------------ */

void FEA_Module_Inertial::compute_moment_gradients(const_host_vec_array design_variables, host_vec_array design_gradients, int moment_component, bool use_initial_coords)
{
    // local number of uniquely assigned elements
    size_t nonoverlap_nelements = element_map->getLocalNumElements();
    // local variable for host view in the dual view
    const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    const_host_vec_array all_initial_node_coords;
    if (use_initial_coords)
    {
        all_initial_node_coords = all_initial_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    const_host_vec_array       all_node_densities;
    if (nodal_density_flag)
    {
        all_node_densities = all_node_densities_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    int    nodes_per_elem = elem->num_basis();
    int    z_quad, y_quad, x_quad, direct_product_count;
    size_t local_node_id;
    LO     ielem;
    GO     global_element_index;

    Solver::node_ordering_convention                             active_node_ordering_convention = Solver_Pointer_->active_node_ordering_convention;
    CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_node_order(max_nodes_per_element);
    if ((active_node_ordering_convention == Solver::ENSIGHT && num_dim == 3) || (active_node_ordering_convention == Solver::IJK && num_dim == 2))
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 3;
        convert_node_order(3) = 2;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 7;
            convert_node_order(7) = 6;
        }
    }
    else
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 2;
        convert_node_order(3) = 3;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 6;
            convert_node_order(7) = 7;
        }
    }

    real_t Jacobian, weight_multiply;
    // CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
    // CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
    CArray<real_t>     legendre_nodes_1D(num_gauss_points);
    CArray<real_t>     legendre_weights_1D(num_gauss_points);
    real_t             pointer_quad_coordinate[num_dim];
    real_t             pointer_quad_coordinate_weight[num_dim];
    real_t             pointer_interpolated_point[num_dim];
    real_t             pointer_JT_row1[num_dim];
    real_t             pointer_JT_row2[num_dim];
    real_t             pointer_JT_row3[num_dim];
    ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate, num_dim);
    ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight, num_dim);
    ViewCArray<real_t> interpolated_point(pointer_interpolated_point, num_dim);
    ViewCArray<real_t> JT_row1(pointer_JT_row1, num_dim);
    ViewCArray<real_t> JT_row2(pointer_JT_row2, num_dim);
    ViewCArray<real_t> JT_row3(pointer_JT_row3, num_dim);

    real_t                                                         pointer_basis_values[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s1[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s2[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s3[elem->num_basis()];
    ViewCArray<real_t>                                             basis_values(pointer_basis_values, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s1(pointer_basis_derivative_s1, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s2(pointer_basis_derivative_s2, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s3(pointer_basis_derivative_s3, elem->num_basis());
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(), num_dim);
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> initial_nodal_positions(elem->num_basis(), num_dim);
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_position(num_dim);

    // initialize weights
    elements::legendre_nodes_1D(legendre_nodes_1D, num_gauss_points);
    elements::legendre_weights_1D(legendre_weights_1D, num_gauss_points);

    // initialize design gradients to 0
    for (int init = 0; init < nlocal_nodes; init++)
    {
        design_gradients(init, 0) = 0;
    }

    // loop over elements and use quadrature rule to compute volume from Jacobian determinant
    for (int ielem = 0; ielem < rnum_elem; ielem++)
    {
        // acquire set of nodes for this local element
        for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
        {
            local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, convert_node_order(node_loop)));
            nodal_positions(node_loop, 0) = all_node_coords(local_node_id, 0);
            nodal_positions(node_loop, 1) = all_node_coords(local_node_id, 1);
            nodal_positions(node_loop, 2) = all_node_coords(local_node_id, 2);
            if (use_initial_coords)
            {
                initial_nodal_positions(node_loop, 0) = all_initial_node_coords(local_node_id, 0);
                initial_nodal_positions(node_loop, 1) = all_initial_node_coords(local_node_id, 1);
                initial_nodal_positions(node_loop, 2) = all_initial_node_coords(local_node_id, 2);
            }
            if (nodal_density_flag)
            {
                nodal_density(node_loop) = all_node_densities(local_node_id, 0);
            }
            /*
            if(myrank==1&&nodal_positions(node_loop,2)>10000000){
              std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
              std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
              std::fflush(stdout);
            }
            */
            // std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " " << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) <<std::endl;
        }

        direct_product_count = std::pow(num_gauss_points, num_dim);

        // loop over quadrature points
        for (int iquad = 0; iquad < direct_product_count; iquad++)
        {
            // set current quadrature point
            if (num_dim == 3)
            {
                z_quad = iquad / (num_gauss_points * num_gauss_points);
            }
            y_quad = (iquad % (num_gauss_points * num_gauss_points)) / num_gauss_points;
            x_quad = iquad % num_gauss_points;
            quad_coordinate(0) = legendre_nodes_1D(x_quad);
            quad_coordinate(1) = legendre_nodes_1D(y_quad);
            if (num_dim == 3)
            {
                quad_coordinate(2) = legendre_nodes_1D(z_quad);
            }

            // set current quadrature weight
            quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
            quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
            if (num_dim == 3)
            {
                quad_coordinate_weight(2) = legendre_weights_1D(z_quad);
            }
            else
            {
                quad_coordinate_weight(2) = 1;
            }
            weight_multiply = quad_coordinate_weight(0) * quad_coordinate_weight(1) * quad_coordinate_weight(2);

            // compute shape functions at this point for the element type
            elem->basis(basis_values, quad_coordinate);

            // compute all the necessary coordinates and derivatives at this point

            // compute shape function derivatives
            elem->partial_xi_basis(basis_derivative_s1, quad_coordinate);
            elem->partial_eta_basis(basis_derivative_s2, quad_coordinate);
            elem->partial_mu_basis(basis_derivative_s3, quad_coordinate);

            // compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
            // derivative of x,y,z w.r.t s
            JT_row1(0) = 0;
            JT_row1(1) = 0;
            JT_row1(2) = 0;
            if (use_initial_coords)
            {
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row1(0) += initial_nodal_positions(node_loop, 0) * basis_derivative_s1(node_loop);
                    JT_row1(1) += initial_nodal_positions(node_loop, 1) * basis_derivative_s1(node_loop);
                    JT_row1(2) += initial_nodal_positions(node_loop, 2) * basis_derivative_s1(node_loop);
                }

                // derivative of x,y,z w.r.t t
                JT_row2(0) = 0;
                JT_row2(1) = 0;
                JT_row2(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row2(0) += initial_nodal_positions(node_loop, 0) * basis_derivative_s2(node_loop);
                    JT_row2(1) += initial_nodal_positions(node_loop, 1) * basis_derivative_s2(node_loop);
                    JT_row2(2) += initial_nodal_positions(node_loop, 2) * basis_derivative_s2(node_loop);
                }

                // derivative of x,y,z w.r.t w
                JT_row3(0) = 0;
                JT_row3(1) = 0;
                JT_row3(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row3(0) += initial_nodal_positions(node_loop, 0) * basis_derivative_s3(node_loop);
                    JT_row3(1) += initial_nodal_positions(node_loop, 1) * basis_derivative_s3(node_loop);
                    JT_row3(2) += initial_nodal_positions(node_loop, 2) * basis_derivative_s3(node_loop);
                }
            }
            else
            {
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row1(0) += nodal_positions(node_loop, 0) * basis_derivative_s1(node_loop);
                    JT_row1(1) += nodal_positions(node_loop, 1) * basis_derivative_s1(node_loop);
                    JT_row1(2) += nodal_positions(node_loop, 2) * basis_derivative_s1(node_loop);
                }

                // derivative of x,y,z w.r.t t
                JT_row2(0) = 0;
                JT_row2(1) = 0;
                JT_row2(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row2(0) += nodal_positions(node_loop, 0) * basis_derivative_s2(node_loop);
                    JT_row2(1) += nodal_positions(node_loop, 1) * basis_derivative_s2(node_loop);
                    JT_row2(2) += nodal_positions(node_loop, 2) * basis_derivative_s2(node_loop);
                }

                // derivative of x,y,z w.r.t w
                JT_row3(0) = 0;
                JT_row3(1) = 0;
                JT_row3(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row3(0) += nodal_positions(node_loop, 0) * basis_derivative_s3(node_loop);
                    JT_row3(1) += nodal_positions(node_loop, 1) * basis_derivative_s3(node_loop);
                    JT_row3(2) += nodal_positions(node_loop, 2) * basis_derivative_s3(node_loop);
                }
            }

            // compute the determinant of the Jacobian
            Jacobian = JT_row1(0) * (JT_row2(1) * JT_row3(2) - JT_row3(1) * JT_row2(2)) -
                       JT_row1(1) * (JT_row2(0) * JT_row3(2) - JT_row3(0) * JT_row2(2)) +
                       JT_row1(2) * (JT_row2(0) * JT_row3(1) - JT_row3(0) * JT_row2(1));
            if (Jacobian < 0)
            {
                Jacobian = -Jacobian;
            }

            // compute current position
            current_position(0) = current_position(1) = current_position(2) = 0;
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                current_position(moment_component) += nodal_positions(node_loop, moment_component) * basis_values(node_loop);
            }

            // assign contribution to every local node this element has
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                if (map->isNodeGlobalElement(nodes_in_elem(ielem, node_loop)))
                {
                    local_node_id = map->getLocalElement(nodes_in_elem(ielem, node_loop));
                    design_gradients(local_node_id, 0) += weight_multiply * basis_values(node_loop) * current_position(moment_component) * Jacobian;
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------------
   Compute the moment of inertia of each element for a specified component of the inertia tensor; estimated with quadrature
--------------------------------------------------------------------------------------------------------------------------- */

void FEA_Module_Inertial::compute_element_moments_of_inertia(const_host_vec_array design_densities, bool max_flag, int inertia_component, bool use_initial_coords)
{
    // local number of uniquely assigned elements
    size_t nonoverlap_nelements = element_map->getLocalNumElements();
    // initialize memory for volume storage
    host_vec_array Element_Masses = Global_Element_Masses->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    host_vec_array Element_Moments_of_Inertia;

    if (inertia_component == 0)
    {
        Element_Moments_of_Inertia = Global_Element_Moments_of_Inertia_xx->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    }
    if (inertia_component == 1)
    {
        Element_Moments_of_Inertia = Global_Element_Moments_of_Inertia_yy->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    }
    if (inertia_component == 2)
    {
        Element_Moments_of_Inertia = Global_Element_Moments_of_Inertia_zz->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    }
    if (inertia_component == 3)
    {
        Element_Moments_of_Inertia = Global_Element_Moments_of_Inertia_xy->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    }
    if (inertia_component == 4)
    {
        Element_Moments_of_Inertia = Global_Element_Moments_of_Inertia_xz->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    }
    if (inertia_component == 5)
    {
        Element_Moments_of_Inertia = Global_Element_Moments_of_Inertia_yz->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    }

    const_host_vec_array Element_Volumes = Global_Element_Volumes->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    // local variable for host view in the dual view
    const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    const_host_vec_array all_initial_node_coords;
    if (use_initial_coords)
    {
        all_initial_node_coords = all_initial_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    const_host_vec_array all_design_densities;
    if (nodal_density_flag)
    {
        all_design_densities = all_node_densities_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    double                     inertia_center[3];
    if (module_params->enable_inertia_center[0])
    {
        inertia_center[0] = module_params->moment_of_inertia_center[0];
    }
    else
    {
        inertia_center[0] = center_of_mass[0];
    }
    if (module_params->enable_inertia_center[1])
    {
        inertia_center[1] = module_params->moment_of_inertia_center[1];
    }
    else
    {
        inertia_center[1] = center_of_mass[1];
    }
    if (num_dim == 3)
    {
        if (module_params->enable_inertia_center[2])
        {
            inertia_center[2] = module_params->moment_of_inertia_center[2];
        }
        else
        {
            inertia_center[2] = center_of_mass[2];
        }
    }
    int    nodes_per_elem = elem->num_basis();
    int    z_quad, y_quad, x_quad, direct_product_count;
    size_t local_node_id;
    LO     ielem;
    GO     global_element_index;
    real_t delx1, delx2;

    real_t Jacobian, current_density, weight_multiply;
    // CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
    // CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
    CArray<real_t>     legendre_nodes_1D(num_gauss_points);
    CArray<real_t>     legendre_weights_1D(num_gauss_points);
    real_t             pointer_quad_coordinate[num_dim];
    real_t             pointer_quad_coordinate_weight[num_dim];
    real_t             pointer_interpolated_point[num_dim];
    real_t             pointer_JT_row1[num_dim];
    real_t             pointer_JT_row2[num_dim];
    real_t             pointer_JT_row3[num_dim];
    ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate, num_dim);
    ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight, num_dim);
    ViewCArray<real_t> interpolated_point(pointer_interpolated_point, num_dim);
    ViewCArray<real_t> JT_row1(pointer_JT_row1, num_dim);
    ViewCArray<real_t> JT_row2(pointer_JT_row2, num_dim);
    ViewCArray<real_t> JT_row3(pointer_JT_row3, num_dim);

    real_t                                                         pointer_basis_values[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s1[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s2[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s3[elem->num_basis()];
    ViewCArray<real_t>                                             basis_values(pointer_basis_values, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s1(pointer_basis_derivative_s1, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s2(pointer_basis_derivative_s2, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s3(pointer_basis_derivative_s3, elem->num_basis());
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(), num_dim);
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> initial_nodal_positions(elem->num_basis(), num_dim);
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_position(num_dim);

    // initialize weights
    elements::legendre_nodes_1D(legendre_nodes_1D, num_gauss_points);
    elements::legendre_weights_1D(legendre_weights_1D, num_gauss_points);

    Solver::node_ordering_convention                             active_node_ordering_convention = Solver_Pointer_->active_node_ordering_convention;
    CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_node_order(max_nodes_per_element);
    if ((active_node_ordering_convention == Solver::ENSIGHT && num_dim == 3) || (active_node_ordering_convention == Solver::IJK && num_dim == 2))
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 3;
        convert_node_order(3) = 2;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 7;
            convert_node_order(7) = 6;
        }
    }
    else
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 2;
        convert_node_order(3) = 3;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 6;
            convert_node_order(7) = 7;
        }
    }

    // loop over elements and use quadrature rule to compute volume from Jacobian determinant
    for (int nonoverlapping_ielem = 0; nonoverlapping_ielem < nonoverlap_nelements; nonoverlapping_ielem++)
    {
        global_element_index = element_map->getGlobalElement(nonoverlapping_ielem);
        ielem = all_element_map->getLocalElement(global_element_index);
        // acquire set of nodes for this local element
        for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
        {
            local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, convert_node_order(node_loop)));
            nodal_positions(node_loop, 0) = all_node_coords(local_node_id, 0);
            nodal_positions(node_loop, 1) = all_node_coords(local_node_id, 1);
            nodal_positions(node_loop, 2) = all_node_coords(local_node_id, 2);
            if (use_initial_coords)
            {
                initial_nodal_positions(node_loop, 0) = all_initial_node_coords(local_node_id, 0);
                initial_nodal_positions(node_loop, 1) = all_initial_node_coords(local_node_id, 1);
                initial_nodal_positions(node_loop, 2) = all_initial_node_coords(local_node_id, 2);
            }
            if (nodal_density_flag)
            {
                nodal_density(node_loop) = all_design_densities(local_node_id, 0);
            }
            /*
            if(myrank==1&&nodal_positions(node_loop,2)>10000000){
              std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
              std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
              std::fflush(stdout);
            }
            */
            // std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " "
            // << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) << " " << nodal_density(node_loop) <<std::endl;
        }

        // debug print of index
        // std::cout << "nonoverlap element id on TASK " << myrank << " is " << nonoverlapping_ielem << std::endl;
        // std::fflush(stdout);

        // initialize element mass
        Element_Moments_of_Inertia(nonoverlapping_ielem, 0) = 0;

        direct_product_count = std::pow(num_gauss_points, num_dim);

        // loop over quadrature points
        for (int iquad = 0; iquad < direct_product_count; iquad++)
        {
            // set current quadrature point
            if (num_dim == 3)
            {
                z_quad = iquad / (num_gauss_points * num_gauss_points);
            }
            y_quad = (iquad % (num_gauss_points * num_gauss_points)) / num_gauss_points;
            x_quad = iquad % num_gauss_points;
            quad_coordinate(0) = legendre_nodes_1D(x_quad);
            quad_coordinate(1) = legendre_nodes_1D(y_quad);
            if (num_dim == 3)
            {
                quad_coordinate(2) = legendre_nodes_1D(z_quad);
            }

            // set current quadrature weight
            quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
            quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
            if (num_dim == 3)
            {
                quad_coordinate_weight(2) = legendre_weights_1D(z_quad);
            }
            else
            {
                quad_coordinate_weight(2) = 1;
            }
            weight_multiply = quad_coordinate_weight(0) * quad_coordinate_weight(1) * quad_coordinate_weight(2);

            // compute shape functions at this point for the element type
            elem->basis(basis_values, quad_coordinate);

            // compute all the necessary coordinates and derivatives at this point

            // compute shape function derivatives
            elem->partial_xi_basis(basis_derivative_s1, quad_coordinate);
            elem->partial_eta_basis(basis_derivative_s2, quad_coordinate);
            elem->partial_mu_basis(basis_derivative_s3, quad_coordinate);

            // compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
            // derivative of x,y,z w.r.t s
            JT_row1(0) = 0;
            JT_row1(1) = 0;
            JT_row1(2) = 0;
            if (use_initial_coords)
            {
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row1(0) += initial_nodal_positions(node_loop, 0) * basis_derivative_s1(node_loop);
                    JT_row1(1) += initial_nodal_positions(node_loop, 1) * basis_derivative_s1(node_loop);
                    JT_row1(2) += initial_nodal_positions(node_loop, 2) * basis_derivative_s1(node_loop);
                }

                // derivative of x,y,z w.r.t t
                JT_row2(0) = 0;
                JT_row2(1) = 0;
                JT_row2(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row2(0) += initial_nodal_positions(node_loop, 0) * basis_derivative_s2(node_loop);
                    JT_row2(1) += initial_nodal_positions(node_loop, 1) * basis_derivative_s2(node_loop);
                    JT_row2(2) += initial_nodal_positions(node_loop, 2) * basis_derivative_s2(node_loop);
                }

                // derivative of x,y,z w.r.t w
                JT_row3(0) = 0;
                JT_row3(1) = 0;
                JT_row3(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row3(0) += initial_nodal_positions(node_loop, 0) * basis_derivative_s3(node_loop);
                    JT_row3(1) += initial_nodal_positions(node_loop, 1) * basis_derivative_s3(node_loop);
                    JT_row3(2) += initial_nodal_positions(node_loop, 2) * basis_derivative_s3(node_loop);
                }
            }
            else
            {
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row1(0) += nodal_positions(node_loop, 0) * basis_derivative_s1(node_loop);
                    JT_row1(1) += nodal_positions(node_loop, 1) * basis_derivative_s1(node_loop);
                    JT_row1(2) += nodal_positions(node_loop, 2) * basis_derivative_s1(node_loop);
                }

                // derivative of x,y,z w.r.t t
                JT_row2(0) = 0;
                JT_row2(1) = 0;
                JT_row2(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row2(0) += nodal_positions(node_loop, 0) * basis_derivative_s2(node_loop);
                    JT_row2(1) += nodal_positions(node_loop, 1) * basis_derivative_s2(node_loop);
                    JT_row2(2) += nodal_positions(node_loop, 2) * basis_derivative_s2(node_loop);
                }

                // derivative of x,y,z w.r.t w
                JT_row3(0) = 0;
                JT_row3(1) = 0;
                JT_row3(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row3(0) += nodal_positions(node_loop, 0) * basis_derivative_s3(node_loop);
                    JT_row3(1) += nodal_positions(node_loop, 1) * basis_derivative_s3(node_loop);
                    JT_row3(2) += nodal_positions(node_loop, 2) * basis_derivative_s3(node_loop);
                }
            }

            // compute the determinant of the Jacobian
            Jacobian = JT_row1(0) * (JT_row2(1) * JT_row3(2) - JT_row3(1) * JT_row2(2)) -
                       JT_row1(1) * (JT_row2(0) * JT_row3(2) - JT_row3(0) * JT_row2(2)) +
                       JT_row1(2) * (JT_row2(0) * JT_row3(1) - JT_row3(0) * JT_row2(1));
            if (Jacobian < 0)
            {
                Jacobian = -Jacobian;
            }

            // compute density
            if (max_flag)
            {
                current_density = 1;
            }
            else
            {
                if (nodal_density_flag)
                {
                    current_density = 0;
                    for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                    {
                        current_density += nodal_density(node_loop) * basis_values(node_loop);
                    }
                } // if
                else
                {
                    current_density = design_densities(nonoverlapping_ielem, 0);
                }
            }

            // compute current position
            current_position(0) = current_position(1) = current_position(2) = 0;
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                current_position(0) += nodal_positions(node_loop, 0) * basis_values(node_loop);
                current_position(1) += nodal_positions(node_loop, 1) * basis_values(node_loop);
                current_position(2) += nodal_positions(node_loop, 2) * basis_values(node_loop);
            }

            if (inertia_component == 0)
            {
                delx1 = current_position(1) - inertia_center[1];
                delx2 = current_position(2) - inertia_center[2];
                Element_Moments_of_Inertia(nonoverlapping_ielem, 0) += current_density * (delx1 * delx1 + delx2 * delx2) * weight_multiply * Jacobian;
            }
            if (inertia_component == 1)
            {
                delx1 = current_position(0) - inertia_center[0];
                delx2 = current_position(2) - inertia_center[2];
                Element_Moments_of_Inertia(nonoverlapping_ielem, 0) += current_density * (delx1 * delx1 + delx2 * delx2) * weight_multiply * Jacobian;
            }
            if (inertia_component == 2)
            {
                delx1 = current_position(0) - inertia_center[0];
                delx2 = current_position(1) - inertia_center[1];
                Element_Moments_of_Inertia(nonoverlapping_ielem, 0) += current_density * (delx1 * delx1 + delx2 * delx2) * weight_multiply * Jacobian;
            }
            if (inertia_component == 3)
            {
                delx1 = current_position(0) - inertia_center[0];
                delx2 = current_position(1) - inertia_center[1];
                Element_Moments_of_Inertia(nonoverlapping_ielem, 0) -= current_density * (delx1 * delx2) * weight_multiply * Jacobian;
            }
            if (inertia_component == 4)
            {
                delx1 = current_position(0) - inertia_center[0];
                delx2 = current_position(2) - inertia_center[2];
                Element_Moments_of_Inertia(nonoverlapping_ielem, 0) -= current_density * (delx1 * delx2) * weight_multiply * Jacobian;
            }
            if (inertia_component == 5)
            {
                delx1 = current_position(1) - inertia_center[1];
                delx2 = current_position(2) - inertia_center[2];
                Element_Moments_of_Inertia(nonoverlapping_ielem, 0) -= current_density * (delx1 * delx2) * weight_multiply * Jacobian;
            }
        }
    }

    // std::ostream &out = std::cout;
    // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    // if(myrank==0)
    // *fos << "Global Element Masses:" << std::endl;
    // Global_Element_Masses->describe(*fos,Teuchos::VERB_EXTREME);
    // *fos << std::endl;
}

/* ---------------------------------------------------------------------------------------------------
   Compute the gradients of the specified moment of inertia component with respect to design densities
------------------------------------------------------------------------------------------------------ */

void FEA_Module_Inertial::compute_moment_of_inertia_gradients(const_host_vec_array design_variables, host_vec_array design_gradients, int inertia_component, bool use_initial_coords)
{
    // local number of uniquely assigned elements
    size_t nonoverlap_nelements = element_map->getLocalNumElements();
    // local variable for host view in the dual view
    const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    const_host_vec_array all_initial_node_coords;
    if (use_initial_coords)
    {
        all_initial_node_coords = all_initial_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    double                     inertia_center[3];

    Solver::node_ordering_convention                             active_node_ordering_convention = Solver_Pointer_->active_node_ordering_convention;
    CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_node_order(max_nodes_per_element);
    if ((active_node_ordering_convention == Solver::ENSIGHT && num_dim == 3) || (active_node_ordering_convention == Solver::IJK && num_dim == 2))
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 3;
        convert_node_order(3) = 2;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 7;
            convert_node_order(7) = 6;
        }
    }
    else
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 2;
        convert_node_order(3) = 3;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 6;
            convert_node_order(7) = 7;
        }
    }

    if (module_params->enable_inertia_center[0])
    {
        inertia_center[0] = module_params->moment_of_inertia_center[0];
    }
    else
    {
        inertia_center[0] = center_of_mass[0];
    }
    if (module_params->enable_inertia_center[1])
    {
        inertia_center[1] = module_params->moment_of_inertia_center[1];
    }
    else
    {
        inertia_center[1] = center_of_mass[1];
    }
    if (num_dim == 3)
    {
        if (module_params->enable_inertia_center[2])
        {
            inertia_center[2] = module_params->moment_of_inertia_center[2];
        }
        else
        {
            inertia_center[2] = center_of_mass[2];
        }
    }
    const_host_vec_array all_node_densities;
    if (nodal_density_flag)
    {
        all_node_densities = all_node_densities_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    }
    int    nodes_per_elem = elem->num_basis();
    int    z_quad, y_quad, x_quad, direct_product_count;
    size_t local_node_id;
    LO     ielem;
    GO     global_element_index;
    real_t delx1, delx2;

    real_t Jacobian, weight_multiply;
    // CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
    // CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
    CArray<real_t>     legendre_nodes_1D(num_gauss_points);
    CArray<real_t>     legendre_weights_1D(num_gauss_points);
    real_t             pointer_quad_coordinate[num_dim];
    real_t             pointer_quad_coordinate_weight[num_dim];
    real_t             pointer_interpolated_point[num_dim];
    real_t             pointer_JT_row1[num_dim];
    real_t             pointer_JT_row2[num_dim];
    real_t             pointer_JT_row3[num_dim];
    ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate, num_dim);
    ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight, num_dim);
    ViewCArray<real_t> interpolated_point(pointer_interpolated_point, num_dim);
    ViewCArray<real_t> JT_row1(pointer_JT_row1, num_dim);
    ViewCArray<real_t> JT_row2(pointer_JT_row2, num_dim);
    ViewCArray<real_t> JT_row3(pointer_JT_row3, num_dim);

    real_t                                                         pointer_basis_values[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s1[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s2[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s3[elem->num_basis()];
    ViewCArray<real_t>                                             basis_values(pointer_basis_values, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s1(pointer_basis_derivative_s1, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s2(pointer_basis_derivative_s2, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s3(pointer_basis_derivative_s3, elem->num_basis());
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(), num_dim);
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> initial_nodal_positions(elem->num_basis(), num_dim);
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_position(num_dim);

    // initialize weights
    elements::legendre_nodes_1D(legendre_nodes_1D, num_gauss_points);
    elements::legendre_weights_1D(legendre_weights_1D, num_gauss_points);

    // initialize design gradients to 0
    for (int init = 0; init < nlocal_nodes; init++)
    {
        design_gradients(init, 0) = 0;
    }

    // loop over elements and use quadrature rule to compute volume from Jacobian determinant
    for (int ielem = 0; ielem < rnum_elem; ielem++)
    {
        // acquire set of nodes for this local element
        for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
        {
            local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, convert_node_order(node_loop)));
            nodal_positions(node_loop, 0) = all_node_coords(local_node_id, 0);
            nodal_positions(node_loop, 1) = all_node_coords(local_node_id, 1);
            nodal_positions(node_loop, 2) = all_node_coords(local_node_id, 2);
            if (use_initial_coords)
            {
                initial_nodal_positions(node_loop, 0) = all_initial_node_coords(local_node_id, 0);
                initial_nodal_positions(node_loop, 1) = all_initial_node_coords(local_node_id, 1);
                initial_nodal_positions(node_loop, 2) = all_initial_node_coords(local_node_id, 2);
            }
            if (nodal_density_flag)
            {
                nodal_density(node_loop) = all_node_densities(local_node_id, 0);
            }
            /*
            if(myrank==1&&nodal_positions(node_loop,2)>10000000){
              std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
              std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
              std::fflush(stdout);
            }
            */
            // std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " " << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) <<std::endl;
        }

        direct_product_count = std::pow(num_gauss_points, num_dim);

        // loop over quadrature points
        for (int iquad = 0; iquad < direct_product_count; iquad++)
        {
            // set current quadrature point
            if (num_dim == 3)
            {
                z_quad = iquad / (num_gauss_points * num_gauss_points);
            }
            y_quad = (iquad % (num_gauss_points * num_gauss_points)) / num_gauss_points;
            x_quad = iquad % num_gauss_points;
            quad_coordinate(0) = legendre_nodes_1D(x_quad);
            quad_coordinate(1) = legendre_nodes_1D(y_quad);
            if (num_dim == 3)
            {
                quad_coordinate(2) = legendre_nodes_1D(z_quad);
            }

            // set current quadrature weight
            quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
            quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
            if (num_dim == 3)
            {
                quad_coordinate_weight(2) = legendre_weights_1D(z_quad);
            }
            else
            {
                quad_coordinate_weight(2) = 1;
            }
            weight_multiply = quad_coordinate_weight(0) * quad_coordinate_weight(1) * quad_coordinate_weight(2);

            // compute shape functions at this point for the element type
            elem->basis(basis_values, quad_coordinate);

            // compute all the necessary coordinates and derivatives at this point

            // compute shape function derivatives
            elem->partial_xi_basis(basis_derivative_s1, quad_coordinate);
            elem->partial_eta_basis(basis_derivative_s2, quad_coordinate);
            elem->partial_mu_basis(basis_derivative_s3, quad_coordinate);

            // compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
            // derivative of x,y,z w.r.t s
            JT_row1(0) = 0;
            JT_row1(1) = 0;
            JT_row1(2) = 0;
            if (use_initial_coords)
            {
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row1(0) += initial_nodal_positions(node_loop, 0) * basis_derivative_s1(node_loop);
                    JT_row1(1) += initial_nodal_positions(node_loop, 1) * basis_derivative_s1(node_loop);
                    JT_row1(2) += initial_nodal_positions(node_loop, 2) * basis_derivative_s1(node_loop);
                }

                // derivative of x,y,z w.r.t t
                JT_row2(0) = 0;
                JT_row2(1) = 0;
                JT_row2(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row2(0) += initial_nodal_positions(node_loop, 0) * basis_derivative_s2(node_loop);
                    JT_row2(1) += initial_nodal_positions(node_loop, 1) * basis_derivative_s2(node_loop);
                    JT_row2(2) += initial_nodal_positions(node_loop, 2) * basis_derivative_s2(node_loop);
                }

                // derivative of x,y,z w.r.t w
                JT_row3(0) = 0;
                JT_row3(1) = 0;
                JT_row3(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row3(0) += initial_nodal_positions(node_loop, 0) * basis_derivative_s3(node_loop);
                    JT_row3(1) += initial_nodal_positions(node_loop, 1) * basis_derivative_s3(node_loop);
                    JT_row3(2) += initial_nodal_positions(node_loop, 2) * basis_derivative_s3(node_loop);
                }
            }
            else
            {
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row1(0) += nodal_positions(node_loop, 0) * basis_derivative_s1(node_loop);
                    JT_row1(1) += nodal_positions(node_loop, 1) * basis_derivative_s1(node_loop);
                    JT_row1(2) += nodal_positions(node_loop, 2) * basis_derivative_s1(node_loop);
                }

                // derivative of x,y,z w.r.t t
                JT_row2(0) = 0;
                JT_row2(1) = 0;
                JT_row2(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row2(0) += nodal_positions(node_loop, 0) * basis_derivative_s2(node_loop);
                    JT_row2(1) += nodal_positions(node_loop, 1) * basis_derivative_s2(node_loop);
                    JT_row2(2) += nodal_positions(node_loop, 2) * basis_derivative_s2(node_loop);
                }

                // derivative of x,y,z w.r.t w
                JT_row3(0) = 0;
                JT_row3(1) = 0;
                JT_row3(2) = 0;
                for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
                {
                    JT_row3(0) += nodal_positions(node_loop, 0) * basis_derivative_s3(node_loop);
                    JT_row3(1) += nodal_positions(node_loop, 1) * basis_derivative_s3(node_loop);
                    JT_row3(2) += nodal_positions(node_loop, 2) * basis_derivative_s3(node_loop);
                }
            }

            // compute the determinant of the Jacobian
            Jacobian = JT_row1(0) * (JT_row2(1) * JT_row3(2) - JT_row3(1) * JT_row2(2)) -
                       JT_row1(1) * (JT_row2(0) * JT_row3(2) - JT_row3(0) * JT_row2(2)) +
                       JT_row1(2) * (JT_row2(0) * JT_row3(1) - JT_row3(0) * JT_row2(1));
            if (Jacobian < 0)
            {
                Jacobian = -Jacobian;
            }

            // compute current position
            current_position(0) = current_position(1) = current_position(2) = 0;
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                current_position(0) += nodal_positions(node_loop, 0) * basis_values(node_loop);
                current_position(1) += nodal_positions(node_loop, 1) * basis_values(node_loop);
                current_position(2) += nodal_positions(node_loop, 2) * basis_values(node_loop);
            }

            // assign contribution to every local node this element has
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                if (map->isNodeGlobalElement(nodes_in_elem(ielem, node_loop)))
                {
                    local_node_id = map->getLocalElement(nodes_in_elem(ielem, node_loop));
                    if (inertia_component == 0)
                    {
                        delx1 = current_position(1) - inertia_center[1];
                        delx2 = current_position(2) - inertia_center[2];
                        design_gradients(local_node_id, 0) += weight_multiply * basis_values(node_loop) * (delx1 * delx1 + delx2 * delx2) * Jacobian;
                    }
                    if (inertia_component == 1)
                    {
                        delx1 = current_position(0) - inertia_center[0];
                        delx2 = current_position(2) - inertia_center[2];
                        design_gradients(local_node_id, 0) += weight_multiply * basis_values(node_loop) * (delx1 * delx1 + delx2 * delx2) * Jacobian;
                    }
                    if (inertia_component == 2)
                    {
                        delx1 = current_position(0) - inertia_center[0];
                        delx2 = current_position(1) - inertia_center[1];
                        design_gradients(local_node_id, 0) += weight_multiply * basis_values(node_loop) * (delx1 * delx1 + delx2 * delx2) * Jacobian;
                    }
                    if (inertia_component == 3)
                    {
                        delx1 = current_position(0) - inertia_center[0];
                        delx2 = current_position(1) - inertia_center[1];
                        design_gradients(local_node_id, 0) -= weight_multiply * basis_values(node_loop) * (delx1 * delx2) * Jacobian;
                    }
                    if (inertia_component == 4)
                    {
                        delx1 = current_position(0) - inertia_center[0];
                        delx2 = current_position(2) - inertia_center[2];
                        design_gradients(local_node_id, 0) -= weight_multiply * basis_values(node_loop) * (delx1 * delx2) * Jacobian;
                    }
                    if (inertia_component == 5)
                    {
                        delx1 = current_position(1) - inertia_center[1];
                        delx2 = current_position(2) - inertia_center[2];
                        design_gradients(local_node_id, 0) -= weight_multiply * basis_values(node_loop) * (delx1 * delx2) * Jacobian;
                    }
                }
            }
        }
    }
}

/* ----------------------------------------------------------------------
   Compute the volume of each element; estimated with quadrature
------------------------------------------------------------------------- */

void FEA_Module_Inertial::compute_element_volumes()
{
    // local number of uniquely assigned elements
    size_t nonoverlap_nelements = element_map->getLocalNumElements();
    // local variable for host view in the dual view
    const_host_vec_array       all_node_coords = all_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    const_host_elem_conn_array nodes_in_elem   = global_nodes_in_elem_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    host_vec_array             Element_Volumes = Global_Element_Volumes->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
    int                        nodes_per_elem  = elem->num_basis();
    int                        z_quad, y_quad, x_quad, direct_product_count;
    size_t                     local_node_id;
    LO                         ielem;
    GO                         global_element_index;

    real_t Jacobian, weight_multiply;
    // CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
    // CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
    CArray<real_t>     legendre_nodes_1D(num_gauss_points);
    CArray<real_t>     legendre_weights_1D(num_gauss_points);
    real_t             pointer_quad_coordinate[num_dim];
    real_t             pointer_quad_coordinate_weight[num_dim];
    real_t             pointer_interpolated_point[num_dim];
    real_t             pointer_JT_row1[num_dim];
    real_t             pointer_JT_row2[num_dim];
    real_t             pointer_JT_row3[num_dim];
    ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate, num_dim);
    ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight, num_dim);
    ViewCArray<real_t> interpolated_point(pointer_interpolated_point, num_dim);
    ViewCArray<real_t> JT_row1(pointer_JT_row1, num_dim);
    ViewCArray<real_t> JT_row2(pointer_JT_row2, num_dim);
    ViewCArray<real_t> JT_row3(pointer_JT_row3, num_dim);

    real_t                                                         pointer_basis_values[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s1[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s2[elem->num_basis()];
    real_t                                                         pointer_basis_derivative_s3[elem->num_basis()];
    ViewCArray<real_t>                                             basis_values(pointer_basis_values, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s1(pointer_basis_derivative_s1, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s2(pointer_basis_derivative_s2, elem->num_basis());
    ViewCArray<real_t>                                             basis_derivative_s3(pointer_basis_derivative_s3, elem->num_basis());
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(), num_dim);

    // initialize weights
    elements::legendre_nodes_1D(legendre_nodes_1D, num_gauss_points);
    elements::legendre_weights_1D(legendre_weights_1D, num_gauss_points);

    Solver::node_ordering_convention                             active_node_ordering_convention = Solver_Pointer_->active_node_ordering_convention;
    CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_node_order(max_nodes_per_element);
    if ((active_node_ordering_convention == Solver::ENSIGHT && num_dim == 3) || (active_node_ordering_convention == Solver::IJK && num_dim == 2))
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 3;
        convert_node_order(3) = 2;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 7;
            convert_node_order(7) = 6;
        }
    }
    else
    {
        convert_node_order(0) = 0;
        convert_node_order(1) = 1;
        convert_node_order(2) = 2;
        convert_node_order(3) = 3;
        if (num_dim == 3)
        {
            convert_node_order(4) = 4;
            convert_node_order(5) = 5;
            convert_node_order(6) = 6;
            convert_node_order(7) = 7;
        }
    }

    // loop over elements and use quadrature rule to compute volume from Jacobian determinant
    for (int nonoverlapping_ielem = 0; nonoverlapping_ielem < nonoverlap_nelements; nonoverlapping_ielem++)
    {
        global_element_index = element_map->getGlobalElement(nonoverlapping_ielem);
        ielem = all_element_map->getLocalElement(global_element_index);
        // debug print
        // std::cout << "ELEMENT INDEX IS: " << ielem << " " <<global_element_index << std::endl;

        // acquire set of nodes for this local element
        for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
        {
            local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, convert_node_order(node_loop)));
            nodal_positions(node_loop, 0) = all_node_coords(local_node_id, 0);
            nodal_positions(node_loop, 1) = all_node_coords(local_node_id, 1);
            nodal_positions(node_loop, 2) = all_node_coords(local_node_id, 2);
            /*
            if(myrank==1&&nodal_positions(node_loop,2)>10000000){
              std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
              std::cout << node_loop+1 <<" " << local_node_id <<" "<< host_elem_conn_(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
              std::fflush(stdout);
            }
            */
            // std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " " << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) <<std::endl;
        }

        // initialize element volume
        Element_Volumes(nonoverlapping_ielem, 0) = 0;

        direct_product_count = std::pow(num_gauss_points, num_dim);

        // loop over quadrature points
        for (int iquad = 0; iquad < direct_product_count; iquad++)
        {
            // set current quadrature point
            if (num_dim == 3)
            {
                z_quad = iquad / (num_gauss_points * num_gauss_points);
            }
            y_quad = (iquad % (num_gauss_points * num_gauss_points)) / num_gauss_points;
            x_quad = iquad % num_gauss_points;
            quad_coordinate(0) = legendre_nodes_1D(x_quad);
            quad_coordinate(1) = legendre_nodes_1D(y_quad);
            if (num_dim == 3)
            {
                quad_coordinate(2) = legendre_nodes_1D(z_quad);
            }

            // set current quadrature weight
            quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
            quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
            if (num_dim == 3)
            {
                quad_coordinate_weight(2) = legendre_weights_1D(z_quad);
            }
            else
            {
                quad_coordinate_weight(2) = 1;
            }
            weight_multiply = quad_coordinate_weight(0) * quad_coordinate_weight(1) * quad_coordinate_weight(2);

            // compute shape functions at this point for the element type
            elem->basis(basis_values, quad_coordinate);

            // compute all the necessary coordinates and derivatives at this point

            // compute shape function derivatives
            elem->partial_xi_basis(basis_derivative_s1, quad_coordinate);
            elem->partial_eta_basis(basis_derivative_s2, quad_coordinate);
            elem->partial_mu_basis(basis_derivative_s3, quad_coordinate);

            // compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
            // derivative of x,y,z w.r.t s
            JT_row1(0) = 0;
            JT_row1(1) = 0;
            JT_row1(2) = 0;
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                JT_row1(0) += nodal_positions(node_loop, 0) * basis_derivative_s1(node_loop);
                JT_row1(1) += nodal_positions(node_loop, 1) * basis_derivative_s1(node_loop);
                JT_row1(2) += nodal_positions(node_loop, 2) * basis_derivative_s1(node_loop);
            }

            // derivative of x,y,z w.r.t t
            JT_row2(0) = 0;
            JT_row2(1) = 0;
            JT_row2(2) = 0;
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                JT_row2(0) += nodal_positions(node_loop, 0) * basis_derivative_s2(node_loop);
                JT_row2(1) += nodal_positions(node_loop, 1) * basis_derivative_s2(node_loop);
                JT_row2(2) += nodal_positions(node_loop, 2) * basis_derivative_s2(node_loop);
            }

            // derivative of x,y,z w.r.t w
            JT_row3(0) = 0;
            JT_row3(1) = 0;
            JT_row3(2) = 0;
            for (int node_loop = 0; node_loop < elem->num_basis(); node_loop++)
            {
                JT_row3(0) += nodal_positions(node_loop, 0) * basis_derivative_s3(node_loop);
                JT_row3(1) += nodal_positions(node_loop, 1) * basis_derivative_s3(node_loop);
                JT_row3(2) += nodal_positions(node_loop, 2) * basis_derivative_s3(node_loop);
                // debug print
                /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
                std::cout << " ELEMENT VOLUME JACOBIAN DEBUG ON TASK " << myrank << std::endl;
                std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
                std::fflush(stdout);
                }*/
            }

            // compute the determinant of the Jacobian
            Jacobian = JT_row1(0) * (JT_row2(1) * JT_row3(2) - JT_row3(1) * JT_row2(2)) -
                       JT_row1(1) * (JT_row2(0) * JT_row3(2) - JT_row3(0) * JT_row2(2)) +
                       JT_row1(2) * (JT_row2(0) * JT_row3(1) - JT_row3(0) * JT_row2(1));
            if (Jacobian < 0)
            {
                Jacobian = -Jacobian;
            }

            Element_Volumes(nonoverlapping_ielem, 0) += weight_multiply * Jacobian;
        }
    }

    std::ostream&                       out = std::cout;
    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    // if(myrank==0)
    // *fos << "Global Element Volumes:" << std::endl;
    // Global_Element_Volumes->describe(*fos,Teuchos::VERB_EXTREME);
    // *fos << std::endl;
}

/* -------------------------------------------------------------------------------------------
   Communicate ghosts using the current optimization design data
---------------------------------------------------------------------------------------------- */

void FEA_Module_Inertial::comm_variables(Teuchos::RCP<const MV> zp)
{
    comm_densities(zp);
}
