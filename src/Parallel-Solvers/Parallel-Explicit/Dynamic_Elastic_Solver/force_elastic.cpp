
#include <math.h>  // fmin, fmax, abs note: fminl is long
#include "mesh.h"
#include "state.h"
#include "FEA_Module_Dynamic_Elasticity.h"
#include "Tpetra_Import.hpp"
#include "Tpetra_Import_Util2.hpp"

// -----------------------------------------------------------------------------
// This function calculates the corner forces and the evolves stress (hypo)
//------------------------------------------------------------------------------
void FEA_Module_Dynamic_Elasticity::get_force_elastic(const DCArrayKokkos <material_t> &material,
                   const mesh_t &mesh,
                   const DViewCArrayKokkos <double> &node_coords,
                   const DViewCArrayKokkos <double> &node_vel,
                   const DViewCArrayKokkos <double> &node_mass,
                   const DViewCArrayKokkos <double> &elem_den,
                   const DViewCArrayKokkos <double> &elem_vol,
                   const DViewCArrayKokkos <double> &elem_div,
                   const DViewCArrayKokkos <size_t> &elem_mat_id,
                   DViewCArrayKokkos <double> &corner_force,
                   const double rk_alpha,
                   const size_t cycle
                   ){

    const size_t rk_level = simparam.dynamic_options.rk_num_bins - 1;    
    const size_t num_dim = mesh.num_dims;
    const real_t damping_constant = fea_params.damping_constant;
    const_vec_array all_initial_node_coords = all_initial_node_coords_distributed->getLocalView<device_type> (Tpetra::Access::ReadOnly);

    // walk over the nodes to update the velocity
    FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
        
        size_t dof_id;
        double node_force[3];
        for (size_t dim = 0; dim < num_dim; dim++){
            node_force[dim] = 0.0;
        } // end for dim
        
        // loop over all corners around the node and calculate the nodal force
        for (size_t corner_lid=0; corner_lid<num_corners_in_node(node_gid); corner_lid++){
        
            // Get corner gid
            size_t corner_gid = corners_in_node(node_gid, corner_lid);
            
            // loop over dimension
            for (size_t dim = 0; dim < num_dim; dim++){
                node_force[dim] += -damping_constant*node_vel(rk_level, node_gid, dim);
            } // end for dim
            
        } // end for corner_lid

        // loop over dimension
        for (size_t idim = 0; idim < num_dim; idim++){
            for(int idof = 0; idof < Stiffness_Matrix_Strides(node_gid*num_dim+idim%num_dim); idof++){
              dof_id = DOF_Graph_Matrix(node_gid*num_dim+idim%num_dim,idof);
              node_force[idim] += -(node_coords(rk_level, dof_id/num_dim, dof_id%num_dim)-all_initial_node_coords(dof_id/num_dim, dof_id%num_dim))*Stiffness_Matrix(node_gid*num_dim+idim%num_dim,idof);
            }
          
          //node_force[idim] += -0.001*(node_coords(rk_level, dof_id/num_dim, dof_id%num_dim)-all_initial_node_coords(dof_id/num_dim, dof_id%num_dim));
        } // end for dim
        
        // update the velocity
        for (int dim = 0; dim < num_dim; dim++){
            node_vel(rk_level, node_gid, dim) = node_vel(0, node_gid, dim) +
                                         rk_alpha * dt*node_force[dim]/node_mass(node_gid);
        } // end for dim
        
    }); // end for parallel for over nodes

    return;
    
} // end of routine

// -----------------------------------------------------------------------------
// This function calculates the corner forces and the evolves stress (hypo)
//------------------------------------------------------------------------------
void FEA_Module_Dynamic_Elasticity::applied_forces(const DCArrayKokkos <material_t> &material,
                   const mesh_t &mesh,
                   const DViewCArrayKokkos <double> &node_coords,
                   const DViewCArrayKokkos <double> &node_vel,
                   const DViewCArrayKokkos <double> &node_mass,
                   const DViewCArrayKokkos <double> &elem_den,
                   const DViewCArrayKokkos <double> &elem_vol,
                   const DViewCArrayKokkos <double> &elem_div,
                   const DViewCArrayKokkos <size_t> &elem_mat_id,
                   DViewCArrayKokkos <double> &corner_force,
                   const double rk_alpha,
                   const size_t cycle
                   ){

    const size_t rk_level = simparam.dynamic_options.rk_num_bins - 1;    
    const size_t num_dim = mesh.num_dims;
    const_vec_array all_initial_node_coords = all_initial_node_coords_distributed->getLocalView<device_type> (Tpetra::Access::ReadOnly);
    const size_t num_lcs = fea_params.loading.size();
    
    const DCArrayKokkos <mat_fill_t> mat_fill = simparam.mat_fill;
    const DCArrayKokkos <loading_t> loading = fea_params.loading;

    //debug check
    //std::cout << "NUMBER OF LOADING CONDITIONS: " << num_lcs << std::endl;

    // walk over the nodes to update the velocity
    FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
        double current_node_coords[3];
        size_t dof_id;
        double node_force[3];
        double applied_force[3];
        double radius;
        for (size_t dim = 0; dim < num_dim; dim++){
            node_force[dim] = 0.0;
            current_node_coords[dim] = all_initial_node_coords(node_gid, dim);
        } // end for dim
        radius = sqrt(current_node_coords[0]*current_node_coords[0]+current_node_coords[1]*current_node_coords[1]+current_node_coords[2]*current_node_coords[2]);
        for(size_t ilc=0; ilc < num_lcs; ilc++){
          //debug check
          //std::cout << "LOADING CONDITION VOLUME TYPE: " << to_string(loading(ilc).volume) << std::endl;

          bool fill_this = loading(ilc).contains(current_node_coords);
          if(fill_this){
            // loop over all corners around the node and calculate the nodal force
            for (size_t corner_lid=0; corner_lid<num_corners_in_node(node_gid); corner_lid++){
            
                // Get corner gid
                size_t corner_gid = corners_in_node(node_gid, corner_lid);
                applied_force[0] = loading(ilc).x;
                applied_force[1] = loading(ilc).y;
                applied_force[2] = loading(ilc).z;
                // loop over dimension
                for (size_t dim = 0; dim < num_dim; dim++){
                    node_force[dim] += applied_force[dim]*(all_initial_node_coords(node_gid, 0) + all_initial_node_coords(node_gid, 1) + all_initial_node_coords(node_gid, 2))/radius;
                } // end for dim
                
            } // end for corner_lid
            
            // update the velocity
            for (int dim = 0; dim < num_dim; dim++){
                node_vel(rk_level, node_gid, dim) +=
                                            rk_alpha * dt*node_force[dim]/node_mass(node_gid);
            } // end for dim
          }
        }
    }); // end for parallel for over nodes
    
    return;
    
} // end of routine

/* ----------------------------------------------------------------------
   Assemble the Sparse Stiffness Matrix
------------------------------------------------------------------------- */

void FEA_Module_Dynamic_Elasticity::assemble_matrix(){
  int num_dim = simparam.num_dims;
  int nodes_per_elem;
  int current_row_n_nodes_scanned;
  int local_dof_index, local_node_index, current_row, current_column;
  int max_stride = 0;
  
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Local_Stiffness_Matrix(num_dim*max_nodes_per_element,num_dim*max_nodes_per_element);

  //initialize stiffness Matrix entries to 0
  //debug print
    //std::cout << "DOF GRAPH MATRIX ENTRIES ON TASK " << myrank << std::endl;
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
    for (int istride = 0; istride < Stiffness_Matrix_Strides(idof); istride++){
      Stiffness_Matrix(idof,istride) = 0;
      //debug print
      //std::cout << "{" <<istride + 1 << "," << DOF_Graph_Matrix(idof,istride) << "} ";
    }
    //debug print
    //std::cout << std::endl;
  }

  //reset unsorted DOF Graph corresponding to assembly mapped values
  //debug print
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
    for (int istride = 0; istride < Stiffness_Matrix_Strides(idof); istride++){
      DOF_Graph_Matrix(idof,istride) = Graph_Matrix(idof/num_dim,istride/num_dim)*num_dim + istride%num_dim;
      //debug print
      //std::cout << "{" <<istride + 1 << "," << DOF_Graph_Matrix(idof,istride) << "} ";
    }
    //debug print
    //std::cout << std::endl;
  }

  //assemble the global stiffness matrix
  if(num_dim==2)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
    nodes_per_elem = elem2D->num_nodes();
    //construct local stiffness matrix for this element
    local_matrix_multiply(ielem, Local_Stiffness_Matrix);
    //assign entries of this local matrix to the sparse global matrix storage;
    for (int inode = 0; inode < nodes_per_elem; inode++){
      //see if this node is local
      local_node_index = nodes_in_elem.host(ielem,inode);
      if(!map->isNodeLocalElement(local_node_index)) continue;
      //set dof row start index
      current_row = num_dim*local_node_index;
      for(int jnode = 0; jnode < nodes_per_elem; jnode++){
        
        current_column = num_dim*Global_Stiffness_Matrix_Assembly_Map(ielem,inode,jnode);
        for (int idim = 0; idim < num_dim; idim++){
          for (int jdim = 0; jdim < num_dim; jdim++){

            //debug print
            //if(current_row + idim==15&&current_column + jdim==4)
            //std::cout << " Local stiffness matrix contribution for row " << current_row + idim +1 << " and column " << current_column + jdim + 1 << " : " <<
            //Local_Stiffness_Matrix(num_dim*inode + idim,num_dim*jnode + jdim) << " from " << ielem +1 << " i: " << num_dim*inode+idim+1 << " j: " << num_dim*jnode + jdim +1 << std::endl << std::endl;
            //end debug

            Stiffness_Matrix(current_row + idim, current_column + jdim) += Local_Stiffness_Matrix(num_dim*inode + idim,num_dim*jnode + jdim);
          }
        }
      }
    }
  }

  if(num_dim==3)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_3Delem_type(Element_Types(ielem), elem);
    nodes_per_elem = elem->num_nodes();
    //construct local stiffness matrix for this element
    local_matrix_multiply(ielem, Local_Stiffness_Matrix);
    //assign entries of this local matrix to the sparse global matrix storage;
    for (int inode = 0; inode < nodes_per_elem; inode++){
      //see if this node is local
      local_node_index = nodes_in_elem.host(ielem,inode);
      if(!map->isNodeLocalElement(local_node_index)) continue;
      //set dof row start index
      current_row = num_dim*local_node_index;
      for(int jnode = 0; jnode < nodes_per_elem; jnode++){
        
        current_column = num_dim*Global_Stiffness_Matrix_Assembly_Map(ielem,inode,jnode);
        for (int idim = 0; idim < num_dim; idim++){
          for (int jdim = 0; jdim < num_dim; jdim++){

            //debug print
            //if(current_row + idim==15&&current_column + jdim==4)
            //std::cout << " Local stiffness matrix contribution for row " << current_row + idim +1 << " and column " << current_column + jdim + 1 << " : " <<
            //Local_Stiffness_Matrix(num_dim*inode + idim,num_dim*jnode + jdim) << " from " << ielem +1 << " i: " << num_dim*inode+idim+1 << " j: " << num_dim*jnode + jdim +1 << std::endl << std::endl;
            //end debug

            Stiffness_Matrix(current_row + idim, current_column + jdim) += Local_Stiffness_Matrix(num_dim*inode + idim,num_dim*jnode + jdim);
          }
        }
      }
    }
  }


  //debug print of A matrix
  //*fos << "Global Stiffness Matrix :" << std::endl;
  //Global_Stiffness_Matrix->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;

  //filter small negative numbers (that should have been 0 from cancellation) from floating point error
  /*
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
    for (int istride = 0; istride < Stiffness_Matrix_Strides(idof); istride++){
      if(Stiffness_Matrix(idof,istride)<0.000000001*simparam.Elastic_Modulus*density_epsilon||Stiffness_Matrix(idof,istride)>-0.000000001*simparam.Elastic_Modulus*density_epsilon)
      Stiffness_Matrix(idof,istride) = 0;
      //debug print
      //std::cout << "{" <<istride + 1 << "," << DOF_Graph_Matrix(idof,istride) << "} ";
    }
    //debug print
    //std::cout << std::endl;
  }
  */
  

}

/* ----------------------------------------------------------------------
   Retrieve material properties associated with a finite element
------------------------------------------------------------------------- */

void FEA_Module_Dynamic_Elasticity::Element_Material_Properties(size_t ielem, real_t &Element_Modulus, real_t &Poisson_Ratio, real_t density){
  real_t unit_scaling = simparam.input_options.unit_scaling;
  real_t penalty_product = 1;
  real_t density_epsilon = simparam.optimization_options.value().density_epsilon;
  if(density < 0) density = 0;
  for(int i = 0; i < penalty_power; i++)
    penalty_product *= density;
  //relationship between density and stiffness
  Element_Modulus = (density_epsilon + (1 - density_epsilon)*penalty_product)*fea_params.Elastic_Modulus/unit_scaling/unit_scaling;
  //Element_Modulus = density*simparam.Elastic_Modulus/unit_scaling/unit_scaling;
  Poisson_Ratio = fea_params.Poisson_Ratio;
}

/* ----------------------------------------------------------------------
   Retrieve derivative of material properties with respect to local density
------------------------------------------------------------------------- */

void FEA_Module_Dynamic_Elasticity::Gradient_Element_Material_Properties(size_t ielem, real_t &Element_Modulus_Derivative, real_t &Poisson_Ratio, real_t density){
  real_t unit_scaling = simparam.input_options.unit_scaling;
  real_t penalty_product = 1;
  real_t density_epsilon = simparam.optimization_options.value().density_epsilon;
  Element_Modulus_Derivative = 0;
  if(density < 0) density = 0;
  for(int i = 0; i < penalty_power - 1; i++)
    penalty_product *= density;
  //relationship between density and stiffness
  Element_Modulus_Derivative = penalty_power*(1 - density_epsilon)*penalty_product*fea_params.Elastic_Modulus/unit_scaling/unit_scaling;
  //Element_Modulus_Derivative = simparam.Elastic_Modulus/unit_scaling/unit_scaling;
  Poisson_Ratio = fea_params.Poisson_Ratio;
}

/* ----------------------------------------------------------------------
   Construct the local stiffness matrix
------------------------------------------------------------------------- */

void FEA_Module_Dynamic_Elasticity::local_matrix_multiply(int ielem, CArrayKokkos<real_t, array_layout, device_type, memory_traits> &Local_Matrix){
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array Element_Densities;
  //local variable for host view of densities from the dual view
  
  const_host_vec_array all_node_densities;
  if(nodal_density_flag){
    if(simparam.optimization_options.value().helmholtz_filter)
      all_node_densities = all_filtered_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    else
      all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  }
  else{
    Element_Densities = Global_Element_Densities->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
  }
  int num_dim = simparam.num_dims;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam.num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id;

  direct_product_count = std::pow(num_gauss_points,num_dim);
  real_t Elastic_Constant, Shear_Term, Pressure_Term, matrix_term;
  real_t matrix_subterm1, matrix_subterm2, matrix_subterm3, invJacobian, Jacobian, weight_multiply;
  real_t Element_Modulus, Poisson_Ratio;
  //CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  //CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(),num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());

  size_t Brows;
  if(num_dim==2) Brows = 3;
  if(num_dim==3) Brows = 6;
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix(Brows,num_dim*elem->num_basis());
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> C_matrix(Brows,Brows);

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

  real_t current_density = 1;

  //acquire set of nodes for this local element
  for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
    local_node_id = nodes_in_elem.host(ielem, node_loop);
    nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
    nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
    nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
    if(nodal_density_flag) nodal_density(node_loop) = all_node_densities(local_node_id,0);
    /*
    if(myrank==1&&nodal_positions(node_loop,2)>10000000){
      std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
      std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
      std::fflush(stdout);
    }
    */
    //std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " " << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) <<std::endl;
  }

  //initialize C matrix
  for(int irow = 0; irow < Brows; irow++)
    for(int icol = 0; icol < Brows; icol++)
      C_matrix(irow,icol) = 0;

  //initialize local stiffness matrix storage
  for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++)
      for(int jfill=0; jfill < num_dim*nodes_per_elem; jfill++)
      Local_Matrix(ifill,jfill) = 0;

  //B matrix initialization
  for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++){
        CB_matrix(irow,icol) = B_matrix(irow,icol) = 0;
      }

  //loop over quadrature points
  for(int iquad=0; iquad < direct_product_count; iquad++){

    //set current quadrature point
    if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
    y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
    x_quad = iquad % num_gauss_points;
    quad_coordinate(0) = legendre_nodes_1D(x_quad);
    quad_coordinate(1) = legendre_nodes_1D(y_quad);
    if(num_dim==3)
    quad_coordinate(2) = legendre_nodes_1D(z_quad);

    //set current quadrature weight
    quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
    quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
    if(num_dim==3)
    quad_coordinate_weight(2) = legendre_weights_1D(z_quad);
    else
    quad_coordinate_weight(2) = 1;
    weight_multiply = quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2);

    //compute shape functions at this point for the element type
    elem->basis(basis_values,quad_coordinate);
    
    //compute density
    current_density = relative_element_densities.host(ielem);

    //debug print
    //std::cout << "Current Density " << current_density << std::endl;

    //look up element material properties at this point as a function of density
    Element_Material_Properties((size_t) ielem,Element_Modulus,Poisson_Ratio, current_density);
    Elastic_Constant = Element_Modulus/((1 + Poisson_Ratio)*(1 - 2*Poisson_Ratio));
    Shear_Term = 0.5-Poisson_Ratio;
    Pressure_Term = 1 - Poisson_Ratio;

    //debug print
    //std::cout << "Element Material Params " << Elastic_Constant << std::endl;

    //compute Elastic (C) matrix
    if(num_dim==2){
      C_matrix(0,0) = Pressure_Term;
      C_matrix(1,1) = Pressure_Term;
      C_matrix(0,1) = Poisson_Ratio;
      C_matrix(1,0) = Poisson_Ratio;
      C_matrix(2,2) = Shear_Term;
    }
    if(num_dim==3){
      C_matrix(0,0) = Pressure_Term;
      C_matrix(1,1) = Pressure_Term;
      C_matrix(2,2) = Pressure_Term;
      C_matrix(0,1) = Poisson_Ratio;
      C_matrix(0,2) = Poisson_Ratio;
      C_matrix(1,0) = Poisson_Ratio;
      C_matrix(1,2) = Poisson_Ratio;
      C_matrix(2,0) = Poisson_Ratio;
      C_matrix(2,1) = Poisson_Ratio;
      C_matrix(3,3) = Shear_Term;
      C_matrix(4,4) = Shear_Term;
      C_matrix(5,5) = Shear_Term;
    }
  
  /*
  //debug print of elasticity matrix
  std::cout << " ------------ELASTICITY MATRIX "<< ielem + 1 <<"--------------"<<std::endl;
  for (int idof = 0; idof < Brows; idof++){
    std::cout << "row: " << idof + 1 << " { ";
    for (int istride = 0; istride < Brows; istride++){
      std::cout << istride + 1 << " = " << C_matrix(idof,istride) << " , " ;
    }
    std::cout << " }"<< std::endl;
  }
  //end debug block
  */

    //compute all the necessary coordinates and derivatives at this point
    //compute shape function derivatives
    elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
    elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
    elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

    //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
    //derivative of x,y,z w.r.t s
    JT_row1(0) = 0;
    JT_row1(1) = 0;
    JT_row1(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
      JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
      JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
    }

    //derivative of x,y,z w.r.t t
    JT_row2(0) = 0;
    JT_row2(1) = 0;
    JT_row2(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
      JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
      JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
    }

    //derivative of x,y,z w.r.t w
    JT_row3(0) = 0;
    JT_row3(1) = 0;
    JT_row3(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
      JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
      JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
      //debug print
    /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
      std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
      std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
      std::fflush(stdout);
    }*/
    }
    
    
    //compute the determinant of the Jacobian
    Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
               JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
               JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
    if(Jacobian<0) Jacobian = -Jacobian;
    invJacobian = 1/Jacobian;
    //compute the contributions of this quadrature point to the B matrix
    if(num_dim==2)
    for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape*num_dim) = 0;
      B_matrix_contribution(2,ishape*num_dim) = 0;
      B_matrix_contribution(3,ishape*num_dim) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(5,ishape*num_dim) = 0;
      B_matrix_contribution(0,ishape*num_dim+1) = 0;
      B_matrix_contribution(1,ishape*num_dim+1) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(2,ishape*num_dim+1) = 0;
      B_matrix_contribution(3,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim+1) = 0;
      B_matrix_contribution(5,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(0,ishape*num_dim+2) = 0;
      B_matrix_contribution(1,ishape*num_dim+2) = 0;
      B_matrix_contribution(2,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(3,ishape*num_dim+2) = 0;
      B_matrix_contribution(4,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(5,ishape*num_dim+2) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
    }
    if(num_dim==3)
    for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape*num_dim) = 0;
      B_matrix_contribution(2,ishape*num_dim) = 0;
      B_matrix_contribution(3,ishape*num_dim) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(5,ishape*num_dim) = 0;
      B_matrix_contribution(0,ishape*num_dim+1) = 0;
      B_matrix_contribution(1,ishape*num_dim+1) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(2,ishape*num_dim+1) = 0;
      B_matrix_contribution(3,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim+1) = 0;
      B_matrix_contribution(5,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(0,ishape*num_dim+2) = 0;
      B_matrix_contribution(1,ishape*num_dim+2) = 0;
      B_matrix_contribution(2,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(3,ishape*num_dim+2) = 0;
      B_matrix_contribution(4,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(5,ishape*num_dim+2) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
    }
    /*
    //debug print of B matrix per quadrature point
    std::cout << " ------------B MATRIX QUADRATURE CONTRIBUTION"<< ielem + 1 <<"--------------"<<std::endl;
    for (int idof = 0; idof < Brows; idof++){
      std::cout << "row: " << idof + 1 << " { ";
      for (int istride = 0; istride < nodes_per_elem*num_dim; istride++){
        std::cout << istride + 1 << " = " << B_matrix_contribution(idof,istride) << " , " ;
      }
      std::cout << " }"<< std::endl;
    }
    //end debug block
    */
    //accumulate B matrix
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++)
      B_matrix(irow,icol) += B_matrix_contribution(irow,icol);

    //compute the previous multiplied by the Elastic (C) Matrix
    for(int irow=0; irow < Brows; irow++){
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++){
        CB_matrix_contribution(irow,icol) = 0;
        for(int span=0; span < Brows; span++){
          CB_matrix_contribution(irow,icol) += C_matrix(irow,span)*B_matrix_contribution(span,icol);
        }
      }
    }

    //accumulate CB matrix
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++)
      CB_matrix(irow,icol) += CB_matrix_contribution(irow,icol);

    //compute the contributions of this quadrature point to all the local stiffness matrix elements
    for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++)
      for(int jfill=ifill; jfill < num_dim*nodes_per_elem; jfill++){
        matrix_term = 0;
        for(int span = 0; span < Brows; span++){
          matrix_term += B_matrix_contribution(span,ifill)*CB_matrix_contribution(span,jfill);
        }
        Local_Matrix(ifill,jfill) += Elastic_Constant*weight_multiply*matrix_term*invJacobian;
        if(ifill!=jfill)
          Local_Matrix(jfill,ifill) = Local_Matrix(ifill,jfill);
      }
    
    }
    
    /*
    //debug print of B matrix per quadrature point
    std::cout << " ------------B MATRIX"<< ielem + 1 <<"--------------"<<std::endl;
    for (int idof = 0; idof < Brows; idof++){
      std::cout << "row: " << idof + 1 << " { ";
      for (int istride = 0; istride < nodes_per_elem*num_dim; istride++){
        std::cout << istride + 1 << " = " << B_matrix(idof,istride) << " , " ;
      }
      std::cout << " }"<< std::endl;
    }
    //end debug block

    //debug print of B matrix per quadrature point
    std::cout << " ------------CB MATRIX"<< ielem + 1 <<"--------------"<<std::endl;
    for (int idof = 0; idof < Brows; idof++){
      std::cout << "row: " << idof + 1 << " { ";
      for (int istride = 0; istride < nodes_per_elem*num_dim; istride++){
        std::cout << istride + 1 << " = " << CB_matrix(idof,istride) << " , " ;
      }
      std::cout << " }"<< std::endl;
    }
    //end debug block
    */

    //debug print of local stiffness matrix
      /*
      std::cout << " ------------LOCAL STIFFNESS MATRIX "<< ielem + 1 <<"--------------"<<std::endl;
      for (int idof = 0; idof < num_dim*nodes_per_elem; idof++){
        std::cout << "row: " << idof + 1 << " { ";
        for (int istride = 0; istride < num_dim*nodes_per_elem; istride++){
          std::cout << istride + 1 << " = " << Local_Matrix(idof,istride) << " , " ;
        }
        std::cout << " }"<< std::endl;
        }
      */
}

/* ----------------------------------------------------------------------
   Compute the gradient of strain energy with respect to nodal densities
------------------------------------------------------------------------- */

void FEA_Module_Dynamic_Elasticity::compute_stiffness_gradients(const_host_vec_array design_variables, host_vec_array design_gradients){
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  
  const_vec_array all_initial_node_coords = all_initial_node_coords_distributed->getLocalView<device_type> (Tpetra::Access::ReadOnly);
  const_host_vec_array Element_Densities;
  //local variable for host view of densities from the dual view
  const_host_vec_array all_node_densities;
  if(nodal_density_flag)
  all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  else
  Element_Densities = Global_Element_Densities->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
  int num_dim = simparam.num_dims;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam.num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id, local_dof_idx, local_dof_idy, local_dof_idz;
  const size_t rk_level = simparam.dynamic_options.rk_num_bins - 1;  
  GO current_global_index;
  
  real_t global_dt;
  bool element_constant_density = true;
  size_t current_data_index, next_data_index;
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_element_adjoint = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(num_nodes_in_elem*num_dim);

  direct_product_count = std::pow(num_gauss_points,num_dim);
  real_t Element_Modulus_Gradient, Poisson_Ratio, gradient_force_density[3];
  real_t Elastic_Constant, Shear_Term, Pressure_Term;
  real_t inner_product, matrix_term, Jacobian, invJacobian, weight_multiply;
  //CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  //CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(),num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_nodal_displacements(elem->num_basis()*num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());

  size_t Brows;
  if(num_dim==2) Brows = 3;
  if(num_dim==3) Brows = 6;
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix(Brows,num_dim*elem->num_basis());
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> C_matrix(Brows,Brows);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Local_Matrix_Contribution(num_dim*nodes_per_elem,num_dim*nodes_per_elem);

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);
  
  real_t current_density = 1;

  //loop through each element and assign the contribution to compliance gradient for each of its local nodes
  for (int cycle = 0; cycle < last_time_step+1; cycle++) {
    //compute timestep from time data
    global_dt = time_data[cycle+1] - time_data[cycle];

    //print
    if(simparam.dynamic_options.output_time_sequence_level==TIME_OUTPUT_LEVEL::extreme){
      if (cycle==0){
        if(myrank==0)
          printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_data[cycle], global_dt);
      }
          // print time step every 10 cycles
      else if (cycle%20==0){
        if(myrank==0)
          printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_data[cycle], global_dt);
      } // end if
    }

    //view scope
    {
      //const_vec_array current_velocity_vector = (*forward_solve_velocity_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      const_vec_array current_adjoint_vector = (*adjoint_vector_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      const_vec_array next_adjoint_vector = (*adjoint_vector_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);

      //set velocity and position for this timestep
      const_vec_array next_velocity_vector = (*forward_solve_velocity_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      const_vec_array current_velocity_vector = (*forward_solve_velocity_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      
      const_vec_array next_coordinate_vector = (*forward_solve_coordinate_data)[cycle+1]->getLocalView<device_type> (Tpetra::Access::ReadOnly);
      const_vec_array current_coordinate_vector = (*forward_solve_coordinate_data)[cycle]->getLocalView<device_type> (Tpetra::Access::ReadOnly);

      //interface of arrays for current implementation of force calculation

      FOR_ALL_CLASS(node_gid, 0, nlocal_nodes+nghost_nodes, {
        for (int idim = 0; idim < num_dim; idim++){
          node_vel(rk_level,node_gid,idim) = current_velocity_vector(node_gid,idim);
          node_coords(0,node_gid,idim) = current_coordinate_vector(node_gid,idim);
          node_coords(rk_level,node_gid,idim) = next_coordinate_vector(node_gid,idim);
        }
      });
      Kokkos::fence();

      for(size_t ielem = 0; ielem < rnum_elem; ielem++){
        nodes_per_elem = elem->num_basis();

        //initialize C matrix
        for(int irow = 0; irow < Brows; irow++)
          for(int icol = 0; icol < Brows; icol++)
            C_matrix(irow,icol) = 0;

        //B matrix initialization
        for(int irow=0; irow < Brows; irow++)
          for(int icol=0; icol < num_dim*nodes_per_elem; icol++){
            CB_matrix(irow,icol) = 0;
          }

        //acquire set of nodes and nodal displacements for this local element
        for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
          local_node_id = nodes_in_elem.host(ielem, node_loop);
          local_dof_idx = nodes_in_elem.host(ielem, node_loop)*num_dim;
          local_dof_idy = local_dof_idx + 1;
          local_dof_idz = local_dof_idx + 2;

          nodal_positions(node_loop,0) = all_initial_node_coords(local_node_id,0);
          nodal_positions(node_loop,1) = all_initial_node_coords(local_node_id,1);
          current_nodal_displacements(node_loop*num_dim) = (node_coords(rk_level, local_node_id, 0)-all_initial_node_coords(local_node_id,0)+node_coords(0, local_node_id, 0)-all_initial_node_coords(local_node_id,0))/2;
          current_nodal_displacements(node_loop*num_dim+1) = (node_coords(rk_level, local_node_id, 1)-all_initial_node_coords(local_node_id,1)+node_coords(0, local_node_id, 1)-all_initial_node_coords(local_node_id,1))/2;
          current_element_adjoint(node_loop*num_dim) = (current_adjoint_vector(local_node_id,0)+next_adjoint_vector(local_node_id,0))/2;
          current_element_adjoint(node_loop*num_dim+1) = (current_adjoint_vector(local_node_id,1)+next_adjoint_vector(local_node_id,1))/2;

          if(num_dim==3){
          nodal_positions(node_loop,2) = all_initial_node_coords(local_node_id,2);
          current_nodal_displacements(node_loop*num_dim+2) = (node_coords(rk_level, local_node_id, 2)-all_initial_node_coords(local_node_id,2)+node_coords(0, local_node_id, 2)-all_initial_node_coords(local_node_id,2))/2;
          current_element_adjoint(node_loop*num_dim+2) = (current_adjoint_vector(local_node_id,2)+next_adjoint_vector(local_node_id,2))/2;
          }

          if(nodal_density_flag) nodal_density(node_loop) = all_node_densities(local_node_id,0);
          //debug print
          /*
          std::cout << "node index access x "<< local_node_id << std::endl;
          std::cout << "local index access x "<< local_dof_idx << " displacement x " << current_nodal_displacements(node_loop*num_dim) <<std::endl;
          std::cout << "local index access y "<< local_dof_idy << " displacement y " << current_nodal_displacements(node_loop*num_dim + 1) << std::endl;
          std::cout << "local index access z "<< local_dof_idz << " displacement z " << current_nodal_displacements(node_loop*num_dim + 2) << std::endl; 
          */
        }

        //debug print of current_nodal_displacements
        /*
        std::cout << " ------------nodal displacements for Element "<< ielem + 1 <<"--------------"<<std::endl;
        std::cout << " { ";
        for (int idof = 0; idof < num_dim*nodes_per_elem; idof++){
          std::cout << idof + 1 << " = " << current_nodal_displacements(idof) << " , " ;
        }
        std::cout << " }"<< std::endl;
        */

        //loop over quadrature points
        for(int iquad=0; iquad < direct_product_count; iquad++){

          //set current quadrature point
          if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
          y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
          x_quad = iquad % num_gauss_points;
          quad_coordinate(0) = legendre_nodes_1D(x_quad);
          quad_coordinate(1) = legendre_nodes_1D(y_quad);
          if(num_dim==3)
          quad_coordinate(2) = legendre_nodes_1D(z_quad);

          //set current quadrature weight
          quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
          quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
          if(num_dim==3)
          quad_coordinate_weight(2) = legendre_weights_1D(z_quad);
          else
          quad_coordinate_weight(2) = 1;
          weight_multiply = quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2);

          //compute shape functions at this point for the element type
          elem->basis(basis_values,quad_coordinate);

          //compute all the necessary coordinates and derivatives at this point
          //compute shape function derivatives
          elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
          elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
          elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

          //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
          //derivative of x,y,z w.r.t s
          JT_row1(0) = 0;
          JT_row1(1) = 0;
          JT_row1(2) = 0;
          for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
            JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
            JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
            JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
          }

          //derivative of x,y,z w.r.t t
          JT_row2(0) = 0;
          JT_row2(1) = 0;
          JT_row2(2) = 0;
          for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
            JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
            JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
            JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
          }

          //derivative of x,y,z w.r.t w
          JT_row3(0) = 0;
          JT_row3(1) = 0;
          JT_row3(2) = 0;
          for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
            JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
            JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
            JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
          }
          
          //compute the determinant of the Jacobian
          Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                    JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                    JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
          if(Jacobian<0) Jacobian = -Jacobian;
          invJacobian = 1/Jacobian;

          //compute density
          current_density = relative_element_densities.host(ielem);

          //debug print
          //std::cout << "Current Density " << current_density << std::endl;

          //compute the contributions of this quadrature point to the B matrix
          if(num_dim==2)
          for(int ishape=0; ishape < nodes_per_elem; ishape++){
            B_matrix_contribution(0,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
                basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
            B_matrix_contribution(1,ishape*num_dim) = 0;
            B_matrix_contribution(2,ishape*num_dim) = 0;
            B_matrix_contribution(3,ishape*num_dim) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
                basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
            B_matrix_contribution(4,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
                basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
                basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
            B_matrix_contribution(5,ishape*num_dim) = 0;
            B_matrix_contribution(0,ishape*num_dim+1) = 0;
            B_matrix_contribution(1,ishape*num_dim+1) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
                basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
            B_matrix_contribution(2,ishape*num_dim+1) = 0;
            B_matrix_contribution(3,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
                basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
            B_matrix_contribution(4,ishape*num_dim+1) = 0;
            B_matrix_contribution(5,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
                basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
                basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
            B_matrix_contribution(0,ishape*num_dim+2) = 0;
            B_matrix_contribution(1,ishape*num_dim+2) = 0;
            B_matrix_contribution(2,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
                basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
                basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
            B_matrix_contribution(3,ishape*num_dim+2) = 0;
            B_matrix_contribution(4,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
                basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
            B_matrix_contribution(5,ishape*num_dim+2) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
                basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
          }
          if(num_dim==3)
          for(int ishape=0; ishape < nodes_per_elem; ishape++){
            B_matrix_contribution(0,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
                basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
            B_matrix_contribution(1,ishape*num_dim) = 0;
            B_matrix_contribution(2,ishape*num_dim) = 0;
            B_matrix_contribution(3,ishape*num_dim) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
                basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
            B_matrix_contribution(4,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
                basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
                basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
            B_matrix_contribution(5,ishape*num_dim) = 0;
            B_matrix_contribution(0,ishape*num_dim+1) = 0;
            B_matrix_contribution(1,ishape*num_dim+1) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
                basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
            B_matrix_contribution(2,ishape*num_dim+1) = 0;
            B_matrix_contribution(3,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
                basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
            B_matrix_contribution(4,ishape*num_dim+1) = 0;
            B_matrix_contribution(5,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
                basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
                basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
            B_matrix_contribution(0,ishape*num_dim+2) = 0;
            B_matrix_contribution(1,ishape*num_dim+2) = 0;
            B_matrix_contribution(2,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
                basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
                basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
            B_matrix_contribution(3,ishape*num_dim+2) = 0;
            B_matrix_contribution(4,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
                basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
            B_matrix_contribution(5,ishape*num_dim+2) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
                basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
          }
          
          //look up element material properties at this point as a function of density
          Gradient_Element_Material_Properties(ielem, Element_Modulus_Gradient, Poisson_Ratio, current_density);
          Elastic_Constant = Element_Modulus_Gradient/((1 + Poisson_Ratio)*(1 - 2*Poisson_Ratio));
          Shear_Term = 0.5 - Poisson_Ratio;
          Pressure_Term = 1 - Poisson_Ratio;

          //debug print
          //std::cout << "Element Material Params " << Elastic_Constant << std::endl;

          //compute Elastic (C) matrix
          if(num_dim==2){
            C_matrix(0,0) = Pressure_Term;
            C_matrix(1,1) = Pressure_Term;
            C_matrix(0,1) = Poisson_Ratio;
            C_matrix(1,0) = Poisson_Ratio;
            C_matrix(2,2) = Shear_Term;
          }
          if(num_dim==3){
            C_matrix(0,0) = Pressure_Term;
            C_matrix(1,1) = Pressure_Term;
            C_matrix(2,2) = Pressure_Term;
            C_matrix(0,1) = Poisson_Ratio;
            C_matrix(0,2) = Poisson_Ratio;
            C_matrix(1,0) = Poisson_Ratio;
            C_matrix(1,2) = Poisson_Ratio;
            C_matrix(2,0) = Poisson_Ratio;
            C_matrix(2,1) = Poisson_Ratio;
            C_matrix(3,3) = Shear_Term;
            C_matrix(4,4) = Shear_Term;
            C_matrix(5,5) = Shear_Term;
          }

          //compute the previous multiplied by the Elastic (C) Matrix
          for(int irow=0; irow < Brows; irow++){
            for(int icol=0; icol < num_dim*nodes_per_elem; icol++){
              CB_matrix_contribution(irow,icol) = 0;
              for(int span=0; span < Brows; span++){
                CB_matrix_contribution(irow,icol) += C_matrix(irow,span)*B_matrix_contribution(span,icol);
              }
            }
          }
          
          //compute the contributions of this quadrature point to all the local stiffness matrix elements
          for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++){
            for(int jfill=ifill; jfill < num_dim*nodes_per_elem; jfill++){
              matrix_term = 0;
              for(int span = 0; span < Brows; span++){
                matrix_term += B_matrix_contribution(span,ifill)*CB_matrix_contribution(span,jfill);
              }
              Local_Matrix_Contribution(ifill,jfill) = matrix_term;
              if(ifill!=jfill)
                Local_Matrix_Contribution(jfill,ifill) = Local_Matrix_Contribution(ifill,jfill);
            }
          }

          //compute inner product for this quadrature point contribution
          inner_product = 0;
          for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++){
            for(int jfill=ifill; jfill < num_dim*nodes_per_elem; jfill++){
              if(ifill==jfill)
                inner_product += Local_Matrix_Contribution(ifill, jfill)*current_element_adjoint(ifill)*current_nodal_displacements(jfill);
              else
                inner_product += 2*Local_Matrix_Contribution(ifill, jfill)*current_element_adjoint(ifill)*current_nodal_displacements(jfill);
              //debug
              //if(Local_Matrix_Contribution(ifill, jfill)<0) Local_Matrix_Contribution(ifill, jfill) = - Local_Matrix_Contribution(ifill, jfill);
              //inner_product += Local_Matrix_Contribution(ifill, jfill);
            }
          }

          //evaluate local stiffness matrix gradient with respect to igradient
          for(int igradient=0; igradient < nodes_per_elem; igradient++){
            if(!map->isNodeLocalElement(nodes_in_elem.host(ielem, igradient))) continue;
            local_node_id = nodes_in_elem.host(ielem, igradient);
            
            //debug print
            //std::cout << "contribution for " << igradient + 1 << " is " << inner_product << std::endl;
            design_gradients(local_node_id,0) -= -inner_product*Elastic_Constant*weight_multiply*invJacobian*global_dt/nodes_per_elem;
          }

          //evaluate gradient of body force (such as gravity which depends on density) with respect to igradient
          if(body_term_flag){
            //look up element material properties at this point as a function of density
            Gradient_Body_Term(ielem, current_density, gradient_force_density);
            for(int igradient=0; igradient < nodes_per_elem; igradient++){
              if(!map->isNodeLocalElement(nodes_in_elem.host(ielem, igradient))) continue;
              local_node_id = nodes_in_elem.host(ielem, igradient);
              
              //compute inner product for this quadrature point contribution
              inner_product = 0;
              for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++){
                inner_product += gradient_force_density[ifill%num_dim]*current_nodal_displacements(ifill)*basis_values(ifill/num_dim);
              }
              
              //debug print
              //std::cout << "contribution for " << igradient + 1 << " is " << inner_product << std::endl;
              design_gradients(local_node_id,0) += inner_product*weight_multiply*Jacobian*global_dt/nodes_per_elem;
            }
          }
        }
      }
    }
    //end view scope
  }
}
