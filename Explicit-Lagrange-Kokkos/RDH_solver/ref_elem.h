#ifndef REFELEM_H
#define REFELEM_H


#include <cmath>
#include "matar.h"

#define EPSILON 1.0e-12
using namespace mtr;

struct ref_elem_t{
    
    size_t num_dim;
    
    // Dofs
    size_t num_ref_dofs_1d;
    size_t num_ref_dofs_in_elem;
    
    // Gauss Points
    size_t num_gauss_lob_1d;
    size_t num_gauss_lob_in_elem;
   
    size_t num_gauss_leg_1d;
    size_t num_gauss_leg_in_elem;
    
    // Zones
    size_t num_zones_1d;
    size_t num_zones_in_elem;
    
    // Num basis functions
    size_t num_basis;
    
    // Basis evaluation at nodes
    CArrayKokkos <double> ref_gauss_lob_basis;
    CArrayKokkos <double> ref_gauss_leg_basis;
     
    // Gradient of basis
    CArrayKokkos <double> ref_gauss_lob_grad_basis;
    CArrayKokkos <double> ref_gauss_leg_grad_basis;
   
    // Gauss and DOF positions 
    CArrayKokkos <double> lob_nodes_1D;
    CArrayKokkos <double> leg_nodes_1D;
    CArrayKokkos <double> ref_gauss_lob_positions;
    CArrayKokkos <double> ref_gauss_leg_positions;
    CArrayKokkos <double> ref_dof_positions;
    CArrayKokkos <double> ref_dof_positions_1d;
    
    // Quadrature Weights
    CArrayKokkos <double> lob_weights_1D;
    CArrayKokkos <double> leg_weights_1D;
    CArrayKokkos <double> ref_gauss_lob_weights;
    CArrayKokkos <double> ref_gauss_leg_weights;
    

void init(int p_order, int num_dim_inp){ 
    
    num_dim = num_dim_inp;

    if(p_order == 0){       
        
        num_gauss_lob_1d = 2; // num gauss lobatto points in 1d
	num_gauss_leg_1d = 1;
        num_ref_dofs_1d = 2;
        num_zones_1d = 1;
        num_zones_in_elem = num_zones_1d*num_zones_1d*num_zones_1d;

    }

    else{
        
        num_gauss_lob_1d = 2 * p_order + 1; // num gauss lobatto points in 1d
	num_gauss_leg_1d = 2*p_order;

        num_ref_dofs_1d = p_order+1;
        num_zones_1d = (num_ref_dofs_1d - 1) / 2;
        num_zones_in_elem = num_zones_1d*num_zones_1d*num_zones_1d;

    }

    num_gauss_lob_in_elem = 1;
   
    num_gauss_leg_in_elem = 1;

    num_ref_dofs_in_elem = 1;
    
    for (int dim = 0; dim < num_dim; dim++){
    
      num_gauss_lob_in_elem *= num_gauss_lob_1d;    
      num_gauss_leg_in_elem *= num_gauss_leg_1d;

      num_ref_dofs_in_elem *= num_ref_dofs_1d; 
    }
    


    // TO DO get rid of one of these !?!?!?!
    num_basis = num_ref_dofs_in_elem;

    // allocate memory
    ref_dof_positions = CArrayKokkos <double> (num_ref_dofs_in_elem, num_dim);
    ref_dof_positions_1d = CArrayKokkos <double> (num_ref_dofs_1d);    
    ref_gauss_lob_weights = CArrayKokkos <double> (num_gauss_lob_in_elem);
    ref_gauss_leg_weights = CArrayKokkos <double> (num_gauss_leg_in_elem);

    // Memory for gradients
    ref_gauss_lob_grad_basis = CArrayKokkos <double> (num_gauss_lob_in_elem, num_basis, num_dim);
    ref_gauss_leg_grad_basis = CArrayKokkos  <double> (num_gauss_leg_in_elem, num_basis, num_dim);

    // Basis evaluation at the nodes
    ref_gauss_lob_basis = CArrayKokkos <double> (num_gauss_lob_in_elem, num_basis);
    ref_gauss_leg_basis = CArrayKokkos <double> (num_gauss_leg_in_elem, num_basis);
    
    ref_gauss_lob_positions = CArrayKokkos <double> (num_gauss_lob_in_elem, num_dim); 
    ref_gauss_leg_positions = CArrayKokkos <double> (num_gauss_leg_in_elem, num_dim);

    // --- build reference index spaces for 3D ---
    if(num_dim == 3){
        
        // --- build gauss nodal positions and weights ---
        lob_nodes_1D = CArrayKokkos <double> (num_gauss_lob_1d);
        lobatto_nodes_1D(lob_nodes_1D, num_gauss_lob_1d);
    
        lob_weights_1D = CArrayKokkos <double> (num_gauss_lob_1d);
        lobatto_weights_1D(lob_weights_1D, num_gauss_lob_1d);
    
        leg_nodes_1D = CArrayKokkos <double> (num_gauss_leg_1d);
        legendre_nodes_1D(leg_nodes_1D, num_gauss_leg_1d);
    
        leg_weights_1D = CArrayKokkos <double> (num_gauss_leg_1d);
        legendre_weights_1D(leg_weights_1D, num_gauss_leg_1d);


        FOR_ALL( k, 0, num_gauss_lob_1d, 
                 j, 0, num_gauss_lob_1d,
                 i, 0, num_gauss_lob_1d, { 
                    
                    int lob_rid = lobatto_rid(i,j,k);
                    
                    ref_gauss_lob_positions(lob_rid,0) = lob_nodes_1D(i);
                    ref_gauss_lob_positions(lob_rid,1) = lob_nodes_1D(j);
                    ref_gauss_lob_positions(lob_rid,2) = lob_nodes_1D(k);
                    
                    ref_gauss_lob_weights(lob_rid) = lob_weights_1D(i)*lob_weights_1D(j)*lob_weights_1D(k);
        });
    
        FOR_ALL( k, 0, num_gauss_leg_1d, 
                 j, 0, num_gauss_leg_1d,
                 i, 0, num_gauss_leg_1d, { 
        
                    int leg_rid = legendre_rid(i,j,k);
                    
                    ref_gauss_leg_positions(leg_rid,0) = leg_nodes_1D(i);
                    ref_gauss_leg_positions(leg_rid,1) = leg_nodes_1D(j);
                    ref_gauss_leg_positions(leg_rid,2) = leg_nodes_1D(k);
                    printf(" leg_weight: %f \n", leg_weights_1D(i)); 
                    ref_gauss_leg_weights(leg_rid) = leg_weights_1D(i)*leg_weights_1D(j)*leg_weights_1D(k);
        });

        // Saving vertex positions in 1D
        if( p_order == 0){
            // dofs same as lobatto quadrature points 
            FOR_ALL(i,  0, num_gauss_lob_1d,{
                ref_dof_positions_1d(i) = lob_nodes_1D(i);
            });
        }

        else{
            RUN({
                int dof_id = 0;
                for(int i = 0; i < num_gauss_lob_1d; i=i+2){

                    ref_dof_positions_1d(dof_id) = lob_nodes_1D(i);

                    dof_id++;
                }
            });  
        }

        FOR_ALL( num_k, 0, num_ref_dofs_1d, 
                 num_j, 0, num_ref_dofs_1d,
                 num_i, 0, num_ref_dofs_1d, { 
        
                    int dof_rlid = dof_rid(num_i, num_j, num_k);

                    ref_dof_positions(dof_rlid, 0) = ref_dof_positions_1d(num_i);
                    ref_dof_positions(dof_rlid, 1) = ref_dof_positions_1d(num_j);
                    ref_dof_positions(dof_rlid, 2) = ref_dof_positions_1d(num_k);
        });

        // basis and grad basis evaluations done at points //
        
        // temp variables hold evaluations at a single point for each dof //
        CArrayKokkos <double> temp_nodal_basis(num_ref_dofs_in_elem);
        
        CArrayKokkos <double> val_1d(num_ref_dofs_in_elem);
        CArrayKokkos <double> val_3d(num_ref_dofs_in_elem, 3);

        CArrayKokkos <double> point(3);
	
        // --- evaluate the basis at the lobatto positions
        FOR_ALL(gauss_lob_rid, 0, num_gauss_lob_in_elem, {

            // Get the nodal coordinates
            for(int dim = 0; dim < 3; dim++){
              point(dim) = ref_gauss_lob_positions(gauss_lob_rid, dim);
            }

            get_basis(temp_nodal_basis, val_1d, val_3d, point);
            
            //double check_basis = 0.0;

            for(int basis_id = 0; basis_id < num_ref_dofs_in_elem; basis_id++){

                ref_gauss_lob_basis(gauss_lob_rid, basis_id) = temp_nodal_basis(basis_id);
                //check_basis += temp_nodal_basis(basis_id);
            	temp_nodal_basis(basis_id) = 0.0;
            }
            //printf(" basis tally = %f \n", check_basis );

        });

	// --- evaluate the basis at the legendre points
        FOR_ALL(gauss_leg_rid,  0, num_gauss_leg_in_elem, {

            // Get the nodal coordinates
            for(int dim = 0; dim < 3; dim++){
                point(dim) = ref_gauss_leg_positions(gauss_leg_rid, dim);
            }

            get_basis(temp_nodal_basis, val_1d, val_3d, point);
            
            //double check_basis = 0.0;

            for(int basis_id = 0; basis_id < num_ref_dofs_in_elem; basis_id++){

                ref_gauss_leg_basis(gauss_leg_rid, basis_id) = temp_nodal_basis(basis_id);
              //  check_basis += temp_nodal_basis(basis_id);
                temp_nodal_basis(basis_id) = 0.0;
	    }

            //printf(" basis tally = %f \n", check_basis );
        });

        // --- evaluate grad_basis functions at the lobatto points ---

        CArrayKokkos <double> temp_partial_xi(num_ref_dofs_in_elem);
        CArrayKokkos <double> temp_partial_eta(num_ref_dofs_in_elem);
        CArrayKokkos <double> temp_partial_mu(num_ref_dofs_in_elem);
        
        CArrayKokkos <double> Dval_1d(num_ref_dofs_in_elem);
        CArrayKokkos <double> Dval_3d(num_ref_dofs_in_elem,3);
        
        FOR_ALL(gauss_lob_rid, 0, num_gauss_lob_in_elem,{

            // Get the lobatto coordinates
            for(int dim = 0; dim < 3; dim++){
                point(dim) = ref_gauss_lob_positions(gauss_lob_rid, dim);
            }

            //double check[3];
            //for (int i = 0; i < 3; i++) check[i] = 0.0;

            partial_xi_basis(temp_partial_xi, val_1d, val_3d, Dval_1d, Dval_3d, point);
            partial_eta_basis(temp_partial_eta, val_1d, val_3d, Dval_1d, Dval_3d, point);
            partial_mu_basis(temp_partial_mu, val_1d, val_3d, Dval_1d, Dval_3d, point);

            for(int basis_id = 0; basis_id < num_ref_dofs_in_elem; basis_id++){


                ref_gauss_lob_grad_basis(gauss_lob_rid, basis_id, 0) = temp_partial_xi(basis_id);
                ref_gauss_lob_grad_basis(gauss_lob_rid, basis_id, 1) = temp_partial_eta(basis_id);
                ref_gauss_lob_grad_basis(gauss_lob_rid, basis_id, 2) = temp_partial_mu(basis_id);
                
              //  check[0] += temp_partial_xi(basis_id);
              //  check[1] += temp_partial_eta(basis_id);
              //  check[2] += temp_partial_mu(basis_id);

                temp_partial_xi(basis_id)  = 0.0;
                temp_partial_eta(basis_id) = 0.0;
                temp_partial_mu(basis_id)  = 0.0;
            }
            
            //printf(" grad_basis tally = %f, %f, %f \n", check[0], check[1], check[2]);
        });


        FOR_ALL(gauss_leg_rid,  0, num_gauss_leg_in_elem, {

            // Get the nodal coordinates
            for(int dim = 0; dim < 3; dim++){
                point(dim) = ref_gauss_leg_positions(gauss_leg_rid, dim);
            }

            partial_xi_basis(temp_partial_xi, val_1d, val_3d, Dval_1d, Dval_3d, point);
            partial_eta_basis(temp_partial_eta, val_1d, val_3d, Dval_1d, Dval_3d, point);
            partial_mu_basis(temp_partial_mu, val_1d, val_3d, Dval_1d, Dval_3d, point);
            
            //double check[3];
            //for (int i = 0; i < 3; i++) check[i] = 0.0;

            for(int basis_id = 0; basis_id < num_ref_dofs_in_elem; basis_id++){


                ref_gauss_leg_grad_basis(gauss_leg_rid, basis_id, 0) = temp_partial_xi(basis_id);
                //printf(" grad basis value : %f \n ", ref_gauss_leg_grad_basis(gauss_leg_rid, basis_id, 0) );
                ref_gauss_leg_grad_basis(gauss_leg_rid, basis_id, 1) = temp_partial_eta(basis_id);
                //printf(" grad basis value : %f \n ", ref_gauss_leg_grad_basis(gauss_leg_rid, basis_id, 1) );
                ref_gauss_leg_grad_basis(gauss_leg_rid, basis_id, 2) = temp_partial_mu(basis_id);
                //printf(" grad basis value : %f \n ", ref_gauss_leg_grad_basis(gauss_leg_rid, basis_id, 2) );
              
              //  check[0] += temp_partial_xi(basis_id);
              //  check[1] += temp_partial_eta(basis_id);
              //  check[2] += temp_partial_mu(basis_id);
                
                temp_partial_xi(basis_id)  = 0.0;
                temp_partial_eta(basis_id) = 0.0;
                temp_partial_mu(basis_id)  = 0.0;
            }
            
            //printf(" grad_basis tally = %f, %f, %f \n", check[0], check[1], check[2]);
        });
        
        
    }// end 3d scope    

}; // end of member function

    
KOKKOS_FUNCTION
void lobatto_nodes_1D(
                      CArrayKokkos <double> &lob_nodes_1D,
                      const int &num){
    if (num == 1){
        lob_nodes_1D(0) = 0.0;
    }
    else if (num == 2){
        lob_nodes_1D(0) = -1.0;
        lob_nodes_1D(1) =  1.0;
    }
    else if (num == 3){
        lob_nodes_1D(0) = -1.0;
        lob_nodes_1D(1) =  0.0;
        lob_nodes_1D(2) =  1.0;
    }
    else if (num == 4){
        lob_nodes_1D(0) = -1.0;
        lob_nodes_1D(1) = -1.0/5.0*sqrt(5.0);
        lob_nodes_1D(2) =  1.0/5.0*sqrt(5.0);
        lob_nodes_1D(3) =  1.0;
    }
    else if (num == 5){
        lob_nodes_1D(0) = -1.0;
        lob_nodes_1D(1) = -1.0/7.0*sqrt(21.0);
        lob_nodes_1D(2) =  0.0;
        lob_nodes_1D(3) =  1.0/7.0*sqrt(21.0);
        lob_nodes_1D(4) =  1.0;
    }
    else if (num == 6){
        lob_nodes_1D(0) = -1.0;
        lob_nodes_1D(1) = -sqrt(1.0/21.0*(7.0 + 2.0*sqrt(7.0)));
        lob_nodes_1D(2) = -sqrt(1.0/21.0*(7.0 - 2.0*sqrt(7.0)));
        lob_nodes_1D(3) =  sqrt(1.0/21.0*(7.0 - 2.0*sqrt(7.0)));
        lob_nodes_1D(4) =  sqrt(1.0/21.0*(7.0 +2.0*sqrt(7.0)));
        lob_nodes_1D(5) =  1.0;
    }
    else if (num == 7){
        lob_nodes_1D(0) =  - 1.0E+00;
        lob_nodes_1D(1) =  - 0.830223896278566929872032213967E+00;
        lob_nodes_1D(2) =  - 0.468848793470714213803771881909E+00;
        lob_nodes_1D(3) =    0.0E+00;
        lob_nodes_1D(4) =    0.468848793470714213803771881909E+00;
        lob_nodes_1D(5) =    0.830223896278566929872032213967E+00;
        lob_nodes_1D(6) =    1.0E+00;
    }
    else if (num == 8){
        lob_nodes_1D(0) =  - 1.0E+00;
        lob_nodes_1D(1) =  - 0.871740148509606615337445761221E+00;
        lob_nodes_1D(2) =  - 0.591700181433142302144510731398E+00;
        lob_nodes_1D(3) =  - 0.209299217902478868768657260345E+00;
        lob_nodes_1D(4) =    0.209299217902478868768657260345E+00;
        lob_nodes_1D(5) =    0.591700181433142302144510731398E+00;
        lob_nodes_1D(6) =    0.871740148509606615337445761221E+00;
        lob_nodes_1D(7) =    1.0E+00;
    }
    else if (num == 9){
        lob_nodes_1D(0) =  - 1.0E+00;
        lob_nodes_1D(1) =  - 0.899757995411460157312345244418E+00;
        lob_nodes_1D(2) =  - 0.677186279510737753445885427091E+00;
        lob_nodes_1D(3) =  - 0.363117463826178158710752068709E+00;
        lob_nodes_1D(4) =    0.0E+00;
        lob_nodes_1D(5) =    0.363117463826178158710752068709E+00;
        lob_nodes_1D(6) =    0.677186279510737753445885427091E+00;
        lob_nodes_1D(7) =    0.899757995411460157312345244418E+00;
        lob_nodes_1D(8) =    1.0E+00;
    }
    else if (num == 10){
        lob_nodes_1D(0) =  - 1.0E+00;
        lob_nodes_1D(1) =  - 0.919533908166458813828932660822E+00;
        lob_nodes_1D(2) =  - 0.738773865105505075003106174860E+00;
        lob_nodes_1D(3) =  - 0.477924949810444495661175092731E+00;
        lob_nodes_1D(4) =  - 0.165278957666387024626219765958E+00;
        lob_nodes_1D(5) =    0.165278957666387024626219765958E+00;
        lob_nodes_1D(6) =    0.477924949810444495661175092731E+00;
        lob_nodes_1D(7) =    0.738773865105505075003106174860E+00;
        lob_nodes_1D(8) =    0.919533908166458813828932660822E+00;
        lob_nodes_1D(9) =   1.0E+00;

    }
    else if (num == 11){
        lob_nodes_1D(0) =  - 1.0E+00;
        lob_nodes_1D(1) =  - 0.934001430408059134332274136099E+00;
        lob_nodes_1D(2) =  - 0.784483473663144418622417816108E+00;
        lob_nodes_1D(3) =  - 0.565235326996205006470963969478E+00;
        lob_nodes_1D(4) =  - 0.295758135586939391431911515559E+00;
        lob_nodes_1D(5) =    0.0E+00;
        lob_nodes_1D(6) =    0.295758135586939391431911515559E+00;
        lob_nodes_1D(7) =    0.565235326996205006470963969478E+00;
        lob_nodes_1D(8) =    0.784483473663144418622417816108E+00;
        lob_nodes_1D(9) =   0.934001430408059134332274136099E+00;
        lob_nodes_1D(10) =   1.0E+00;
    }

    else if (num == 12){
        lob_nodes_1D(0) =  - 1.0E+00;
        lob_nodes_1D(1) =  - 0.944899272222882223407580138303E+00;
        lob_nodes_1D(2) =  - 0.819279321644006678348641581717E+00;
        lob_nodes_1D(3) =  - 0.632876153031869677662404854444E+00;
        lob_nodes_1D(4) =  - 0.399530940965348932264349791567E+00;
        lob_nodes_1D(5) =  - 0.136552932854927554864061855740E+00;
        lob_nodes_1D(6) =    0.136552932854927554864061855740E+00;
        lob_nodes_1D(7) =    0.399530940965348932264349791567E+00;
        lob_nodes_1D(8) =    0.632876153031869677662404854444E+00;
        lob_nodes_1D(9) =   0.819279321644006678348641581717E+00;
        lob_nodes_1D(10) =   0.944899272222882223407580138303E+00;
        lob_nodes_1D(11) =   1.0E+00;
    }

    else if (num == 13){
        lob_nodes_1D(0) =  - 1.0E+00;
        lob_nodes_1D(1) =  - 0.953309846642163911896905464755E+00;
        lob_nodes_1D(2) =  - 0.846347564651872316865925607099E+00;
        lob_nodes_1D(3) =  - 0.686188469081757426072759039566E+00;
        lob_nodes_1D(4) =  - 0.482909821091336201746937233637E+00;
        lob_nodes_1D(5) =  - 0.249286930106239992568673700374E+00;
        lob_nodes_1D(6) =    0.0E+00;
        lob_nodes_1D(7) =    0.249286930106239992568673700374E+00;
        lob_nodes_1D(8) =    0.482909821091336201746937233637E+00;
        lob_nodes_1D(9) =   0.686188469081757426072759039566E+00;
        lob_nodes_1D(10) =   0.846347564651872316865925607099E+00;
        lob_nodes_1D(11) =   0.953309846642163911896905464755E+00;
        lob_nodes_1D(12) =   1.0E+00;
    }

    else if (num == 14){
        lob_nodes_1D(0) =  - 1.0E+00;
        lob_nodes_1D(1) =  - 0.959935045267260901355100162015E+00;
        lob_nodes_1D(2) =  - 0.867801053830347251000220202908E+00;
        lob_nodes_1D(3) =  - 0.728868599091326140584672400521E+00;
        lob_nodes_1D(4) =  - 0.550639402928647055316622705859E+00;
        lob_nodes_1D(5) =  - 0.342724013342712845043903403642E+00;
        lob_nodes_1D(6) =  - 0.116331868883703867658776709736E+00;
        lob_nodes_1D(7) =    0.116331868883703867658776709736E+00;
        lob_nodes_1D(8) =    0.342724013342712845043903403642E+00;
        lob_nodes_1D(9) =   0.550639402928647055316622705859E+00;
        lob_nodes_1D(10) =   0.728868599091326140584672400521E+00;
        lob_nodes_1D(11) =   0.867801053830347251000220202908E+00;
        lob_nodes_1D(12) =   0.959935045267260901355100162015E+00;
        lob_nodes_1D(13) =   1.0E+00;
    }

    else if (num == 15){
        lob_nodes_1D(0) =  - 1.0E+00;
        lob_nodes_1D(1) =  - 0.965245926503838572795851392070E+00;
        lob_nodes_1D(2) =  - 0.885082044222976298825401631482E+00;
        lob_nodes_1D(3) =  - 0.763519689951815200704118475976E+00;
        lob_nodes_1D(4) =  - 0.606253205469845711123529938637E+00;
        lob_nodes_1D(5) =  - 0.420638054713672480921896938739E+00;
        lob_nodes_1D(6) =  - 0.215353955363794238225679446273E+00;
        lob_nodes_1D(7) =    0.0E+00;
        lob_nodes_1D(8) =    0.215353955363794238225679446273E+00;
        lob_nodes_1D(9) =   0.420638054713672480921896938739E+00;
        lob_nodes_1D(10) =   0.606253205469845711123529938637E+00;
        lob_nodes_1D(11) =   0.763519689951815200704118475976E+00;
        lob_nodes_1D(12) =   0.885082044222976298825401631482E+00;
        lob_nodes_1D(13) =   0.965245926503838572795851392070E+00;
        lob_nodes_1D(14) =   1.0E+00;
    }

    else if (num == 16){
        lob_nodes_1D(0) =  - 1.0E+00;
        lob_nodes_1D(1) =  - 0.969568046270217932952242738367E+00;
        lob_nodes_1D(2) =  - 0.899200533093472092994628261520E+00;
        lob_nodes_1D(3) =  - 0.792008291861815063931088270963E+00;
        lob_nodes_1D(4) =  - 0.652388702882493089467883219641E+00;
        lob_nodes_1D(5) =  - 0.486059421887137611781890785847E+00;
        lob_nodes_1D(6) =  - 0.299830468900763208098353454722E+00;
        lob_nodes_1D(7) =  - 0.101326273521949447843033005046E+00;
        lob_nodes_1D(8) =    0.101326273521949447843033005046E+00;
        lob_nodes_1D(9) =   0.299830468900763208098353454722E+00;
        lob_nodes_1D(10) =   0.486059421887137611781890785847E+00;
        lob_nodes_1D(11) =   0.652388702882493089467883219641E+00;
        lob_nodes_1D(12) =   0.792008291861815063931088270963E+00;
        lob_nodes_1D(13) =   0.899200533093472092994628261520E+00;
        lob_nodes_1D(14) =   0.969568046270217932952242738367E+00;
        lob_nodes_1D(15) =   1.0E+00;
    }

    else if (num == 17){
        lob_nodes_1D(0) =  - 1.0E+00;
        lob_nodes_1D(1) =  - 0.973132176631418314156979501874E+00;
        lob_nodes_1D(2) =  - 0.910879995915573595623802506398E+00;
        lob_nodes_1D(3) =  - 0.815696251221770307106750553238E+00;
        lob_nodes_1D(4) =  - 0.691028980627684705394919357372E+00;
        lob_nodes_1D(5) =  - 0.541385399330101539123733407504E+00;
        lob_nodes_1D(6) =  - 0.372174433565477041907234680735E+00;
        lob_nodes_1D(7) =  - 0.189511973518317388304263014753E+00;
        lob_nodes_1D(8) =    0.0E+00;
        lob_nodes_1D(9) =    0.189511973518317388304263014753E+00;
        lob_nodes_1D(10) =   0.372174433565477041907234680735E+00;
        lob_nodes_1D(11) =   0.541385399330101539123733407504E+00;
        lob_nodes_1D(12) =   0.691028980627684705394919357372E+00;
        lob_nodes_1D(13) =   0.815696251221770307106750553238E+00;
        lob_nodes_1D(14) =   0.910879995915573595623802506398E+00;
        lob_nodes_1D(15) =   0.973132176631418314156979501874E+00;
        lob_nodes_1D(16) =   1.0E+00;
    }

    else if (num == 18){
        lob_nodes_1D(0) =  - 1.0E+00;
        lob_nodes_1D(1) =  - 0.976105557412198542864518924342E+00;
        lob_nodes_1D(2) =  - 0.920649185347533873837854625431E+00;
        lob_nodes_1D(3) =  - 0.835593535218090213713646362328E+00;
        lob_nodes_1D(4) =  - 0.723679329283242681306210365302E+00;
        lob_nodes_1D(5) =  - 0.588504834318661761173535893194E+00;
        lob_nodes_1D(6) =  - 0.434415036912123975342287136741E+00;
        lob_nodes_1D(7) =  - 0.266362652878280984167665332026E+00;
        lob_nodes_1D(8) =  - 0.897490934846521110226450100886E-01;
        lob_nodes_1D(9) =    0.897490934846521110226450100886E-01;
        lob_nodes_1D(10) =   0.266362652878280984167665332026E+00;
        lob_nodes_1D(11) =   0.434415036912123975342287136741E+00;
        lob_nodes_1D(12) =   0.588504834318661761173535893194E+00;
        lob_nodes_1D(13) =   0.723679329283242681306210365302E+00;
        lob_nodes_1D(14) =   0.835593535218090213713646362328E+00;
        lob_nodes_1D(15) =   0.920649185347533873837854625431E+00;
        lob_nodes_1D(16) =   0.976105557412198542864518924342E+00;
        lob_nodes_1D(17) =   1.0E+00;
    }

    else if (num == 19) {
        lob_nodes_1D(0)=   - 1.0E+00;
        lob_nodes_1D(1)=   - 0.978611766222080095152634063110E+00;
        lob_nodes_1D(2)=   - 0.928901528152586243717940258797E+00;
        lob_nodes_1D(3)=   - 0.852460577796646093085955970041E+00;
        lob_nodes_1D(4)=   - 0.751494202552613014163637489634E+00;
        lob_nodes_1D(5)=   - 0.628908137265220497766832306229E+00;
        lob_nodes_1D(6)=   - 0.488229285680713502777909637625E+00;
        lob_nodes_1D(7)=   - 0.333504847824498610298500103845E+00;
        lob_nodes_1D(8)=   - 0.169186023409281571375154153445E+00;
        lob_nodes_1D(9)=     0.0E+00;
        lob_nodes_1D(10) =   0.169186023409281571375154153445E+00;
        lob_nodes_1D(11) =   0.333504847824498610298500103845E+00;
        lob_nodes_1D(12) =   0.488229285680713502777909637625E+00;
        lob_nodes_1D(13) =   0.628908137265220497766832306229E+00;
        lob_nodes_1D(14) =   0.751494202552613014163637489634E+00;
        lob_nodes_1D(15) =   0.852460577796646093085955970041E+00;
        lob_nodes_1D(16) =   0.928901528152586243717940258797E+00;
        lob_nodes_1D(17) =   0.978611766222080095152634063110E+00;
        lob_nodes_1D(18) =   1.0E+00;

    } // end if

    else if (num == 20){
        lob_nodes_1D(0) =  - 1.0E+00;
        lob_nodes_1D(1) =  - 0.980743704893914171925446438584E+00;
        lob_nodes_1D(2) =  - 0.935934498812665435716181584931E+00;
        lob_nodes_1D(3) =  - 0.866877978089950141309847214616E+00;
        lob_nodes_1D(4) =  - 0.775368260952055870414317527595E+00;
        lob_nodes_1D(5) =  - 0.663776402290311289846403322971E+00;
        lob_nodes_1D(6) =  - 0.534992864031886261648135961829E+00;
        lob_nodes_1D(7) =  - 0.392353183713909299386474703816E+00;
        lob_nodes_1D(8) =  - 0.239551705922986495182401356927E+00;
        lob_nodes_1D(9) =  - 0.805459372388218379759445181596E-01;
        lob_nodes_1D(10) =   0.805459372388218379759445181596E-01;
        lob_nodes_1D(11) =   0.239551705922986495182401356927E+00;
        lob_nodes_1D(12) =   0.392353183713909299386474703816E+00;
        lob_nodes_1D(13) =   0.534992864031886261648135961829E+00;
        lob_nodes_1D(14) =   0.663776402290311289846403322971E+00;
        lob_nodes_1D(15) =   0.775368260952055870414317527595E+00;
        lob_nodes_1D(16) =   0.866877978089950141309847214616E+00;
        lob_nodes_1D(17) =   0.935934498812665435716181584931E+00;
        lob_nodes_1D(18) =   0.980743704893914171925446438584E+00;
        lob_nodes_1D(19) =   1.0E+00;
    }
} // end of lobbato_nodes_1D function


/**************************************************************************************//**
*  lobatto_weights_1D creates quadrature weights corresponding the nodal positions defined 
*  on [-1,1] using the Gauss-Lobatto quadrature points. The CArray lob_weights_1D is passed 
*  in my reference and modified in place. The integer num is the numer of points being 
*  defined in 1D.
*****************************************************************************************/

KOKKOS_FUNCTION
void lobatto_weights_1D(
                        CArrayKokkos <double> &lob_weights_1D,  // Lobbatto weights
                        const int &num){                     // Interpolation order
    if (num == 1){
        lob_weights_1D(0) = 2.0;
    }
    else if (num == 2){
        lob_weights_1D(0) = 1.0;
        lob_weights_1D(1) = 1.0;
    }
    else if (num == 3){
        lob_weights_1D(0) = 1.0/3.0;
        lob_weights_1D(1) = 4.0/3.0;
        lob_weights_1D(2) = 1.0/3.0;
    }
    else if (num == 4){
        lob_weights_1D(0) = 1.0/6.0;
        lob_weights_1D(1) = 5.0/6.0;
        lob_weights_1D(2) = 5.0/6.0;
        lob_weights_1D(3) = 1.0/6.0;
    }
    else if (num == 5){
        lob_weights_1D(0) = 1.0/10.0;
        lob_weights_1D(1) = 49.0/90.0;
        lob_weights_1D(2) = 32.0/45.0;
        lob_weights_1D(3) = 49.0/90.0;
        lob_weights_1D(4) = 1.0/10.0;
    }
    else if (num == 6){
        lob_weights_1D(0) = 1.0/15.0;
        lob_weights_1D(1) = 1.0/30.0*(14.0 - sqrt(7.0));
        lob_weights_1D(2) = 1.0/30.0*(14.0 + sqrt(7.0));
        lob_weights_1D(3) = 1.0/30.0*(14.0 + sqrt(7.0));
        lob_weights_1D(4) = 1.0/30.0*(14.0 - sqrt(7.0));
        lob_weights_1D(5) = 1.0/15.0;
    }
    else if (num == 7){
        lob_weights_1D(0) =  0.476190476190476190476190476190E-01;
        lob_weights_1D(1) =  0.276826047361565948010700406290E+00;
        lob_weights_1D(2) =  0.431745381209862623417871022281E+00;
        lob_weights_1D(3) =  0.487619047619047619047619047619E+00;
        lob_weights_1D(4) =  0.431745381209862623417871022281E+00;
        lob_weights_1D(5) =  0.276826047361565948010700406290E+00;
        lob_weights_1D(6) =  0.476190476190476190476190476190E-01;
    }
    else if (num == 8){
        lob_weights_1D(0) =  0.357142857142857142857142857143E-01;
        lob_weights_1D(1) =  0.210704227143506039382991065776E+00;
        lob_weights_1D(2) =  0.341122692483504364764240677108E+00;
        lob_weights_1D(3) =  0.412458794658703881567052971402E+00;
        lob_weights_1D(4) =  0.412458794658703881567052971402E+00;
        lob_weights_1D(5) =  0.341122692483504364764240677108E+00;
        lob_weights_1D(6) =  0.210704227143506039382991065776E+00;
        lob_weights_1D(7) =  0.357142857142857142857142857143E-01;
    }
    else if (num == 9){
        lob_weights_1D(0) =  0.277777777777777777777777777778E-01;
        lob_weights_1D(1) =  0.165495361560805525046339720029E+00;
        lob_weights_1D(2) =  0.274538712500161735280705618579E+00;
        lob_weights_1D(3) =  0.346428510973046345115131532140E+00;
        lob_weights_1D(4) =  0.371519274376417233560090702948E+00;
        lob_weights_1D(5) =  0.346428510973046345115131532140E+00;
        lob_weights_1D(6) =  0.274538712500161735280705618579E+00;
        lob_weights_1D(7) =  0.165495361560805525046339720029E+00;
        lob_weights_1D(8) =  0.277777777777777777777777777778E-01;
    }
    else if (num == 10){
        lob_weights_1D(0) =  0.222222222222222222222222222222E-01;
        lob_weights_1D(1) =  0.133305990851070111126227170755E+00;
        lob_weights_1D(2) =  0.224889342063126452119457821731E+00;
        lob_weights_1D(3) =  0.292042683679683757875582257374E+00;
        lob_weights_1D(4) =  0.327539761183897456656510527917E+00;
        lob_weights_1D(5) =  0.327539761183897456656510527917E+00;
        lob_weights_1D(6) =  0.292042683679683757875582257374E+00;
        lob_weights_1D(7) =  0.224889342063126452119457821731E+00;
        lob_weights_1D(8) =  0.133305990851070111126227170755E+00;
        lob_weights_1D(9) =  0.222222222222222222222222222222E-01;
    }
    else if (num == 11){
        lob_weights_1D(0) =  0.181818181818181818181818181818E-01;
        lob_weights_1D(1) =  0.109612273266994864461403449580E+00;
        lob_weights_1D(2) =  0.187169881780305204108141521899E+00;
        lob_weights_1D(3) =  0.248048104264028314040084866422E+00;
        lob_weights_1D(4) =  0.286879124779008088679222403332E+00;
        lob_weights_1D(5) =  0.300217595455690693785931881170E+00;
        lob_weights_1D(6) =  0.286879124779008088679222403332E+00;
        lob_weights_1D(7) =  0.248048104264028314040084866422E+00;
        lob_weights_1D(8) =  0.187169881780305204108141521899E+00;
        lob_weights_1D(9) =  0.109612273266994864461403449580E+00;
        lob_weights_1D(10)=  0.181818181818181818181818181818E-01;
    }

    else if (num == 12){
        lob_weights_1D(0) =  0.151515151515151515151515151515E-01;
        lob_weights_1D(1) =  0.916845174131961306683425941341E-01;
        lob_weights_1D(2) =  0.157974705564370115164671062700E+00;
        lob_weights_1D(3) =  0.212508417761021145358302077367E+00;
        lob_weights_1D(4) =  0.251275603199201280293244412148E+00;
        lob_weights_1D(5) =  0.271405240910696177000288338500E+00;
        lob_weights_1D(6) =  0.271405240910696177000288338500E+00;
        lob_weights_1D(7) =  0.251275603199201280293244412148E+00;
        lob_weights_1D(8) =  0.212508417761021145358302077367E+00;
        lob_weights_1D(9) =  0.157974705564370115164671062700E+00;
        lob_weights_1D(10) = 0.916845174131961306683425941341E-01;
        lob_weights_1D(11) = 0.151515151515151515151515151515E-01;
    }

    else if (num == 13){
        lob_weights_1D(0) =  0.128205128205128205128205128205E-01;
        lob_weights_1D(1) =  0.778016867468189277935889883331E-01;
        lob_weights_1D(2) =  0.134981926689608349119914762589E+00;
        lob_weights_1D(3) =  0.183646865203550092007494258747E+00;
        lob_weights_1D(4) =  0.220767793566110086085534008379E+00;
        lob_weights_1D(5) =  0.244015790306676356458578148360E+00;
        lob_weights_1D(6) =  0.251930849333446736044138641541E+00;
        lob_weights_1D(7) =  0.244015790306676356458578148360E+00;
        lob_weights_1D(8) =  0.220767793566110086085534008379E+00;
        lob_weights_1D(9) =  0.183646865203550092007494258747E+00;
        lob_weights_1D(10) = 0.134981926689608349119914762589E+00;
        lob_weights_1D(11) = 0.778016867468189277935889883331E-01;
        lob_weights_1D(12) = 0.128205128205128205128205128205E-01;
    }

    else if (num == 14){
        lob_weights_1D(0) =  0.109890109890109890109890109890E-01;
        lob_weights_1D(1) =  0.668372844976812846340706607461E-01;
        lob_weights_1D(2) =  0.116586655898711651540996670655E+00;
        lob_weights_1D(3) =  0.160021851762952142412820997988E+00;
        lob_weights_1D(4) =  0.194826149373416118640331778376E+00;
        lob_weights_1D(5) =  0.219126253009770754871162523954E+00;
        lob_weights_1D(6) =  0.231612794468457058889628357293E+00;
        lob_weights_1D(7) =  0.231612794468457058889628357293E+00;
        lob_weights_1D(8) =  0.219126253009770754871162523954E+00;
        lob_weights_1D(9) =  0.194826149373416118640331778376E+00;
        lob_weights_1D(10) = 0.160021851762952142412820997988E+00;
        lob_weights_1D(11) = 0.116586655898711651540996670655E+00;
        lob_weights_1D(12) = 0.668372844976812846340706607461E-01;
        lob_weights_1D(13) = 0.109890109890109890109890109890E-01;
    }


    else if (num == 15){
        lob_weights_1D(0) =  0.952380952380952380952380952381E-02;
        lob_weights_1D(1) =  0.580298930286012490968805840253E-01;
        lob_weights_1D(2) =  0.101660070325718067603666170789E+00;
        lob_weights_1D(3) =  0.140511699802428109460446805644E+00;
        lob_weights_1D(4) =  0.172789647253600949052077099408E+00;
        lob_weights_1D(5) =  0.196987235964613356092500346507E+00;
        lob_weights_1D(6) =  0.211973585926820920127430076977E+00;
        lob_weights_1D(7) =  0.217048116348815649514950214251E+00;
        lob_weights_1D(8) =  0.211973585926820920127430076977E+00;
        lob_weights_1D(9) =  0.196987235964613356092500346507E+00;
        lob_weights_1D(10) = 0.172789647253600949052077099408E+00;
        lob_weights_1D(11) = 0.140511699802428109460446805644E+00;
        lob_weights_1D(12) = 0.101660070325718067603666170789E+00;
        lob_weights_1D(13) = 0.580298930286012490968805840253E-01;
        lob_weights_1D(14) = 0.952380952380952380952380952381E-02;
    }


    else if (num == 16){
        lob_weights_1D(0) =  0.833333333333333333333333333333E-02;
        lob_weights_1D(1) =  0.508503610059199054032449195655E-01;
        lob_weights_1D(2) =  0.893936973259308009910520801661E-01;
        lob_weights_1D(3) =  0.124255382132514098349536332657E+00;
        lob_weights_1D(4) =  0.154026980807164280815644940485E+00;
        lob_weights_1D(5) =  0.177491913391704125301075669528E+00;
        lob_weights_1D(6) =  0.193690023825203584316913598854E+00;
        lob_weights_1D(7) =  0.201958308178229871489199125411E+00;
        lob_weights_1D(8) =  0.201958308178229871489199125411E+00;
        lob_weights_1D(9) =  0.193690023825203584316913598854E+00;
        lob_weights_1D(10) = 0.177491913391704125301075669528E+00;
        lob_weights_1D(11) = 0.154026980807164280815644940485E+00;
        lob_weights_1D(12) = 0.124255382132514098349536332657E+00;
        lob_weights_1D(13) = 0.893936973259308009910520801661E-01;
        lob_weights_1D(14) = 0.508503610059199054032449195655E-01;
        lob_weights_1D(15) = 0.833333333333333333333333333333E-02;
    }


    else if (num == 17){
        lob_weights_1D(0) =  0.735294117647058823529411764706E-02;
        lob_weights_1D(1) =  0.449219405432542096474009546232E-01;
        lob_weights_1D(2) =  0.791982705036871191902644299528E-01;
        lob_weights_1D(3) =  0.110592909007028161375772705220E+00;
        lob_weights_1D(4) =  0.137987746201926559056201574954E+00;
        lob_weights_1D(5) =  0.160394661997621539516328365865E+00;
        lob_weights_1D(6) =  0.177004253515657870436945745363E+00;
        lob_weights_1D(7) =  0.187216339677619235892088482861E+00;
        lob_weights_1D(8) =  0.190661874753469433299407247028E+00;
        lob_weights_1D(9) =  0.187216339677619235892088482861E+00;
        lob_weights_1D(10) = 0.177004253515657870436945745363E+00;
        lob_weights_1D(11) = 0.160394661997621539516328365865E+00;
        lob_weights_1D(12) = 0.137987746201926559056201574954E+00;
        lob_weights_1D(13) = 0.110592909007028161375772705220E+00;
        lob_weights_1D(14) = 0.791982705036871191902644299528E-01;
        lob_weights_1D(15) = 0.449219405432542096474009546232E-01;
        lob_weights_1D(16) = 0.735294117647058823529411764706E-02;
    }

    else if (num == 18){
        lob_weights_1D(0) =  0.653594771241830065359477124183E-02;
        lob_weights_1D(1) =  0.399706288109140661375991764101E-01;
        lob_weights_1D(2) =  0.706371668856336649992229601678E-01;
        lob_weights_1D(3) =  0.990162717175028023944236053187E-01;
        lob_weights_1D(4) =  0.124210533132967100263396358897E+00;
        lob_weights_1D(5) =  0.145411961573802267983003210494E+00;
        lob_weights_1D(6) =  0.161939517237602489264326706700E+00;
        lob_weights_1D(7) =  0.173262109489456226010614403827E+00;
        lob_weights_1D(8) =  0.179015863439703082293818806944E+00;
        lob_weights_1D(9) =  0.179015863439703082293818806944E+00;
        lob_weights_1D(10) = 0.173262109489456226010614403827E+00;
        lob_weights_1D(11) = 0.161939517237602489264326706700E+00;
        lob_weights_1D(12) = 0.145411961573802267983003210494E+00;
        lob_weights_1D(13) = 0.124210533132967100263396358897E+00;
        lob_weights_1D(14) = 0.990162717175028023944236053187E-01;
        lob_weights_1D(15) = 0.706371668856336649992229601678E-01;
        lob_weights_1D(16) = 0.399706288109140661375991764101E-01;
        lob_weights_1D(17) = 0.653594771241830065359477124183E-02;
    }

    else if (num == 19) {
        lob_weights_1D(0) =  0.584795321637426900584795321637E-02;
        lob_weights_1D(1) =  0.357933651861764771154255690351E-01;
        lob_weights_1D(2) =  0.633818917626297368516956904183E-01;
        lob_weights_1D(3) =  0.891317570992070844480087905562E-01;
        lob_weights_1D(4) =  0.112315341477305044070910015464E+00;
        lob_weights_1D(5) =  0.132267280448750776926046733910E+00;
        lob_weights_1D(6) =  0.148413942595938885009680643668E+00;
        lob_weights_1D(7) =  0.160290924044061241979910968184E+00;
        lob_weights_1D(8) =  0.167556584527142867270137277740E+00;
        lob_weights_1D(9) =  0.170001919284827234644672715617E+00;
        lob_weights_1D(10) = 0.167556584527142867270137277740E+00;
        lob_weights_1D(11) = 0.160290924044061241979910968184E+00;
        lob_weights_1D(12) = 0.148413942595938885009680643668E+00;
        lob_weights_1D(13) = 0.132267280448750776926046733910E+00;
        lob_weights_1D(14) = 0.112315341477305044070910015464E+00;
        lob_weights_1D(15) = 0.891317570992070844480087905562E-01;
        lob_weights_1D(16) = 0.633818917626297368516956904183E-01;
        lob_weights_1D(17) = 0.357933651861764771154255690351E-01;
        lob_weights_1D(18) = 0.584795321637426900584795321637E-02;
    } // end if

    else if (num == 20) {
        lob_weights_1D(0) =  0.526315789473684210526315789474E-02;
        lob_weights_1D(1) =  0.322371231884889414916050281173E-01;
        lob_weights_1D(2) =  0.571818021275668260047536271732E-01;
        lob_weights_1D(3) =  0.806317639961196031447768461137E-01;
        lob_weights_1D(4) =  0.101991499699450815683781205733E+00;
        lob_weights_1D(5) =  0.120709227628674725099429705002E+00;
        lob_weights_1D(6) =  0.136300482358724184489780792989E+00;
        lob_weights_1D(7) =  0.148361554070916825814713013734E+00;
        lob_weights_1D(8) =  0.156580102647475487158169896794E+00;
        lob_weights_1D(9) =  0.160743286387845749007726726449E+00;
        lob_weights_1D(10) = 0.160743286387845749007726726449E+00;
        lob_weights_1D(11) = 0.156580102647475487158169896794E+00;
        lob_weights_1D(12) = 0.148361554070916825814713013734E+00;
        lob_weights_1D(13) = 0.136300482358724184489780792989E+00;
        lob_weights_1D(14) = 0.120709227628674725099429705002E+00;
        lob_weights_1D(15) = 0.101991499699450815683781205733E+00;
        lob_weights_1D(16) = 0.806317639961196031447768461137E-01;
        lob_weights_1D(17) = 0.571818021275668260047536271732E-01;
        lob_weights_1D(18) = 0.322371231884889414916050281173E-01;
        lob_weights_1D(19) = 0.526315789473684210526315789474E-02;
    } // end if
} // end of lobatto_weights_1D function

KOKKOS_FUNCTION
void legendre_nodes_1D(
                       CArrayKokkos <double> &leg_nodes_1D,
                       const int &num){

    if (num == 1){
        leg_nodes_1D(0) = 0.0;
    }
    else if (num == 2){
        leg_nodes_1D(0) = -0.577350269189625764509148780501;
        leg_nodes_1D(1) =  0.577350269189625764509148780501;
    }
    else if (num == 3){
        leg_nodes_1D(0) = -0.774596669241483377035853079956;
        leg_nodes_1D(1) =  0.0;
        leg_nodes_1D(2) =  0.774596669241483377035853079956;
    }
    else if (num == 4){
        leg_nodes_1D(0) = -0.861136311594052575223946488892;
        leg_nodes_1D(1) = -0.339981043584856264802665759103;
        leg_nodes_1D(2) =  0.339981043584856264802665759103;
        leg_nodes_1D(3) =  0.861136311594052575223946488892;
    }
    else if (num == 5){
        leg_nodes_1D(0) = -0.906179845938663992797626878299;
        leg_nodes_1D(1) = -0.538469310105683091036314420700;
        leg_nodes_1D(2) =  0.0;
        leg_nodes_1D(3) =  0.538469310105683091036314420700;
        leg_nodes_1D(4) =  0.906179845938663992797626878299;
    }
    else if (num == 6){
        leg_nodes_1D(0) = -0.932469514203152027812301554493;
        leg_nodes_1D(1) = -0.661209386466264513661399595019;
        leg_nodes_1D(2) = -0.238619186083196908630501721680;
        leg_nodes_1D(3) =  0.238619186083196908630501721680;
        leg_nodes_1D(4) =  0.661209386466264513661399595019;
        leg_nodes_1D(5) =  0.932469514203152027812301554493;
    }
    else if (num == 7){
        leg_nodes_1D(0) = -0.949107912342758524526189684047;
        leg_nodes_1D(1) = -0.741531185599394439863864773280;
        leg_nodes_1D(2) = -0.405845151377397166906606412076;
        leg_nodes_1D(3) =  0.0E+00;
        leg_nodes_1D(4) =  0.405845151377397166906606412076;
        leg_nodes_1D(5) =  0.741531185599394439863864773280;
        leg_nodes_1D(6) =  0.949107912342758524526189684047;
    }
    else if (num == 8){
        leg_nodes_1D(0) = -0.960289856497536231683560868569;
        leg_nodes_1D(1) = -0.796666477413626739591553936475;
        leg_nodes_1D(2) = -0.525532409916328985817739049189;
        leg_nodes_1D(3) = -0.183434642495649804939476142360;
        leg_nodes_1D(4) =  0.183434642495649804939476142360;
        leg_nodes_1D(5) =  0.525532409916328985817739049189;
        leg_nodes_1D(6) =  0.796666477413626739591553936475;
        leg_nodes_1D(7) =  0.960289856497536231683560868569;
    }
    else if (num == 9){
        leg_nodes_1D(0) = -0.968160239507626089835576202903;
        leg_nodes_1D(1) = -0.836031107326635794299429788069;
        leg_nodes_1D(2) = -0.613371432700590397308702039341;
        leg_nodes_1D(3) = -0.324253423403808929038538014643;
        leg_nodes_1D(4) =  0.0E+00;
        leg_nodes_1D(5) =  0.324253423403808929038538014643;
        leg_nodes_1D(6) =  0.613371432700590397308702039341;
        leg_nodes_1D(7) =  0.836031107326635794299429788069;
        leg_nodes_1D(8) =  0.968160239507626089835576202903;
    }
    else if (num == 10){
        leg_nodes_1D(0) = -0.9739065285171717200779640120844;
        leg_nodes_1D(1) = -0.8650633666889845107320966884234;
        leg_nodes_1D(2) = -0.6794095682990244062343273651148;
        leg_nodes_1D(3) = -0.4333953941292471907992659431657;
        leg_nodes_1D(4) = -0.1488743389816312108848260011297;
        leg_nodes_1D(5) =  0.1488743389816312108848260011297;
        leg_nodes_1D(6) =  0.4333953941292471907992659431657;
        leg_nodes_1D(7) =  0.6794095682990244062343273651148;
        leg_nodes_1D(8) =  0.8650633666889845107320966884234;
        leg_nodes_1D(9) =  0.9739065285171717200779640120844;
        
    }
    else if (num == 11){
        leg_nodes_1D(0) = -0.9782286581460569928039380011228;
        leg_nodes_1D(1) = -0.8870625997680952990751577693039;
        leg_nodes_1D(2) = -0.7301520055740493240934162520311;
        leg_nodes_1D(3) = -0.5190961292068118159257256694586;
        leg_nodes_1D(4) = -0.2695431559523449723315319854008;
        leg_nodes_1D(5) =  0.0E+00;
        leg_nodes_1D(6) =  0.2695431559523449723315319854008;
        leg_nodes_1D(7) =  0.5190961292068118159257256694586;
        leg_nodes_1D(8) =  0.7301520055740493240934162520311;
        leg_nodes_1D(9) =  0.8870625997680952990751577693039;
        leg_nodes_1D(10) = 0.9782286581460569928039380011228;
    }
    
    else if (num == 12){
        leg_nodes_1D(0) = -0.9815606342467192506905490901492;
        leg_nodes_1D(1) = -0.9041172563704748566784658661190;
        leg_nodes_1D(2) = -0.7699026741943046870368938332128;
        leg_nodes_1D(3) = -0.5873179542866174472967024189405;
        leg_nodes_1D(4) = -0.3678314989981801937526915366437;
        leg_nodes_1D(5) = -0.1252334085114689154724413694638;
        leg_nodes_1D(6) =  0.1252334085114689154724413694638;
        leg_nodes_1D(7) =  0.3678314989981801937526915366437;
        leg_nodes_1D(8) =  0.5873179542866174472967024189405;
        leg_nodes_1D(9) =  0.7699026741943046870368938332128;
        leg_nodes_1D(10) = 0.9041172563704748566784658661190;
        leg_nodes_1D(11) = 0.9815606342467192506905490901492;
    }
    
    else if (num == 13){
        leg_nodes_1D(0) = -0.98418305471858814947282944880710;
        leg_nodes_1D(1) = -0.91759839922297796520654783650071;
        leg_nodes_1D(2) = -0.80157809073330991279420648958285;
        leg_nodes_1D(3) = -0.64234933944034022064398460699551;
        leg_nodes_1D(4) = -0.44849275103644685287791285212763;
        leg_nodes_1D(5) = -0.23045831595513479406552812109798;
        leg_nodes_1D(6) =  0.0E+00;
        leg_nodes_1D(7) =  0.23045831595513479406552812109798;
        leg_nodes_1D(8) =  0.44849275103644685287791285212763;
        leg_nodes_1D(9) =  0.64234933944034022064398460699551;
        leg_nodes_1D(10) = 0.80157809073330991279420648958285;
        leg_nodes_1D(11) = 0.91759839922297796520654783650071;
        leg_nodes_1D(12) = 0.98418305471858814947282944880710;
    }
    
    else if (num == 14){
        leg_nodes_1D(0) = -0.986283808696812338841597266704052;
        leg_nodes_1D(1) = -0.928434883663573517336391139377874;
        leg_nodes_1D(2) = -0.827201315069764993189794742650394;
        leg_nodes_1D(3) = -0.687292904811685470148019803019334;
        leg_nodes_1D(4) = -0.515248636358154091965290718551188;
        leg_nodes_1D(5) = -0.319112368927889760435671824168475;
        leg_nodes_1D(6) = -0.108054948707343662066244650219834;
        leg_nodes_1D(7) =  0.108054948707343662066244650219834;
        leg_nodes_1D(8) =  0.319112368927889760435671824168475;
        leg_nodes_1D(9) =  0.515248636358154091965290718551188;
        leg_nodes_1D(10) = 0.687292904811685470148019803019334;
        leg_nodes_1D(11) = 0.827201315069764993189794742650394;
        leg_nodes_1D(12) = 0.928434883663573517336391139377874;
        leg_nodes_1D(13) = 0.986283808696812338841597266704052;
    }
    
    else if (num == 15){
        leg_nodes_1D(0) = -0.987992518020485428489565718586612;
        leg_nodes_1D(1) = -0.937273392400705904307758947710209;
        leg_nodes_1D(2) = -0.848206583410427216200648320774216;
        leg_nodes_1D(3) = -0.724417731360170047416186054613938;
        leg_nodes_1D(4) = -0.570972172608538847537226737253910;
        leg_nodes_1D(5) = -0.394151347077563369897207370981045;
        leg_nodes_1D(6) = -0.201194093997434522300628303394596;
        leg_nodes_1D(7) =  0.0E+00;
        leg_nodes_1D(8) =  0.201194093997434522300628303394596;
        leg_nodes_1D(9) =  0.394151347077563369897207370981045;
        leg_nodes_1D(10) = 0.570972172608538847537226737253910;
        leg_nodes_1D(11) = 0.724417731360170047416186054613938;
        leg_nodes_1D(12) = 0.848206583410427216200648320774216;
        leg_nodes_1D(13) = 0.937273392400705904307758947710209;
        leg_nodes_1D(14) = 0.987992518020485428489565718586612;
    }
    
    else if (num == 16){
        leg_nodes_1D(0) = -0.989400934991649932596154173450332;
        leg_nodes_1D(1) = -0.944575023073232576077988415534608;
        leg_nodes_1D(2) = -0.865631202387831743880467897712393;
        leg_nodes_1D(3) = -0.755404408355003033895101194847442;
        leg_nodes_1D(4) = -0.617876244402643748446671764048791;
        leg_nodes_1D(5) = -0.458016777657227386342419442983577;
        leg_nodes_1D(6) = -0.281603550779258913230460501460496;
        leg_nodes_1D(7) = -0.095012509837637440185319335424958;
        leg_nodes_1D(8) =  0.095012509837637440185319335424958;
        leg_nodes_1D(9) =  0.281603550779258913230460501460496;
        leg_nodes_1D(10) = 0.458016777657227386342419442983577;
        leg_nodes_1D(11) = 0.617876244402643748446671764048791;
        leg_nodes_1D(12) = 0.755404408355003033895101194847442;
        leg_nodes_1D(13) = 0.865631202387831743880467897712393;
        leg_nodes_1D(14) = 0.944575023073232576077988415534608;
        leg_nodes_1D(15) = 0.989400934991649932596154173450332;
    }
    
    else if (num == 17){
        leg_nodes_1D(0) = -0.990575475314417335675434019940665;
        leg_nodes_1D(1) = -0.950675521768767761222716957895803;
        leg_nodes_1D(2) = -0.880239153726985902122955694488155;
        leg_nodes_1D(3) = -0.781514003896801406925230055520476;
        leg_nodes_1D(4) = -0.657671159216690765850302216643002;
        leg_nodes_1D(5) = -0.512690537086476967886246568629551;
        leg_nodes_1D(6) = -0.351231763453876315297185517095346;
        leg_nodes_1D(7) = -0.178484181495847855850677493654065;
        leg_nodes_1D(8) =  0.0E+00;
        leg_nodes_1D(9) =  0.178484181495847855850677493654065;
        leg_nodes_1D(10) = 0.351231763453876315297185517095346;
        leg_nodes_1D(11) = 0.512690537086476967886246568629551;
        leg_nodes_1D(12) = 0.657671159216690765850302216643002;
        leg_nodes_1D(13) = 0.781514003896801406925230055520476;
        leg_nodes_1D(14) = 0.880239153726985902122955694488155;
        leg_nodes_1D(15) = 0.950675521768767761222716957895803;
        leg_nodes_1D(16) = 0.990575475314417335675434019940665;
    }
    
    else if (num == 18){
        leg_nodes_1D(0) = -0.991565168420930946730016004706150;
        leg_nodes_1D(1) = -0.955823949571397755181195892929776;
        leg_nodes_1D(2) = -0.892602466497555739206060591127145;
        leg_nodes_1D(3) = -0.803704958972523115682417455014590;
        leg_nodes_1D(4) = -0.691687043060353207874891081288848;
        leg_nodes_1D(5) = -0.559770831073947534607871548525329;
        leg_nodes_1D(6) = -0.411751161462842646035931793833051;
        leg_nodes_1D(7) = -0.251886225691505509588972854877911;
        leg_nodes_1D(8) = -0.084775013041735301242261852935783;
        leg_nodes_1D(9) =  0.084775013041735301242261852935783;
        leg_nodes_1D(10) = 0.251886225691505509588972854877911;
        leg_nodes_1D(11) = 0.411751161462842646035931793833051;
        leg_nodes_1D(12) = 0.559770831073947534607871548525329;
        leg_nodes_1D(13) = 0.691687043060353207874891081288848;
        leg_nodes_1D(14) = 0.803704958972523115682417455014590;
        leg_nodes_1D(15) = 0.892602466497555739206060591127145;
        leg_nodes_1D(16) = 0.955823949571397755181195892929776;
        leg_nodes_1D(17) = 0.991565168420930946730016004706150;
    }
    
    else if (num == 19) {
        leg_nodes_1D(0) = -0.992406843843584403189017670253260;
        leg_nodes_1D(1) = -0.960208152134830030852778840687651;
        leg_nodes_1D(2) = -0.903155903614817901642660928532312;
        leg_nodes_1D(3) = -0.822714656537142824978922486712713;
        leg_nodes_1D(4) = -0.720966177335229378617095860823781;
        leg_nodes_1D(5) = -0.600545304661681023469638164946239;
        leg_nodes_1D(6) = -0.464570741375960945717267148104102;
        leg_nodes_1D(7) = -0.316564099963629831990117328849844;
        leg_nodes_1D(8) = -0.160358645640225375868096115740743;
        leg_nodes_1D(9) =  0.0E+00;
        leg_nodes_1D(10) = 0.160358645640225375868096115740743;
        leg_nodes_1D(11) = 0.316564099963629831990117328849844;
        leg_nodes_1D(12) = 0.464570741375960945717267148104102;
        leg_nodes_1D(13) = 0.600545304661681023469638164946239;
        leg_nodes_1D(14) = 0.720966177335229378617095860823781;
        leg_nodes_1D(15) = 0.822714656537142824978922486712713;
        leg_nodes_1D(16) = 0.903155903614817901642660928532312;
        leg_nodes_1D(17) = 0.960208152134830030852778840687651;
        leg_nodes_1D(18) = 0.992406843843584403189017670253260;
        
    } // end if
    
    
}; // end of legendre_nodes_1D function


KOKKOS_FUNCTION
void legendre_weights_1D(
                         CArrayKokkos <double> &leg_weights_1D,  // Legendre weights
                         const int &num){                  // Interpolation order
    if (num == 1){
        leg_weights_1D(0) = 2.0;
    }
    else if (num == 2){
        leg_weights_1D(0) = 1.0;
        leg_weights_1D(1) = 1.0;
    }
    else if (num == 3){
        leg_weights_1D(0) = 0.555555555555555555555555555555555;
        leg_weights_1D(1) = 0.888888888888888888888888888888888;
        leg_weights_1D(2) = 0.555555555555555555555555555555555;
    }
    else if (num == 4){
        leg_weights_1D(0) = 0.347854845137453857373063949221999;
        leg_weights_1D(1) = 0.652145154862546142626936050778000;
        leg_weights_1D(2) = 0.652145154862546142626936050778000;
        leg_weights_1D(3) = 0.347854845137453857373063949221999;
    }
    else if (num == 5){
        leg_weights_1D(0) = 0.236926885056189087514264040719917;
        leg_weights_1D(1) = 0.478628670499366468041291514835638;
        leg_weights_1D(2) = 0.568888888888888888888888888888888;
        leg_weights_1D(3) = 0.478628670499366468041291514835638;
        leg_weights_1D(4) = 0.236926885056189087514264040719917;
    }
    else if (num == 6){
        leg_weights_1D(0) = 0.171324492379170345040296142172732;
        leg_weights_1D(1) = 0.360761573048138607569833513837716;
        leg_weights_1D(2) = 0.467913934572691047389870343989550;
        leg_weights_1D(3) = 0.467913934572691047389870343989550;
        leg_weights_1D(4) = 0.360761573048138607569833513837716;
        leg_weights_1D(5) = 0.171324492379170345040296142172732;
    }
    else if (num == 7){
        leg_weights_1D(0) = 0.129484966168869693270611432679082;
        leg_weights_1D(1) = 0.279705391489276667901467771423779;
        leg_weights_1D(2) = 0.381830050505118944950369775488975;
        leg_weights_1D(3) = 0.417959183673469387755102040816326;
        leg_weights_1D(4) = 0.381830050505118944950369775488975;
        leg_weights_1D(5) = 0.279705391489276667901467771423779;
        leg_weights_1D(6) = 0.129484966168869693270611432679082;
    }
    else if (num == 8){
        leg_weights_1D(0) = 0.101228536290376259152531354309962;
        leg_weights_1D(1) = 0.222381034453374470544355994426240;
        leg_weights_1D(2) = 0.313706645877887287337962201986601;
        leg_weights_1D(3) = 0.362683783378361982965150449277195;
        leg_weights_1D(4) = 0.362683783378361982965150449277195;
        leg_weights_1D(5) = 0.313706645877887287337962201986601;
        leg_weights_1D(6) = 0.222381034453374470544355994426240;
        leg_weights_1D(7) = 0.101228536290376259152531354309962;
    }
    else if (num == 9){
        leg_weights_1D(0) = 0.081274388361574411971892158110523;
        leg_weights_1D(1) = 0.180648160694857404058472031242912;
        leg_weights_1D(2) = 0.260610696402935462318742869418632;
        leg_weights_1D(3) = 0.312347077040002840068630406584443;
        leg_weights_1D(4) = 0.330239355001259763164525069286974;
        leg_weights_1D(5) = 0.312347077040002840068630406584443;
        leg_weights_1D(6) = 0.260610696402935462318742869418632;
        leg_weights_1D(7) = 0.180648160694857404058472031242912;
        leg_weights_1D(8) = 0.081274388361574411971892158110523;
    }
    else if (num == 10){
        leg_weights_1D(0) = 0.066671344308688137593568809893331;
        leg_weights_1D(1) = 0.149451349150580593145776339657697;
        leg_weights_1D(2) = 0.219086362515982043995534934228163;
        leg_weights_1D(3) = 0.269266719309996355091226921569469;
        leg_weights_1D(4) = 0.295524224714752870173892994651338;
        leg_weights_1D(5) = 0.295524224714752870173892994651338;
        leg_weights_1D(6) = 0.269266719309996355091226921569469;
        leg_weights_1D(7) = 0.219086362515982043995534934228163;
        leg_weights_1D(8) = 0.149451349150580593145776339657697;
        leg_weights_1D(9) = 0.066671344308688137593568809893331;
    }
    else if (num == 11){
        leg_weights_1D(0) = 0.055668567116173666482753720442548;
        leg_weights_1D(1) = 0.125580369464904624634694299223940;
        leg_weights_1D(2) = 0.186290210927734251426097641431655;
        leg_weights_1D(3) = 0.233193764591990479918523704843175;
        leg_weights_1D(4) = 0.262804544510246662180688869890509;
        leg_weights_1D(5) = 0.272925086777900630714483528336342;
        leg_weights_1D(6) = 0.262804544510246662180688869890509;
        leg_weights_1D(7) = 0.233193764591990479918523704843175;
        leg_weights_1D(8) = 0.186290210927734251426097641431655;
        leg_weights_1D(9) = 0.125580369464904624634694299223940;
        leg_weights_1D(10)= 0.055668567116173666482753720442548;
    }
    
    else if (num == 12){
        leg_weights_1D(0) =  0.04717533638651182719461596148501;
        leg_weights_1D(1) =  0.10693932599531843096025471819399;
        leg_weights_1D(2) =  0.16007832854334622633465252954335;
        leg_weights_1D(3) =  0.20316742672306592174906445580979;
        leg_weights_1D(4) =  0.23349253653835480876084989892487;
        leg_weights_1D(5) =  0.24914704581340278500056243604295;
        leg_weights_1D(6) =  0.24914704581340278500056243604295;
        leg_weights_1D(7) =  0.23349253653835480876084989892487;
        leg_weights_1D(8) =  0.20316742672306592174906445580979;
        leg_weights_1D(9) =  0.16007832854334622633465252954335;
        leg_weights_1D(10) = 0.10693932599531843096025471819399;
        leg_weights_1D(11) = 0.04717533638651182719461596148501;
    }
    
    else if (num == 13){
        leg_weights_1D(0) =  0.04048400476531587952002159220098;
        leg_weights_1D(1) =  0.09212149983772844791442177595379;
        leg_weights_1D(2) =  0.13887351021978723846360177686887;
        leg_weights_1D(3) =  0.17814598076194573828004669199609;
        leg_weights_1D(4) =  0.20781604753688850231252321930605;
        leg_weights_1D(5) =  0.22628318026289723841209018603977;
        leg_weights_1D(6) =  0.23255155323087391019458951526883;
        leg_weights_1D(7) =  0.22628318026289723841209018603977;
        leg_weights_1D(8) =  0.20781604753688850231252321930605;
        leg_weights_1D(9) =  0.17814598076194573828004669199609;
        leg_weights_1D(10) = 0.13887351021978723846360177686887;
        leg_weights_1D(11) = 0.09212149983772844791442177595379;
        leg_weights_1D(12) = 0.04048400476531587952002159220098;
    }
    
    else if (num == 14){
        leg_weights_1D(0) =  0.03511946033175186303183287613819;
        leg_weights_1D(1) =  0.08015808715976020980563327706285;
        leg_weights_1D(2) =  0.12151857068790318468941480907247;
        leg_weights_1D(3) =  0.15720316715819353456960193862384;
        leg_weights_1D(4) =  0.18553839747793781374171659012515;
        leg_weights_1D(5) =  0.20519846372129560396592406566121;
        leg_weights_1D(6) =  0.21526385346315779019587644331626;
        leg_weights_1D(7) =  0.21526385346315779019587644331626;
        leg_weights_1D(8) =  0.20519846372129560396592406566121;
        leg_weights_1D(9) =  0.18553839747793781374171659012515;
        leg_weights_1D(10) = 0.15720316715819353456960193862384;
        leg_weights_1D(11) = 0.12151857068790318468941480907247;
        leg_weights_1D(12) = 0.08015808715976020980563327706285;
        leg_weights_1D(13) = 0.03511946033175186303183287613819;
    }
    
    
    else if (num == 15){
        leg_weights_1D(0) =  0.03075324199611726835462839357720;
        leg_weights_1D(1) =  0.07036604748810812470926741645066;
        leg_weights_1D(2) =  0.10715922046717193501186954668586;
        leg_weights_1D(3) =  0.13957067792615431444780479451102;
        leg_weights_1D(4) =  0.16626920581699393355320086048120;
        leg_weights_1D(5) =  0.18616100001556221102680056186642;
        leg_weights_1D(6) =  0.19843148532711157645611832644383;
        leg_weights_1D(7) =  0.20257824192556127288062019996751;
        leg_weights_1D(8) =  0.19843148532711157645611832644383;
        leg_weights_1D(9) =  0.18616100001556221102680056186642;
        leg_weights_1D(10) = 0.16626920581699393355320086048120;
        leg_weights_1D(11) = 0.13957067792615431444780479451102;
        leg_weights_1D(12) = 0.10715922046717193501186954668586;
        leg_weights_1D(13) = 0.07036604748810812470926741645066;
        leg_weights_1D(14) = 0.03075324199611726835462839357720;
    }
    
    
    else if (num == 16){
        leg_weights_1D(0) =  0.02715245941175409485178057245601;
        leg_weights_1D(1) =  0.06225352393864789286284383699437;
        leg_weights_1D(2) =  0.09515851168249278480992510760224;
        leg_weights_1D(3) =  0.12462897125553387205247628219201;
        leg_weights_1D(4) =  0.14959598881657673208150173054747;
        leg_weights_1D(5) =  0.16915651939500253818931207903035;
        leg_weights_1D(6) =  0.18260341504492358886676366796921;
        leg_weights_1D(7) =  0.18945061045506849628539672320828;
        leg_weights_1D(8) =  0.18945061045506849628539672320828;
        leg_weights_1D(9) =  0.18260341504492358886676366796921;
        leg_weights_1D(10) = 0.16915651939500253818931207903035;
        leg_weights_1D(11) = 0.14959598881657673208150173054747;
        leg_weights_1D(12) = 0.12462897125553387205247628219201;
        leg_weights_1D(13) = 0.09515851168249278480992510760224;
        leg_weights_1D(14) = 0.06225352393864789286284383699437;
        leg_weights_1D(15) = 0.02715245941175409485178057245601;
    }
    
    
    else if (num == 17){
        leg_weights_1D(0) =  0.02414830286854793196011002628756;
        leg_weights_1D(1) =  0.05545952937398720112944016535824;
        leg_weights_1D(2) =  0.08503614831717918088353537019106;
        leg_weights_1D(3) =  0.11188384719340397109478838562635;
        leg_weights_1D(4) =  0.13513636846852547328631998170235;
        leg_weights_1D(5) =  0.15404576107681028808143159480195;
        leg_weights_1D(6) =  0.16800410215645004450997066378832;
        leg_weights_1D(7) =  0.17656270536699264632527099011319;
        leg_weights_1D(8) =  0.17944647035620652545826564426188;
        leg_weights_1D(9) =  0.17656270536699264632527099011319;
        leg_weights_1D(10) = 0.16800410215645004450997066378832;
        leg_weights_1D(11) = 0.15404576107681028808143159480195;
        leg_weights_1D(12) = 0.13513636846852547328631998170235;
        leg_weights_1D(13) = 0.11188384719340397109478838562635;
        leg_weights_1D(14) = 0.08503614831717918088353537019106;
        leg_weights_1D(15) = 0.05545952937398720112944016535824;
        leg_weights_1D(16) = 0.02414830286854793196011002628756;
    }
    
    else if (num == 18){
        leg_weights_1D(0) =  0.02161601352648331031334271026645;
        leg_weights_1D(1) =  0.04971454889496979645333494620263;
        leg_weights_1D(2) =  0.07642573025488905652912967761663;
        leg_weights_1D(3) =  0.10094204410628716556281398492483;
        leg_weights_1D(4) =  0.12255520671147846018451912680020;
        leg_weights_1D(5) =  0.14064291467065065120473130375194;
        leg_weights_1D(6) =  0.15468467512626524492541800383637;
        leg_weights_1D(7) =  0.16427648374583272298605377646592;
        leg_weights_1D(8) =  0.16914238296314359184065647013498;
        leg_weights_1D(9) =  0.16914238296314359184065647013498;
        leg_weights_1D(10) = 0.16427648374583272298605377646592;
        leg_weights_1D(11) = 0.15468467512626524492541800383637;
        leg_weights_1D(12) = 0.14064291467065065120473130375194;
        leg_weights_1D(13) = 0.12255520671147846018451912680020;
        leg_weights_1D(14) = 0.10094204410628716556281398492483;
        leg_weights_1D(15) = 0.07642573025488905652912967761663;
        leg_weights_1D(16) = 0.04971454889496979645333494620263;
        leg_weights_1D(17) = 0.02161601352648331031334271026645;
    }
    
    else if (num == 19) {
        leg_weights_1D(0) =  0.01946178822972647703631204146443;
        leg_weights_1D(1) =  0.04481422676569960033283815740199;
        leg_weights_1D(2) =  0.06904454273764122658070825800601;
        leg_weights_1D(3) =  0.09149002162244999946446209412383;
        leg_weights_1D(4) =  0.11156664554733399471602390168176;
        leg_weights_1D(5) =  0.12875396253933622767551578485687;
        leg_weights_1D(6) =  0.14260670217360661177574610944190;
        leg_weights_1D(7) =  0.15276604206585966677885540089766;
        leg_weights_1D(8) =  0.15896884339395434764995643946504;
        leg_weights_1D(9) =  0.16105444984878369597916362532091;
        leg_weights_1D(10) = 0.15896884339395434764995643946504;
        leg_weights_1D(11) = 0.15276604206585966677885540089766;
        leg_weights_1D(12) = 0.14260670217360661177574610944190;
        leg_weights_1D(13) = 0.12875396253933622767551578485687;
        leg_weights_1D(14) = 0.11156664554733399471602390168176;
        leg_weights_1D(15) = 0.09149002162244999946446209412383;
        leg_weights_1D(16) = 0.06904454273764122658070825800601;
        leg_weights_1D(17) = 0.04481422676569960033283815740199;
        leg_weights_1D(18) = 0.01946178822972647703631204146443;
    } // end if
    
    
} // end of legendre_weights_1D function
 

// --- ref index access member functions ---
KOKKOS_INLINE_FUNCTION
int dof_rid(int i, int j, int k) const 
{
    return i + j*num_ref_dofs_1d + k*num_ref_dofs_1d*num_ref_dofs_1d;
};
KOKKOS_INLINE_FUNCTION
int lobatto_rid(int i, int j, int k) const 
{
    return i + j*num_gauss_lob_1d + k*num_gauss_lob_1d*num_gauss_lob_1d;
};
KOKKOS_INLINE_FUNCTION
int legendre_rid(int i, int j, int k) const 
{
    return i + j*num_gauss_leg_1d + k*num_gauss_leg_1d*num_gauss_leg_1d;
};

KOKKOS_INLINE_FUNCTION
void get_basis(const CArrayKokkos <double> &basis,
               const CArrayKokkos <double> &val_1d,
               const CArrayKokkos <double> &val_3d,
               const CArrayKokkos <double> &point){

        
        // initialize to zero //
        for (int i =0; i< num_ref_dofs_1d; i++){
          val_1d(i) = 0.0;
        }
        
        // Calculate 1D basis for the X coordinate of the point
        lagrange_basis_1D(val_1d, point(0));
        
        // Save the basis value at the point to a temp array and zero out the temp array
        for(int i = 0; i < num_ref_dofs_1d; i++){
            val_3d(i,0) = val_1d(i);
            val_1d(i) = 0.0;
        }

        // Calculate 1D basis for the Y coordinate of the point
        lagrange_basis_1D(val_1d, point(1));
        
        // Save the basis value at the point to a temp array and zero out the temp array
        for(int i = 0; i < num_ref_dofs_1d; i++){
            val_3d(i,1) = val_1d(i);
            val_1d(i) = 0.0;
        }

        // Calculate 1D basis for the Z coordinate of the point
        lagrange_basis_1D(val_1d, point(2));
        
        // Save the basis value at the point to a temp array and zero out the temp array
        for(int i = 0; i < num_ref_dofs_1d; i++){
            val_3d(i,2) = val_1d(i);
            val_1d(i) = 0.0;
        }
        
        // Multiply the i, j, k components of the basis from each node
        // to get the tensor product basis for the node
        for(int k = 0; k < num_ref_dofs_1d; k++){
            for(int j = 0; j < num_ref_dofs_1d; j++){
                for(int i = 0; i < num_ref_dofs_1d; i++){

                    int dof_rlid = dof_rid(i,j,k);
                    basis(dof_rlid) = val_3d(i,0)*val_3d(j,1)*val_3d(k,2);
                }
            }
        }

        for (int i =0; i< num_ref_dofs_1d; i++){
          val_1d(i) = 0.0;
          val_3d(i,0) = 0.0;
          val_3d(i,1) = 0.0;
          val_3d(i,2) = 0.0;
        }
};

KOKKOS_INLINE_FUNCTION
void partial_xi_basis(const CArrayKokkos <double> &partial_xi,
                      const CArrayKokkos <double> &val_1d,
                      const CArrayKokkos <double> &val_3d,
                      const CArrayKokkos <double> &Dval_1d,
                      const CArrayKokkos <double> &Dval_3d,
                      const CArrayKokkos <double> &point)
{


        
        //initialize//
        for (int i = 0; i < num_ref_dofs_1d; i++){
           val_1d(i) = 0.0;
           Dval_1d(i) = 0.0;
        }


        // Calculate 1D partial w.r.t. xi for the X coordinate of the point
        lagrange_derivative_1D(Dval_1d, point(0));


        // Save the basis value at the point to a temp array and zero out the temp array
        for(int i = 0; i < num_ref_dofs_1d; i++){
            
            Dval_3d(i,0) = Dval_1d(i);
            Dval_1d(i) = 0.0;
        }


        // Calculate 1D basis for the Y coordinate of the point
        lagrange_basis_1D(val_1d, point(1));
        
        // Save the basis value at the point to a temp array and zero out the temp array
        for(int i = 0; i < num_ref_dofs_1d; i++){
            
            val_3d(i,1) = val_1d(i);
            val_1d(i) = 0.0;
        }


        // Calculate 1D basis for the Z coordinate of the point
        lagrange_basis_1D(val_1d, point(2));
        
        // Save the basis value at the point to a temp array and zero out the temp array
        for(int i = 0; i < num_ref_dofs_1d; i++){
            
            val_3d(i,2) = val_1d(i);
            val_1d(i) = 0.0;
        }

        // Multiply the i, j, k components of the basis and partial_xi from each node
        // to get the tensor product partial derivatives of the basis at each node
        for(int k = 0; k < num_ref_dofs_1d; k++){
            for(int j = 0; j < num_ref_dofs_1d; j++){
                for(int i = 0; i < num_ref_dofs_1d; i++){
                    
                    int dof_rlid = dof_rid(i,j,k);

                    // Partial w.r.t xi
                    partial_xi(dof_rlid) = Dval_3d(i, 0)*val_3d(j, 1)*val_3d(k, 2);

                }
            }
        }

        for (int i =0; i< num_ref_dofs_1d; i++){
          val_1d(i) = 0.0;
          val_3d(i,0) = 0.0;
          val_3d(i,1) = 0.0;
          val_3d(i,2) = 0.0;
          Dval_1d(i) = 0.0;
          Dval_3d(i,0) = 0.0;
          Dval_3d(i,1) = 0.0;
          Dval_3d(i,2) = 0.0;
        }
    };

KOKKOS_INLINE_FUNCTION
void partial_eta_basis(const CArrayKokkos <double> &partial_eta,
                       const CArrayKokkos <double> &val_1d,
                       const CArrayKokkos <double> &val_3d,
                       const CArrayKokkos <double> &Dval_1d,
                       const CArrayKokkos <double> &Dval_3d,
                       const CArrayKokkos <double> &point){   

        //initialize//
        for (int i = 0; i < num_ref_dofs_1d; i++){
           val_1d(i) = 0.0;
           Dval_1d(i) = 0.0;
        }

        // Calculate 1D basis for the Y coordinate of the point
        lagrange_basis_1D(val_1d, point(0));
        
        // Save the basis value at the point to a temp array and zero out the temp array
        for(int i = 0; i < num_ref_dofs_1d; i++){            
            val_3d(i,0) = val_1d(i);
            val_1d(i) = 0.0;
        }

        // Calculate 1D partial w.r.t. eta for the Y coordinate of the point
        lagrange_derivative_1D(Dval_1d, point(1));
        
        // Save the basis value at the point to a temp array and zero out the temp array
        for(int i = 0; i < num_ref_dofs_1d; i++){
            
            Dval_3d(i,1) = Dval_1d(i);

            Dval_1d(i) = 0.0;
        }


        // Calculate 1D basis for the Z coordinate of the point
        lagrange_basis_1D(val_1d, point(2));
        
        // Save the basis value at the point to a temp array and zero out the temp array
        for(int i = 0; i < num_ref_dofs_1d; i++){
            
            val_3d(i,2) = val_1d(i);
            val_1d(i) = 0.0;
        }

        // Multiply the i, j, k components of the basis and partial_eta from each node
        // to get the tensor product partial derivatives of the basis at each node
        for(int k = 0; k < num_ref_dofs_1d; k++){
            for(int j = 0; j < num_ref_dofs_1d; j++){
                for(int i = 0; i < num_ref_dofs_1d; i++){
                    
                    int dof_rlid = dof_rid(i,j,k);

                    // Partial w.r.t xi
                    partial_eta(dof_rlid) = val_3d(i, 0)*Dval_3d(j, 1)*val_3d(k, 2);

                }
            }
        }

        for (int i =0; i< num_ref_dofs_1d; i++){
          val_1d(i) = 0.0;
          val_3d(i,0) = 0.0;
          val_3d(i,1) = 0.0;
          val_3d(i,2) = 0.0;
          Dval_1d(i) = 0.0;
          Dval_3d(i,0) = 0.0;
          Dval_3d(i,1) = 0.0;
          Dval_3d(i,2) = 0.0;
        }
    };


KOKKOS_INLINE_FUNCTION
void partial_mu_basis(const CArrayKokkos <double> &partial_mu, 
                      const CArrayKokkos <double> &val_1d,
                      const CArrayKokkos <double> &val_3d,
                      const CArrayKokkos <double> &Dval_1d,
                      const CArrayKokkos <double> &Dval_3d,
                      const CArrayKokkos <double> &point){

        //initialize//
        for (int i = 0; i < num_ref_dofs_1d; i++){
           val_1d(i) = 0.0;
           Dval_1d(i) = 0.0;
        }

        // Calculate 1D basis for the X coordinate of the point
        lagrange_basis_1D(val_1d, point(0));
        
        // Save the basis value at the point to a temp array and zero out the temp array
        for(int i = 0; i < num_ref_dofs_1d; i++){
            
            val_3d(i,0) = val_1d(i);
            val_1d(i) = 0.0;
        }


        // Calculate 1D basis for the Y coordinate of the point
        lagrange_basis_1D(val_1d, point(1));
        
        // Save the basis value at the point to a temp array and zero out the temp array
        for(int i = 0; i < num_ref_dofs_1d; i++){
            
            val_3d(i,1) = val_1d(i);
            val_1d(i) = 0.0;
        }


        // Calculate 1D partial w.r.t. mu for the Z coordinate of the point
        lagrange_derivative_1D(Dval_1d, point(2));
        
        // Save the basis value at the point to a temp array and zero out the temp array
        for(int i = 0; i < num_ref_dofs_1d; i++){
            
            Dval_3d(i,2) = Dval_1d(i);
            val_1d(i) = 0.0;
        }

        // Multiply the i, j, k components of the basis and partial_xi from each node
        // to get the tensor product partial derivatives of the basis at each node
        for(int k = 0; k < num_ref_dofs_1d; k++){
            for(int j = 0; j < num_ref_dofs_1d; j++){
                for(int i = 0; i < num_ref_dofs_1d; i++){
                    
                    int dof_rlid = dof_rid(i,j,k);

                    // Partial w.r.t mu
                    partial_mu(dof_rlid) = val_3d(i, 0)*val_3d(j, 1)*Dval_3d(k, 2);

                }
            }
        }

        for (int i =0; i< num_ref_dofs_1d; i++){
          val_1d(i) = 0.0;
          val_3d(i,0) = 0.0;
          val_3d(i,1) = 0.0;
          val_3d(i,2) = 0.0;
          Dval_1d(i) = 0.0;
          Dval_3d(i,0) = 0.0;
          Dval_3d(i,1) = 0.0;
          Dval_3d(i,2) = 0.0;
        }
};


KOKKOS_INLINE_FUNCTION
void lagrange_basis_1D(
        const CArrayKokkos <double> &interp,    // interpolant from each basis
        const double x_point){     // point of interest in element
                         
        
        // calculate the basis value associated with each node_i
        for(int vert_i = 0; vert_i < num_ref_dofs_1d; vert_i++){ 
            
            double numerator = 1.0;         // placeholder numerator
            double denominator = 1.0;       // placeholder denominator
            double interpolant = 1.0;       // placeholder value of numerator/denominator
            

            for(int vert_j = 0; vert_j < num_ref_dofs_1d; vert_j++){  // looping over the verts !=vert_i
                if (vert_j != vert_i ){
                    
                    // Calculate the numerator
                    numerator = numerator*(x_point - ref_dof_positions_1d(vert_j));
                    
                    // Calculate the denominator 
                    denominator = denominator*(ref_dof_positions_1d(vert_i) - ref_dof_positions_1d(vert_j));
                
                }//end if
                
                interpolant = numerator/denominator; // storing a single value for interpolation for node vert_i
                
            } // end looping over nodes != vert_i

            // writing value to vectors for later use
            interp(vert_i)   = interpolant;           // Interpolant value at given point

        } // end loop over all nodes
} // end of Lagrange_1D function


KOKKOS_INLINE_FUNCTION
void lagrange_derivative_1D(
                            const CArrayKokkos <double> &derivative,    // derivative
                            const double x_point){         // point of interest in element

        for(int vert_i = 0; vert_i < num_ref_dofs_1d; vert_i++){ // looping over the nodes
        

            double denominator = 1.0;       // placeholder denominator
            double num_gradient = 0.0;      // placeholder for numerator of the gradient
            double gradient = 0.0;

            for(int vert_j = 0; vert_j < num_ref_dofs_1d; vert_j++){  // looping over the nodes !=vert_i
                if (vert_j != vert_i ){

                    // Calculate the denominator that is the same for 
                    // both the basis and the gradient of the basis
                    denominator = denominator*(ref_dof_positions_1d(vert_i) - ref_dof_positions_1d(vert_j));
                    
                    double product_gradient = 1.0;
                    
                    // Calculate the numerator of the gradient
                    for(int N = 0; N < num_ref_dofs_1d; N++){  // looping over the nodes !=vert_i
                        
                        if (N != vert_j && N != vert_i ){
                            product_gradient = product_gradient * (x_point - ref_dof_positions_1d(N));
                        }// end if
                    }//end for
                    
                    // Sum over the product of the numerator 
                    // contributions from each node
                    num_gradient += product_gradient; 
                
                }//end if
                
                gradient = (num_gradient/denominator); // storing the derivative of the interpolating function
            
            } // end looping over nodes != vert_i

            // writing value to vectors for later use
            derivative(vert_i)  = gradient;    // derivative of each function

        } // end loop over all nodes
} // end of Lagrange_1D function


};

#endif 
