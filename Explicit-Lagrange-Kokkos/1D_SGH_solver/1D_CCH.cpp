// -----------------------------------------------------------------------------
//
//    This is a 1D c++ cell-centered finite volume code for material dynamics
//    written using MATAR+Kokkos for performance portabilityover CPUs and GPUs.
//    The code structure is designed to immulate the 3D high-order DG code
//
//    Written by: Nathaniel Morgan
//      June 27, 2023
//
//    To run the code with e.g., 4 threads, type
//      ./hydro --kokkos-threads=4
//
//    openMP settings:
//      setenv OMP_PROC_BIND true
//      setenv OMP_PROC_BIND spread
//      setenv OMP_PLACES threads
//      setenv OMP_NUM_THREADS 2
// -----------------------------------------------------------------------------

#include "matar.h"
#include <stdio.h>
#include <math.h>  // c math lib
#include <chrono>

using namespace mtr;

// -----------------------------------------------------------------------------
//    Global variables
// -----------------------------------------------------------------------------
const double fuzz = 1.0E-15;
const double huge = 1.0E15;


// -----------------------------------------------------------------------------
//    A region
// -----------------------------------------------------------------------------
struct region_t{
   double x_min;
   double x_max;
   double den;
   double sie;
   double vel;
   double gamma;
};


// -----------------------------------------------------------------------------
//    Functions
// -----------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION
int get_corners_in_cell(int,int);

KOKKOS_INLINE_FUNCTION
int get_corners_in_node(int,int);

KOKKOS_INLINE_FUNCTION
double limiter(double, double, double, double, double);

KOKKOS_FUNCTION
void Riemann_solve(DCArrayKokkos <double> &node_coords,
                   DCArrayKokkos <double> &cell_coords,
                   DCArrayKokkos <double> &cell_vel,
                   DCArrayKokkos <double> &cell_vel_n,
                   DCArrayKokkos <double> &node_vel,
                   DCArrayKokkos <double> &cell_pres,
                   DCArrayKokkos <double> &cell_den,
                   DCArrayKokkos <double> &cell_sspd,
                   DCArrayKokkos <double> &corner_vel,
                   DCArrayKokkos <double> &corner_pres,
                   DCArrayKokkos <double> &corner_imp,
                   DCArrayKokkos <double> &corner_force,
                   DCArrayKokkos <double> &cell_dpres_dx,
                   DCArrayKokkos <double> &cell_dvel_dx,
                   const double integral_grad_basis[2],
                   const int num_cells,
                   const int num_nodes);



// -----------------------------------------------------------------------------
//    The Main function
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]){
    
    // -------------------------------
    //    User settable variables
    // -------------------------------

    // time step settings
    const double time_max = 20.0;
    double       dt       = 0.01;
    const double dt_max   = 100;
    const double dt_cfl   = 0.3;
    const int    num_rk_stages = 2;
    const int    max_cycles = 2000000;

    // mesh information
    const double  x_min = 0.0;
    const double  x_max = 100.0;
    const int num_cells = 1000;

    // global model parameters
    const double sspd_min = 1.0E-3;

    // intial conditions for each region
    
    // Sod
    const int num_regions = 2;
    region_t ics[num_regions];
    
    ics[0].x_min = 0.0;
    ics[0].x_max = 50.0;
    ics[0].den   = 1.0;
    ics[0].sie   = 2.5;
    ics[0].vel   = 0.0;
    ics[0].gamma = 1.4;
    
    ics[1].x_min = 50.0;
    ics[1].x_max = 100.0;
    ics[1].den   = 0.125;
    ics[1].sie   = 2.0;
    ics[1].vel   = 0.0;
    ics[1].gamma = 1.4;

    // -------------------------------
    
    printf("\nstarting FE code\n");
    
    FILE * myfile;
    
    // This is the meat in the code, it must be inside Kokkos scope
    Kokkos::initialize(argc, argv);
    {
        
        // 1D linear element int( -Grad(phi) dot jJ^{-1} )
        const double integral_grad_basis[2] = {1.0, -1.0};
        
        // calculate mesh information based on inputs
        const int num_nodes = num_cells+1;
        const int num_corners = 2*num_cells;
        double dx = (x_max-x_min)/num_cells;
        
        // initialize the time to zero
        double time = 0.0;
    
        
        // --- setup variables based on user inputs ---
        
        // cell variables
        DCArrayKokkos <double> cell_den(num_cells);   // density
        DCArrayKokkos <double> cell_pres(num_cells);  // pressure
        DCArrayKokkos <double> cell_sspd(num_cells);  // sound speed
        DCArrayKokkos <double> cell_sie(num_cells);   // specific internal energy
        DCArrayKokkos <double> cell_tau(num_cells);   // specific total energy
        DCArrayKokkos <double> cell_tau_n(num_cells); // specific total energy at t_n
        DCArrayKokkos <double> cell_mass(num_cells);  // mass in the cell
        DCArrayKokkos <double> cell_vel(num_cells);   // velocity
        DCArrayKokkos <double> cell_vel_n(num_cells); // velocity at t_n
        
        DCArrayKokkos <double> cell_dpres_dx(num_cells);  // pressure gradient
        DCArrayKokkos <double> cell_dvel_dx(num_cells);   // velocity gradient
        
        DCArrayKokkos <double> cell_gamma(num_cells);  // gamma law gas
        
        // nodal variables
        DCArrayKokkos <double> node_vel(num_nodes);    // velocity
        DCArrayKokkos <double> node_vel_n(num_nodes);  // the velocity at t_n
        
        // corner variables
        DCArrayKokkos <double> corner_force(num_corners); // force from cell to node
        DCArrayKokkos <double> corner_pres(num_corners);  // pres in cell corner
        DCArrayKokkos <double> corner_vel(num_corners);   // velocity in cell corner
        DCArrayKokkos <double> corner_imp(num_corners);   // accoustic impedance in cell corner
        
        // mesh variables
        DCArrayKokkos <double> cell_coords(num_cells);   // coordinates of cell
        DCArrayKokkos <double> cell_vol(num_cells);      // volume of the cell
        
        DCArrayKokkos <double> node_coords(num_nodes);   // coordinates of nodes
        DCArrayKokkos <double> node_coords_n(num_nodes); // coordinates at t_n
        
        
        
        // --- build the mesh ---
        
        // calculate nodal coordinates of the mesh
        FOR_ALL (node_id, 0, num_nodes, {
           node_coords(node_id) = double(node_id) * dx;
        }); // end parallel for on device

        
        
        // calculate cell center coordinates of the mesh
        FOR_ALL (cell_id, 0, num_cells, {
            cell_coords(cell_id) =
                           0.5*(node_coords(cell_id) + node_coords(cell_id+1));
            
            cell_vol(cell_id)  = node_coords(cell_id+1) - node_coords(cell_id);
        }); // end parallel for on device

        
        
        // --- initial state on the mesh ---
        
        // initial cell state
        FOR_ALL (cell_id, 0, num_cells, {
            
            // loop over the regions
            for (int reg=0; reg<num_regions; reg++){

               if (cell_coords(cell_id) >= ics[reg].x_min &&
                   cell_coords(cell_id) <= ics[reg].x_max){
                     
                   cell_den(cell_id) = ics[reg].den;
                   cell_sie(cell_id) = ics[reg].sie;
                   cell_gamma(cell_id) = ics[reg].gamma;
                   
                   cell_pres(cell_id)  = cell_den(cell_id)*cell_sie(cell_id)*
                                         (cell_gamma(cell_id) - 1.0);
                   
                   cell_sspd(cell_id)  = sqrt( cell_gamma(cell_id)*
                                         cell_pres(cell_id)/cell_den(cell_id) );
                   
                   cell_mass(cell_id)  = ics[reg].den*cell_vol(cell_id);
                   
                   cell_vel(cell_id) = ics[reg].vel;
                   double ke         = 0.5*cell_vel(cell_id)*cell_vel(cell_id);
                   cell_tau(cell_id) = cell_sie(cell_id) + ke;
                   
                   
               } // end if
                  
            } // end for
            
        }); // end parallel for on device
        
        
       
        
        // calculate intial Riemann velocity and cooresponding corner forces
        Riemann_solve(node_coords,
                      cell_coords,
                      cell_vel,
                      cell_vel_n,
                      node_vel,
                      cell_pres,
                      cell_den,
                      cell_sspd,
                      corner_vel,
                      corner_pres,
                      corner_imp,
                      corner_force,
                      cell_dpres_dx,
                      cell_dvel_dx,
                      integral_grad_basis,
                      num_cells,
                      num_nodes);
        
        
        // -------------------------------
        //    Print intiial state to file
        // -------------------------------
        
        // update the host side to print (i.e., copy from device to host)
        cell_coords.update_host();
        cell_den.update_host();
        cell_pres.update_host();
        cell_sie.update_host();
        cell_tau.update_host();
        cell_vel.update_host();
        node_vel.update_host();
        
        // write out the intial conditions to a file on the host
        myfile=fopen("time0.txt","w");
        fprintf(myfile,"# x  den  pres  sie vel \n");
        
        // write data on the host side
        for (int cell_id=0; cell_id<num_cells; cell_id++){
        fprintf( myfile,"%f\t%f\t%f\t%f\t%f\n",
                 cell_coords.host(cell_id),
                 cell_den.host(cell_id),
                 cell_pres.host(cell_id),
                 cell_sie.host(cell_id),
                 cell_vel.host(cell_id) );
        }
        fclose(myfile);
        
        // total energy check
        double total_e = 0.0;
        double e_lcl = 0.0;
        REDUCE_SUM(cell_id, 0, num_cells, e_lcl, {
            double ke = 0.5*cell_vel(cell_id)*cell_vel(cell_id);
            e_lcl += cell_mass(cell_id)*(cell_sie(cell_id) + ke);
        }, total_e);
        
        
        auto time_1 = std::chrono::high_resolution_clock::now();
        
        // -------------------------------------
        // Solve equations until time=time_max
        // -------------------------------------
        for (int cycle = 0; cycle<max_cycles; cycle++){
            
            printf("time %f \n", time);
           
            
            // get the new time step
            double dt_ceiling = dt*1.1;
            
            // parallel reduction with min
            double dt_lcl;
            double min_dt_calc;
            REDUCE_MIN(cell_id, 0, num_cells, dt_lcl, {
                // mesh size
                double dx_lcl = node_coords(cell_id+1) - node_coords(cell_id);
                
                // local dt calc
                double dt_lcl_ = dt_cfl*dx_lcl/(cell_sspd(cell_id) + fuzz);
                
                // make dt be in bounds
                dt_lcl_ = fmin(dt_lcl_, dt_max);
                dt_lcl_ = fmin(dt_lcl_, time_max-time);
                dt_lcl_ = fmin(dt_lcl_, dt_ceiling);
        
                if (dt_lcl_ < dt_lcl) dt_lcl = dt_lcl_;
                        
            }, min_dt_calc); // end parallel reduction on min
            Kokkos::fence();
            
            // save the min dt
            if(min_dt_calc < dt) dt = min_dt_calc;
            
            
            //printf("time = %f, dt = %f \n", time, dt);
            if (dt<=fuzz) break;
            
            
            
            // --- integrate forward in time ---
            
            // Runge-Kutta loop
            for (int rk_stage=0; rk_stage<num_rk_stages; rk_stage++ ){
               
                // rk coefficient on dt
                double rk_alpha = 1.0/(double(num_rk_stages) - double(rk_stage));
            
                
                // save the state at t_n
                if (rk_stage==0){
                    
                    // nodal state
                    FOR_ALL (node_id, 0, num_nodes, {
                          node_vel_n(node_id)    = node_vel(node_id);
                          node_coords_n(node_id) = node_coords(node_id);
                    }); // end parallel for on device
                    
                    
                    // cell state
                    FOR_ALL (cell_id, 0, num_cells, {
                        cell_vel_n(cell_id) = cell_vel(cell_id);
                        cell_tau_n(cell_id) = cell_tau(cell_id);
                    }); // end parallel for on device
                    Kokkos::fence();
                    
                } // end if
                
                
                // --- Evolve mesh positions ---
                
                // x_new = x_n + alpha*dt*vel
                FOR_ALL (node_id, 0, num_nodes, {
                    
                    // update node coords
                    node_coords(node_id) = node_coords_n(node_id) +
                                           rk_alpha*dt*node_vel(node_id);
                    
                }); // end parallel for on device
                

                
                // --- Evolve specific total energy and velocity ---
                
                FOR_ALL (cell_id, 0, num_cells, {
                    
                    int node_id_0 = cell_id;     // left node
                    int node_id_1 = cell_id+1;   // right node
                    
                    // get the global indices for the corners of this cell
                    int corner_id_0 = get_corners_in_cell(cell_id, 0); // left
                    int corner_id_1 = get_corners_in_cell(cell_id, 1); // right

                    double force_tally = corner_force(corner_id_0) +
                                         corner_force(corner_id_1);
                    
                    // update velocity
                    // v_new = v_n + alpha*dt/mass*Sum(forces)
                    cell_vel(cell_id) = cell_vel_n(cell_id) +
                                        rk_alpha*dt/cell_mass(cell_id)*force_tally;
                    
                    
                    double power_tally = corner_force(corner_id_0)*node_vel(node_id_0) +
                                         corner_force(corner_id_1)*node_vel(node_id_1);
                    
                    // update specific total energy
                    // tau_new = tau_n + alpha*dt/mass*Sum(forces*vel^*)
                    cell_tau(cell_id) = cell_tau_n(cell_id) +
                                        rk_alpha*dt/cell_mass(cell_id)*power_tally;
                    
                    
                }); // end parallel for on device
                Kokkos::fence();
                
                
                
                
                // --- Update the state on the mesh ---
                
                // Calculate the ke and ie at the mat_pnts using the limited fields
                FOR_ALL (cell_id, 0, num_cells, {
                    
                    int node_id_0 = cell_id;     // left node
                    int node_id_1 = cell_id+1;   // right node
                    
                    // specific ke
                    double ke = 0.5*cell_vel(cell_id)*cell_vel(cell_id);
                    
                    // specific interal energy
                    cell_sie(cell_id) = cell_tau(cell_id) - ke;
                    
                    
                    // update vol
                    cell_vol(cell_id)  = node_coords(node_id_1) - node_coords(node_id_0);
                    
                    // update coordinates
                    cell_coords(cell_id) = 0.5*(node_coords(node_id_1) + node_coords(node_id_0));
                    
                    // update density
                    cell_den(cell_id) = cell_mass(cell_id)/cell_vol(cell_id);
                    
                    // update pressure via EOS
                    cell_pres(cell_id) = cell_den(cell_id)*cell_sie(cell_id)*
                                         (cell_gamma(cell_id) - 1.0);
                    
                    // update sound speed
                    cell_sspd(cell_id) = sqrt( cell_gamma(cell_id)*
                                         cell_pres(cell_id)/cell_den(cell_id) );
                    cell_sspd(cell_id) = fmax(cell_sspd(cell_id), sspd_min);
                    
                }); // end parallel for on device
                Kokkos::fence();
                
                
                
                
                // calculate the Riemann velocity and cooresponding corner forces
                Riemann_solve(node_coords,
                              cell_coords,
                              cell_vel,
                              cell_vel_n,
                              node_vel,
                              cell_pres,
                              cell_den,
                              cell_sspd,
                              corner_vel,
                              corner_pres,
                              corner_imp,
                              corner_force,
                              cell_dpres_dx,
                              cell_dvel_dx,
                              integral_grad_basis,
                              num_cells,
                              num_nodes);
                
            } // end rk loop

            
            
            // update the time
            time += dt;
            if (abs(time-time_max)<=fuzz) time=time_max;
            
        } // end for cycles in calculation
        //------------- Done with calculation ------------------
        
        auto time_2 = std::chrono::high_resolution_clock::now();
        
        auto calc_time = std::chrono::duration_cast
                           <std::chrono::nanoseconds>(time_2 - time_1).count();
        printf("\nCalculation time in seconds: %f \n", calc_time * 1e-9);
        
        // -------------------------------
        //    Print final state to a file
        // -------------------------------
        
        // update the host side to print (i.e., copy from device to host)
        cell_coords.update_host();
        cell_den.update_host();
        cell_pres.update_host();
        cell_sie.update_host();
        cell_tau.update_host();
        cell_vel.update_host();
        node_vel.update_host();
        
        // write out the intial conditions to a file on the host
        myfile=fopen("timeEnd.txt","w");
        fprintf(myfile,"# x  den  pres  sie vel \n");
        
        // write data on the host side
        for (int cell_id=0; cell_id<num_cells; cell_id++){
           fprintf( myfile,"%f\t%f\t%f\t%f\t%f\n",
                    cell_coords.host(cell_id),
                    cell_den.host(cell_id),
                    cell_pres.host(cell_id),
                    cell_sie.host(cell_id),
                    cell_vel.host(cell_id) );
        }
        fclose(myfile);
        
        
        // total energy check
        double total_e_final = 0.0;
        double ef_lcl = 0.0;
        REDUCE_SUM(cell_id, 0, num_cells, ef_lcl, {
            double ke = 0.5*cell_vel(cell_id)*cell_vel(cell_id);
            ef_lcl += cell_mass(cell_id)*(cell_sie(cell_id) + ke);
        }, total_e_final);
        Kokkos::fence();
        
        printf("total energy, TE(t=0) = %f", total_e);
        printf(" , TE(t=end) = %f" ,  total_e_final);
        printf(" , TE error = %f \n", total_e_final-total_e);
        
        
        // ======== Done using Kokkos ============
        
    } // end of kokkos scope
    Kokkos::finalize();
    
    
    printf("\nfinished\n\n");
    return 0;
  
} // end main function



KOKKOS_INLINE_FUNCTION
int get_corners_in_cell(int cell_gid, int corner_lid){
    // corner_lid is 0 to 1
    return (2*cell_gid + corner_lid);
}

KOKKOS_INLINE_FUNCTION
int get_corners_in_node(int node_gid, int corner_lid){
    // corner_lid is 0 to 1
    return (2*node_gid - 1 + corner_lid);
}

KOKKOS_INLINE_FUNCTION
double limiter(double var, double var_avg, double var_n, double max_var, double min_var){
    
    double ratio;
    double half_vel = 0.5*(var + var_n);
    
    if (var > var_avg+fuzz){
        ratio = (max_var - half_vel)/(0.5*(var - var_avg));
    }
    else if (var < var_avg-fuzz){
        ratio = (min_var - half_vel)/(0.5*(var - var_avg));
    }
    else {
        ratio = 1.0;
    }
    
    return fmax(0.0, fmin(1.0, ratio));
    
}



KOKKOS_FUNCTION
void Riemann_solve(DCArrayKokkos <double> &node_coords,
                   DCArrayKokkos <double> &cell_coords,
                   DCArrayKokkos <double> &cell_vel,
                   DCArrayKokkos <double> &cell_vel_n,
                   DCArrayKokkos <double> &node_vel,
                   DCArrayKokkos <double> &cell_pres,
                   DCArrayKokkos <double> &cell_den,
                   DCArrayKokkos <double> &cell_sspd,
                   DCArrayKokkos <double> &corner_vel,
                   DCArrayKokkos <double> &corner_pres,
                   DCArrayKokkos <double> &corner_imp,
                   DCArrayKokkos <double> &corner_force,
                   DCArrayKokkos <double> &cell_dpres_dx,
                   DCArrayKokkos <double> &cell_dvel_dx,
                   const double integral_grad_basis[2],
                   const int num_cells,
                   const int num_nodes){
    
    
    // --- calculate the gradients ---
    
    FOR_ALL (cell_id, 1, num_cells-1, {
        
        double dx = cell_coords(cell_id+1) - cell_coords(cell_id-1);
        
        cell_dpres_dx(cell_id) = (cell_pres(cell_id+1) - cell_pres(cell_id-1))/dx;
        cell_dvel_dx(cell_id)  = (cell_vel(cell_id+1) - cell_vel(cell_id-1))/dx;
        
        //cell_dpres_dx(cell_id) = 0.0;  // WARNING WARNING WARING P0 reconstructions
        //cell_dvel_dx(cell_id)  = 0.0;
        
    });
    RUN({
        cell_dpres_dx(0) = 0.;
        cell_dpres_dx(num_cells-1) = 0.;
        
        cell_dvel_dx(0) = 0.;
        cell_dvel_dx(num_cells-1) = 0.;
        
    });
    Kokkos::fence();
    
    
    
    // --- Calculate corner velocity, forces, and impedance ---
    
    FOR_ALL (cell_id, 0, num_cells, {
        
        int node_id_0 = cell_id;   // global id for node on the left
        int node_id_1 = cell_id+1; // global id for node on the right
        
        // calculate the velocity and pressure in the cell corners
        
        // get the corners id's
        int corner_id_0 = get_corners_in_cell(cell_id, 0); // left
        int corner_id_1 = get_corners_in_cell(cell_id, 1); // right
        
        // calculate the distance to the nodes
        double dx_0 = node_coords(node_id_0) - cell_coords(cell_id);
        double dx_1 = node_coords(node_id_1) - cell_coords(cell_id);
        
        
        // calculate the corner velocities on the left and right side
        corner_vel(corner_id_0) = cell_vel(cell_id) + dx_0*cell_dvel_dx(cell_id);
        corner_vel(corner_id_1) = cell_vel(cell_id) + dx_1*cell_dvel_dx(cell_id);
        
        corner_pres(corner_id_0) = cell_pres(cell_id) + dx_0*cell_dpres_dx(cell_id);
        corner_pres(corner_id_1) = cell_pres(cell_id) + dx_1*cell_dpres_dx(cell_id);
        
        
        // limit the velocity to be in bounds of the Riemann velocity
        
        // velocity bounds based on node_vel, the Riemann velocity
        double max_vel = fmax(node_vel(cell_id), node_vel(cell_id+1));
        double min_vel = fmin(node_vel(cell_id), node_vel(cell_id+1));
        
        // limiter value for the left side
        double phi_0 = limiter(corner_vel(corner_id_0),
                               cell_vel(cell_id),
                               cell_vel_n(cell_id),
                               max_vel,
                               min_vel);
        
        // limiter value for the right side
        double phi_1 = limiter(corner_vel(corner_id_1),
                               cell_vel(cell_id),
                               cell_vel_n(cell_id),
                               max_vel,
                               min_vel);
        
        
        // limit the velocity reconstruction to left and right side
        double phi = fmin(phi_0,phi_1);
        corner_vel(corner_id_0) = cell_vel(cell_id) + phi*(corner_vel(corner_id_0) - cell_vel(cell_id));
        corner_vel(corner_id_1) = cell_vel(cell_id) + phi*(corner_vel(corner_id_1) - cell_vel(cell_id));
        
        
        // left and right corner impedance
        corner_imp(corner_id_0) = cell_den(cell_id)*cell_sspd(cell_id);
        corner_imp(corner_id_1) = cell_den(cell_id)*cell_sspd(cell_id);
        
        // left outward corner force, must add mu*(u^* - u)
        corner_force(corner_id_0) =
                 -integral_grad_basis[0]*(-corner_pres(corner_id_0));
        
        // right outward corner force, must add mu*(u^* - u)
        corner_force(corner_id_1) =
                 -integral_grad_basis[1]*(-corner_pres(corner_id_1));
        
    }); // end parallel for on device
    Kokkos::fence();
    
    
    
    // --- solve for Riemann velocity, the node vel ---
    FOR_ALL (node_id, 1, num_nodes-1, {
        
        
        double top = 0.0;
        double bottom = 0.0;
        
        for (int corner_lid=0; corner_lid<2; corner_lid++){
            
            // get the global index for the corner
            int corner_gid = get_corners_in_node(node_id, corner_lid);
            
            top += corner_imp(corner_gid)*corner_vel(corner_gid) - corner_force(corner_gid);
            bottom += corner_imp(corner_gid);
        } // end for
        
        node_vel(node_id) = top/bottom;
        
    }); // end parallel for on device
    
    // applying a wall BC on velocity
    RUN ({
        node_vel(0) = 0.0;
        node_vel(num_nodes-1) = 0.0;
    });  // end run once on the device
    // Warning: free surface BC requires including node_vel calc,
    Kokkos::fence();
    
    // calculate corner forces
    FOR_ALL (node_id, 0, num_nodes, {
        
        double conservation_check = 0;
        for (int corner_lid=0; corner_lid<2; corner_lid++){
            
            // get the global index for the corner
            int corner_gid = get_corners_in_node(node_id, corner_lid);
            
            corner_force(corner_gid) += corner_imp(corner_gid)*
                                        (node_vel(node_id)-corner_vel(corner_gid));
            
            conservation_check += corner_force(corner_gid);

        } // end for
        
        if(node_id>0 && node_id<num_nodes-1){
            if (fabs(conservation_check)>1e-13) printf("conservation error = %f\n", conservation_check);
        }
    }); // end parallel for on device
    Kokkos::fence();
}
