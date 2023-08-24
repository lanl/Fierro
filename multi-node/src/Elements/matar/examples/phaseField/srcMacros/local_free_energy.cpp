#include "local_free_energy.h"


double calculate_total_free_energy(int* nn, double* delta, double kappa, DCArrayKokkos<double> &comp)
{
    // this function calculates the total free energy of the system.

    // unpack simimulation parameters needed 
    // for calculations in this function
    int nx    = nn[0];
    int ny    = nn[1];
    int nz    = nn[2];
    double dx = delta[0];
    double dy = delta[1];
    double dz = delta[2];

    // 
    double total_energy = 0.0;
    double loc_sum = 0.0;
    REDUCE_SUM(i, 1, nx-1,
               j, 1, ny-1,
               k, 1, nz-1,
               loc_sum, {
           // central difference spatial derivative of comp 
           double dcdx = (comp(i+1,j,k) - comp(i-1,j,k)) / (2.0 * dx);
           double dcdy = (comp(i,j+1,k) - comp(i,j-1,k)) / (2.0 * dy);
           double dcdz = (comp(i,j,k+1) - comp(i,j,k-1)) / (2.0 * dz);
           loc_sum +=   comp(i,j,k) * comp(i,j,k) * (1.0 - comp(i,j,k)) * (1.0 - comp(i,j,k))
                      + 0.5 * kappa * (dcdx * dcdx + dcdy * dcdy + dcdz * dcdz);
               }, total_energy);

    return total_energy;
}

void calculate_dfdc(int* nn, DCArrayKokkos<double> &comp, CArrayKokkos<double> &dfdc)
{
    // this function calculates the derivitive of local free energy density (f) 
    // with respect to composition (c) (df/dc).

    // unpack simimulation parameters needed 
    // for calculations in this function
    int nx = nn[0];
    int ny = nn[1];
    int nz = nn[2];
    
    FOR_ALL(i, 0, nx, 
            j, 0, ny,
            k, 0, nz,{
        dfdc(i,j,k) =   4.0 * comp(i,j,k) * comp(i,j,k) * comp(i,j,k)
                      - 6.0 * comp(i,j,k) * comp(i,j,k)
                      + 2.0 * comp(i,j,k);
    });
}
