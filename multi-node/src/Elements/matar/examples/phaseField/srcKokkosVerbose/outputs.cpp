#include <stdio.h>

#include "outputs.h"
#include "local_free_energy.h"

void track_progress(int iter, int* nn, DCArrayKokkos<double> &comp)
{
    // unpack simimulation parameters needed 
    // for calculations in this function
    int nx = nn[0];
    int ny = nn[1];
    int nz = nn[2];

    // sum of comp field
    double sum_comp = 0.0;
    double loc_sum = 0.0;
    Kokkos::parallel_reduce(
        Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0}, {nx,ny,nz}),
        KOKKOS_LAMBDA(const int i, const int j, const int k, double& loc_sum){
            loc_sum += comp(i,j,k);
        }, sum_comp);


    // max of comp field
    double max_comp;
    double loc_max;
    Kokkos::parallel_reduce(
        Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0}, {nx,ny,nz}),
        KOKKOS_LAMBDA(const int i, const int j, const int k, double& loc_max){
            if(loc_max < comp(i,j,k)) loc_max = comp(i,j,k);
        },
        Kokkos::Max<double>(max_comp)
    );


    // min of comp field
    double min_comp;
    double loc_min;
    Kokkos::parallel_reduce(
        Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0}, {nx,ny,nz}),
        KOKKOS_LAMBDA(const int i, const int j, const int k, double& loc_min){
            if(loc_min > comp(i,j,k)) loc_min = comp(i,j,k);
        },
        Kokkos::Min<double>(min_comp)
    );


    printf("\n----------------------------------------------------\n");
    printf("Iteration : %d\n", iter);
    printf("Conservation of comp : %f\n", sum_comp);
    printf("Max comp : %f\n", max_comp);
    printf("Min comp : %f\n", min_comp);

}



void write_vtk(int iter, int* nn, double* delta, DCArrayKokkos<double> &comp)
{
    // unpack simimulation parameters needed 
    // for calculations in this function
    int nx    = nn[0];
    int ny    = nn[1];
    int nz    = nn[2];
    double dx = delta[0];
    double dy = delta[1];
    double dz = delta[2];

    // update host copy of comp
    comp.update_host();

    // output file management
    FILE* output_file;
    char filename[50];

    // create name of output vtk file
    sprintf(filename, "outputComp_%d.vtk", iter);

    // open output vtk file
    output_file = fopen(filename, "w");

    // write vtk file heading
    fprintf(output_file, "%s\n", "# vtk DataFile Version 2.0");
    fprintf(output_file, "%s\n", filename);
    fprintf(output_file, "%s\n", "ASCII");
    fprintf(output_file, "%s\n", "DATASET STRUCTURED_GRID");
    fprintf(output_file, "%s %i  %i  %i\n", "DIMENSIONS", nx, ny, nz);
    fprintf(output_file, "%s   %i   %s\n", "POINTS", nx*ny*nz, "double");

    // write grid point values 
    // Note: order of for loop is important (k,j,i)
    double x, y, z;
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {           
                x = double(i)*dx;
                y = double(j)*dy;
                z = double(k)*dz;
                fprintf(output_file, "  %12.6E  %12.6E  %12.6E\n", x, y, z);
            }
        }
    }

    // write data values
    // Note: order of for loop is important (k,j,i)
    fprintf(output_file, "%s %i\n", "POINT_DATA", nx*ny*nz);
    fprintf(output_file, "%s\n", "SCALARS data double");
    fprintf(output_file, "%s\n", "LOOKUP_TABLE default");
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {         
                fprintf(output_file, " %12.6E\n", comp.host(i,j,k));
            }
        }
    }
    
    // close file
    fclose(output_file);
}



void output_total_free_energy(int iter, int print_rate, int num_steps, int* nn, 
                              double* delta, double kappa, DCArrayKokkos<double> &comp)
{
    // get total_free_energy
    double total_free_energy = calculate_total_free_energy(nn, delta, kappa, comp);

    // output file management
    static FILE* output_file;
    static char filename[50];

    // open output vtk file
    if (iter == print_rate)
    {
        // create name of output vtk file
        sprintf(filename, "total_free_energy.csv");
        output_file = fopen(filename, "w");
    }

    // write total_free_energy to file
    fprintf(output_file, "%i,%12.6E\n", iter, total_free_energy);

    // close file
    if (iter == num_steps)
    {
        fclose(output_file);
    }
}

