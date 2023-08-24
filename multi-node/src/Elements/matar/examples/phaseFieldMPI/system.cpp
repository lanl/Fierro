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


#include "system.h"
#include "stdlib.h"
#include "string"

System::System(MPI_Comm comm_, const SimParameters & sp_) :
comm(comm_),
my_rank(heffte::mpi::comm_rank(comm)),
num_ranks(heffte::mpi::comm_size(comm)),
sp(sp_),
fft(comm, sp.nn),
ga(sp.nn, fft.localRealBoxSizes[my_rank]),
ca(sp, fft.localComplexBoxSizes[my_rank], fft.myComplexBox.low),
total_free_energy_file(NULL),
vtk_writer(comm, fft.globalRealBoxSize, fft.localRealBoxSizes[my_rank], fft.localRealBoxes[my_rank].low, "%12.6E\n")
{
    // print simulation parameters
    if (root == my_rank) 
    {
        sp.print();
        total_free_energy_file = fopen("total_free_energy.csv", "w");
    }
}

System::~System()
{
    if (root == my_rank)
    {
        fclose(total_free_energy_file);
    }
}

void System::initialize_comp()
{
    // start a non-blocking recieve
    MPI_Request request;
    MPI_Status status;
    MPI_Irecv(ga.comp.host_pointer(), ga.comp.size(), fft.mpiType, root, 999, comm, &request);
 
    if (root == my_rank)
    {
        // seed random number generator
        srand(sp.iseed);

        for (int k = 0; k < ga.comp_all.dims(0); ++k) {
        for (int j = 0; j < ga.comp_all.dims(1); ++j) {
        for (int i = 0; i < ga.comp_all.dims(2); ++i) {
            // random number between 0.0 and 1.0
            double r = (double) rand()/RAND_MAX;

            // initialize "comp" with stochastic thermal fluctuations
            ga.comp_all(k,j,i) = sp.c0 + (2.0*r - 1.0)*sp.noise;
        }
        }
        }

        // send subarrays to ranks
        MPI_Request requests[num_ranks];
        MPI_Status statuses[num_ranks];
        for (size_t i = 0; i < num_ranks; i++)
        {
            MPI_Isend(ga.comp_all.pointer(), 1, fft.mpiSubarrayTypes[i], i, 999, comm, &requests[i]);
        }

        // wait for all other messages to be sent
        MPI_Waitall(num_ranks, requests, statuses);
    }

    // wait for non-blocking recieves to complete
    MPI_Wait(&request, &status);

    // update device
    ga.comp.update_device();
}

void System::calculate_dfdc()
{
    // this function calculates the derivitive of local free energy density (f) 
    // with respect to composition (c) (df/dc).

    FOR_ALL(k, 0, ga.dfdc.dims(0),
            j, 0, ga.dfdc.dims(1),
            i, 0, ga.dfdc.dims(2), {
        ga.dfdc(k,j,i) =   4.0 * ga.comp(k,j,i) * ga.comp(k,j,i) * ga.comp(k,j,i)
                         - 6.0 * ga.comp(k,j,i) * ga.comp(k,j,i)
                         + 2.0 * ga.comp(k,j,i);
    });
    Kokkos::fence();
}

double System::calculate_total_free_energy()
{
    // this function calculates the total free energy of the system.

    // unpack simimulation parameters needed 
    // for calculations in this function
    double dx = sp.delta[0];
    double dy = sp.delta[1];
    double dz = sp.delta[2];
    double kappa = sp.kappa;

    // 
    double total_energy = 0.0;
    double loc_sum = 0.0;

#if 0
    // bulk free energy + interfacial energy
    REDUCE_SUM(k, 1, ga.comp.dims(0)-1,
               j, 1, ga.comp.dims(1)-1,
               i, 1, ga.comp.dims(2)-1,
               loc_sum, {
           // central difference spatial derivative of comp 
           double dcdz = (ga.comp(k+1,j,i) - ga.comp(k-1,j,i)) / (2.0 * dz);
           double dcdy = (ga.comp(i,j+1,k) - ga.comp(i,j-1,k)) / (2.0 * dy);
           double dcdx = (ga.comp(k,j,i+1) - ga.comp(k,j,i-1)) / (2.0 * dx);
           loc_sum += ga.comp(k,j,i) * ga.comp(k,j,i)
                    * (1.0 - ga.comp(k,j,i)) * (1.0 - ga.comp(k,j,i))
                    + 0.5 * kappa * (dcdx * dcdx + dcdy * dcdy + dcdz * dcdz);
    }, total_energy);
#endif

    // bulk free energy only
    REDUCE_SUM(k, 0, ga.comp.dims(0),
               j, 0, ga.comp.dims(1),
               i, 0, ga.comp.dims(2),
               loc_sum, {
           loc_sum += ga.comp(k,j,i) * ga.comp(k,j,i) * (1.0 - ga.comp(k,j,i)) * (1.0 - ga.comp(k,j,i));
    }, total_energy);

    return total_energy;
}

void System::time_march()
{
    // get foward fft of comp
    Profile::start_barrier(Profile::fft_forward);
    fft.forward(ga.comp.device_pointer(), ca.comp_img.device_pointer());
    Profile::stop_barrier(Profile::fft_forward);

    // get foward fft of dfdc
    Profile::start_barrier(Profile::fft_forward);
    fft.forward(ga.dfdc.device_pointer(), ca.dfdc_img.device_pointer());
    Profile::stop_barrier(Profile::fft_forward);
    Kokkos::fence();

    // solve Cahn Hilliard equation in fourier space
    FOR_ALL(k, 0, ca.comp_img.dims(0),
            j, 0, ca.comp_img.dims(1),
            i, 0, ca.comp_img.dims(2), {
        ca.comp_img(k,j,i,0) =   (ca.comp_img(k,j,i,0) - (sp.dt * sp.M * ca.kpow2(k,j,i)) * ca.dfdc_img(k,j,i,0))
                               / (ca.denominator(k,j,i));
    
        ca.comp_img(k,j,i,1) =   (ca.comp_img(k,j,i,1) - (sp.dt * sp.M * ca.kpow2(k,j,i)) * ca.dfdc_img(k,j,i,1))
                               / (ca.denominator(k,j,i));
    });
    Kokkos::fence();
    
    // get backward fft of comp_img (note fft.backward was set to scale the result already.
    // you can chnage if needed in FFT3D_R2C class)
    Profile::start_barrier(Profile::fft_backward);
    fft.backward(ca.comp_img.device_pointer(), ga.comp.device_pointer());
    Profile::stop_barrier(Profile::fft_backward);
    Kokkos::fence();
}

void System::track_progress(int iter)
{
    // sum of comp field
    double sum_comp = 0.0;
    double loc_sum = 0.0;
    REDUCE_SUM(k, 0, ga.comp.dims(0),
               j, 0, ga.comp.dims(1),
               i, 0, ga.comp.dims(2),
               loc_sum, {
                   loc_sum += ga.comp(k,j,i);
    }, sum_comp);

    // max of comp field
    double max_comp;
    double loc_max;
    REDUCE_MAX(k, 0, ga.comp.dims(0),
               j, 0, ga.comp.dims(1),
               i, 0, ga.comp.dims(2),
               loc_max, {
                   if(loc_max < ga.comp(k,j,i)){
                       loc_max = ga.comp(k,j,i);
                   }    
    }, max_comp);

    // min of comp field
    double min_comp;
    double loc_min;
    REDUCE_MIN(k, 0, ga.comp.dims(0),
               j, 0, ga.comp.dims(1),
               i, 0, ga.comp.dims(2),
               loc_min, {                  
                   if(loc_min > ga.comp(k,j,i)){
                       loc_min = ga.comp(k,j,i);
                   }                   
    }, min_comp);

    double glob_sum_comp = 0.0;
    double glob_max_comp;
    double glob_min_comp;
    MPI_Reduce(&sum_comp, &glob_sum_comp, 1, MPI_DOUBLE, MPI_SUM, root, comm);
    MPI_Reduce(&max_comp, &glob_max_comp, 1, MPI_DOUBLE, MPI_MAX, root, comm);
    MPI_Reduce(&min_comp, &glob_min_comp, 1, MPI_DOUBLE, MPI_MIN, root, comm);

    if (root == my_rank)
    {
        printf("\n----------------------------------------------------\n");
        printf("Iteration : %d\n", iter);
        printf("Conservation of comp : %E\n", glob_sum_comp);
        printf("Max comp : %E\n", glob_max_comp);
        printf("Min comp : %E\n", glob_min_comp);
    }
}

void System::output_total_free_energy(int iter)
{
    // get total_free_energy
    double total_free_energy = calculate_total_free_energy();

    double glob_total_free_energy = 0.0;
    MPI_Reduce(&total_free_energy, &glob_total_free_energy, 1, MPI_DOUBLE, MPI_SUM, root, comm);

    if (root == my_rank)
        fprintf(total_free_energy_file, "%i,%12.6E\n", iter, glob_total_free_energy);
}

void System::solve()
{
    Profile::start_barrier(Profile::total);

    initialize_comp();

    // time stepping loop
    for (int iter = 1; iter <= sp.num_steps; iter++) {
        // calculate df/dc
        calculate_dfdc();

        // Cahn Hilliard equation time step
        time_march();

        // report simulation progress and output vtk files
        if (iter % sp.print_rate == 0) {

            track_progress(iter);

            output_total_free_energy(iter);

            ga.comp.update_host();
            vtk_writer.write(iter, ga.comp.host_pointer());
        }
    }

    Profile::stop_barrier(Profile::total);
    if (root == my_rank) {
        Profile::print();
    }
}
