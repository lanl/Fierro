#include <mpi.h>
#include <matar.h>
#include <stdio.h>
#include <iostream>
#include <chrono>
#include <math.h>
#include <string>

// Dont change ROOT
#define ROOT  0
//----------------

// Change to 0 or 1 as needed
#define TRACK_PROGRESS  0

#if defined HAVE_CUDA || defined HAVE_HIP
  #define GPU 1
#else
  #define GPU 0
#endif

using namespace mtr; // matar namespace

int width = 1000;
int height = 1000;
int max_num_iterations = 1000;
double temp_tolerance = 0.01;

void initialize(DCArrayKokkos<double> &temperature_previous, int height, int width);
void track_progress(int iteration, DCArrayKokkos<double> &temperature);
void parse_command_line(int argc, char *argv[]);


int main(int argc, char *argv[])
{

  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  { // kokkos scope

  // Parse command line options
  parse_command_line(argc, argv);

  // start timing total code 
  double begin_time_total = MPI_Wtime();

  int world_size,
      rank,
      width_loc,
      height_loc,
      size_loc;

  CArray <int> all_size_loc;
  CArray <int> offsets;

  // get world_size and rank
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  // divide work along the height 
  // Note: +2 is added for boundary
  width_loc = width+2;
  height_loc = (height+2) / world_size;
  if (rank < ((height+2) % world_size)) {
    height_loc++;
  }
  size_loc = width_loc * height_loc;

  // root should keep an array of size_loc and offset
  // for all processes
  if (rank == ROOT) {
    all_size_loc = CArray <int> (world_size);
    offsets = CArray <int> (world_size);

    all_size_loc(ROOT) = size_loc;
    for (int i = 1; i < world_size; i++) {
      MPI_Recv(&all_size_loc(i), 1, MPI_INT, i,
               MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  
    offsets(0) = 0;
    for (int i = 1; i < world_size; i++) {
      offsets(i) = offsets(i-1) + all_size_loc(i-1);
    }

  } else {
    MPI_Send(&size_loc, 1, MPI_INT, ROOT, 999, MPI_COMM_WORLD);
  }

  // declare arrays
  CArrayKokkos <double> temperature_loc;
  DCArrayKokkos <double> temperature_previous_loc;
  DCArrayKokkos <double> temperature_previous_glob;

  // define and allocate arrays
  temperature_loc = CArrayKokkos <double> (height_loc, width_loc);
  temperature_previous_loc = DCArrayKokkos <double> (height_loc, width_loc);
  if (rank == ROOT) {
    temperature_previous_glob = DCArrayKokkos <double> (height+2, width+2);
  }

  // initialize temperature field.
  if (rank == ROOT) {
    initialize(temperature_previous_glob, height, width);
  }

  // distribut work to all processes
  if (rank == ROOT) {
    temperature_previous_glob.update_host();
    MPI_Scatterv(&temperature_previous_glob.host(0,0), &all_size_loc(0), &offsets(0),
                 MPI_DOUBLE, &temperature_previous_loc.host(0,0), size_loc, MPI_DOUBLE,
                 ROOT, MPI_COMM_WORLD);
  } else {
    MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, &temperature_previous_loc.host(0,0), size_loc,
                 MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  }
  //
  temperature_previous_loc.update_device();

  // define neighbours
  int up = rank-1;
  int down = rank+1;
  int up_tag = 888;
  int down_tag = 999;
  int nrequests;
  MPI_Request requests_send[4];
  MPI_Request requests_recv[4];

  DCArrayKokkos <double> halo_up, halo_down;
  DCArrayKokkos <double> halo_up_out, halo_down_out;
  if (up != -1) {
    halo_up = DCArrayKokkos <double> (width_loc);
    halo_up_out = DCArrayKokkos <double> (width_loc);
  }
  if (down != world_size) {
    halo_down = DCArrayKokkos <double> (width_loc);
    halo_down_out = DCArrayKokkos <double> (width_loc);
  }

  int height_index_start, height_index_end;
  if (rank == 0) {
    height_index_start = 1;    
  } else {
    height_index_start = 0;
  }
  if (rank == world_size-1) { 
    height_index_end = height_loc-1;
  } else {
    height_index_end = height_loc;
  }

  //
  int iteration = 1;
  double worst_dt = 100.0;
  double worst_dt_loc;

  double begin_time_main_loop = MPI_Wtime();
  // main loop
  while (worst_dt > temp_tolerance && iteration <= max_num_iterations) {
#if !defined GPU
    // communicate halo nodes
    nrequests = 0;
    if (up != -1) {
      MPI_Irecv(halo_up.device_pointer(), halo_up.size(), MPI_DOUBLE, 
                up, up_tag, MPI_COMM_WORLD, &requests_recv[nrequests]);
      MPI_Isend(temperature_previous_loc.device_pointer()+(0+(0*width_loc)), 
                halo_up_out.size(), MPI_DOUBLE, up, down_tag, MPI_COMM_WORLD, 
                &requests_send[nrequests]);
      nrequests++;
    }

    if (down != world_size) {
      MPI_Irecv(halo_down.device_pointer(), halo_down.size(), MPI_DOUBLE,
                down, down_tag, MPI_COMM_WORLD, &requests_recv[nrequests]);
      MPI_Isend(temperature_previous_loc.device_pointer()+(0+((height_loc-1)*width_loc)),
                halo_down_out.size(), MPI_DOUBLE, down, up_tag, MPI_COMM_WORLD, 
                &requests_send[nrequests]);
      nrequests++;
    }

#else
    // fill halo with data
    if (up != -1) {
      FOR_ALL(j, 0, width_loc, {
        halo_up_out(j) = temperature_previous_loc(0,j);
      });
      halo_up_out.update_host();
    }

    if (down != world_size) {
      FOR_ALL(j, 0, width_loc, {
        halo_down_out(j) = temperature_previous_loc(height_loc-1, j);
      });
      halo_up_out.update_host();
    }
    //
    Kokkos::fence();
    
    // communicate halo nodes
    nrequests = 0;
    if (up != -1) {
      MPI_Irecv(halo_up.host_pointer(), halo_up.size(), MPI_DOUBLE, 
                up, up_tag, MPI_COMM_WORLD, &requests_recv[nrequests]);
      MPI_Isend(halo_up_out.host_pointer(), halo_up_out.size(), MPI_DOUBLE, 
                up, down_tag, MPI_COMM_WORLD, &requests_send[nrequests]);
      nrequests++;
    }

    if (down != world_size) {
      MPI_Irecv(halo_down.host_pointer(), halo_down.size(), MPI_DOUBLE,
                down, down_tag, MPI_COMM_WORLD, &requests_recv[nrequests]);
      MPI_Isend(halo_down_out.host_pointer(), halo_down_out.size(), MPI_DOUBLE,
                down, up_tag, MPI_COMM_WORLD, &requests_send[nrequests]);
      nrequests++;
    }
#endif 

    // finite difference for internal nodes
    FOR_ALL(i, 1, height_loc-1,
            j, 1, width_loc-1, {
        temperature_loc(i,j) = 0.25 * (temperature_previous_loc(i+1,j)
                                    + temperature_previous_loc(i-1,j)
                                    + temperature_previous_loc(i,j+1)
                                    + temperature_previous_loc(i,j-1));
    });

    // Wait for all halo exchange to complete
    if (nrequests > 0) {
      MPI_Waitall(nrequests, requests_send, MPI_STATUSES_IGNORE);
      MPI_Waitall(nrequests, requests_recv, MPI_STATUSES_IGNORE);
    }

    // finite difference on surface nodes
    if (up != -1) {
#if defined GPU
      halo_up.update_device();
#endif
      int i = 0;
      FOR_ALL(j, 1, width_loc-1, {
        temperature_loc(i,j) = 0.25 * (temperature_previous_loc(i+1,j)
                                    + halo_up(j)
                                    + temperature_previous_loc(i,j+1)
                                    + temperature_previous_loc(i,j-1));
      });
    } // end if (up != -1)

    if (down != world_size) {
#if defined GPU
      halo_down.update_device();
#endif
      int i = height_loc-1;
      FOR_ALL(j, 1, width_loc-1, {
        temperature_loc(i,j) = 0.25 * (halo_down(j)
                                    + temperature_previous_loc(i-1,j)
                                    + temperature_previous_loc(i,j+1)
                                    + temperature_previous_loc(i,j-1));
      });
    } // end if (down != world_size)

    // calculate max difference between temperature and temperature_previous
    double loc_max_value = 100.0;
    REDUCE_MAX(i, height_index_start, height_index_end,
               j, 1, width_loc-1,
               loc_max_value, {
      double value = fabs(temperature_loc(i,j) - temperature_previous_loc(i,j));
      if (value > loc_max_value) loc_max_value = value;
    }, worst_dt_loc);

    // update temperature_previous
    FOR_ALL(i, height_index_start, height_index_end,
            j, 1, width_loc-1, {
        temperature_previous_loc(i,j) = temperature_loc(i,j);
    });

    // wait for all kokkos kernals to complete
    Kokkos::fence();

    // all reduce for worst_dt
    MPI_Allreduce(&worst_dt_loc, &worst_dt, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);


#if TRACK_PROGRESS
    // track progress
    if (iteration % 100 == 0) {

      temperature_previous_loc.update_host();

      if (rank == ROOT) {
        MPI_Gatherv(&temperature_previous_loc.host(0,0), size_loc, MPI_DOUBLE, 
                    &temperature_previous_glob.host(0,0), &all_size_loc(0), &offsets(0), 
                    MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

      } else {

        MPI_Gatherv(&temperature_previous_loc.host(0,0), size_loc, MPI_DOUBLE, 
                    NULL, NULL, NULL, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
      }
      
      if (rank == ROOT) {
        track_progress(iteration, temperature_previous_glob);
      }

    } // end if (iteration % 100 == 0)
#endif

    iteration++;
  } // end while loop

  // stop timing
  double end_time = MPI_Wtime();

  if (rank == ROOT) {
    printf("\n");
    printf("Number of MPI processes = %d\n", world_size);
    printf("height = %d; width = %d\n", height, width);
    printf("Total code time was %10.6e seconds.\n", end_time-begin_time_total);
    printf("Main loop time was %10.6e seconds.\n", end_time-begin_time_main_loop);
    printf("Max error at iteration %d was %10.6e\n", iteration-1, worst_dt);
  }


  } // end kokkos scope
  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}


void initialize(DCArrayKokkos<double> &temperature_previous, int height, int width) {
  // initialize temperature_previous to 0.0
  FOR_ALL(i, 0, height+2,
          j, 0, width+2, {
      temperature_previous(i,j) = 0.0;    
  });

  // setting the left and right boundary conditions
  FOR_ALL(i, 0, height+2, {
      temperature_previous(i,0) = 0.0;
      temperature_previous(i,width+1) = (100.0/height)*i;
  });

  // setting the top and bottom boundary condition
  FOR_ALL(j, 0, width+2, {
      temperature_previous(0,j) = 0.0;
      temperature_previous(height+1,j) = (100.0/width)*j; 

  });
}

void track_progress(int iteration, DCArrayKokkos<double> &temperature) { 

  printf("---------- Iteration number: %d ----------\n", iteration);
  for (int i = height-5; i <= height; i++) {
    printf("[%d,%d]: %5.2f  ", i,i, temperature.host(i,i));
  }
  printf("\n");
}

void parse_command_line(int argc, char *argv[])
{
  std::string opt;
  int i = 1;
  while (i < argc && argv[i][0] == '-') 
  {
    opt = std::string(argv[i]);

    if(opt == "-height") 
      height = atoi(argv[++i]);

    if(opt == "-width")
      width = atoi(argv[++i]);

    ++i;
  }
}
