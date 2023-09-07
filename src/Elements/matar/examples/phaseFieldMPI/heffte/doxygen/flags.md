HEFFTE Flags
============

Benchmarks are available in the *benchmarks* folder. The best choice of parameters are tunned for each backend and available by default. Parameters can also be manually defined by the user via the following flags:

- `-reorder`: reorder the elements of the arrays so that each 1-D FFT will use contiguous data.

- `-no-reorder`: some of the 1-D will be strided (non contiguous).

- `-a2a`: use MPI_Alltoallv() communication method.

- `-p2p`: use MPI_Send() and MPI_Irecv() communication methods.

- `-pencils`: use pencil reshape logic.

- `-slabs`: use slab reshape logic.

- `-io_pencils`: if input and output proc grids are pencils, useful for comparison with other libraries.

- `-mps`: for the CUFFT backend and multiple GPUs, it associates the mpi ranks with different cuda devices, using CudaSetDevice(my_rank%device_count). It is deactivated by default.
