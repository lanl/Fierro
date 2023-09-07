This folder shows the correct and incorrect use of MATAR arrays and the consequences on performance on both the CPU and GPU. The numerical solution of the Laplace equation for steady-state temperature distribution using Jacobi iteration is used as a case study. The results of the tests are detailed in "report.pdf".

To run the Laplace solver routine (that is based on the c_indexing convention) using the Kokkos thread parallelization backend with a different number threads, the syntax is
./carraykokkos_c_indexing --kokkos-threads=4
where the number 4 is the number of threads.
