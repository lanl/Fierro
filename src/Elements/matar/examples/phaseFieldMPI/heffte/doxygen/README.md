Documentation
=============

We provide a detailed [Online Documentation](https://mkstoyanov.bitbucket.io/heffte/) of all HeFFTe classes and functions.


Installing heFFTe
=================

* [HeFFTe Installation](https://mkstoyanov.bitbucket.io/heffte/md_doxygen_installation.html): Using CMake or standard GNU Make.

Testing installation
====================

To ensure HeFFTe was properly built, we provide several tests for all kernels. Using Cake functionality, simply do as follows:

~~~
cd heffte/build
    ctests -V
~~~

Performance benchmarks
======================

Once HeFFTe is built, several tests are available in folder `heffte/build/test/`. These tests allow to evaluate correctness and performance, and they should be used to validate new developments. 

To evaluate scalability and make performance comparison with other parallel FFT libraries, refer to folder `heffte/build/benchmarks/`, where you will find two executables: `speed3d_c2c` for complex-complex transforms, and `speed3d_r2c` for real-to-complex transforms. To run these tests on an MPI supported cluster, follow the examples:

~~~
mpirun -n 12 ./speed3d_r2c fftw double 512 256 512 -p2p -pencils -no-reorder
mpirun -n 5 --map-by node  ./speed3d_c2c mkl float 1024 256 512  -a2a -slabs -reorder
mpirun -n 2 ./speed3d_c2c cufft double 512 256 512  -mps -a2a
~~~

Should you have questions about the use of flags, please refer to `flags.md` for detailed information. For systems, such as Summit supercomputer, which support execution with `jsrun` by default, follow the examples:

~~~
jsrun  -n1280 -a1 -c1 -r40 ./speed3d_r2c fftw double 1024 256 512 -pencils 
jsrun --smpiargs="-gpu" -n192 -a1 -c1 -g1 -r6 ./speed3d_c2c cufft double 1024 1024 1024 -p2p -reorder
~~~

For comparison to other libraries, make sure to use equivalent flags. Some libraries only provide benchmarks for evaluating FFT performance starting and ending at a pencils-shaped FFT grids. For such cases, use the flag `-io_pencils`.

 Refer to [HeFFTe Papers](https://www.icl.utk.edu/publications/fft-ecp-fast-fourier-transform) for scalability and performance speedup results.
