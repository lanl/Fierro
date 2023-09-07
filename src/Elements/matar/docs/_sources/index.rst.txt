Introduction to MATAR 
======================================================

.. toctree::

  about 
  api/library_root
  

Indices and tables
==================
  * :ref:`genindex`
  * :ref:`modindex`
  * :ref:`search`


Table of Contents
^^^^^^^^^^^^^^^^^

.. image:: https://github.com/lanl/MATAR/blob/main/MATAR_Logo.png
   :width: 350

MATAR is a C++ library that addresses the need for simple, fast, and memory-efficient multi-dimensional data representations for dense and sparse storage that arise with numerical methods and in software applications. The data representations are designed to perform well across multiple computer architectures, including CPUs and GPUs. MATAR allows users to easily create and use intricate data representations that are also portable across disparate architectures using Kokkos. The performance aspect is achieved by forcing contiguous memory layout (or as close to contiguous as possible) for multi-dimensional and multi-size dense or sparse MATrix and ARray (hence, MATAR) types. Results show that MATAR has the capability to improve memory utilization, performance, and programmer productivity in scientific computing. This is achieved by fitting more work into the available memory, minimizing memory loads required, and by loading memory in the most efficient order. 


## Examples
* [ELEMENTS](https://github.com/lanl/ELEMENTS/):   MATAR is a part of the ELEMENTS Library (LANL C# C20058) and it underpins the routines implemented in ELEMENTS.  MATAR is available in a stand-alone directory outside of the ELEMENTS directory because it can aid many code applications.  The dense and sparse storage types in MATAR are the foundation for the ELEMENTS library, which contains mathematical functions to support a very broad range of element types including: linear, quadratic, and cubic serendipity elements in 2D and 3D; high-order spectral elements; and a linear 4D element. An unstructured high-order mesh class is available in ELEMENTS and it takes advantage of MATAR for efficient access of various mesh entities. 

* [Fierro](https://github.com/lanl/Fierro): The MATAR library underpins the Fierro code that is designed to simulate quasi-static solid mechanics problems and material dynamics problems.  
    
* Simple examples are in the /test folder

## Descriptions

* All Array MATAR types (e.g., CArray, ViewCArray, FArray, RaggedRightArray, etc.) start with an index of 0 and stop at an index of N-1, where N is the number of entries.  

* All Matrix MATAR types  (e.g., CMatrix, ViewCMatrix, FMatrix, etc.)  start with an index of 1 and stop at an index of N, where N is the number of entries. 

* The MATAR View types (e.g., ViewCArray, ViewCMatrix, ViewFArray, etc. ) are designed to accept a pointer to an existing 1D array and then access that 1D data as a multi-dimensional array.  The MATAR View types can also be used to slice an existing View.  

* The C dense storage and View types (e.g., CArray, ViewCArray, CMatrix, etc.) access the data following the C/C++ language convection of having the last index in a multi-dimensional array vary the quickest.  In a 2D CArray A, the index j in A(i,j) varies first followed by the index i, so the optimal performance is achieved using the following loop ordering.

  
 ```
  // Optimal use of CArray
  for (i=0,i<N,i++){
      for (j=0,j<N,j++){
          A(i,j) = 0.0;
      }
   }
  ```

* The F dense storage and View types (e.g., FArray, ViewFArray, FMatrix, etc.) access the data following the Fortran language convection of having the first index in a multi-dimensional array vary the quickest.  In a 2D FMatrix M, the index i in M(i,j) varies first followed by the index j, so the optimal performance is achieved using the following loop ordering.

  ```
// Optimal use of FMatrix
for (j=1,j<=N,j++){
    for (i=1,i<=N,i++){
        M(i,j) = 0.0;
    }
}
  ```

* The ragged data types (e.g., RaggedRightArray, RaggedDownArray, etc) in MATAR are special sparse storage types.  The Right access types are for R(i,j) where the number of column entries varies in width across the array.  The Down access types are for D(i,j) where the number of row entries vary in length across the array.

* The SparseRowArray MATAR type is the idetical to the Compressed Sparse Row (CSR) or Compressed Row Storage (CSR) respresentation.

* The SparseColumnArray MATAR type is identical to the Compressed Sparse Column (CSC) or Compressed Column Storage (CCS) respresentation.


## Usage
```
// create a 1D array of integers and then access as a 2D array
int A[9];
auto A_array = ViewCArray <int> (A, 3, 3); // access as A(i,j) 

// create a 3D array of doubles
auto B = CArray <double> (3,3,3); // access as B(i,j,k)

// create a slice of the 3D array at index 1
auto C = ViewCArray <double> (&B(1,0,0),3,3); // access as C(j,k)


// create a 4D matrix of doubles, indices start at 1 
auto D = CMatrix <double> (10,9,8,7); // access as D(i,j,k,l)


// create a 2D view of a standard array
std::array<int, 9> E1d;
auto E = ViewCArray<int> (&E1d[0], 3, 3);
E(0,0) = 1;  // and so on


// create a ragged-right array of integers
//
// [1, 2, 3]
// [4, 5]
// [6]
// [7, 8, 9, 10]
//
size_t my_strides[4] = {3, 2, 1, 4};
RaggedRightArray <int> ragged(my_strides, 4);
    
int value = 1;
for (int i=0; i<4; i++){
    for (int j=0; j<my_ragged.stride(i); j++){
        ragged(i,j) = value;
        value++;
    }
}


```
    
## Cloning the code
If your SSH keys are set in github, then from the terminal type:
```
git clone --recursive ssh://git@github.com/lanl/MATAR.git    
```

## Basic build
The basic build is for users only interested in the serial CPU only MATAR data types.  For this build, we recommend making a folder perhaps called build then go into the build folder and type
```
cmake ..
make
```
The compiled code will be in the build folder.

## Debug basic build 

To build serial CPU only MATAR data types in the debug mode, please use
```
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```
The debug flag includes checks on array and matrix dimensions and index bounds.


## Building MATAR with Kokkos
A suite of build scripts are provided to build MATAR with Kokkos for performance portability across computer architectures (CPUs and GPUs).  The scripts for various Kokkos backends (e.g., CUDA, HIP, OpenMP, and pthreads) are located within the scripts folder.  The provided scripts are configured for particular hardware, the user will likely need to alter the inputs to reflect their hardware.  There are three scripts in each folder that are sourced to build MATAR with Kokkos.  The scripts are
```
sourceme-env.sh
kokkos-install.sh
backend-cmake-build.sh
```
The word backend denotes cuda, hip, openMP, and so forth.  Scripts are also provided to build MATAR without Kokkos, and in that case there is no backend listed since it doesn't use Kokko.  The backend-cmake-build.sh script will run cmake and make for the project.  Afterwords, the user can just runs make inside the respective build directory to compile the project.  For clarity, running all the scripts is only necessary to set up and compile the code the first time, afterwards, the use can compile the code using make in the build directory.  The environment variables will need to be set when logging into a compute node or when changing to a different kokkos backend. For all builds, a single script is provided in each script folder to automate the entire build process, it runs the three aforementioned scripts sequentially. 
```
build-it.sh
```
Before using the build-it.sh script, the user must verify that the settings in the other scripts that build MATAR with a Kokkos backend are correctly set.  After running the build-it.sh script, the entire project is compiled and stored in a directory that is named with the respective Kokkos backend e.g., build-kokkos-cuda.  Further details are provided on the three scripts to configure and build MATAR with a Kokkos backend.


### Environment configuration script
To start, the environment variables and modules must be configured by sourcing the following script
```
source sourceme-env.sh
```
This script is where the user will load the necessary module files for their given machine/architecture combination.  This script also creates the build directory for the project e.g., build-kokkos-cuda, build-kokkos-hip, build-kokkos-openmp, etc.


### Install Kokkos script
The next step is to install Kokkos, using the version that was cloned recursively within MATAR, and configure the Kokkos build for specific hardware and a backend.  
```
source kokkos-install.sh
```
Within this script, the user will need to set any Kokkos specific variables for their project. The architecture variables will need to be modified based on the architecture being used. The provided scripts are set for a particular hardware that might differ from what a user might be using.  CPU architecture information needs to be listed if running with the Kokkos serial or OpenMP backends; GPU architecture information must be listed if using a Kokkos GPU backend.  We refer the user to Kokkos compiling page to see the large list of compilation options,
```
https://github.com/kokkos/kokkos/wiki/Compiling
```



### CUDA compilation script
To build the project with cuda, the last step is to type
```
source cuda-cmake-build.sh
```


### HIP compilation script
To build the project with hip, the last step is to type
```
source hip-cmake-build.sh
```


### openMP compilation script
To build the project with openMP, the last step is to type

```
source openmp-cmake-build.sh
```
The sourceme-env.sh script (the first step) sets the number of threads to 16 by default.  Changing the number of threads used with openMP requires manually setting the environment variable OMP_NUM_THREADS.  


    
### pthreads compilation script
To build the project with ptheads, the last step is to type
```
source pthreads-cmake-build.sh
```    
To specify number of threads when running a code with the Kokkos pthread backend, add the following command line arguments
```
--kokkos-threads=4
```


### Automate build process
A build-it.sh script is provided that runs all scripts sequentially for the user.  The build-it.sh script obviates the need to manually source each script.  The user must verify the settings are correct in each script prior to using the build-it.sh script.  If the build-it.sh script fails to build the project correctly, the user should carefully look at the loaded modules and settings for building Kokkos.   



## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.


## License
This program is open source under the BSD-3 License.

## Citation
```
@article{MATAR,
title = "{MATAR: A Performance Portability and Productivity Implementation of Data-Oriented Design with Kokkos}",
journal = {Journal of Parallel and Distributed Computing},
pages = {86-104},
volume = {157},
year = {2021},
author = {Daniel J. Dunning and Nathaniel R. Morgan and Jacob L. Moore and Eappen Nelluvelil and Tanya V. Tafolla and Robert W. Robey},
keywords = {Performance, Portability, Productivity, Memory Efficiency, GPUs, dense and sparse storage}
```


