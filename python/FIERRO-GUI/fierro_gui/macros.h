#ifndef MACROS_H
#define MACROS_H
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

/**********************************************************************************************
 This file has suite of MACROS to build serial and parallel loops that are more readable and
 are written with the same syntax. The parallel loops use kokkos (i.e., the MACROS hide the
 complexity) and the serial loops are done using functions located in this file. The goal is to
 help users add kokkos to their code projects for performance portability across architectures.

 The loop order with the MACRO enforces the inner loop varies the fastest and the outer most
 loop varies the slowest.  Optiminal performance will be achieved by ensureing the loop indices
 align with the access pattern of the MATAR datatype.
 
 1.  The syntax to use the FOR_ALL MACRO is as follows:

 // parallelization over a single loop
 FOR_ALL(k, 0, 10,
        { loop contents is here });

 // parallellization over two loops
 FOR_ALL(m, 0, 3,
         n, 0, 3,
        { loop contents is here });

 // parallellization over two loops
 FOR_ALL(i, 0, 3,
         j, 0, 3,
         k, 0, 3,
        { loop contents is here });

 2.  The syntax to use the FOR_REDUCE is as follows:

 // reduce over a single loop
 REDUCE_SUM(i, 0, 100,
            local_answer,
            { loop contents is here }, answer);

 REDUCE_SUM(i, 0, 100,
            j, 0, 100,
            local_answer,
           { loop contents is here }, answer);
 
 REDUCE_SUM(i, 0, 100,
            j, 0, 100,
            k, 0, 100,
            local_answer,
           { loop contents is here }, answer);
 
 // other reduces are: RDUCE_MAX and REDUCE_MIN
 **********************************************************************************************/


#include <stdio.h>
#include <iostream>





// -----------------------------------------
// MACROS used with both Kokkos and non-kokkos versions
// -----------------------------------------
// a macro to select the name of a macro based on the number of inputs
#define \
    GET_MACRO(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, NAME,...) NAME


// -----------------------------------------
// MACROS for kokkos
// -----------------------------------------

#ifdef HAVE_KOKKOS

// CArray nested loop convention use Right, use Left for outermost loop first
#define LOOP_ORDER Kokkos::Iterate::Right

// FArray nested loop convention use Right
#define F_LOOP_ORDER Kokkos::Iterate::Right


// run once on the device
#define \
    RUN(fcn) \
    Kokkos::parallel_for( Kokkos::RangePolicy<> ( 0, 1), \
                          KOKKOS_LAMBDA(const int ijkabc){fcn} )

// run once on the device inside a class
#define \
    RUN_CLASS(fcn) \
    Kokkos::parallel_for( Kokkos::RangePolicy<> ( 0, 1), \
                          KOKKOS_CLASS_LAMBDA(const int ijkabc){fcn} )
              

// the FOR_ALL loop
#define \
    FOR1D(i, x0, x1,fcn) \
    Kokkos::parallel_for( Kokkos::RangePolicy<> ( (x0), (x1)), \
                          KOKKOS_LAMBDA( const int (i) ){fcn} )

#define \
    FOR2D(i, x0, x1, j, y0, y1,fcn) \
    Kokkos::parallel_for( \
        Kokkos::MDRangePolicy< Kokkos::Rank<2,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0)}, {(x1), (y1)} ), \
        KOKKOS_LAMBDA( const int (i), const int (j) ){fcn} )

#define \
    FOR3D(i, x0, x1, j, y0, y1, k, z0, z1, fcn) \
    Kokkos::parallel_for( \
         Kokkos::MDRangePolicy< Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0), (z0)}, {(x1), (y1), (z1)} ), \
         KOKKOS_LAMBDA( const int (i), const int (j), const int (k) ) {fcn} )

#define \
    FOR_ALL(...) \
    GET_MACRO(__VA_ARGS__, _13, _12, _11, FOR3D, _9, _8, FOR2D, _6, _5, FOR1D)(__VA_ARGS__)


// the DO_ALL loop
#define \
    DO1D(i, x0, x1,fcn) \
    Kokkos::parallel_for( Kokkos::RangePolicy<> ( (x0), (x1)+1), \
                          KOKKOS_LAMBDA( const int (i) ){fcn} )

#define \
    DO2D(i, x0, x1, j, y0, y1,fcn) \
    Kokkos::parallel_for( \
        Kokkos::MDRangePolicy< Kokkos::Rank<2,F_LOOP_ORDER, F_LOOP_ORDER> > ( {(x0), (y0)}, {(x1)+1, (y1)+1} ), \
        KOKKOS_LAMBDA( const int (i), const int (j) ){fcn} )

#define \
    DO3D(i, x0, x1, j, y0, y1, k, z0, z1, fcn) \
    Kokkos::parallel_for( \
         Kokkos::MDRangePolicy< Kokkos::Rank<3,F_LOOP_ORDER,F_LOOP_ORDER> > ( {(x0), (y0), (z0)}, {(x1)+1, (y1)+1, (z1)+1} ), \
         KOKKOS_LAMBDA( const int (i), const int (j), const int (k) ) {fcn} )

#define \
    DO_ALL(...) \
    GET_MACRO(__VA_ARGS__, _13, _12, _11, DO3D, _9, _8, DO2D, _6, _5, DO1D)(__VA_ARGS__)


// the REDUCE SUM loop
#define \
    RSUM1D(i, x0, x1, var, fcn, result) \
    Kokkos::parallel_reduce( Kokkos::RangePolicy<> ( (x0), (x1) ),  \
                             KOKKOS_LAMBDA(const int (i), decltype(var) &(var)){fcn}, (result))

#define \
    RSUM2D(i, x0, x1, j, y0, y1, var, fcn, result) \
    Kokkos::parallel_reduce( \
        Kokkos::MDRangePolicy< Kokkos::Rank<2,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0)}, {(x1), (y1)} ), \
        KOKKOS_LAMBDA( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
           (result) )

#define \
    RSUM3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
    Kokkos::parallel_reduce( \
        Kokkos::MDRangePolicy< Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0), (z0)}, {(x1), (y1), (z1)} ), \
        KOKKOS_LAMBDA( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
            (result) )

#define \
    REDUCE_SUM(...) \
    GET_MACRO(__VA_ARGS__, _13, RSUM3D, _11, _10, RSUM2D, _8, _7, RSUM1D)(__VA_ARGS__)


// the DO_REDUCE_SUM loop
#define \
    DO_RSUM1D(i, x0, x1, var, fcn, result) \
    Kokkos::parallel_reduce( Kokkos::RangePolicy<> ( (x0), (x1)+1 ),  \
                             KOKKOS_LAMBDA(const int (i), decltype(var) &(var)){fcn}, (result))

#define \
    DO_RSUM2D(i, x0, x1, j, y0, y1, var, fcn, result) \
    Kokkos::parallel_reduce( \
        Kokkos::MDRangePolicy< Kokkos::Rank<2,F_LOOP_ORDER,F_LOOP_ORDER> > ( {(x0), (y0)}, {(x1)+1, (y1)+1} ), \
        KOKKOS_LAMBDA( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
           (result) )

#define \
    DO_RSUM3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
    Kokkos::parallel_reduce( \
        Kokkos::MDRangePolicy< Kokkos::Rank<3,F_LOOP_ORDER,F_LOOP_ORDER> > ( {(x0), (y0), (z0)}, {(x1)+1, (y1)+1, (z1)+1} ), \
        KOKKOS_LAMBDA( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
            (result) )

#define \
    DO_REDUCE_SUM(...) \
    GET_MACRO(__VA_ARGS__, _13, DO_RSUM3D, _11, _10, DO_RSUM2D, _8, _7, DO_RSUM1D)(__VA_ARGS__)


// the REDUCE MAX loop
#define \
    RMAX1D(i, x0, x1, var, fcn, result) \
    Kokkos::parallel_reduce( \
                        Kokkos::RangePolicy<> ( (x0), (x1) ),  \
                        KOKKOS_LAMBDA(const int (i), decltype(var) &(var)){fcn}, \
                        Kokkos::Max< decltype(result) > ( (result) ) )

#define \
    RMAX2D(i, x0, x1, j, y0, y1, var, fcn, result) \
    Kokkos::parallel_reduce( \
                        Kokkos::MDRangePolicy< Kokkos::Rank<2,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0)}, {(x1), (y1)} ), \
                        KOKKOS_LAMBDA( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
                        Kokkos::Max< decltype(result) > ( (result) ) )

#define \
    RMAX3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
    Kokkos::parallel_reduce( \
                        Kokkos::MDRangePolicy< Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0), (z0)}, {(x1), (y1), (z1)} ), \
                        KOKKOS_LAMBDA( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
                        Kokkos::Max< decltype(result) > ( (result) ) )

#define \
    REDUCE_MAX(...) \
    GET_MACRO(__VA_ARGS__, _13, RMAX3D, _11, _10, RMAX2D, _8, _7, RMAX1D)(__VA_ARGS__)


// the DO_REDUCE_MAX loop
#define \
    DO_RMAX1D(i, x0, x1, var, fcn, result) \
    Kokkos::parallel_reduce( \
                        Kokkos::RangePolicy<> ( (x0), (x1)+1 ),  \
                        KOKKOS_LAMBDA(const int (i), decltype(var) &(var)){fcn}, \
                        Kokkos::Max< decltype(result) > ( (result) ) )

#define \
    DO_RMAX2D(i, x0, x1, j, y0, y1, var, fcn, result) \
    Kokkos::parallel_reduce( \
                        Kokkos::MDRangePolicy< Kokkos::Rank<2,F_LOOP_ORDER,F_LOOP_ORDER> > ( {(x0), (y0)}, {(x1)+1, (y1)+1} ), \
                        KOKKOS_LAMBDA( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
                        Kokkos::Max< decltype(result) > ( (result) ) )

#define \
    DO_RMAX3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
    Kokkos::parallel_reduce( \
                        Kokkos::MDRangePolicy< Kokkos::Rank<3,F_LOOP_ORDER,F_LOOP_ORDER> > ( {(x0), (y0), (z0)}, {(x1)+1, (y1)+1, (z1)+1} ), \
                        KOKKOS_LAMBDA( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
                        Kokkos::Max< decltype(result) > ( (result) ) )

#define \
    DO_REDUCE_MAX(...) \
    GET_MACRO(__VA_ARGS__, _13, DO_RMAX3D, _11, _10, DO_RMAX2D, _8, _7, DO_RMAX1D)(__VA_ARGS__)



// the REDUCE MIN loop
#define \
    RMIN1D(i, x0, x1, var, fcn, result) \
    Kokkos::parallel_reduce( \
                        Kokkos::RangePolicy<> ( (x0), (x1) ),  \
                        KOKKOS_LAMBDA( const int (i), decltype(var) &(var) ){fcn}, \
                        Kokkos::Min< decltype(result) >(result))

#define \
    RMIN2D(i, x0, x1, j, y0, y1, var, fcn, result) \
    Kokkos::parallel_reduce( \
                        Kokkos::MDRangePolicy< Kokkos::Rank<2,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0)}, {(x1), (y1)} ), \
                        KOKKOS_LAMBDA( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
                        Kokkos::Min< decltype(result) >(result) )

#define \
    RMIN3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
    Kokkos::parallel_reduce( \
                        Kokkos::MDRangePolicy< Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0), (z0)}, {(x1), (y1), (z1)} ), \
                        KOKKOS_LAMBDA( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
                        Kokkos::Min< decltype(result) >(result) )

#define \
    REDUCE_MIN(...) \
    GET_MACRO(__VA_ARGS__, _13, RMIN3D, _11, _10, RMIN2D, _8, _7, RMIN1D)(__VA_ARGS__)


// the DO_REDUCE MIN loop
#define \
    DO_RMIN1D(i, x0, x1, var, fcn, result) \
    Kokkos::parallel_reduce( \
                        Kokkos::RangePolicy<> ( (x0), (x1)+1 ),  \
                        KOKKOS_LAMBDA( const int (i), decltype(var) &(var) ){fcn}, \
                        Kokkos::Min< decltype(result) >(result))

#define \
    DO_RMIN2D(i, x0, x1, j, y0, y1, var, fcn, result) \
    Kokkos::parallel_reduce( \
                        Kokkos::MDRangePolicy< Kokkos::Rank<2,F_LOOP_ORDER,F_LOOP_ORDER> > ( {(x0), (y0)}, {(x1)+1, (y1)+1} ), \
                        KOKKOS_LAMBDA( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
                        Kokkos::Min< decltype(result) >(result) )

#define \
    DO_RMIN3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
    Kokkos::parallel_reduce( \
                        Kokkos::MDRangePolicy< Kokkos::Rank<3,F_LOOP_ORDER,F_LOOP_ORDER> > ( {(x0), (y0), (z0)}, {(x1)+1, (y1)+1, (z1)+1} ), \
                        KOKKOS_LAMBDA( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
                        Kokkos::Min< decltype(result) >(result) )

#define \
    DO_REDUCE_MIN(...) \
    GET_MACRO(__VA_ARGS__, _13, DO_RMIN3D, _11, _10, DO_RMIN2D, _8, _7, DO_RMIN1D)(__VA_ARGS__)



// the FOR_ALL loop with variables in a class
#define \
FORCLASS1D(i, x0, x1,fcn) \
Kokkos::parallel_for( Kokkos::RangePolicy<> ( (x0), (x1)), \
                     KOKKOS_CLASS_LAMBDA( const int (i) ){fcn} )

#define \
FORCLASS2D(i, x0, x1, j, y0, y1,fcn) \
Kokkos::parallel_for( \
                     Kokkos::MDRangePolicy< Kokkos::Rank<2,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0)}, {(x1), (y1)} ), \
                     KOKKOS_CLASS_LAMBDA( const int (i), const int (j) ){fcn} )

#define \
FORCLASS3D(i, x0, x1, j, y0, y1, k, z0, z1, fcn) \
Kokkos::parallel_for( \
                     Kokkos::MDRangePolicy< Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0), (z0)}, {(x1), (y1), (z1)} ), \
                     KOKKOS_CLASS_LAMBDA( const int (i), const int (j), const int (k) ) {fcn} )

#define \
FOR_ALL_CLASS(...) \
GET_MACRO(__VA_ARGS__, _13, _12, _11, FORCLASS3D, _9, _8, FORCLASS2D, _6, _5, FORCLASS1D)(__VA_ARGS__)


// the REDUCE SUM loop
#define \
RSUMCLASS1D(i, x0, x1, var, fcn, result) \
Kokkos::parallel_reduce( Kokkos::RangePolicy<> ( (x0), (x1) ),  \
                        KOKKOS_CLASS_LAMBDA(const int (i), decltype(var) &(var)){fcn}, (result))

#define \
RSUMCLASS2D(i, x0, x1, j, y0, y1, var, fcn, result) \
Kokkos::parallel_reduce( \
                        Kokkos::MDRangePolicy< Kokkos::Rank<2,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0)}, {(x1), (y1)} ), \
                        KOKKOS_CLASS_LAMBDA( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
                        (result) )

#define \
RSUMCLASS3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
Kokkos::parallel_reduce( \
                        Kokkos::MDRangePolicy< Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0), (z0)}, {(x1), (y1), (z1)} ), \
                        KOKKOS_CLASS_LAMBDA( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
                        (result) )

#define \
REDUCE_SUM_CLASS(...) \
GET_MACRO(__VA_ARGS__, _13, RSUMCLASS3D, _11, _10, RSUMCLASS2D, _8, _7, RSUMCLASS1D)(__VA_ARGS__)



// the REDUCE MAX loop with variables in a class

#define \
RMAXCLASS1D(i, x0, x1, var, fcn, result) \
Kokkos::parallel_reduce( \
                        Kokkos::RangePolicy<> ( (x0), (x1) ),  \
                        KOKKOS_CLASS_LAMBDA(const int (i), decltype(var) &(var)){fcn}, \
                        Kokkos::Max< decltype(result) > ( (result) ) )

#define \
RMAXCLASS2D(i, x0, x1, j, y0, y1, var, fcn, result) \
Kokkos::parallel_reduce( \
                        Kokkos::MDRangePolicy< Kokkos::Rank<2,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0)}, {(x1), (y1)} ), \
                        KOKKOS_CLASS_LAMBDA( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
                        Kokkos::Max< decltype(result) > ( (result) ) )

#define \
RMAXCLASS3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
Kokkos::parallel_reduce( \
                        Kokkos::MDRangePolicy< Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0), (z0)}, {(x1), (y1), (z1)} ), \
                        KOKKOS_CLASS_LAMBDA( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
                        Kokkos::Max< decltype(result) > ( (result) ) )

#define \
REDUCE_MAX_CLASS(...) \
GET_MACRO(__VA_ARGS__, _13, RMAXCLASS3D, _11, _10, RMAXCLASS2D, _8, _7, RMAXCLASS1D)(__VA_ARGS__)


// the REDUCE MIN loop with variables in a class
#define \
RMINCLASS1D(i, x0, x1, var, fcn, result) \
Kokkos::parallel_reduce( \
                        Kokkos::RangePolicy<> ( (x0), (x1) ),  \
                        KOKKOS_CLASS_LAMBDA( const int (i), decltype(var) &(var) ){fcn}, \
                        Kokkos::Min< decltype(result) >(result))

#define \
RMINCLASS2D(i, x0, x1, j, y0, y1, var, fcn, result) \
Kokkos::parallel_reduce( \
                        Kokkos::MDRangePolicy< Kokkos::Rank<2,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0)}, {(x1), (y1)} ), \
                        KOKKOS_CLASS_LAMBDA( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
                        Kokkos::Min< decltype(result) >(result) )

#define \
RMINCLASS3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
Kokkos::parallel_reduce( \
                        Kokkos::MDRangePolicy< Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER> > ( {(x0), (y0), (z0)}, {(x1), (y1), (z1)} ), \
                        KOKKOS_CLASS_LAMBDA( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
                        Kokkos::Min< decltype(result) >(result) )

#define \
REDUCE_MIN_CLASS(...) \
GET_MACRO(__VA_ARGS__, _13, RMINCLASS3D, _11, _10, RMINCLASS2D, _8, _7, RMINCLASS1D)(__VA_ARGS__)

#endif


// end of KOKKOS routines




// -----------------------------------------
// The for_all is used for serial loops and
// with the non-kokkos MACROS
// -----------------------------------------

template <typename F>
void for_all (int i_start, int i_end,
              const F &lambda_fcn){
    
    for (int i=i_start; i<i_end; i++){
        lambda_fcn(i);
    }
    
}; // end for_all


template <typename F>
void for_all (int i_start, int i_end,
              int j_start, int j_end,
              const F &lambda_fcn){
    
    for (int i=i_start; i<i_end; i++){
        for (int j=j_start; j<j_end; j++){
            lambda_fcn(i,j);
        }
    }
    
}; // end for_all


template <typename F>
void for_all (int i_start, int i_end,
              int j_start, int j_end,
              int k_start, int k_end,
              const F &lambda_fcn){
    
    for (int i=i_start; i<i_end; i++){
        for (int j=j_start; j<j_end; j++){
            for (int k=k_start; k<k_end; k++){
                lambda_fcn(i,j,k);
            }
        }
    }
    
}; // end for_all


template <typename F>
void for_all_delta (int i_start, int i_end, int i_delta,
                    const F &lambda_fcn){
    
    for (int i=i_start; i<i_end; i+=i_delta){
        lambda_fcn(i);
    }
    
}; // end for_all

template <typename F>
void for_all_delta (int i_start, int i_end, int i_delta,
                    int j_start, int j_end, int j_delta,
                    const F &lambda_fcn){
    
    for (int i=i_start; i<i_end; i+=i_delta){
        for (int j=j_start; j<j_end; j+=j_delta){
            lambda_fcn(i,j);
        }
    }
    
}; // end for_all


template <typename F>
void for_all_delta (int i_start, int i_end, int i_delta,
                    int j_start, int j_end, int j_delta,
                    int k_start, int k_end, int k_delta,
                    const F &lambda_fcn){
    
    for (int i=i_start; i<i_end; i+=i_delta){
        for (int j=j_start; j<j_end; j+=j_delta){
            for (int k=k_start; k<k_end; k+=k_delta){
                lambda_fcn(i,j,k);
            }
        }
    }
    
}; // end for_all



// the FOR_LOOP
// 1D FOR loop has 4 inputs
#define \
    FOR1DLOOP(i, x0, x1, fcn) \
    for_all( (x0), (x1), \
             [&]( const int (i) ){fcn} )

// 1D FOR loop with increment has 5 inputs
#define \
    FOR1DLOOPDELTA(i, x0, x1, i_delta, fcn) \
    for_all_delta( (x0), (x1), (i_delta), \
             [&]( const int (i) ){fcn} )

// 2D FOR loop has 7 inputs
#define \
    FOR2DLOOP(i, x0, x1, j, y0, y1, fcn)  \
    for_all( (x0), (x1), (y0), (y1), \
             [&]( const int (i), const int (j) ){fcn} )

// 2D FOR loop with increments has 9 inputs
#define \
    FOR2DLOOPDELTA(i, x0, x1, i_delta, j, y0, y1, j_delta, fcn)  \
    for_all_delta( (x0), (x1), (i_delta), (y0), (y1), (j_delta), \
                [&]( const int (i), const int (j) ){fcn} )

// 3D FOR loop has 10 inputs
#define \
    FOR3DLOOP(i, x0, x1, j, y0, y1, k, z0, z1, fcn) \
    for_all( (x0), (x1), (y0), (y1), (z0), (z1), \
             [&]( const int (i), const int (j), const int (k) ) {fcn} )

// 3D FOR loop with increments has 13 inputs
#define \
    FOR3DLOOPDELTA(i, x0, x1, i_delta, j, y0, y1, j_delta, k, z0, z1, k_delta, fcn) \
    for_all_delta( (x0), (x1), (i_delta), (y0), (y1), (j_delta), (z0), (z1), (k_delta), \
             [&]( const int (i), const int (j), const int (k) ) {fcn} )

#define \
    FOR_LOOP(...) \
    GET_MACRO(__VA_ARGS__, FOR3DLOOPDELTA, _12, _11, FOR3DLOOP, FOR2DLOOPDELTA, _8, FOR2DLOOP, _6, FOR1DLOOPDELTA, FOR1DLOOP)(__VA_ARGS__)


// the DO_ALL loop
// 1D DOloop has 4 inputs
#define \
    DO1DLOOP(i, x0, x1, fcn) \
    for_all( (x0), (x1)+1, \
             [&]( const int (i) ){fcn} )
// 1D FOR loop with increment has 5 inputs
#define \
    DO1DLOOPDELTA(i, x0, x1, i_delta, fcn) \
    for_all_delta( (x0), (x1)+1, (i_delta), \
             [&]( const int (i) ){fcn} )
// 2D DO loop has 7 inputs
#define \
    DO2DLOOP(i, x0, x1, j, y0, y1, fcn)  \
    for_all( (x0), (x1)+1, (y0), (y1)+1, \
             [&]( const int (i), const int (j) ){fcn} )
// 2D FOR loop with increments has 9 inputs
#define \
    DO2DLOOPDELTA(i, x0, x1, i_delta, j, y0, y1, j_delta, fcn)  \
    for_all_delta( (x0), (x1)+1, (i_delta), (y0), (y1)+1, (j_delta), \
                [&]( const int (i), const int (j) ){fcn} )
// 3D DO loop has 10 inputs
#define \
    DO3DLOOP(i, x0, x1, j, y0, y1, k, z0, z1, fcn) \
    for_all( (x0), (x1)+1, (y0), (y1)+1, (z0), (z1)+1, \
             [&]( const int (i), const int (j), const int (k) ) {fcn} )
// 3D FOR loop with increments has 13 inputs
#define \
    DO3DLOOPDELTA(i, x0, x1, i_delta, j, y0, y1, j_delta, k, z0, z1, k_delta, fcn) \
    for_all_delta( (x0), (x1)+1, (i_delta), (y0), (y1)+1, (j_delta), (z0), (z1)+1, (k_delta), \
             [&]( const int (i), const int (j), const int (k) ) {fcn} )
#define \
    DO_LOOP(...) \
    GET_MACRO(__VA_ARGS__, DO3DLOOPDELTA, _12, _11, DO3DLOOP, DO2DLOOPDELTA, _8, DO2DLOOP, _6, DO1DLOOPDELTA, DO1DLOOP)(__VA_ARGS__)




// -----------------------------------------
// The for_all and for_reduce functions that
// are used with the non-kokkos MACROS
// -----------------------------------------

#ifndef HAVE_KOKKOS
#include <limits>  // for the max and min values of a int, double, etc.

// SUM
template <typename T, typename F>
void reduce_sum (int i_start, int i_end,
                 T var,
                 const F &lambda_fcn, T &result){
    var = 0;
    for (int i=i_start; i<i_end; i++){
        lambda_fcn(i, var);
    }
    result = var;
};  // end for_reduce


template <typename T, typename F>
void reduce_sum (int i_start, int i_end,
                 int j_start, int j_end,
                 T var,
                 const F &lambda_fcn, T &result){
    var = 0;
    for (int i=i_start; i<i_end; i++){
        for (int j=j_start; j<j_end; j++){
            lambda_fcn(i,j,var);
        }
    }
    
    result = var;
};  // end for_reduce


template <typename T, typename F>
void reduce_sum (int i_start, int i_end,
                 int j_start, int j_end,
                 int k_start, int k_end,
                 T  var,
                 const F &lambda_fcn,  T &result){
    var = 0;
    for (int i=i_start; i<i_end; i++){
        for (int j=j_start; j<j_end; j++){
            for (int k=k_start; k<k_end; k++){
                lambda_fcn(i,j,k,var);
            }
        }
    }
    
    result = var;
};  // end for_reduce


// MIN
template <typename T, typename F>
void reduce_min (int i_start, int i_end,
                 T var,
                 const F &lambda_fcn, T &result){
    var = std::numeric_limits<T>::max(); //2147483647;
    for (int i=i_start; i<i_end; i++){
        lambda_fcn(i, var);
    }
    result = var;
};  // end for_reduce


template <typename T, typename F>
void reduce_min (int i_start, int i_end,
                 int j_start, int j_end,
                 T var,
                 const F &lambda_fcn, T &result){
    var = std::numeric_limits<T>::max(); //2147483647;
    for (int i=i_start; i<i_end; i++){
        for (int j=j_start; j<j_end; j++){
            lambda_fcn(i,j,var);
        }
    }
    
    result = var;
};  // end for_reduce


template <typename T, typename F>
void reduce_min (int i_start, int i_end,
                 int j_start, int j_end,
                 int k_start, int k_end,
                 T  var,
                 const F &lambda_fcn,  T &result){
    var = std::numeric_limits<T>::max(); //2147483647;
    for (int i=i_start; i<i_end; i++){
        for (int j=j_start; j<j_end; j++){
            for (int k=k_start; k<k_end; k++){
                lambda_fcn(i,j,k,var);
            }
        }
    }
    
    result = var;
};  // end for_reduce

// MAX
template <typename T, typename F>
void reduce_max (int i_start, int i_end,
                 T var,
                 const F &lambda_fcn, T &result){
    var = std::numeric_limits<T>::min(); // -2147483647 - 1;
    for (int i=i_start; i<i_end; i++){
        lambda_fcn(i, var);
    }
    result = var;
};  // end for_reduce


template <typename T, typename F>
void reduce_max (int i_start, int i_end,
                 int j_start, int j_end,
                 T var,
                 const F &lambda_fcn, T &result){
    var = std::numeric_limits<T>::min(); //-2147483647 - 1;
    for (int i=i_start; i<i_end; i++){
        for (int j=j_start; j<j_end; j++){
            lambda_fcn(i,j,var);
        }
    }
    
    result = var;
};  // end for_reduce


template <typename T, typename F>
void reduce_max (int i_start, int i_end,
                 int j_start, int j_end,
                 int k_start, int k_end,
                 T  var,
                 const F &lambda_fcn,  T &result){
    var = std::numeric_limits<T>::min(); // -2147483647 - 1;
    for (int i=i_start; i<i_end; i++){
        for (int j=j_start; j<j_end; j++){
            for (int k=k_start; k<k_end; k++){
                lambda_fcn(i,j,k,var);
            }
        }
    }
    
    result = var;
};  // end for_reduce

#endif  // if not kokkos


// -----------------------------------------
// MACROS for none kokkos loops
// -----------------------------------------

#ifndef HAVE_KOKKOS

// replace the CLASS loops to be the nominal loops
#define FOR_ALL_CLASS FOR_ALL
#define REDUCE_SUM_CLASS REDUCE_SUM
#define REDUCE_MAX_CLASS REDUCE_MAX
#define REDUCE_MIN_CLASS REDUCE_MIN

// the FOR_ALL loop is chosen based on the number of inputs

// the FOR_ALL loop
// 1D FOR loop has 4 inputs
#define \
    FOR1D(i, x0, x1, fcn) \
    for_all( (x0), (x1), \
             [&]( const int (i) ){fcn} )
// 2D FOR loop has 7 inputs
#define \
    FOR2D(i, x0, x1, j, y0, y1, fcn)  \
    for_all( (x0), (x1), (y0), (y1), \
             [&]( const int (i), const int (j) ){fcn} )
// 3D FOR loop has 10 inputs
#define \
    FOR3D(i, x0, x1, j, y0, y1, k, z0, z1, fcn) \
    for_all( (x0), (x1), (y0), (y1), (z0), (z1), \
             [&]( const int (i), const int (j), const int (k) ) {fcn} )
#define \
    FOR_ALL(...) \
    GET_MACRO(__VA_ARGS__, _13, _12, _11, FOR3D, _9, _8, FOR2D, _6, _5, FOR1D)(__VA_ARGS__)


// the DO_ALL loop
// 1D DOloop has 4 inputs
#define \
    DO1D(i, x0, x1, fcn) \
    for_all( (x0), (x1)+1, \
             [&]( const int (i) ){fcn} )
// 2D DO loop has 7 inputs
#define \
    DO2D(i, x0, x1, j, y0, y1, fcn)  \
    for_all( (x0), (x1)+1, (y0), (y1)+1, \
             [&]( const int (i), const int (j) ){fcn} )
// 3D DO loop has 10 inputs
#define \
    DO3D(i, x0, x1, j, y0, y1, k, z0, z1, fcn) \
    for_all( (x0), (x1)+1, (y0), (y1)+1, (z0), (z1)+1, \
             [&]( const int (i), const int (j), const int (k) ) {fcn} )
#define \
    DO_ALL(...) \
    GET_MACRO(__VA_ARGS__, _13, _12, _11, DO3D, _9, _8, DO2D, _6, _5, DO1D)(__VA_ARGS__)


// the REDUCE loops, no kokkos
#define \
    RSUM1D(i, x0, x1, var, fcn, result) \
    reduce_sum( (x0), (x1), (var),  \
                [=]( const int (i), decltype(var) &(var) ){fcn}, \
                (result) )
#define \
    RSUM2D(i, x0, x1, j, y0, y1, var, fcn, result) \
    reduce_sum( (x0), (x1), (y0), (y1), (var),  \
                [=]( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
                (result) )
#define \
    RSUM3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
    reduce_sum( (x0), (x1), (y0), (y1), (z0), (z1), (var),  \
                [=]( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
                (result) )

#define \
    REDUCE_SUM(...) \
    GET_MACRO(__VA_ARGS__, _13, RSUM3D, _11, _10, RSUM2D, _8, _7, RSUM1D)(__VA_ARGS__)


// DO_REDUCE_SUM
#define \
    DO_RSUM1D(i, x0, x1, var, fcn, result) \
    reduce_sum( (x0), (x1)+1, (var),  \
                [=]( const int (i), decltype(var) &(var) ){fcn}, \
                (result) )
#define \
    DO_RSUM2D(i, x0, x1, j, y0, y1, var, fcn, result) \
    reduce_sum( (x0), (x1)+1, (y0), (y1)+1, (var),  \
                [=]( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
                (result) )
#define \
    DO_RSUM3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
    reduce_sum( (x0), (x1)+1, (y0), (y1)+1, (z0), (z1)+1, (var),  \
                [=]( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
                (result) )

#define \
    DO_REDUCE_SUM(...) \
    GET_MACRO(__VA_ARGS__, _13, DO_RSUM3D, _11, _10, DO_RSUM2D, _8, _7, DO_RSUM1D)(__VA_ARGS__)


// Reduce max
#define \
    RMAX1D(i, x0, x1, var, fcn, result) \
    reduce_max( (x0), (x1), (var),  \
                [=]( const int (i), decltype(var) &(var) ){fcn}, \
                (result) )
#define \
    RMAX2D(i, x0, x1, j, y0, y1, var, fcn, result) \
    reduce_max( (x0), (x1), (y0), (y1), (var),  \
                [=]( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
                (result) )
#define \
    RMAX3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
    reduce_max( (x0), (x1), (y0), (y1), (z0), (z1), (var),  \
                [=]( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
                (result) )

#define \
    REDUCE_MAX(...) \
    GET_MACRO(__VA_ARGS__, _13, RMAX3D, _11, _10, RMAX2D, _8, _7, RMAX1D)(__VA_ARGS__)


// DO_REDUCE_MAX
#define \
    DO_RMAX1D(i, x0, x1, var, fcn, result) \
    reduce_max( (x0), (x1)+1, (var),  \
                [=]( const int (i), decltype(var) &(var) ){fcn}, \
                (result) )
#define \
    DO_RMAX2D(i, x0, x1, j, y0, y1, var, fcn, result) \
    reduce_max( (x0), (x1)+1, (y0), (y1)+1, (var),  \
                [=]( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
                (result) )
#define \
    DO_RMAX3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
    reduce_max( (x0), (x1)+1, (y0), (y1)+1, (z0), (z1)+1, (var),  \
                [=]( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
                (result) )

#define \
    DO_REDUCE_MAX(...) \
    GET_MACRO(__VA_ARGS__, _13, DO_RMAX3D, _11, _10, DO_RMAX2D, _8, _7, DO_RMAX1D)(__VA_ARGS__)


// reduce min
#define \
    RMIN1D(i, x0, x1, var, fcn, result) \
    reduce_min( (x0), (x1), (var),  \
                [=]( const int (i), decltype(var) &(var) ){fcn}, \
                (result) )
#define \
    RMIN2D(i, x0, x1, j, y0, y1, var, fcn, result) \
    reduce_min( (x0), (x1), (y0), (y1), (var),  \
                [=]( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
                (result) )
#define \
    RMIN3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
    reduce_min( (x0), (x1), (y0), (y1), (z0), (z1), (var),  \
                [=]( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
                (result) )

#define \
    REDUCE_MIN(...) \
    GET_MACRO(__VA_ARGS__, _13, RMIN3D, _11, _10, RMIN2D, _8, _7, RMIN1D)(__VA_ARGS__)


// DO_REDUCE_MIN
#define \
    DO_RMIN1D(i, x0, x1, var, fcn, result) \
    reduce_min( (x0), (x1)+1, (var),  \
                [=]( const int (i), decltype(var) &(var) ){fcn}, \
                (result) )
#define \
    DO_RMIN2D(i, x0, x1, j, y0, y1, var, fcn, result) \
    reduce_min( (x0), (x1)+1, (y0), (y1)+1, (var),  \
                [=]( const int (i),const int (j), decltype(var) &(var) ){fcn}, \
                (result) )
#define \
    DO_RMIN3D(i, x0, x1, j, y0, y1, k, z0, z1, var, fcn, result) \
    reduce_min( (x0), (x1)+1, (y0), (y1)+1, (z0), (z1)+1, (var),  \
                [=]( const int (i), const int (j), const int (k), decltype(var) &(var) ){fcn}, \
                (result) )

#define \
    DO_REDUCE_MIN(...) \
    GET_MACRO(__VA_ARGS__, _13, DO_RMIN3D, _11, _10, DO_RMIN2D, _8, _7, DO_RMIN1D)(__VA_ARGS__)


#endif  // if not kokkos



#endif // MACROS_H

