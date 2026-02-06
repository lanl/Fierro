#pragma once

#include <string>
#include "matar.h"

namespace utils {
// SINGLE determines whether single or double precision is built
#ifdef SINGLE
  #define MPI_REAL_T MPI_FLOAT
  typedef float real_t;  // define native type for EVPFFT as single precision
  #define FMT1 "%24.14f" //"%g"      // format argument for floats    
  #define EMT1 "%24.14E" // format argument for eng floats
#else
  #define MPI_REAL_T MPI_DOUBLE
  typedef double real_t; // define native type for EVPFFT as double precision
  #define FMT1 "%24.14f" //"%lg"     // format argument for doubles   
  #define EMT1 "%24.14E" // format argument for eng doubles
#endif
} // end of namespace utils scope
//
using namespace utils;
using namespace mtr;

// math functions
#define POW2(x) ( (x)*(x) )
#define POW3(x) ( (x)*(x)*(x) )

#ifdef SINGLE
  #define SIN(x)          sinf(x)
  #define COS(x)          cosf(x)
  #define TAN(x)          tanf(x)

  #define ASIN(x)         asinf(x)
  #define ACOS(x)         acosf(x)
  #define ATAN(x)         atanf(x)
  #define ATAN2(x,y)      atan2f(x,y)

  #define ABS(x)          fabsf(x)
  #define SQRT(x)         sqrtf(x)
  #define EXP(x)          expf(x)
  #define LOG(x)          logf(x)
  #define POW(x,y)        powf(x,y)
  #define COPYSIGN(x,y)   copysignf(x,y)

#else
  #define SIN(x)          sin(x)
  #define COS(x)          cos(x)
  #define TAN(x)          tan(x)

  #define ASIN(x)         asin(x)
  #define ACOS(x)         acos(x)
  #define ATAN(x)         atan(x)
  #define ATAN2(x,y)      atan2(x,y)

  #define ABS(x)          fabs(x)
  #define SQRT(x)         sqrt(x)
  #define EXP(x)          exp(x)
  #define LOG(x)          log(x)
  #define POW(x,y)        pow(x,y)
  #define COPYSIGN(x,y)   copysign(x,y)
#endif


// aliases
using MatrixTypeIntDevice   = FMatrixKokkos <int>;
using MatrixTypeRealDevice  = FMatrixKokkos <real_t>;
using MatrixTypeRealDual    = DFMatrixKokkos <real_t>;
using MatrixTypeIntDual     = DFMatrixKokkos <int>;
using MatrixTypeRealHost    = FMatrix <real_t>;
using MatrixTypeIntHost     = FMatrix <int>;
using MatrixTypeStringHost  = FMatrix <std::string>;
// Views dont need Device or Host designations
// Views can be created on both device or host.
using ViewMatrixTypeInt     = ViewFMatrixKokkos <int>;
using ViewMatrixTypeReal    = ViewFMatrixKokkos <real_t>;

// CArray nested loop convention use Right, FArray use Left
#define LOOP_ORDER Kokkos::Iterate::Right

// Definitions for printing and file manipulation
#define screenOut stdout
#define CLEAR_LINE(ifstream) ( ifstream.ignore(std::numeric_limits<std::streamsize>::max(), '\n') );
#define FULL_INPUT_PATH(filename) ( _INPUT_DIR filename )
#define FULL_OUTPUT_PATH(filename) ( _OUTPUT_DIR filename )

// *********************** For debuging ***********************
#define PRINT_LINE Kokkos::fence(); printf("HERE LINE %d IN FILE %s\n", __LINE__, __FILE__); fflush(stdout);

#define PRINT_ARRAY_DIMS(x) \
        for (int i = 0; i < x.order(); i++) { \
            printf("size of dim %d = %d\n", i, x.dims(i)); \
            fflush(stdout); \
        }

#define PRINT_MATRIX_DIMS(x) \
        for (int i = 1; i <= x.order(); i++) { \
            printf("size of dim %d = %d\n", i, x.dims(i)); \
            fflush(stdout); \
        }
        
#define PAUSE Kokkos::fence(); std::cout << "\nC++ PAUSE: enter <return> or <ctrl>d to continue>\n"; std::cin.get();

#define PRINT_ARRAY(arr) Kokkos::fence(); print_array(arr.pointer(), arr.size());
        
#define PRINT_DARRAY(arr) Kokkos::fence(); print_array(arr.host_pointer(), arr.size());
