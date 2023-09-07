#ifndef KOKKOS_ALIAS_H
#define KOKKOS_ALIAS_H

#include <stdlib.h> 
#include "parents.h"
#include <Kokkos_Core.hpp>
#include "matar.h"

//MACROS to make the code less scary
//#define kmalloc(size) ( Kokkos::kokkos_malloc<DefaultMemSpace>(size) )
//#define kfree(pnt)        (  Kokkos::kokkos_free(pnt) ) 
//#define ProfileRegionStart  ( Kokkos::Profiling::pushRegion )
//#define ProfileRegionEnd  ( Kokkos::Profiling::popRegion )

using real_t = double;
using u_int  = unsigned int;

/*
#ifdef HAVE_CUDA
//using UVMMemSpace     = Kokkos::CudaUVMSpace;
using DefaultMemSpace  = Kokkos::CudaSpace;
using DefaultExecSpace = Kokkos::Cuda;
using DefaultLayout    = Kokkos::LayoutLeft;
#elif HAVE_OPENMP
using DefaultMemSpace  = Kokkos::HostSpace;
using DefaultExecSpace = Kokkos::OpenMP;
using DefaultLayout    = Kokkos::LayoutRight;
#elif TRILINOS_INTERFACE
using DefaultMemSpace  = void;
using DefaultExecSpace = void;
using DefaultLayout    = void;
#elif HAVE_HIP
using DefaultMemSpace  = Kokkos::HipSpace;
using DefaultExecSpace = Kokkos::Hip;
using DefaultLayout    = Kokkos::LayoutLeft;
#endif
*/

using TeamPolicy = Kokkos::TeamPolicy<DefaultExecSpace>;
using mdrange_policy2 = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
using mdrange_policy3 = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;

using RMatrix1D    = Kokkos::View<real_t *,DefaultLayout,DefaultExecSpace>;
using RMatrix2D    = Kokkos::View<real_t **,DefaultLayout,DefaultExecSpace>;
using RMatrix3D    = Kokkos::View<real_t ***,DefaultLayout,DefaultExecSpace>;
using RMatrix4D    = Kokkos::View<real_t ****,DefaultLayout,DefaultExecSpace>;
using RMatrix5D    = Kokkos::View<real_t *****,DefaultLayout,DefaultExecSpace>;
using IMatrix1D    = Kokkos::View<int *,DefaultLayout,DefaultExecSpace>;
using IMatrix2D    = Kokkos::View<int **,DefaultLayout,DefaultExecSpace>;
using IMatrix3D    = Kokkos::View<int ***,DefaultLayout,DefaultExecSpace>;
using IMatrix4D    = Kokkos::View<int ****,DefaultLayout,DefaultExecSpace>;
using IMatrix5D    = Kokkos::View<int *****,DefaultLayout,DefaultExecSpace>;
using SVar         = Kokkos::View<size_t,DefaultLayout,DefaultExecSpace>;
using SArray1D     = Kokkos::View<size_t *,DefaultLayout,DefaultExecSpace>;
using SArray2D     = Kokkos::View<size_t **,DefaultLayout,DefaultExecSpace>;
using SArray3D     = Kokkos::View<size_t ***,DefaultLayout,DefaultExecSpace>;
using SArray4D     = Kokkos::View<size_t ****,DefaultLayout,DefaultExecSpace>;
using SArray5D     = Kokkos::View<size_t *****,DefaultLayout,DefaultExecSpace>;

using SHArray1D     = Kokkos::View<size_t *,DefaultLayout,Kokkos::HostSpace>;

using Parent1D     = Kokkos::View<parent_models*,DefaultLayout,DefaultExecSpace>;
using ParentHost1D = Kokkos::View<parent_models*,DefaultLayout,Kokkos::HostSpace>;

#endif
