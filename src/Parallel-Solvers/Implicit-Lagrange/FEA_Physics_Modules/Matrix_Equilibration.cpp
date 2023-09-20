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

//Derived from examples in the Trilinos Ifpack2 package

// Teuchos
#include <Teuchos_ScalarTraits.hpp>

// Xpetra
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <KokkosBlas1_abs.hpp>
#include <Tpetra_leftAndOrRightScaleCrsMatrix.hpp>
#include <Tpetra_computeRowAndColumnOneNorms.hpp>
#include "KokkosBlas1_abs_impl.hpp"
#include <MueLu_Utilities.hpp>
#include "Implicit_Solver.h"

using namespace utils;

/* ----------------------------------------------------------------------
   Equilibrate Matrix
------------------------------------------------------------------------- */
void Implicit_Solver::equilibrateMatrix(Teuchos::RCP<Xpetra::Matrix<real_t,LO,GO,node_type> > &Axpetra, std::string equilibrate) {
  using Tpetra::computeRowAndColumnOneNorms;
  using Tpetra::leftAndOrRightScaleCrsMatrix;
  bool equilibrate_1norm = (equilibrate == "1-norm");
  bool equilibrate_diag  = (equilibrate == "diag");
  bool assumeSymmetric = false;

  Teuchos::RCP<Tpetra::CrsMatrix<real_t,LO,GO,node_type> > A = MueLu::Utilities<real_t,LO,GO,node_type>::Op2NonConstTpetraCrs(Axpetra);

     equibResult = computeRowAndColumnOneNorms (*A, assumeSymmetric);
      if(equilibrate_diag) {

        scaling_view_type rowDiagAbsVals ("rowDiagAbsVals",
                                  equibResult.rowDiagonalEntries.extent (0));
        KokkosBlas::abs (rowDiagAbsVals, equibResult.rowDiagonalEntries);
        scaling_view_type colDiagAbsVals ("colDiagAbsVals",
                                  equibResult.colDiagonalEntries.extent (0));
        KokkosBlas::abs (colDiagAbsVals, equibResult.colDiagonalEntries);

        leftAndOrRightScaleCrsMatrix (*A, rowDiagAbsVals, colDiagAbsVals,
                                      true, true, equibResult.assumeSymmetric,
                                      Tpetra::SCALING_DIVIDE);
      }
      else{
        auto colScalingFactors = equibResult.assumeSymmetric ?
          equibResult.colNorms :
          equibResult.rowScaledColNorms;
        leftAndOrRightScaleCrsMatrix (*A, equibResult.rowNorms,
                                      colScalingFactors, true, true,
                                      equibResult.assumeSymmetric,
                                      Tpetra::SCALING_DIVIDE);
      }
}

/* ----------------------------------------------------------------------
   Scale RHS vector by equilibration row scaling
------------------------------------------------------------------------- */

void Implicit_Solver::preScaleRightHandSides (MV& B, std::string equilibrate)
  {
    
    bool equilibrate_1norm = (equilibrate == "1-norm");
    bool equilibrate_diag  = (equilibrate == "diag");
      if(equilibrate_diag) {
        bool takeSquareRootsOfScalingFactors = false; // just use the diagonal entries
        elementWiseDivideMultiVector (B, equibResult.rowDiagonalEntries,
                                      takeSquareRootsOfScalingFactors);
      }
      else{
        bool takeSquareRootsOfScalingFactors = equibResult.assumeSymmetric;
        elementWiseDivideMultiVector (B, equibResult.rowNorms,
                                      takeSquareRootsOfScalingFactors);
      }
      
    
  }

  /* ----------------------------------------------------------------------
   Scale Solution vector initial guess by equilibration column scaling
------------------------------------------------------------------------- */

  void Implicit_Solver::preScaleInitialGuesses (MV& X, std::string equilibrate)
  {
    bool equilibrate_1norm = (equilibrate == "1-norm");
    bool equilibrate_diag  = (equilibrate == "diag");
      if(equilibrate_diag) {
        bool takeSquareRootsOfScalingFactors = false; // just use the diagonal entries
        elementWiseMultiplyMultiVector (X, equibResult.colDiagonalEntries,
                                        takeSquareRootsOfScalingFactors);
      }
      else {
        auto colScalingFactors = equibResult.assumeSymmetric ?
          equibResult.colNorms :
          equibResult.rowScaledColNorms;
        bool takeSquareRootsOfScalingFactors =
          equibResult.assumeSymmetric;
        elementWiseMultiplyMultiVector (X, colScalingFactors,
                                        takeSquareRootsOfScalingFactors);
      }
  }

  /* ----------------------------------------------------------------------
   Scale Solution vector by equilibration column scaling
------------------------------------------------------------------------- */

  void Implicit_Solver::postScaleSolutionVectors (MV& X, std::string equilibrate)
  {
    bool equilibrate_1norm = (equilibrate == "1-norm");
    bool equilibrate_diag  = (equilibrate == "diag");
      if(equilibrate_diag) {
        bool takeSquareRootsOfScalingFactors = false; // just use the diagonal entries
        elementWiseDivideMultiVector (X, equibResult.colDiagonalEntries,
                                      takeSquareRootsOfScalingFactors);
      }
      else {
        auto colScalingFactors = equibResult.assumeSymmetric ?
          equibResult.colNorms :
          equibResult.rowScaledColNorms;
        bool takeSquareRootsOfScalingFactors =
          equibResult.assumeSymmetric;
        elementWiseDivideMultiVector (X, colScalingFactors,
                                      takeSquareRootsOfScalingFactors);
      }
    
  }

/* ----------------------------------------------------------------------
   Scale vector through multiplication
------------------------------------------------------------------------- */

void Implicit_Solver::elementWiseMultiply ( vec_array& X,
                      scaling_view_type& scalingFactors,
                      LO numRows,
                      bool takeSquareRootsOfScalingFactors)
{
  using execution_space = typename vec_array::device_type::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, LO>;
  bool takeAbsoluteValueOfScalingFactors = true;
  if (takeAbsoluteValueOfScalingFactors) {
    constexpr bool takeAbsVal = true;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseMultiply;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors,takeSquareRootsOfScalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseMultiply;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors,takeSquareRootsOfScalingFactors));
    }
  }
  else {
    constexpr bool takeAbsVal = false;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseMultiply;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors,takeSquareRootsOfScalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseMultiply;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors,takeSquareRootsOfScalingFactors));
    }
  }
}

void Implicit_Solver::elementWiseMultiplyMultiVector (MV& X,
                                scaling_view_type& scalingFactors,
                                bool takeSquareRootsOfScalingFactors)
{
  bool takeAbsoluteValueOfScalingFactors = true;
  using index_type = typename MV::local_ordinal_type;
  const index_type lclNumRows = static_cast<index_type> (X.getLocalLength ());

  auto X_lcl = X.getLocalViewDevice (Tpetra::Access::ReadWrite);
  /*
  if (static_cast<std::size_t> (X.getNumVectors ()) == std::size_t (1)) {
    using pair_type = Kokkos::pair<index_type, index_type>;
    auto X_lcl_1d = Kokkos::subview (X_lcl, pair_type (0, lclNumRows), 0);
    elementWiseMultiply (X_lcl_1d, scalingFactors, lclNumRows,
                         takeSquareRootsOfScalingFactors);
  }
  */
  //else {
    elementWiseMultiply (X_lcl, scalingFactors, lclNumRows,
                         takeSquareRootsOfScalingFactors);
  //}
}

/* ----------------------------------------------------------------------
   Scale vector through division
------------------------------------------------------------------------- */
void Implicit_Solver::elementWiseDivide (vec_array& X,
                   scaling_view_type& scalingFactors,
                  LO numRows,
                   bool takeSquareRootsOfScalingFactors)
{
  bool takeAbsoluteValueOfScalingFactors = true;
  using execution_space = typename vec_array::device_type::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, LO>;

  if (takeAbsoluteValueOfScalingFactors) {
    constexpr bool takeAbsVal = true;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseDivide;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors,takeSquareRootsOfScalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseDivide;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors,takeSquareRootsOfScalingFactors));
    }
  }
  else {
    constexpr bool takeAbsVal = false;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseDivide;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors,takeSquareRootsOfScalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseDivide;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors,takeSquareRootsOfScalingFactors));
    }
  }
}

void Implicit_Solver::elementWiseDivideMultiVector (MV& X,
                              scaling_view_type& scalingFactors,
                              bool takeSquareRootsOfScalingFactors)
{
  bool takeAbsoluteValueOfScalingFactors = true;
  using index_type = typename MV::local_ordinal_type;
  const index_type lclNumRows = static_cast<index_type> (X.getLocalLength ());

  auto X_lcl = X.getLocalViewDevice (Tpetra::Access::ReadWrite);
  /*
  if (static_cast<std::size_t> (X.getNumVectors ()) == std::size_t (1)) {
    using pair_type = Kokkos::pair<index_type, index_type>;
    auto X_lcl_1d = Kokkos::subview (X_lcl, pair_type (0, lclNumRows), 0);
    elementWiseDivide (X_lcl_1d, scalingFactors, lclNumRows,
                       takeSquareRootsOfScalingFactors);
  }
  */
  //else {
    elementWiseDivide (X_lcl, scalingFactors, lclNumRows,
                       takeSquareRootsOfScalingFactors);
  //}
}

