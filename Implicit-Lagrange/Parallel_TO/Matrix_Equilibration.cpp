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
#include "Parallel_Nonlinear_Solver.h"

using namespace utils;

/* ----------------------------------------------------------------------
   Equilibrate Matrix
------------------------------------------------------------------------- */
void Parallel_Nonlinear_Solver::equilibrateMatrix(Teuchos::RCP<Xpetra::Matrix<real_t,LO,GO,node_type> > &Axpetra, std::string equilibrate) {
  using Tpetra::computeRowAndColumnOneNorms;
  using Tpetra::leftAndOrRightScaleCrsMatrix;
  bool equilibrate_1norm = (equilibrate == "1-norm");
  bool equilibrate_diag  = (equilibrate == "diag");
  bool assumeSymmetric = false;

  Teuchos::RCP<Tpetra::CrsMatrix<real_t,LO,GO,node_type> > A = MueLu::Utilities<real_t,LO,GO,node_type>::Op2NonConstTpetraCrs(Axpetra);

     equibResult = computeRowAndColumnOneNorms (*A, assumeSymmetric);
      if(equilibrate_diag) {
        using device_type = typename node_type::device_type;
        using mag_type = typename Kokkos::ArithTraits<real_t>::mag_type;
        using view_type = Kokkos::View<mag_type*, device_type>;

        view_type rowDiagAbsVals ("rowDiagAbsVals",
                                  equibResult.rowDiagonalEntries.extent (0));
        KokkosBlas::abs (rowDiagAbsVals, equibResult.rowDiagonalEntries);
        view_type colDiagAbsVals ("colDiagAbsVals",
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

void Parallel_Nonlinear_Solver::preScaleRightHandSides (Tpetra::MultiVector<real_t,LO,GO,node_type>& B, std::string equilibrate) const
  {
    
    bool equilibrate_1norm = (equilibrate == "1-norm");
    bool equilibrate_diag  = (equilibrate == "diag");
      if(equilibrate_diag) {
        const bool takeSquareRootsOfScalingFactors = false; // just use the diagonal entries
        elementWiseDivideMultiVector (B, equibResult.rowDiagonalEntries,
                                      takeSquareRootsOfScalingFactors);
      }
      else{
        const bool takeSquareRootsOfScalingFactors = equibResult.assumeSymmetric;
        elementWiseDivideMultiVector (B, equibResult.rowNorms,
                                      takeSquareRootsOfScalingFactors);
      }
      
    
  }

  /* ----------------------------------------------------------------------
   Scale Solution vector initial guess by equilibration column scaling
------------------------------------------------------------------------- */

  void Parallel_Nonlinear_Solver::preScaleInitialGuesses (Tpetra::MultiVector<real_t,LO,GO,node_type>& X, std::string equilibrate) const
  {
    bool equilibrate_1norm = (equilibrate == "1-norm");
    bool equilibrate_diag  = (equilibrate == "diag");
      if(equilibrate_diag) {
        const bool takeSquareRootsOfScalingFactors = false; // just use the diagonal entries
        elementWiseMultiplyMultiVector (X, equibResult.colDiagonalEntries,
                                        takeSquareRootsOfScalingFactors);
      }
      else {
        auto colScalingFactors = equibResult.assumeSymmetric ?
          equibResult.colNorms :
          equibResult.rowScaledColNorms;
        const bool takeSquareRootsOfScalingFactors =
          equibResult.assumeSymmetric;
        elementWiseMultiplyMultiVector (X, colScalingFactors,
                                        takeSquareRootsOfScalingFactors);
      }
  }

  /* ----------------------------------------------------------------------
   Scale Solution vector by equilibration column scaling
------------------------------------------------------------------------- */

  void Parallel_Nonlinear_Solver::postScaleSolutionVectors (Tpetra::MultiVector<real_t,LO,GO,node_type>& X, std::string equilibrate) const
  {
    bool equilibrate_1norm = (equilibrate == "1-norm");
    bool equilibrate_diag  = (equilibrate == "diag");
      if(equilibrate_diag) {
        const bool takeSquareRootsOfScalingFactors = false; // just use the diagonal entries
        elementWiseDivideMultiVector (X, equibResult.colDiagonalEntries,
                                      takeSquareRootsOfScalingFactors);
      }
      else {
        auto colScalingFactors = equibResult.assumeSymmetric ?
          equibResult.colNorms :
          equibResult.rowScaledColNorms;
        const bool takeSquareRootsOfScalingFactors =
          equibResult.assumeSymmetric;
        elementWiseDivideMultiVector (X, colScalingFactors,
                                      takeSquareRootsOfScalingFactors);
      }
    
  }

/* ----------------------------------------------------------------------
   Scale vector through multiplication
------------------------------------------------------------------------- */

void
elementWiseMultiply (const MV& X,
                     const ScalingFactorsViewType& scalingFactors,
                     const IndexType numRows,
                     const bool takeSquareRootsOfScalingFactors,
                     const bool takeAbsoluteValueOfScalingFactors =
                       ! std::is_same<
                           typename Kokkos::ArithTraits<
                             typename MV::non_const_value_type
                           >::mag_type,
                           typename ScalingFactorsViewType::non_const_value_type
                         >::value)
{
  using execution_space = typename MV::device_type::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, IndexType>;

  if (takeAbsoluteValueOfScalingFactors) {
    constexpr bool takeAbsVal = true;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseMultiply<MV,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseMultiply<MV,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
  }
  else {
    constexpr bool takeAbsVal = false;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseMultiply<MV,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseMultiply<MV,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseMultiply",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
  }
}

void
elementWiseMultiplyMultiVector (MV& X,
                                const ScalingFactorsViewType& scalingFactors,
                                const bool takeSquareRootsOfScalingFactors,
                                const bool takeAbsoluteValueOfScalingFactors =
                                  ! std::is_same<
                                      typename Kokkos::ArithTraits<
                                        typename MV::scalar_type
                                      >::mag_type,
                                      typename ScalingFactorsViewType::non_const_value_type
                                    >::value)
{
  using index_type = typename MV::local_ordinal_type;
  const index_type lclNumRows = static_cast<index_type> (X.getLocalLength ());

  auto X_lcl = X.getLocalViewDevice (Tpetra::Access::ReadWrite);
  if (static_cast<std::size_t> (X.getNumVectors ()) == std::size_t (1)) {
    using pair_type = Kokkos::pair<index_type, index_type>;
    auto X_lcl_1d = Kokkos::subview (X_lcl, pair_type (0, lclNumRows), 0);
    elementWiseMultiply (X_lcl_1d, scalingFactors, lclNumRows,
                         takeSquareRootsOfScalingFactors,
                         takeAbsoluteValueOfScalingFactors);
  }
  else {
    elementWiseMultiply (X_lcl, scalingFactors, lclNumRows,
                         takeSquareRootsOfScalingFactors,
                         takeAbsoluteValueOfScalingFactors);
  }
}

/* ----------------------------------------------------------------------
   Scale vector through division
------------------------------------------------------------------------- */
void
elementWiseDivide (const MV& X,
                   const ScalingFactorsViewType& scalingFactors,
                   const IndexType numRows,
                   const bool takeSquareRootsOfScalingFactors,
                   const bool takeAbsoluteValueOfScalingFactors =
                     ! std::is_same<
                         typename Kokkos::ArithTraits<
                           typename MV::non_const_value_type
                         >::mag_type,
                         typename ScalingFactorsViewType::non_const_value_type
                       >::value)
{
  using execution_space = typename MV::device_type::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, IndexType>;

  if (takeAbsoluteValueOfScalingFactors) {
    constexpr bool takeAbsVal = true;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseDivide<MV,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseDivide<MV,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
  }
  else {
    constexpr bool takeAbsVal = false;
    if (takeSquareRootsOfScalingFactors) {
      constexpr bool takeSquareRoots = true;
      using functor_type = ElementWiseDivide<MV,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
    else {
      constexpr bool takeSquareRoots = false;
      using functor_type = ElementWiseDivide<MV,
        ScalingFactorsViewType, IndexType, takeSquareRoots, takeAbsVal>;
      Kokkos::parallel_for ("elementWiseDivide",
                            range_type (0, numRows),
                            functor_type (X, scalingFactors));
    }
  }
}

void
elementWiseDivideMultiVector (MV& X,
                              const ScalingFactorsViewType& scalingFactors,
                              const bool takeSquareRootsOfScalingFactors,
                              const bool takeAbsoluteValueOfScalingFactors =
                                ! std::is_same<
                                    typename Kokkos::ArithTraits<
                                      typename MV::scalar_type
                                    >::mag_type,
                                    typename ScalingFactorsViewType::non_const_value_type
                                  >::value)
{
  using index_type = typename MV::local_ordinal_type;
  const index_type lclNumRows = static_cast<index_type> (X.getLocalLength ());

  auto X_lcl = X.getLocalViewDevice (Tpetra::Access::ReadWrite);
  if (static_cast<std::size_t> (X.getNumVectors ()) == std::size_t (1)) {
    using pair_type = Kokkos::pair<index_type, index_type>;
    auto X_lcl_1d = Kokkos::subview (X_lcl, pair_type (0, lclNumRows), 0);
    elementWiseDivide (X_lcl_1d, scalingFactors, lclNumRows,
                       takeSquareRootsOfScalingFactors,
                       takeAbsoluteValueOfScalingFactors);
  }
  else {
    elementWiseDivide (X_lcl, scalingFactors, lclNumRows,
                       takeSquareRootsOfScalingFactors,
                       takeAbsoluteValueOfScalingFactors);
  }
}

/* ----------------------------------------------------------------------
   Division functor
------------------------------------------------------------------------- */

class ElementWiseDivide {
public:
  static_assert (ViewType1::Rank == 2, "ViewType1 must be a rank-2 "
                 "Kokkos::View in order to use this specialization.");

  ElementWiseDivide (const ViewType1& X,
                     const ViewType2& scalingFactors) :
    X_ (X),
    scalingFactors_ (scalingFactors)
  {}

  KOKKOS_INLINE_FUNCTION void operator () (const IndexType i) const {
    using val_type = typename ViewType2::non_const_value_type;
    using KAT = Kokkos::ArithTraits<val_type>;
    using mag_type = typename KAT::mag_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    for (IndexType j = 0; j < static_cast<IndexType> (X_.extent (1)); ++j) {
      if (takeAbsoluteValueOfScalingFactors) {
        const mag_type scalFactAbs = KAT::abs (scalingFactors_(i));
        const mag_type scalFinalVal = takeSquareRootsOfScalingFactors ?
          KAM::sqrt (scalFactAbs) : scalFactAbs;
        X_(i,j) = X_(i,j) / scalFinalVal;
      }
      else {
        const val_type scalFact = scalingFactors_(i);
        const val_type scalFinalVal = takeSquareRootsOfScalingFactors ?
          KAT::sqrt (scalFact) : scalFact;
        X_(i,j) = X_(i,j) / scalFinalVal;
      }
    }
  }

private:
  ViewType1 X_;
  typename ViewType2::const_type scalingFactors_;
};

/* ----------------------------------------------------------------------
   Multiplication functor
------------------------------------------------------------------------- */

class ElementWiseMultiply {
public:
  static_assert (ViewType1::Rank == 1, "ViewType1 must be a rank-1 "
                 "Kokkos::View in order to use this specialization.");

  ElementWiseMultiply (const ViewType1& X,
                       const ViewType2& scalingFactors) :
    X_ (X),
    scalingFactors_ (scalingFactors)
  {}

  KOKKOS_INLINE_FUNCTION void operator () (const IndexType i) const {
    using val_type = typename ViewType2::non_const_value_type;
    using KAT = Kokkos::ArithTraits<val_type>;
    using mag_type = typename KAT::mag_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    if (takeAbsoluteValueOfScalingFactors) {
      const mag_type scalFactAbs = KAT::abs (scalingFactors_(i));
      const mag_type scalFinalVal = takeSquareRootsOfScalingFactors ?
        KAM::sqrt (scalFactAbs) : scalFactAbs;
      X_(i) = X_(i) * scalFinalVal;
    }
    else {
      const val_type scalFact = scalingFactors_(i);
      const val_type scalFinalVal = takeSquareRootsOfScalingFactors ?
        KAT::sqrt (scalFact) : scalFact;
      X_(i) = X_(i) * scalFinalVal;
    }
  }

private:
  ViewType1 X_;
  typename ViewType2::const_type scalingFactors_;
};
