#pragma once

#include "common.h"

template <typename NumType>
struct LagrangeElement {
  // Element specifications
  static const SizeType Nd = 3;
  SizeType Np;
  SizeType N;
  SizeType Ne;
  const NumType *Z;
  NumType *w;

  // Arrays for converting from flat to multidimensional indices
  SizeType ijk[Nd];
  SizeType rad[Nd];

  // Work array for intermediate coefficients
  NumType *C;

  LagrangeElement(const SizeType, const NumType *);
  ~LagrangeElement();

  // Basis functions and basis function gradients
  NumType eval_basis(const SizeType, const NumType *);
  void eval_grad_basis(const SizeType, const NumType *, NumType *);

  // Function approximation over element
  NumType eval_approx(const NumType *, const NumType *);
  void eval_grad_approx(const NumType *, const NumType *, NumType *);

  // Jacobian of spatial mapping
  void eval_jac(const NumType *, const NumType *, const NumType *, 
      const NumType *, NumType *);
  NumType eval_det_jac(const NumType *, const NumType *, const NumType *, 
      const NumType *);
  void eval_inv_jac(const NumType *, const NumType *, const NumType *, 
      const NumType *, NumType *);
};
