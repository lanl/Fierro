#include "matar_blas_lapack_interface.h"
#include <iostream>

/*
 * Transpose a 2D MATAR CArray and return the result in B
 *
 * B = A^{T}
 *
 */
void matar2blas::transpose(const CArray<Real> &A, 
    CArray<Real> &B) {
  // Assert compatibility of inputs
  int 
  num_dim_a = A.order(),
  num_dim_b = B.order();

  bool both_are_2d = num_dim_a == 2 and num_dim_b == 2;
  assert(both_are_2d and "Error: arrays are not both 2D");

  int 
  m_a = A.dims(0),
  n_a = A.dims(1),
  m_b = B.dims(0),
  n_b = B.dims(1);

  bool matching_dims = m_a == n_b and n_a == m_b;
  assert(matching_dims and "Error: array dimensions don't match");

  // Put the transpose of A into B
  for (int i = 0; i < m_a; i++)
    for (int j = 0; j < n_a; j++)
      B(j,i) = A(i,j);
}

/*
 * Compute the matrix-vector multiplication of two MATAR CArrays, A (2D) and x
 * (1D), using BLAS's matrix-vector multiplication routine (dgemv) and return
 * the result in y
 *
 * y = A \times x
 *
 */
void matar2blas::matvec(const CArray<Real> &A, 
    const CArray<Real> &x, CArray<Real> &y) {
  // Assert compatibility of inputs
  int num_dim_a = A.order();
  bool a_is_2d = num_dim_a == 2;
  assert(a_is_2d and "Error: matrix has wrong dimensions");

  int num_dim_x = x.order(), num_dim_y = y.order();
  bool x_and_y_are_1d = num_dim_x == 1 and num_dim_y == 1;
  assert(x_and_y_are_1d and "Error: vectors have wrong dimensions");

  int 
  m_a = A.dims(0),
  n_a = A.dims(1),
  n_x = x.dims(0),
  m_y = y.dims(0);
  bool compatible_dims = m_y == m_a and n_a == n_x;
  assert(compatible_dims and "Error: matrix and vector dims incompatible");
  
  int k = n_a;

  // Copy the contents of A into regular, unwrapped arrays in column major
  // order, which what the BLAS routine requires
  Real *a = new Real[m_a*n_a];
  for (int i = 0; i < m_a; i++)
    for (int j = 0; j < n_a; j++)
      a[i+n_a*j] = A(i,j);

  // Use dgemv from BLAS to compute the matrix-vector multiplication. Note that
  // dgemv computes 
  //
  //   y = \alpha A \times x + \beta y. 
  //
  // In this case, we set \alpha = 1 and \beta = 0.
  const char if_transpose = 'N';  // whether to tranpose matrix A

  Real 
  alpha = 1.0,  // scaling of matrix-vector multiplication
  beta  = 0.0;  // Scaling of input/output vector y

  int 
  lda  = m_a,  // leading dimension of matrix A
  incx = 1,    // increment of vector x
  incy = 1;    // increment of vector y

  blas_gemv( 	
    &if_transpose,
		&m_a, &n_a,
		&alpha,
		a, &lda,
		x.pointer(), &incx,
		&beta,
		y.pointer(), &incy 	
	);

  delete[] a;
}

/*
 * Compute the multiplication of two 2D MATAR CArrays, A and B, using BLAS's
 * matrix-matrix multiplication routine (dgemm) and return the result in C
 *
 * C = A \times B
 *
 */
void matar2blas::matmul(const CArray<Real> &A, 
    const CArray<Real> &B, CArray<Real> &C) {
  // Assert compatibility of inputs
  int 
  num_dim_a = A.order(),
  num_dim_b = B.order(),
  num_dim_c = C.order();
  bool all_are_2d = num_dim_a == 2 and num_dim_b == 2 and num_dim_c == 2;
  assert(all_are_2d and "Error: arrays are not all 2D");

  int 
  m_a = A.dims(0),
  n_a = A.dims(1),
  m_b = B.dims(0),
  n_b = B.dims(1);
  bool compatible_dims = n_a == m_b;
  assert(compatible_dims and "Error: dimensions are incompatible");
  
  int k = n_a;

  // Copy the contents of the input MATAR CArrays into regular, unwrapped
  // arrays in column major order, which what the BLAS routine requires
  Real *a = new Real[m_a*n_a];
  for (int i = 0; i < m_a; i++)
    for (int j = 0; j < n_a; j++)
      a[i+n_a*j] = A(i,j);

  Real *b = new Real[m_b*n_b];
  for (int i = 0; i < m_b; i++)
    for (int j = 0; j < n_b; j++)
      b[i+n_b*j] = B(i,j);

  // Use dgemm from BLAS to compute the matrix-matrix multiplication. Note that
  // dgemm computes 
  //
  //   C = \alpha A \times B + \beta C. 
  //
  // In this case, we set \alpha = 1 and \beta = 0.
  const char if_transpose_a = 'N';  // whether to tranpose matrix A
  const char if_transpose_b = 'N';  // whether to tranpose matrix B

  Real 
  alpha = 1.0,  // scaling of matrix-matrix multiplication
  beta  = 0.0;  // Scaling of input/output matrix C

  int 
  lda = m_a,  // leading dimension of matrix A
  ldb = m_b,  // leading dimension of matrix B
  ldc = m_a;  // leading dimension of matrix C

  Real *c = new Real[m_a*n_b];

  blas_gemm(
      &if_transpose_a, &if_transpose_b, 
      &m_a, &n_b, &k, 
      &alpha, 
      a, &lda, 
      b, &ldb, 
      &beta, 
      c, &lda
    );

  // Copy the result into the output MATAR CArray
  for (int j = 0; j < n_b; j++) {
    for (int i = 0; i < m_a; i++) {
      C(i,j) = c[i+n_b*j];
    }
  }

  delete[] a, b, c;
}

/*
 * Compute the inverse of a 2D MATAR CArray A using LAPACK's LU factorization
 * and solution routines (dgetrf/dgetrf) and return the result in B
 *
 * B = A^{-1}
 *
 */
void matar2lapack::invert(const CArray<Real> &A, CArray<Real> &B) {
  // Assert compatibility of inputs
  int num_dim_a = A.order();
  int num_dim_b = B.order();
  bool both_are_2d = num_dim_a == 2 and num_dim_b == 2;
  assert(both_are_2d and "Error: arrays are not both 2D");
  
  int 
  m_a = A.dims(0),
  n_a = A.dims(1),
  m_b = B.dims(0),
  n_b = B.dims(1);
  bool matrices_have_same_dims = m_a == m_b and n_a == n_b;
  assert(matrices_have_same_dims and "Error: non-matching dimensions");
  
  int m = m_a;
  int n = n_a;
  bool matrices_are_square = m == n;
  assert(matrices_are_square and "Error: can't invert non-square matrix");

  // Copy the contents of the input MATAR CArray into a regular, unwrapped
  // array in column major order, which what the LAPACK routine requires
  Real *a = new Real[m*n];
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      a[i+n*j] = A(i,j);

  // Compute LU factorization of matrix A
  int 
  lda = m,      // leading dimension of matrix A
  factor_info;  // error code

  int *ipiv = new int[n];  // pivot indices in LU factorization

  try {
    lapack_getrf(&m, &n, a, &lda, ipiv, &factor_info);
    if (factor_info != 0) throw FactorizationError(
        "Error: LAPACK LU factorization failed");
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
  }

  // Populate right-hand side array (identity since A = L*U, L*U*A^{-1} = I)
  Real *b = new Real[m*n];
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++)
      b[i+n*j] = (i == j) ? 1 : 0;

  // Use the LU factorization to compute the inverse
  const char if_transpose = 'N';  // whether A should be transposed

  int 
  nrhs  = n,   // number of right-hand sides
  ldb   = m,   // leading dimension of matrix B
  solve_info;  // error code 

  try {
    lapack_getrs(&if_transpose, &n, &nrhs, a, &lda, ipiv, b, &ldb, &solve_info);
    if (solve_info != 0) throw SolutionError(
        "Error: LAPACK linear system solution failed");
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
  }

  // Copy the result into the output MATAR CArray
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < m; i++) {
      B(i,j) = b[i+n*j];
    }
  }
  
  delete[] a, ipiv, b;
}

/*
 * Compute eigenvalues and eigenvectors of real symmetric triadiagonal matrix,
 * as defined by its diagonal and subdiagonal values
 */
void matar2lapack::eig_sym_tri(const CArray<Real> &diag, 
    const CArray<Real> &subdiag, CArray<Real> &eigvals, 
    CArray<Real> &eigvecs) {
  // Assert compatibilty of inputs
  int num_dim_diag    = diag.order();
  int num_dim_subdiag = subdiag.order();
  int num_dim_eigvals = eigvals.order();
  int num_dim_eigvecs = eigvecs.order();
  bool correct_input_shape = num_dim_diag == 1 and num_dim_subdiag == 1;
  bool correct_output_shape = num_dim_eigvals == 1 and num_dim_eigvecs == 2;
  assert(correct_input_shape and correct_output_shape 
      and "Error: incorrect input/output shapes");
  
  int n_diag = diag.dims(0);
  int n_subdiag = subdiag.dims(0);
  int n_eigvals = eigvals.dims(0);
  int m_eigvecs = eigvecs.dims(0);
  int n_eigvecs = eigvecs.dims(1);
  bool correct_input_size = n_diag == n_subdiag + 1;
  bool correct_output_size = n_eigvals == n_diag and m_eigvecs == n_eigvals 
      and n_eigvecs == m_eigvecs;
  assert(correct_input_size and correct_output_size 
      and "Error: incorrect input/output sizes");

  // Copy diagonal entries into raw C array (LAPACK will overwrite the entries)
  Real *diag_copy = new Real[n_diag];
  for (int i = 0; i < n_diag; i++) diag_copy[i] = diag(i);

  // Compute eigenvalues and eigenvectors
  const char compute_eigenvectors = 'V';    // request eigenvectors
  Real *work_array = new Real[2*n_diag-2];  // allocate work array
  int info = 0;

  try {
    lapack_stev(&compute_eigenvectors, &n_diag, diag_copy, 
        subdiag.pointer(), eigvecs.pointer(), &n_diag, 
        work_array, &info);
    if (info != 0) throw FactorizationError(
        "Error: LAPACK eigensolution for symmetric tridiagonal matrix failed");
  } catch(std::exception &e) {
    std::cerr << e.what() << std::endl;
  }

  // Extract eigenvalues from overwritten copy of diagonal entry array
  for (int i = 0; i < n_diag; i++) eigvals(i) = diag_copy[i];

  // Deallocate work array
  delete[] work_array;
}
