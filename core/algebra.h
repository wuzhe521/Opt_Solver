#ifndef GONS_ALGEBRA_H
#define GONS_ALGEBRA_H

#include "config.h"
#include "matrix.h"
#include "utilites.h"
#include "vector.h"

namespace gons {
// Algebraic operations and functions
template <typename T, GONS_UINT R, GONS_UINT C>
Vector<T, R> SolveLinearSystem(Matrix<T, R, C> &A, Vector<T, R> &b) {
  // Implement a method to solve the linear system Ax = b
  // For simplicity, we can use Gaussian elimination here
  // Note: This is a placeholder implementation and may not be efficient for
  // large systems

  // Back substitution
  Vector<T, R> x;
  Matrix<T, R, C> A_inv = A.Inverse();
  for (GONS_UINT i = 0; i < R; ++i) {
    x(i) = 0;
    for (GONS_UINT j = 0; j < C; ++j) {
      x(i) += A_inv(i, j) * b(j);
    }
  }
  return x;
}

// LU Decomposition

// QR Decomposition

// Eigenvalue Computation

// Singular Value Decomposition (SVD)

// Cholesky Decomposition

// Identity matrix
template <typename T, GONS_UINT N> Matrix<T, N, N> Identity() {
  Matrix<T, N, N> I;
  for (GONS_UINT i = 0; i < N; ++i) {
    I(i, i) = T(1);
  }
  return I;
}
template <typename T, GONS_UINT N> Matrix<T, N, N> Zero() {
  Matrix<T, N, N> I;
  for (GONS_UINT i = 0; i < N; ++i) {
    I(i, i) = T(0);
  }
  return I;
}

template <typename T, GONS_UINT N> Vector<T, N> Zero() {
  Vector<T, N> I;
  for (GONS_UINT i = 0; i < N; ++i) {
    I(i) = T(0);
  }
  return I;
}
} // namespace gons

#endif