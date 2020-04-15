#ifndef STAN_MATH_PRIM_FUN_COMPLEX_SCHUR_DECOMPOSE_HPP
#define STAN_MATH_PRIM_FUN_COMPLEX_SCHUR_DECOMPOSE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Return the unitary matrix of the complex Schur decomposition of the
 * specified square matrix.
 *
 * The complex Schur decomposition of a square matrix `A` produces a
 * complex unitary matrix `U` and a complex upper-triangular Schur
 * form matrix `T` such that `A = U * T * inv(U)`.  Further, the
 * unitary matrix's inverse is equal to its conjugate transpose,
 * `inv(U) = U*`, where `U*(i, j) = conj(U(j, i))`
 *
 * @tparam value type of matrix
 * @param m real matrix to decompose
 * @return complex unitary matrix of the complex Schur decomposition of the
 * specified matrix
 * @see complex_schur_decompose_t
 */
template <typename T>
Eigen::Matrix<complex_return_t<T>, -1, -1> complex_schur_decompose_u(
    const Eigen::Matrix<T, -1, -1>& m) {
  check_nonzero_size("complex_schur_decompose_u", "m", m);
  check_square("complex_schur_decompose_u", "m", m);
  Eigen::ComplexSchur<Eigen::Matrix<T, -1, -1>> cs(m);
  return cs.matrixU();
}

/**
 * Return the Schur form matrix of the complex Schur decomposition of the
 * specified square matrix.
 *
 * @tparam value type of matrix
 * @param m real matrix to decompose
 * @return Schur form matrix of the complex Schur decomposition of the
 * specified matrix
 * @see complex_schur_decompose_u for a definition of the complex
 * Schur decomposition
 */
template <typename T>
Eigen::Matrix<complex_return_t<T>, -1, -1> complex_schur_decompose_t(
    const Eigen::Matrix<T, -1, -1>& m) {
  check_nonzero_size("complex_schur_decompose_t", "m", m);
  check_square("complex_schur_decompose_t", "m", m);
  Eigen::ComplexSchur<Eigen::Matrix<T, -1, -1>> cs(m);
  return cs.matrixT();
}

}  // namespace math
}  // namespace stan
#endif
