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
 * @tparam M type of matrix
 * @param m real matrix to decompose
 * @return complex unitary matrix of the complex Schur decomposition of the
 * specified matrix
 * @see complex_schur_decompose_t
 */
template <typename M, require_eigen_dense_dynamic_t<M>* = nullptr>
inline Eigen::Matrix<complex_return_t<scalar_type_t<M>>, -1, -1>
complex_schur_decompose_u(const M& m) {
  if (unlikely(m.size() == 0)) {
    return m;
  }
  check_square("complex_schur_decompose_u", "m", m);
  using MatType = Eigen::Matrix<scalar_type_t<M>, -1, -1>;
  // copy because ComplexSchur requires Eigen::Matrix type
  Eigen::ComplexSchur<MatType> cs{MatType(m)};
  return cs.matrixU();
}

/**
 * Return the Schur form matrix of the complex Schur decomposition of the
 * specified square matrix.
 *
 * @tparam M type of matrix
 * @param m real matrix to decompose
 * @return Schur form matrix of the complex Schur decomposition of the
 * specified matrix
 * @see complex_schur_decompose_u for a definition of the complex
 * Schur decomposition
 */
template <typename M, require_eigen_dense_dynamic_t<M>* = nullptr>
inline Eigen::Matrix<complex_return_t<scalar_type_t<M>>, -1, -1>
complex_schur_decompose_t(const M& m) {
  if (unlikely(m.size() == 0)) {
    return m;
  }
  check_square("complex_schur_decompose_t", "m", m);
  using MatType = Eigen::Matrix<scalar_type_t<M>, -1, -1>;
  // copy because ComplexSchur requires Eigen::Matrix type
  Eigen::ComplexSchur<MatType> cs{MatType(m), false};
  return cs.matrixT();
}

/**
 * Return the complex Schur decomposition of the
 * specified square matrix.
 *
 * The complex Schur decomposition of a square matrix `A` produces a
 * complex unitary matrix `U` and a complex upper-triangular Schur
 * form matrix `T` such that `A = U * T * inv(U)`.  Further, the
 * unitary matrix's inverse is equal to its conjugate transpose,
 * `inv(U) = U*`, where `U*(i, j) = conj(U(j, i))`
 *
 * @tparam M type of matrix
 * @param m real matrix to decompose
 * @return a tuple (U,T) where U is the complex unitary matrix of the complex
 * Schur decomposition of `m` and T is the Schur form matrix of
 * the complex Schur decomposition of `m`
 */
template <typename M, require_eigen_dense_dynamic_t<M>* = nullptr>
inline std::tuple<Eigen::Matrix<complex_return_t<scalar_type_t<M>>, -1, -1>,
                  Eigen::Matrix<complex_return_t<scalar_type_t<M>>, -1, -1>>
complex_schur_decompose(const M& m) {
  if (unlikely(m.size() == 0)) {
    return std::make_tuple(m, m);
  }
  check_square("complex_schur_decompose", "m", m);
  using MatType = Eigen::Matrix<scalar_type_t<M>, -1, -1>;
  // copy because ComplexSchur requires Eigen::Matrix type
  Eigen::ComplexSchur<MatType> cs{MatType(m)};
  return std::make_tuple(std::move(cs.matrixU()), std::move(cs.matrixT()));
}

}  // namespace math
}  // namespace stan
#endif
