#ifndef STAN_MATH_PRIM_FUN_READ_CORR_L_HPP
#define STAN_MATH_PRIM_FUN_READ_CORR_L_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <cstddef>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the Cholesky factor of the correlation matrix of the
 * specified dimensionality corresponding to the specified
 * canonical partial correlations.
 *
 * It is generally better to work with the Cholesky factor rather
 * than the correlation matrix itself when the determinant,
 * inverse, etc. of the correlation matrix is needed for some
 * statistical calculation.
 *
 * <p>See <code>read_corr_matrix(Array, size_t, T)</code>
 * for more information.
 *
 * @tparam T type of the array (must be derived from \c Eigen::ArrayBase and
 * have one compile-time dimension equal to 1)
 * @param CPCs The (K choose 2) canonical partial correlations in
 * (-1, 1).
 * @param K Dimensionality of correlation matrix.
 * @return Cholesky factor of correlation matrix for specified
 * canonical partial correlations.
 */
template <typename T, require_eigen_vector_t<T>* = nullptr>
Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic> read_corr_L(
    const T& CPCs,  // on (-1, 1)
    size_t K) {
  using T_scalar = value_type_t<T>;
  if (K == 0) {
    return {};
  }
  if (K == 1) {
    return Eigen::Matrix<T_scalar, Eigen::Dynamic, Eigen::Dynamic>::Identity(1,
                                                                             1);
  }

  using std::sqrt;
  Eigen::Array<T_scalar, Eigen::Dynamic, 1> temp;
  Eigen::Array<T_scalar, Eigen::Dynamic, 1> acc(K - 1);
  acc.setOnes();
  // Cholesky factor of correlation matrix
  Eigen::Matrix<T_scalar, Eigen::Dynamic, Eigen::Dynamic> L(K, K);
  L.setZero();

  size_t position = 0;
  size_t pull = K - 1;

  L(0, 0) = 1.0;
  L.col(0).tail(pull) = temp = CPCs.head(pull);
  acc.tail(pull) = T_scalar(1.0) - temp.square();
  for (size_t i = 1; i < (K - 1); i++) {
    position += pull;
    pull--;
    temp = CPCs.segment(position, pull);
    L(i, i) = sqrt(acc(i - 1));
    L.col(i).tail(pull) = temp * acc.tail(pull).sqrt();
    acc.tail(pull) *= T_scalar(1.0) - temp.square();
  }
  L(K - 1, K - 1) = sqrt(acc(K - 2));
  return L;
}

/**
 * Return the Cholesky factor of the correlation matrix of the
 * specified dimensionality corresponding to the specified
 * canonical partial correlations, incrementing the specified
 * scalar reference with the log absolute determinant of the
 * Jacobian of the transformation.
 *
 * <p>The implementation is Ben Goodrich's Cholesky
 * factor-based approach to the C-vine method of:
 *
 * <ul><li> Daniel Lewandowski, Dorota Kurowicka, and Harry Joe,
 * Generating random correlation matrices based on vines and
 * extended onion method Journal of Multivariate Analysis 100
 * (2009) 1989–2001 </li></ul>
 *
 * @tparam T type of the array (must be derived from \c Eigen::ArrayBase and
 * have one compile-time dimension equal to 1)
 * @param CPCs The (K choose 2) canonical partial correlations in
 * (-1, 1).
 * @param K Dimensionality of correlation matrix.
 * @param log_prob Reference to variable to increment with the log
 * Jacobian determinant.
 * @return Cholesky factor of correlation matrix for specified
 * partial correlations.
 */
template <typename T, require_eigen_vector_t<T>* = nullptr>
Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic> read_corr_L(
    const T& CPCs, size_t K, value_type_t<T>& log_prob) {
  using T_scalar = value_type_t<T>;
  if (K == 0) {
    return {};
  }
  if (K == 1) {
    return Eigen::Matrix<T_scalar, Eigen::Dynamic, Eigen::Dynamic>::Identity(1,
                                                                             1);
  }

  const Eigen::Ref<const plain_type_t<T>>& CPCs_ref = CPCs;
  size_t pos = 0;
  T_scalar acc = 0;
  // no need to abs() because this Jacobian determinant
  // is strictly positive (and triangular)
  // see inverse of Jacobian in equation 11 of LKJ paper
  for (size_t k = 1; k <= (K - 2); k++) {
    for (size_t i = k + 1; i <= K; i++) {
      acc += (K - k - 1) * log1m(square(CPCs_ref(pos)));
      pos++;
    }
  }

  log_prob += 0.5 * acc;
  return read_corr_L(CPCs_ref, K);
}

}  // namespace math
}  // namespace stan

#endif
