#ifndef STAN_MATH_PRIM_FUN_CHOLESKY_FACTOR_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_CHOLESKY_FACTOR_CONSTRAIN_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the Cholesky factor of the specified size read from the
 * specified vector.  A total of (N choose 2) + N + (M - N) * N
 * elements are required to read an M by N Cholesky factor.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase and
 * have one compile-time dimension equal to 1)
 * @param x Vector of unconstrained values
 * @param M number of rows
 * @param N number of columns
 * @return Cholesky factor
 */
template <typename T, require_eigen_col_vector_t<T>* = nullptr>
inline Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>
cholesky_factor_constrain(const T& x, int M, int N) {
  using std::exp;
  using T_scalar = value_type_t<T>;
  check_greater_or_equal("cholesky_factor_constrain",
                         "num rows (must be greater or equal to num cols)", M,
                         N);
  check_size_match("cholesky_factor_constrain", "constrain size", x.size(),
                   "((N * (N + 1)) / 2 + (M - N) * N)",
                   ((N * (N + 1)) / 2 + (M - N) * N));
  Eigen::Matrix<T_scalar, Eigen::Dynamic, Eigen::Dynamic> y(M, N);
  int pos = 0;

  const auto& x_ref = to_ref(x);
  for (int m = 0; m < N; ++m) {
    y.row(m).head(m) = x_ref.segment(pos, m);
    pos += m;
    y.coeffRef(m, m) = exp(x_ref.coeff(pos++));
    y.row(m).tail(N - m - 1).setZero();
  }

  for (int m = N; m < M; ++m) {
    y.row(m) = x_ref.segment(pos, N);
    pos += N;
  }
  return y;
}

/**
 * Return the Cholesky factor of the specified size read from the
 * specified vector and increment the specified log probability
 * reference with the log absolute Jacobian determinant adjustment of the
 * transform.  A total of (N choose 2) + N + N * (M - N) free parameters are
 * required to read an M by N Cholesky factor.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase and
 * have one compile-time dimension equal to 1)
 * @param x Vector of unconstrained values
 * @param M number of rows
 * @param N number of columns
 * @param lp Log probability that is incremented with the log absolute Jacobian
 * determinant
 * @return Cholesky factor
 */
template <typename T, require_eigen_vector_t<T>* = nullptr>
inline Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>
cholesky_factor_constrain(const T& x, int M, int N, return_type_t<T>& lp) {
  check_size_match("cholesky_factor_constrain", "x.size()", x.size(),
                   "((N * (N + 1)) / 2 + (M - N) * N)",
                   ((N * (N + 1)) / 2 + (M - N) * N));
  int pos = 0;
  const auto& x_ref = to_ref(x);
  for (int n = 0; n < N; ++n) {
    pos += n;
    lp += x_ref.coeff(pos++);
  }
  return cholesky_factor_constrain(x_ref, M, N);
}

/**
 * Return the Cholesky factor of the specified size read from the specified
 * vector. A total of (N choose 2) + N + N * (M - N) free parameters are
 * required to read an M by N Cholesky factor. If the `Jacobian` parameter is
 * `true`, the log density accumulator is incremented with the log absolute
 * Jacobian determinant of the transform.  All of the transforms are specified
 * with their Jacobians in the *Stan Reference Manual* chapter Constraint
 * Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A type inheriting from `Eigen::DenseBase` or a `var_value` with
 *  inner type inheriting from `Eigen::DenseBase` with compile time dynamic rows
 *  and 1 column
 * @param x Vector of unconstrained values
 * @param M number of rows
 * @param N number of columns
 * @param[in,out] lp log density accumulator
 * @return Cholesky factor
 */
template <bool Jacobian, typename T, require_not_std_vector_t<T>* = nullptr>
inline auto cholesky_factor_constrain(const T& x, int M, int N,
                                      return_type_t<T>& lp) {
  if (Jacobian) {
    return cholesky_factor_constrain(x, M, N, lp);
  } else {
    return cholesky_factor_constrain(x, M, N);
  }
}

/**
 * Return the Cholesky factor of the specified size read from the specified
 * vector. A total of (N choose 2) + N + N * (M - N) free parameters are
 * required to read an M by N Cholesky factor. If the `Jacobian` parameter is
 * `true`, the log density accumulator is incremented with the log absolute
 * Jacobian determinant of the transform.  All of the transforms are specified
 * with their Jacobians in the *Stan Reference Manual* chapter Constraint
 * Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @param x Vector of unconstrained values
 * @param M number of rows
 * @param N number of columns
 * @param[in,out] lp log density accumulator
 * @return Cholesky factor
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto cholesky_factor_constrain(const T& x, int M, int N,
                                      return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(x, [&lp, M, N](auto&& v) {
    return cholesky_factor_constrain<Jacobian>(v, M, N, lp);
  });
}

}  // namespace math
}  // namespace stan

#endif
