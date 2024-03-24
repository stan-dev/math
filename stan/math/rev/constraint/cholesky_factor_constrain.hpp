#ifndef STAN_MATH_REV_CONSTRAINT_CHOLESKY_FACTOR_CONSTRAIN_HPP
#define STAN_MATH_REV_CONSTRAINT_CHOLESKY_FACTOR_CONSTRAIN_HPP

#include <stan/math/rev/core.hpp>
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
 * @tparam T type of input vector (must be a `var_value<S>` where `S`
 *  inherits from EigenBase)
 * @param x Vector of unconstrained values
 * @param M number of rows
 * @param N number of columns
 * @return Cholesky factor
 */
template <typename T, require_var_vector_t<T>* = nullptr>
var_value<Eigen::MatrixXd> cholesky_factor_constrain(const T& x, int M, int N) {
  using std::exp;
  check_greater_or_equal("cholesky_factor_constrain",
                         "num rows (must be greater or equal to num cols)", M,
                         N);
  check_size_match("cholesky_factor_constrain", "x.size()", x.size(),
                   "((N * (N + 1)) / 2 + (M - N) * N)",
                   ((N * (N + 1)) / 2 + (M - N) * N));
  arena_t<Eigen::MatrixXd> y_val = Eigen::MatrixXd::Zero(M, N);

  int pos = 0;
  for (int m = 0; m < N; ++m) {
    y_val.row(m).head(m) = x.val().segment(pos, m);
    pos += m;
    y_val.coeffRef(m, m) = exp(x.val().coeff(pos));
    pos++;
  }

  for (int m = N; m < M; ++m) {
    y_val.row(m) = x.val().segment(pos, N);
    pos += N;
  }

  var_value<Eigen::MatrixXd> y = y_val;

  reverse_pass_callback([x, M, N, y]() mutable {
    int pos = x.size();
    for (int m = M - 1; m >= N; --m) {
      pos -= N;
      x.adj().segment(pos, N) += y.adj().row(m);
    }

    for (int m = N - 1; m >= 0; --m) {
      pos--;
      x.adj().coeffRef(pos) += y.adj().coeff(m, m) * y.val().coeff(m, m);
      pos -= m;
      x.adj().segment(pos, m) += y.adj().row(m).head(m);
    }
  });

  return y;
}

/**
 * Return the Cholesky factor of the specified size read from the
 * specified vector and increment the specified log probability
 * reference with the log Jacobian adjustment of the transform.  A total
 * of (N choose 2) + N + N * (M - N) free parameters are required to read
 * an M by N Cholesky factor.
 *
 * @tparam T type of input vector (must be a `var_value<S>` where `S`
 *  inherits from EigenBase)
 * @param x Vector of unconstrained values
 * @param M number of rows
 * @param N number of columns
 * @param[out] lp Log density that is incremented with the log Jacobian
 * @return Cholesky factor
 */
template <typename T, require_var_vector_t<T>* = nullptr>
var_value<Eigen::MatrixXd> cholesky_factor_constrain(const T& x, int M, int N,
                                                     scalar_type_t<T>& lp) {
  check_size_match("cholesky_factor_constrain", "x.size()", x.size(),
                   "((N * (N + 1)) / 2 + (M - N) * N)",
                   ((N * (N + 1)) / 2 + (M - N) * N));
  int pos = 0;
  double lp_val = 0.0;
  for (int n = 0; n < N; ++n) {
    pos += n;
    lp_val += x.val().coeff(pos);
    pos++;
  }
  lp += lp_val;

  reverse_pass_callback([x, N, lp]() mutable {
    int pos = 0;
    for (int n = 0; n < N; ++n) {
      pos += n;
      x.adj().coeffRef(pos) += lp.adj();
      pos++;
    }
  });

  return cholesky_factor_constrain(x, M, N);
}

}  // namespace math
}  // namespace stan

#endif
