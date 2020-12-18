#ifndef STAN_MATH_REV_FUN_MATRIX_POWER_HPP
#define STAN_MATH_REV_FUN_MATRIX_POWER_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the nth power of the specific matrix. M^n = M * M * ... * M.
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param[in] M a square matrix
 * @param[in] n exponent
 * @return nth power of M
 * @throw std::domain_error if the matrix contains NaNs or infinities.
 * @throw std::invalid_argument if the exponent is negative or the matrix is not
 * square.
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline plain_type_t<T> matrix_power(const T& M, const int n) {
  check_square("matrix_power", "M", M);
  check_nonnegative("matrix_power", "n", n);

  if (M.size() == 0)
    return M;

  const auto& M_ref = to_ref(M);
  check_finite("matrix_power", "M", M_ref);

  size_t N = M.rows();

  if (n == 0)
    return Eigen::MatrixXd::Identity(N, N);

  if (n == 1)
    return M_ref;

  arena_t<std::vector<Eigen::MatrixXd>> arena_powers(n + 1);
  arena_t<plain_type_t<T>> arena_M = M_ref;

  arena_powers[0] = Eigen::MatrixXd::Identity(N, N);
  arena_powers[1] = M_ref.val();
  for (size_t i = 2; i <= n; ++i) {
    arena_powers[i] = arena_powers[1] * arena_powers[i - 1];
  }
  using ret_type = return_var_matrix_t<T>;
  arena_t<ret_type> res = arena_powers[arena_powers.size() - 1];

  reverse_pass_callback([arena_M, n, res, arena_powers]() mutable {
    const auto& M_val = arena_powers[1];
    Eigen::MatrixXd adj_C = res.adj();
    Eigen::MatrixXd adj_M = Eigen::MatrixXd::Zero(M_val.rows(), M_val.cols());
    for (size_t i = n; i > 1; --i) {
      adj_M += adj_C * arena_powers[i - 1].transpose();
      adj_C = M_val.transpose() * adj_C;
    }
    arena_M.adj() += adj_M + adj_C;
  });

  return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif
