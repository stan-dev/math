#ifndef STAN_MATH_REV_FUN_REP_MATRIX_HPP
#define STAN_MATH_REV_FUN_REP_MATRIX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <>
inline auto rep_matrix<double, var_value<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>(const double& x, int m, int n) {
  check_nonnegative("rep_matrix", "rows", m);
  check_nonnegative("rep_matrix", "cols", n);
  using eig_mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  return var_value<eig_mat>(eig_mat::Constant(m, n, x));
}

template <>
inline auto rep_matrix<var, var_value<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>(const var& x, int m, int n) {
  check_nonnegative("rep_matrix", "rows", m);
  check_nonnegative("rep_matrix", "cols", n);
  using eig_mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  return make_callback_var(eig_mat::Constant(m, n, x.val()), [x](auto& rep) mutable {
    x.adj() += rep.adj().sum();
  });
}

}  // namespace math
}  // namespace stan

#endif
