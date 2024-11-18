#ifndef STAN_MATH_REV_FUN_COLUMNS_DOT_SELF_HPP
#define STAN_MATH_REV_FUN_COLUMNS_DOT_SELF_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/columns_dot_self.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each column of a matrix with itself.
 *
 * @tparam Mat An Eigen matrix with a `var` scalar type.
 * @param x Matrix.
 */
template <typename Mat, require_eigen_vt<is_var, Mat>* = nullptr>
inline Eigen::Matrix<var, 1, Mat::ColsAtCompileTime> columns_dot_self(
    const Mat& x) {
  Eigen::Matrix<var, 1, Mat::ColsAtCompileTime> ret(1, x.cols());
  for (size_type i = 0; i < x.cols(); i++) {
    ret(i) = dot_self(x.col(i));
  }
  return ret;
}

/**
 * Returns the dot product of each column of a matrix with itself.
 *
 * @tparam Mat A `var_value<>` with an inner matrix type.
 * @param x Matrix.
 */
template <typename Mat, require_var_matrix_t<Mat>* = nullptr>
inline auto columns_dot_self(const Mat& x) {
  using ret_type
      = return_var_matrix_t<decltype(x.val().colwise().squaredNorm()), Mat>;
  arena_t<ret_type> res = x.val().colwise().squaredNorm();
  if (x.size() >= 0) {
    reverse_pass_callback([res, x]() mutable {
      x.adj() += x.val() * (2 * res.adj()).asDiagonal();
    });
  }
  return res;
}

}  // namespace math
}  // namespace stan
#endif
