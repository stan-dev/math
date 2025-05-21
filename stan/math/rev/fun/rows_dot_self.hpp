#ifndef STAN_MATH_REV_FUN_ROWS_DOT_SELF_HPP
#define STAN_MATH_REV_FUN_ROWS_DOT_SELF_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/rows_dot_self.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each row of a matrix with itself.
 *
 * @tparam Mat An Eigen matrix with a `var` scalar type.
 * @param x Matrix.
 */
template <typename Mat, require_eigen_vt<is_var, Mat>* = nullptr>
inline Eigen::Matrix<var, Mat::RowsAtCompileTime, 1> rows_dot_self(
    const Mat& x) {
  Eigen::Matrix<var, Mat::RowsAtCompileTime, 1> ret(x.rows());
  for (size_type i = 0; i < x.rows(); i++) {
    ret(i) = dot_self(x.row(i));
  }
  return ret;
}

/**
 * Returns the dot product of row row of a matrix with itself.
 *
 * @tparam Mat A `var_value<>` with an inner matrix type.
 * @param x Matrix.
 */
template <typename Mat, require_var_matrix_t<Mat>* = nullptr>
inline auto rows_dot_self(const Mat& x) {
  using ret_type = var_value<Eigen::VectorXd>;
  arena_t<ret_type> res = x.val().rowwise().squaredNorm();
  if (x.size() >= 0) {
    reverse_pass_callback([res, x]() mutable {
      x.adj() += (2 * res.adj()).asDiagonal() * x.val();
    });
  }
  return res;
}

}  // namespace math
}  // namespace stan
#endif
