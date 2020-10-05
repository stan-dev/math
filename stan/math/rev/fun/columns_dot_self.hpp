#ifndef STAN_MATH_REV_FUN_COLUMNS_DOT_SELF_HPP
#define STAN_MATH_REV_FUN_COLUMNS_DOT_SELF_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each column of a matrix with itself.
 *
 * @tparam Mat Type of the matrix.
 * @param x Matrix.
 */
template <typename Mat, require_eigen_vt<is_var, Mat>* = nullptr>
inline Eigen::Matrix<var, 1, Mat::ColsAtCompileTime> columns_dot_self(
    const Mat& x) {
  const auto& x_ref = to_ref(x);

  Eigen::Matrix<var, 1, Mat::ColsAtCompileTime> ret(1, x_ref.cols());
  for (size_type i = 0; i < x_ref.cols(); i++) {
    ret(i) = dot_self(x_ref.col(i));
  }
  return ret;
}

}  // namespace math
}  // namespace stan
#endif
