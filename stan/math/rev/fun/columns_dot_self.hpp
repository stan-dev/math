#ifndef STAN_MATH_REV_FUN_COLUMNS_DOT_SELF_HPP
#define STAN_MATH_REV_FUN_COLUMNS_DOT_SELF_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each column of a matrix with itself.
 *
 * @tparam Mat Type of the matrix.
 * @param x Matrix.
 */
template <typename Mat, require_rev_matrix_t<Mat>* = nullptr>
inline auto columns_dot_self(const Mat& x) {
  arena_t<Mat> arena_x = x;
  arena_t<Eigen::MatrixXd> arena_x_val = value_of(arena_x);
  using ret_inner_type = decltype(arena_x_val.colwise().squaredNorm());
  using ret_type = promote_var_matrix_t<ret_inner_type, Mat>;
  arena_t<ret_type> res = arena_x_val.colwise().squaredNorm();

  if (x.size() == 0) {
    return ret_type(res);
  } else {
    reverse_pass_callback([res, arena_x, arena_x_val]() mutable {
        arena_x.adj() += arena_x_val * (2 * res.adj()).asDiagonal();
      });
    return ret_type(res);
  }
}

}  // namespace math
}  // namespace stan
#endif
