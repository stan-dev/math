#ifndef STAN_MATH_REV_FUN_TCROSSPROD_HPP
#define STAN_MATH_REV_FUN_TCROSSPROD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of post-multiplying a matrix by its
 * own transpose.
 *
 * @tparam T Type of the matrix
 * @param M Matrix to multiply.
 * @return M times its transpose.
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto tcrossprod(const T& M) {
  using ret_type = return_var_matrix_t<
      Eigen::Matrix<double, T::RowsAtCompileTime, T::RowsAtCompileTime>, T>;
  arena_t<T> arena_M = M;
  arena_t<ret_type> res = arena_M.val_op() * arena_M.val_op().transpose();

  if (likely(M.size() > 0)) {
    reverse_pass_callback([res, arena_M]() mutable {
      arena_M.adj()
          += (res.adj_op() + res.adj_op().transpose()) * arena_M.val_op();
    });
  }

  return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif
