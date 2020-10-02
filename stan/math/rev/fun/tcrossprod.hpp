#ifndef STAN_MATH_REV_FUN_TCROSSPROD_HPP
#define STAN_MATH_REV_FUN_TCROSSPROD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/Eigen_NumTraits.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
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
  using ret_type = promote_var_matrix_t<
      Eigen::Matrix<double, T::RowsAtCompileTime, T::RowsAtCompileTime>, T>;
  arena_t<plain_type_t<T>> arena_M = M;
  auto arena_M_val = to_arena(value_of(arena_M));

  arena_t<ret_type> res = arena_M_val * arena_M_val.transpose();

  if (M.size() > 0) {
    reverse_pass_callback([res, arena_M, arena_M_val]() mutable {
      const auto& adj = to_ref(res.adj());
      arena_M.adj() += (adj + adj.transpose()) * arena_M_val;
    });
  }

  return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif
