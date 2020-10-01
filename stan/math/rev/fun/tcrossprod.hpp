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
inline promote_var_matrix_t<Eigen::MatrixXd, T> tcrossprod(const T& M) {
  if (M.size() == 0) {
    return M;
  }

  arena_t<plain_type_t<T>> arena_M = M;
  arena_t<Eigen::MatrixXd> arena_M_val = value_of(arena_M);

  arena_t<plain_type_t<T>> res = arena_M_val * arena_M_val.transpose();

  reverse_pass_callback([res, arena_M, arena_M_val]() mutable {
    const auto& adj = to_ref(res.adj());
    arena_M.adj() += (adj + adj.transpose()) * arena_M_val;
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
