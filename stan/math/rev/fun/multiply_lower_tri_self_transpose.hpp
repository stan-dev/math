#ifndef STAN_MATH_REV_FUN_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP
#define STAN_MATH_REV_FUN_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline plain_type_t<T> multiply_lower_tri_self_transpose(const T& L) {
  if (L.size() == 0)
    return L;

  arena_t<plain_type_t<T>> arena_L = L;
  arena_t<Eigen::MatrixXd> arena_L_val
      = value_of(arena_L).template triangularView<Eigen::Lower>();

  arena_t<plain_type_t<T>> res = arena_L_val.template triangularView<Eigen::Lower>()
                       * arena_L_val.transpose();

  reverse_pass_callback([res, arena_L, arena_L_val]() mutable {
    const auto& adj = to_ref(res.adj());

    const auto& adjL = (adj.transpose() + adj)
      * arena_L_val.template triangularView<Eigen::Lower>();

    arena_L.adj() += adjL.template triangularView<Eigen::Lower>();
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
