#ifndef STAN_MATH_REV_FUN_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP
#define STAN_MATH_REV_FUN_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/multiply_lower_tri_self_transpose.hpp>

namespace stan {
namespace math {

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto multiply_lower_tri_self_transpose(const T& L) {
  using ret_type = return_var_matrix_t<T>;
  if (L.size() == 0) {
    return ret_type(decltype(multiply_lower_tri_self_transpose(value_of(L)))());
  }

  arena_t<T> arena_L = L;
  arena_t<promote_scalar_t<double, T>> arena_L_val
      = arena_L.val().template triangularView<Eigen::Lower>();

  arena_t<ret_type> res = arena_L_val.template triangularView<Eigen::Lower>()
                          * arena_L_val.transpose();

  reverse_pass_callback([res, arena_L, arena_L_val]() mutable {
    arena_L.adj() += ((res.adj().transpose() + res.adj())
                      * arena_L_val.template triangularView<Eigen::Lower>())
                         .template triangularView<Eigen::Lower>();
  });

  return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif
