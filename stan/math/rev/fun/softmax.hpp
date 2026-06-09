#ifndef STAN_MATH_REV_FUN_SOFTMAX_HPP
#define STAN_MATH_REV_FUN_SOFTMAX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/softmax.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>

namespace stan {
namespace math {

/**
 * Return the softmax of the specified vector or row vector.
 *
 * @tparam T a `var_value` or Eigen vector/row_vector with `var` scalar
 * @param x input
 * @return softmax of the input, or an empty result if the input is empty
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto softmax(T&& x) {
  auto x_arena = to_arena(std::forward<T>(x));
  if (x_arena.size() == 0)
    return x_arena;
  using return_t
      = return_var_matrix_t<plain_type_t<decltype(x_arena.val())>, T>;
  arena_t<return_t> res = softmax(x_arena.val());
  reverse_pass_callback([x_arena, res]() mutable {
    x_arena.adj().array()
        += res.val().array() * (res.adj().array() - res.val().dot(res.adj()));
  });
  return res;
}

/**
 * Return the softmax of each vector in an array.
 *
 * @tparam T `std::vector` whose scalar type is `var`
 * @param x array of vectors to transform
 * @return array of softmax results
 */
template <typename T, require_std_vector_st<is_var, T>* = nullptr>
inline auto softmax(T&& x) {
  return apply_vector_unary<T>::apply(std::forward<T>(x), [](auto&& v) {
    return softmax(std::forward<decltype(v)>(v));
  });
}

}  // namespace math
}  // namespace stan
#endif
