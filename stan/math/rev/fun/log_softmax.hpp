#ifndef STAN_MATH_REV_FUN_LOG_SOFTMAX_HPP
#define STAN_MATH_REV_FUN_LOG_SOFTMAX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/log_softmax.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>

namespace stan {
namespace math {

/**
 * Return the log softmax of the specified vector or row vector.
 *
 * @tparam T a `var_value` or Eigen vector/row_vector with `var` scalar
 * @param x input
 * @return log softmax of the input
 * @throw std::domain_error if the input size is 0
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto log_softmax(T&& x) {
  check_nonzero_size("log_softmax", "x", x);
  auto x_arena = to_arena(std::forward<T>(x));
  using return_t
      = return_var_matrix_t<plain_type_t<decltype(x_arena.val())>, T>;
  arena_t<return_t> res = log_softmax(x_arena.val());
  reverse_pass_callback([x_arena, res]() mutable {
    const auto s = res.val().array().exp().eval();
    const auto& res_adj = to_ref(res.adj());
    x_arena.adj().array() += res_adj.array() - res_adj.sum() * s;
  });
  return return_t(res);
}

/**
 * Return the log softmax of each vector in an array.
 *
 * @tparam T `std::vector` whose scalar type is `var`
 * @param x array of vectors to transform
 * @return array of log softmax results
 * @throw std::domain_error if any element size is 0
 */
template <typename T, require_std_vector_st<is_var, T>* = nullptr>
inline auto log_softmax(T&& x) {
  return apply_vector_unary<T>::apply(std::forward<T>(x), [](auto&& v) {
    return log_softmax(std::forward<decltype(v)>(v));
  });
}

}  // namespace math
}  // namespace stan
#endif
