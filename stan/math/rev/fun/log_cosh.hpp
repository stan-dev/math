#ifndef STAN_MATH_REV_FUN_LOG_COSH_HPP
#define STAN_MATH_REV_FUN_LOG_COSH_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/prim/fun/log_cosh.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the hyperbolic cosine of the argument.
 *
 * @tparam T arithmetic type
 * @param x `var_value` with inner arithmetic type
 * @return `log(cosh(x))`, with derivative `tanh(x)`
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline var log_cosh(const var_value<T>& x) {
  const internal::log_cosh_parts parts = internal::log_cosh_tanh_sech(x.val());
  const double tanh_x = parts.tanh;
  return make_callback_var(parts.log_cosh, [x, tanh_x](auto& vi) mutable {
    x.adj() += vi.adj() * tanh_x;
  });
}

/**
 * Return the elementwise natural logarithm of the hyperbolic cosine of the
 * argument.
 *
 * @tparam T `var_value` with inner Eigen type
 * @param x `var_value` with inner Eigen type
 * @return elementwise `log(cosh(x))`, with derivative `tanh(x)`
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto log_cosh(T&& x) {
  using val_t = plain_type_t<value_type_t<T>>;
  auto arena_x = to_arena(std::forward<T>(x));
  val_t value(arena_x.rows(), arena_x.cols());
  arena_t<val_t> tanh_x(arena_x.rows(), arena_x.cols());
  for (Eigen::Index i = 0; i < arena_x.size(); ++i) {
    const internal::log_cosh_parts parts
        = internal::log_cosh_tanh_sech(arena_x.val().coeff(i));
    value.coeffRef(i) = parts.log_cosh;
    tanh_x.coeffRef(i) = parts.tanh;
  }
  return make_callback_var(
      std::move(value), [arena_x, tanh_x](auto& vi) mutable {
        arena_x.adj().array() += vi.adj().array() * tanh_x.array();
      });
}

}  // namespace math
}  // namespace stan
#endif
