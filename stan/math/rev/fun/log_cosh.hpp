#ifndef STAN_MATH_REV_FUN_LOG_COSH_HPP
#define STAN_MATH_REV_FUN_LOG_COSH_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
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
  const std::pair<double, double> value_and_tanh
      = internal::log_cosh_and_tanh(x.val());
  const double tanh_x = value_and_tanh.second;
  return make_callback_var(value_and_tanh.first, [x, tanh_x](auto& vi) mutable {
    x.adj() += vi.adj() * tanh_x;
  });
}

/**
 * Return the elementwise natural logarithm of the hyperbolic cosine of the
 * argument.
 *
 * @tparam T a type derived from `Eigen::EigenBase`
 * @param x `var_value` with inner Eigen type
 * @return elementwise `log(cosh(x))`, with derivative `tanh(x)`
 */
template <typename T, require_eigen_t<T>* = nullptr>
inline auto log_cosh(const var_value<T>& x) {
  using val_t = plain_type_t<T>;
  val_t value(x.rows(), x.cols());
  arena_t<val_t> tanh_x(x.rows(), x.cols());
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    const std::pair<double, double> value_and_tanh
        = internal::log_cosh_and_tanh(x.val().coeff(i));
    value.coeffRef(i) = value_and_tanh.first;
    tanh_x.coeffRef(i) = value_and_tanh.second;
  }
  return make_callback_var(std::move(value), [x, tanh_x](auto& vi) mutable {
    x.adj().array() += vi.adj().array() * tanh_x.array();
  });
}

}  // namespace math
}  // namespace stan
#endif
