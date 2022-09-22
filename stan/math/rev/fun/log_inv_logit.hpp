#ifndef STAN_MATH_REV_FUN_LOG_INV_LOGIT_HPP
#define STAN_MATH_REV_FUN_LOG_INV_LOGIT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/log_inv_logit.hpp>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the inverse logit of the
 * specified argument.
 *
 * @tparam T An arithmetic type
 * @param u `var_value` with inner arithmetic type
 * @return log inverse logit of the argument
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline auto log_inv_logit(const var_value<T>& u) {
  return make_callback_var(log_inv_logit(u.val()), [u](auto& vi) mutable {
    u.adj() += vi.adj() * inv_logit(-u.val());
  });
}

/**
 * Return the natural logarithm of the inverse logit of the
 * specified argument.
 *
 * @tparam T A type derived from `Eigen::EigenBase`
 * @param u `var_value` with inner Eigen type
 * @return log inverse logit of the argument
 */
template <typename T, require_eigen_t<T>* = nullptr>
inline auto log_inv_logit(const var_value<T>& u) {
  return make_callback_var(log_inv_logit(u.val()), [u](auto& vi) mutable {
    u.adj().array() += vi.adj().array() * inv_logit(-u.val()).array();
  });
}

}  // namespace math
}  // namespace stan
#endif
