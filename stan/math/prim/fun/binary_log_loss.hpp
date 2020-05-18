#ifndef STAN_MATH_PRIM_FUN_BINARY_LOG_LOSS_HPP
#define STAN_MATH_PRIM_FUN_BINARY_LOG_LOSS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/meta/apply_scalar_binary.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the log loss function for binary classification
 * with specified reference and response values.
 *
 * The log loss function for prediction \f$\hat{y} \in [0, 1]\f$
 * given outcome \f$y \in \{ 0, 1 \}\f$ is
 *
 * \f$\mbox{logloss}(1, \hat{y}) = -\log \hat{y} \f$, and
 *
 * \f$\mbox{logloss}(0, \hat{y}) = -\log (1 - \hat{y}) \f$.
 *
 * @tparam T value type
 * @param[in] y reference value, either 0 or 1
 * @param[in] y_hat response value in [0, 1]
 * @return Log loss for response given reference value
 */
template <typename T1, typename T2, require_st_integral<T1>* = nullptr,
          require_st_arithmetic<T2>* = nullptr>
inline auto binary_log_loss(T1 y, const T2& y_hat) {
  return apply_scalar_binary<T1, T2>::apply(
      y, y_hat, [&](const auto& v, const auto& v_hat) {
        return v ? -log(v_hat) : -log1m(v_hat);
      });
}

}  // namespace math
}  // namespace stan

#endif
