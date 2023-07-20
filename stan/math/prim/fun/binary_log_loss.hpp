#ifndef STAN_MATH_PRIM_FUN_BINARY_LOG_LOSS_HPP
#define STAN_MATH_PRIM_FUN_BINARY_LOG_LOSS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
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
template <typename T, require_arithmetic_t<T>* = nullptr>
inline T binary_log_loss(int y, const T& y_hat) {
  using std::log;
  return y ? -log(y_hat) : -log1m(y_hat);
}

/**
 * Enables the vectorized application of the binary log loss function, when
 * the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return Binary log loss function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr,
          require_not_var_matrix_t<T2>* = nullptr>
inline auto binary_log_loss(const T1& a, const T2& b) {
  return apply_scalar_binary(a, b, [&](const auto& c, const auto& d) {
    return binary_log_loss(c, d);
  });
}

}  // namespace math
}  // namespace stan

#endif
