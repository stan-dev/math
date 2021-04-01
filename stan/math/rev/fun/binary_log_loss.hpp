#ifndef STAN_MATH_REV_FUN_BINARY_LOG_LOSS_HPP
#define STAN_MATH_REV_FUN_BINARY_LOG_LOSS_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The log loss function for variables (stan).
 *
 * See binary_log_loss() for the double-based version.
 *
 * The derivative with respect to the variable \f$\hat{y}\f$ is
 *
 * \f$\frac{d}{d\hat{y}} \mbox{logloss}(1, \hat{y}) = - \frac{1}{\hat{y}}\f$,
 and
 *
 * \f$\frac{d}{d\hat{y}} \mbox{logloss}(0, \hat{y}) = \frac{1}{1 - \hat{y}}\f$.
 *
 *
   \f[
   \mbox{binary\_log\_loss}(y, \hat{y}) =
   \begin{cases}
     y \log \hat{y} + (1 - y) \log (1 - \hat{y}) & \mbox{if } 0\leq \hat{y}\leq
 1, y\in\{ 0, 1 \}\\[6pt] \textrm{NaN} & \mbox{if } \hat{y} = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{binary\_log\_loss}(y, \hat{y})}{\partial \hat{y}} =
   \begin{cases}
     \frac{y}{\hat{y}}-\frac{1-y}{1-\hat{y}} & \mbox{if } 0\leq \hat{y}\leq 1,
     y\in\{ 0, 1 \}\\[6pt]
     \textrm{NaN} & \mbox{if } \hat{y} = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param y Reference value.
 * @param y_hat Response variable.
 * @return Log loss of response versus reference value.
 */
inline var binary_log_loss(int y, const var& y_hat) {
  if (y == 0) {
    return make_callback_var(-log1p(-y_hat.val()), [y_hat](auto& vi) mutable {
      y_hat.adj() += vi.adj() / (1.0 - y_hat.val());
    });
  } else {
    return make_callback_var(-std::log(y_hat.val()), [y_hat](auto& vi) mutable {
      y_hat.adj() -= vi.adj() / y_hat.val();
    });
  }
}

}  // namespace math
}  // namespace stan
#endif
