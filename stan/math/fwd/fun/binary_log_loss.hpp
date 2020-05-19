#ifndef STAN_MATH_FWD_FUN_BINARY_LOG_LOSS_HPP
#define STAN_MATH_FWD_FUN_BINARY_LOG_LOSS_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/binary_log_loss.hpp>

namespace stan {
namespace math {

template <typename T1, typename T2, require_st_integral<T1>* = nullptr,
          require_st_fvar<T2>* = nullptr>
inline auto binary_log_loss(T1 y, const T2& y_hat) {
  return apply_scalar_binary<T1, T2>::apply(
      y, y_hat, [&](const auto& v, const auto& v_hat) {
        using T = typename scalar_type_t<T2>::Scalar;
        if (v) {
          return fvar<T>(binary_log_loss(v, v_hat.val_),
                         -v_hat.d_ / v_hat.val_);
        } else {
          return fvar<T>(binary_log_loss(v, v_hat.val_),
                         v_hat.d_ / (1.0 - v_hat.val_));
        }
      });
}
}  // namespace math
}  // namespace stan
#endif
