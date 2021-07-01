#ifndef STAN_MATH_REV_FUN_LOGIT_HPP
#define STAN_MATH_REV_FUN_LOGIT_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/logit.hpp>

namespace stan {
namespace math {

/**
 * Return the log odds of the specified argument.
 *
 * @param u argument
 * @return log odds of argument
 */
inline var logit(const var& u) {
  auto denom = (1.0 / (u.val() - u.val() * u.val()));
  return make_callback_var(logit(u.val()), [u, denom](auto& vi) mutable {
    u.adj() += vi.adj() * denom;
  });
}

template <typename T, require_eigen_t<T>* = nullptr>
inline auto logit(const var_value<T>& u) {
  auto denom
      = to_arena(1.0 / (u.val().array() - u.val().array() * u.val().array()));
  return make_callback_var(logit(u.val()), [u, denom](auto& vi) mutable {
    u.adj().array() += vi.adj().array() * denom;
  });
}

}  // namespace math
}  // namespace stan
#endif
