#ifndef STAN_MATH_REV_FUN_LOGIT_HPP
#define STAN_MATH_REV_FUN_LOGIT_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/logit.hpp>
#include <stan/math/prim/fun/square.hpp>

namespace stan {
namespace math {

/**
 * Return the log odds of the specified argument.
 *
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param u The variable.
 * @return log odds of argument
 */
template <typename T, require_stan_scalar_or_eigen_t<T>* = nullptr>
inline auto logit(const var_value<T>& u) {
  auto denom = to_arena(1.0 / as_array_or_scalar(u.val() - square(u.val())));
  return make_callback_var(logit(u.val()), [u, denom](auto& vi) mutable {
    as_array_or_scalar(u.adj()) += as_array_or_scalar(vi.adj()) * denom;
  });
}

}  // namespace math
}  // namespace stan
#endif
