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
inline auto logit(const var u) {
  auto denom = to_arena(1.0 / (u.val() - square(u.val())));
  return make_callback_var(logit(u.val()), [u, denom](auto& vi) mutable {
    u.adj() += vi.adj() * denom;
  });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto logit(const T& u) {
  auto u_arena = to_arena(u);
  auto denom
      = to_arena((u_arena.val() - square(u_arena.val())).array().inverse());
  return make_callback_rev_matrix<T>(
      logit(u_arena.val()), [u_arena, denom](auto&& vi) mutable {
        u_arena.adj().array() += vi.adj().array() * denom;
      });
}

}  // namespace math
}  // namespace stan
#endif
