#ifndef STAN_MATH_REV_FUN_FALLING_FACTORIAL_HPP
#define STAN_MATH_REV_FUN_FALLING_FACTORIAL_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/falling_factorial.hpp>

namespace stan {
namespace math {

inline var falling_factorial(const var& a, int b) {
  return make_callback_var(
      falling_factorial(a.val(), b), [a, b](auto& vi) mutable {
        a.adj() += vi.adj() * vi.val()
                   * (digamma(a.val() + 1) - digamma(a.val() - b + 1));
      });
}

}  // namespace math
}  // namespace stan
#endif
