#ifndef STAN_MATH_REV_FUN_FALLING_FACTORIAL_HPP
#define STAN_MATH_REV_FUN_FALLING_FACTORIAL_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/digamma.hpp>
#include <stan/math/prim/fun/falling_factorial.hpp>

namespace stan {
namespace math {

inline var falling_factorial(const var& a, int b) {
  auto digamma_ab = digamma(a.val() + 1) - digamma(a.val() - b + 1);
  return make_callback_var(falling_factorial(a.val(), b),
                           [a, digamma_ab](auto& vi) mutable {
                             a.adj() += vi.adj() * vi.val() * digamma_ab;
                           });
}

template <typename T1, typename T2, require_eigen_t<T1>* = nullptr,
          require_st_integral<T2>* = nullptr>
inline auto falling_factorial(const var_value<T1>& a, const T2& b) {
  auto b_map = as_array_or_scalar(b);
  auto digamma_ab = to_arena(digamma(a.val().array() + 1)
                             - digamma(a.val().array() - b_map + 1));
  return make_callback_var(
      falling_factorial(a.val(), b), [a, digamma_ab](auto& vi) mutable {
        a.adj().array() += vi.adj().array() * vi.val().array() * digamma_ab;
      });
}

}  // namespace math
}  // namespace stan
#endif
