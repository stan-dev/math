#ifndef STAN_MATH_REV_FUN_LGAMMA_HPP
#define STAN_MATH_REV_FUN_LGAMMA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>

namespace stan {
namespace math {

/**
 * The log gamma function for variables (C99).
 *
 * The derivative is the digamma function,
 *
 * \f$\frac{d}{dx} \Gamma(x) = \psi^{(0)}(x)\f$.
 *
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param a The variable.
 * @return Log gamma of the variable.
 */
inline auto lgamma(const var a) {
  return make_callback_var(lgamma(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() * digamma(a.val());
  });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto lgamma(const T& x) {
  auto x_arena = to_arena(x);
  return make_callback_rev_matrix<T>(lgamma(x_arena.val()), [x_arena](auto& vi) mutable {
    x_arena.adj().array() += vi.adj().array() * digamma(x_arena.val()).array();
  });
}

}  // namespace math
}  // namespace stan
#endif
