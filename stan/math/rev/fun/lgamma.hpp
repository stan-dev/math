#ifndef STAN_MATH_REV_FUN_LGAMMA_HPP
#define STAN_MATH_REV_FUN_LGAMMA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/digamma.hpp>
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
template <typename T, require_stan_scalar_or_eigen_t<T>* = nullptr>
inline auto lgamma(const var_value<T>& a) {
  return make_callback_var(lgamma(a.val()), [a](auto& vi) mutable {
    as_array_or_scalar(a.adj())
        += as_array_or_scalar(vi.adj()) * as_array_or_scalar(digamma(a.val()));
  });
}

}  // namespace math
}  // namespace stan
#endif
