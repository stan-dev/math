#ifndef STAN_MATH_REV_FUN_LOG1M_HPP
#define STAN_MATH_REV_FUN_LOG1M_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/log1m.hpp>

namespace stan {
namespace math {

/**
 * The log (1 - x) function for variables.
 *
 * The derivative is given by
 *
 * \f$\frac{d}{dx} \log (1 - x) = -\frac{1}{1 - x}\f$.
 *
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param a The variable.
 * @return The variable representing log of 1 minus the variable.
 */
template <typename T, require_stan_scalar_or_eigen_t<T>* = nullptr>
inline auto log1m(const var_value<T>& a) {
  return make_callback_var(log1m(a.val()), [a](auto& vi) mutable {
    as_array_or_scalar(a.adj())
        += as_array_or_scalar(vi.adj()) / (as_array_or_scalar(a.val()) - 1.0);
  });
}

}  // namespace math
}  // namespace stan
#endif
