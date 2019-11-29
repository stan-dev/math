#ifndef STAN_MATH_REV_CORE_STD_ISINF_HPP
#define STAN_MATH_REV_CORE_STD_ISINF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/is_inf.hpp>
#include <stan/math/rev/core/var.hpp>

namespace std {

/**
 * Return <code>1</code> if the specified argument is positive
 * infinity or negative infinity and <code>0</code> otherwise.
 *
 * @param a Argument.
 * @return 1 if argument is infinite and 0 otherwise.
 */
template <typename T, stan::require_var_t<T>...>
inline auto isinf(T&& a) {
  return stan::math::is_inf(a.val());
}

}  // namespace std
#endif
