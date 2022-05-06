#ifndef STAN_MATH_REV_FUN_ABS_HPP
#define STAN_MATH_REV_FUN_ABS_HPP

#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/fabs.hpp>
#include <stan/math/rev/fun/hypot.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the absolute value of the argument.
 *
 * @tparam T floating point or var_mat type
 * @param[in] x argument
 * @return absolute value of argument
 */
template <typename T>
inline auto abs(const var_value<T>& x) {
  return fabs(x);
}
  
/**
 * Return the absolute value of the complex argument.
 *
 * @param[in] z argument
 * @return absolute value of the argument
 */
inline auto abs(const std::complex<var>& z) { return internal::complex_abs(z); }

}  // namespace math
}  // namespace stan
#endif
