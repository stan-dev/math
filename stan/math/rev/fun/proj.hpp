#ifndef STAN_MATH_REV_FUN_PROJ_HPP
#define STAN_MATH_REV_FUN_PROJ_HPP

#include <stan/math/prim/fun/proj.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/is_inf.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the projection of the complex argument onto the Riemann
 * sphere.
 *
 * @param[in] z argument
 * @return projection of the argument onto the Riemann sphere
 */
inline std::complex<var> proj(const std::complex<var>& z) {
  return internal::complex_proj(z);
}

}  // namespace math
}  // namespace stan

#endif
