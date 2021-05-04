#ifndef STAN_MATH_FWD_FUN_PROJ_HPP
#define STAN_MATH_FWD_FUN_PROJ_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/proj.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the projection of the complex argument onto the Riemann
 * sphere.
 *
 * @tparam T value type of autodiff variable
 * @param[in] z argument
 * @return projection of the argument onto the Riemann sphere
 */
template <typename T>
inline std::complex<fvar<T>> proj(const std::complex<fvar<T>>& z) {
  return internal::complex_proj(z);
}

}  // namespace math
}  // namespace stan

#endif
