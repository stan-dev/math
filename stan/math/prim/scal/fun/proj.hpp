#ifndef STAN_MATH_PRIM_SCAL_FUN_PROJ_HPP
#define STAN_MATH_PRIM_SCAL_FUN_PROJ_HPP

#include <stan/math/prim/scal/fun/copysign.hpp>
#include <stan/math/prim/scal/fun/is_inf.hpp>
#include <complex>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return a projection of a complex number onto
 * the Riemann sphere.
 *
 * This overloads an erroneous definition in
 * libstdc++ for general types.
 *
 * @tparam T complex type
 * @param[in] t complex variable input
 * @return projection onto Riemann sphere
 */

template <class T>
inline std::complex<T> proj(std::complex<T> const& t) {
  if (is_inf(t.real()) || is_inf(t.imag()))
    return std::complex<T>(INFINITY, copysign(T(0), t.imag()));
  return t;
}

}  // namespace math
}  // namespace stan
#endif
