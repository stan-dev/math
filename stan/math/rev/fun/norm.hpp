#ifndef STAN_MATH_REV_FUN_NORM_HPP
#define STAN_MATH_REV_FUN_NORM_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/square.hpp>
#include <stan/math/prim/fun/norm.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the squared magnitude of the complex argument.
 *
 * @param[in] z argument
 * @return squared magnitude of the argument
 */
inline var norm(const std::complex<var>& z) {
  return internal::complex_norm(z);
}

}  // namespace math
}  // namespace stan

#endif
