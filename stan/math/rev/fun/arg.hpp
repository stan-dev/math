#ifndef STAN_MATH_REV_FUN_ARG_HPP
#define STAN_MATH_REV_FUN_ARG_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/atan2.hpp>
#include <stan/math/prim/fun/arg.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the phase angle of the complex argument.
 *
 * @param[in] z argument
 * @return phase angle of the argument
 */
inline var arg(const std::complex<var>& z) { return internal::complex_arg(z); }

}  // namespace math
}  // namespace stan

#endif
