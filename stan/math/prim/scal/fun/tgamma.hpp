#ifndef STAN_MATH_PRIM_SCAL_FUN_TGAMMA_HPP
#define STAN_MATH_PRIM_SCAL_FUN_TGAMMA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/scal/fun/is_nonpositive_integer.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the gamma function applied to the specified argument.
 *
 * @param x Argument.
 * @return The gamma function applied to argument.
 */
inline double tgamma(double x) {
  if (x == 0.0 || is_nonpositive_integer(x)) {
    throw_domain_error("tgamma", "x", x, "x == 0 or negative integer");
  }
  return std::tgamma(x);
}

}  // namespace math
}  // namespace stan
#endif
