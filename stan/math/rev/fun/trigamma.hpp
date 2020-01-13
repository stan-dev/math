#ifndef STAN_MATH_REV_FUN_TRIGAMMA_HPP
#define STAN_MATH_REV_FUN_TRIGAMMA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/floor.hpp>
#include <stan/math/rev/fun/sin.hpp>
#include <stan/math/rev/fun/square.hpp>
#include <stan/math/rev/fun/inv_square.hpp>
#include <stan/math/prim/fun/trigamma.hpp>

namespace stan {
namespace math {

/**
 * Return the value of the trigamma function at the specified
 * argument (i.e., the second derivative of the log Gamma function
 * at the specified argument).
 *
 * @param u argument
 * @return trigamma function at argument
 */
inline var trigamma(const var& u) { return trigamma_impl(u); }

}  // namespace math
}  // namespace stan
#endif
