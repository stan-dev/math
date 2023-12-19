#ifndef STAN_MATH_PRIM_ERR_CONSTRAINT_TOLERANCE_HPP
#define STAN_MATH_PRIM_ERR_CONSTRAINT_TOLERANCE_HPP

#ifndef STAN_MATH_CONSTRAINT_TOLERANCE
#define STAN_MATH_CONSTRAINT_TOLERANCE 1E-8
#endif

#include <stan/math/prim/meta.hpp>
namespace stan {
namespace math {

/**
 * The tolerance for checking arithmetic bounds in rank and in
 * simplexes.  The default value is <code>1E-8</code>.
 * This can changed by defining <code>STAN_MATH_CONSTRAINT_TOLERANCE</code>
 * at compile time.
 */
const double CONSTRAINT_TOLERANCE = STAN_MATH_CONSTRAINT_TOLERANCE;

}  // namespace math
}  // namespace stan
#endif
