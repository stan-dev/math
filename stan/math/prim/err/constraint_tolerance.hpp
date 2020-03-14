#ifndef STAN_MATH_PRIM_ERR_CONSTRAINT_TOLERANCE_HPP
#define STAN_MATH_PRIM_ERR_CONSTRAINT_TOLERANCE_HPP

#include <stan/math/prim/meta.hpp>
namespace stan {
namespace math {

/**
 * The tolerance for checking arithmetic bounds in rank and in
 * simplexes.  The default value is <code>1E-8</code>.
 */
const double CONSTRAINT_TOLERANCE = 1E-8;

}  // namespace math
}  // namespace stan
#endif
