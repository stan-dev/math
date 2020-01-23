#ifndef STAN_MATH_REV_FUN_PRIMITIVE_VALUE_HPP
#define STAN_MATH_REV_FUN_PRIMITIVE_VALUE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/primitive_value.hpp>

namespace stan {
namespace math {

/**
 * Return the primitive double value for the specified autodiff
 * variable.
 *
 * @param v input variable.
 * @return value of input.
 */
inline double primitive_value(const var& v) { return v.val(); }

}  // namespace math
}  // namespace stan
#endif
