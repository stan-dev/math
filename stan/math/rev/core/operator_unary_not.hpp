#ifndef STAN_MATH_REV_CORE_OPERATOR_UNARY_NOT_HPP
#define STAN_MATH_REV_CORE_OPERATOR_UNARY_NOT_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return the negation of the value of the argument as defined by
 * <code>!</code>.
 *
 * @param[in] x argument
 * @return negation of argument value
 */
template <typename T>
inline bool operator!(var_value<T> x) { return !x.val(); }

}  // namespace math
}  // namespace stan
#endif
