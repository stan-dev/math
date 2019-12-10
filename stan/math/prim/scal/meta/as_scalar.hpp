#ifndef STAN_MATH_PRIM_SCAL_FUN_AS_SCALAR_HPP
#define STAN_MATH_PRIM_SCAL_FUN_AS_SCALAR_HPP

#include <vector>
#include <stdexcept>

namespace stan {
namespace math {

/** \ingroup type_trait
 * Converts input to a scalar. For scalar arguments this is an identity
 * function.
 * @param a Input value
 * @return Same value
 */
inline double as_scalar(double a) { return a; }

/** \ingroup type_trait
 * Converts input to a scalar. For scalar arguments this is an identity
 * function.
 * @param a Input value
 * @return Same value
 */
inline int as_scalar(int a) { return a; }

/** \ingroup type_trait
 * Converts input to a scalar. As this is not possible for vectors it always
 * throws. This is intended to never be called, only used in templated functions
 * in branches that will be optimized out - to prevent compiler from complaining
 * about expressions with incompatible types.
 * @param a Input expression
 * @throws runtime_error Always throws
 * @return Never returns
 */
template <typename T>
inline double as_scalar(const std::vector<T>& a) {
  throw std::runtime_error("A vector can not be used as a scalar!");
}

}  // namespace math
}  // namespace stan

#endif
