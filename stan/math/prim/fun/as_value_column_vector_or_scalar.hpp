#ifndef STAN_MATH_PRIM_FUN_AS_VALUE_COLUMN_VECTOR_OR_SCALAR
#define STAN_MATH_PRIM_FUN_AS_VALUE_COLUMN_VECTOR_OR_SCALAR

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Extract values from input argument and transform to a column vector. For
 * scalar this returns a scalar. For arithmetic types this is an identity
 * function.
 *
 * @tparam T Type of scalar element.
 * @param a Specified scalar.
 * @return the scalar.
 */
template <typename T>
inline auto as_value_column_vector_or_scalar(T&& a) {
  return value_of(as_column_vector_or_scalar(std::forward<T>(a)));
}

}  // namespace math
}  // namespace stan

#endif
