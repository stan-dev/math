#ifndef STAN_MATH_PRIM_FUN_AS_VALUE_COLUMN_ARRAY_OR_SCALAR
#define STAN_MATH_PRIM_FUN_AS_VALUE_COLUMN_ARRAY_OR_SCALAR

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Extract the value from an object and for eigen vectors and `std::vectors`
 * convert to an eigen column array and for scalars return a scalar.
 * @tparam T A stan scalar, eigen vector, or `std::vector`.
 * @param a Specified scalar.
 * @return the scalar.
 */
template <typename T>
inline auto as_value_column_array_or_scalar(T&& a) {
  return value_of(
      as_array_or_scalar(as_column_vector_or_scalar(std::forward<T>(a))));
}

}  // namespace math
}  // namespace stan

#endif
