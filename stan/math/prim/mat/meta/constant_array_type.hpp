#include <stan/math/prim/mat/meta/is_vector.hpp>

#ifndef STAN_MATH_PRIM_MAT_META_CONSTANT_ARRAY_TYPE_HPP
#define STAN_MATH_PRIM_MAT_META_CONSTANT_ARRAY_TYPE_HPP

namespace stan {
namespace math {
/**
 * constant_array_type provides a uniform wrapper around either a
 * scalar or an array of scalars.
 *
 * @tparam T the container type; will be the scalar type if wrapping a scalar
 */
  template<typename T, typename = void>
  struct constant_array_type {
      typedef typename partials_return_type<T>::type value_type;
  };

/**
 * constant_array_type provides a uniform wrapper around either a
 * scalar or an array of scalars.
 *
 * @tparam T the container type; will be the scalar type if wrapping a scalar
 */
  template<typename T>
  struct constant_array_type<T, typename std::enable_if<is_vector<T>::value>::type> {
      typedef Eigen::Array<typename partials_return_type<T>::type, Eigen::Dynamic, 1> value_type;
  };
}  // namespace math
}  // namespace stan
#endif
