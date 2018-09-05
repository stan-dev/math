#include <stan/math/prim/mat/meta/is_vector.hpp>

#ifndef STAN_MATH_PRIM_MAT_META_CONSTANT_ARRAY_TYPE_HPP
#define STAN_MATH_PRIM_MAT_META_CONSTANT_ARRAY_TYPE_HPP

namespace stan {
namespace math {
/**
 * constant_array_type returns either an array or a scalar type.
 *
 * @tparam T the type that is used to determine the shape of the return type:
 * a scalar type if T is a scalar type and an array type if T is an array type
 * @tparam S the type of the entries of the return type
 */
  template<typename T, typename S, typename = void>
  struct constant_array_type {
      typedef S value_type;
  };

/**
 * constant_array_type returns either an array or a scalar type.
 *
 * @tparam T the type that is used to determine the shape of the return type:
 * a scalar type if T is a scalar type and an array type if T is an array type
 * @tparam S the type of the entries of the return type
 */
  template<typename T, typename S>
  struct constant_array_type<T, S, typename std::enable_if<is_vector<T>::value>::type> {
      typedef Eigen::Array<S, Eigen::Dynamic, 1> value_type;
  };
}  // namespace math
}  // namespace stan
#endif
