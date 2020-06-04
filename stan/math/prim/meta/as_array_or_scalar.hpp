#ifndef STAN_MATH_PRIM_META_AS_ARRAY_OR_SCALAR_HPP
#define STAN_MATH_PRIM_META_AS_ARRAY_OR_SCALAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup type_trait
 * Returns specified input value.
 *
 * @tparam T Type of element.
 * @param v Specified value.
 * @return Same value.
 */
template <typename T, typename = require_stan_scalar_t<T>>
inline T as_array_or_scalar(T&& v) {
  return std::forward<T>(v);
}

/** \ingroup type_trait
 * Converts a matrix type to an array.
 *
 * @tparam T Type of \c Eigen \c Matrix or expression
 * @param v Specified \c Eigen \c Matrix or expression.
 * @return Matrix converted to an array.
 */
template <typename T, typename = require_eigen_t<T>>
inline auto as_array_or_scalar(const T& v) {
  return v.array();
}

/** \ingroup type_trait
 * Converts a std::vector type to an array.
 *
 * @tparam T Type of scalar element.
 * @param v Specified vector.
 * @return Matrix converted to an array.
 */
template <typename T>
inline Eigen::Map<const Eigen::Array<T, Eigen::Dynamic, 1>> as_array_or_scalar(
    const std::vector<T>& v) {
  return Eigen::Map<const Eigen::Array<T, Eigen::Dynamic, 1>>(v.data(),
                                                              v.size());
}

}  // namespace math
}  // namespace stan

#endif
