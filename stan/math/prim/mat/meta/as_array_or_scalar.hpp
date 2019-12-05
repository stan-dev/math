#ifndef STAN_MATH_PRIM_MAT_FUN_AS_ARRAY_OR_SCALAR_HPP
#define STAN_MATH_PRIM_MAT_FUN_AS_ARRAY_OR_SCALAR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup type_trait
 * Converts a matrix type to an array.
 *
 * @tparam T Type of scalar element.
 * @tparam R Row type of input matrix.
 * @tparam C Column type of input matrix.
 * @param v Specified matrix.
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
