#ifndef STAN_MATH_PRIM_META_AS_COLUMN_VECTOR_OR_SCALAR_HPP
#define STAN_MATH_PRIM_META_AS_COLUMN_VECTOR_OR_SCALAR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup type_trait
 * Converts input argument to a column vector or a scalar. For scalar inputs
 * that is an identity function.
 *
 * @tparam T Type of scalar element.
 * @param a Specified scalar.
 * @return 1x1 matrix that contains the value of scalar.
 */
template <typename T>
inline const T& as_column_vector_or_scalar(const T& a) {
  return a;
}

/** \ingroup type_trait
 * Converts input argument to a column vector or a scalar. For column vector
 * inputs this is an identity function.
 *
 * @tparam T Type of scalar element.
 * @param a Specified vector.
 * @return Same vector.
 */
template <typename T>
inline const Eigen::Matrix<T, Eigen::Dynamic, 1>& as_column_vector_or_scalar(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& a) {
  return a;
}

/** \ingroup type_trait
 * Converts input argument to a column vector or a scalar. For a row vector
 * input this is transpose.
 *
 * @tparam T Type of scalar element.
 * @param a Specified vector.
 * @return Transposed vector.
 */
template <typename T>
inline Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>
as_column_vector_or_scalar(const Eigen::Matrix<T, 1, Eigen::Dynamic>& a) {
  return Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>(
      a.data(), a.size());  // uses Eigen::Map instead of .transpose() so that
                            // there are less possible output types
}

/** \ingroup type_trait
 * Converts input argument to a column vector or a scalar. std::vector will be
 * converted to a column vector.
 *
 * @tparam T Type of scalar element.
 * @param a Specified vector.
 * @return intut converted to a column vector.
 */
template <typename T>
inline Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>
as_column_vector_or_scalar(const std::vector<T>& a) {
  return Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>(a.data(),
                                                               a.size());
}

}  // namespace math
}  // namespace stan

#endif
