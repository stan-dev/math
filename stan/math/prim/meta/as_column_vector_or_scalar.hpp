#ifndef STAN_MATH_PRIM_META_AS_COLUMN_VECTOR_OR_SCALAR_HPP
#define STAN_MATH_PRIM_META_AS_COLUMN_VECTOR_OR_SCALAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
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
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline T&& as_column_vector_or_scalar(T&& a) {
  return std::forward<T>(a);
}

/** \ingroup type_trait
 * Converts input argument to a column vector or a scalar. For column vector
 * inputs this is an identity function.
 *
 * @tparam T Type of scalar element.
 * @param a Specified vector.
 * @return Same vector.
 */
template <typename T, require_t<is_eigen_col_vector<T>>* = nullptr>
inline T&& as_column_vector_or_scalar(T&& a) {
  return std::forward<T>(a);
}

/** \ingroup type_trait
 * Converts input argument to a column vector or a scalar. For a row vector
 * input this is transpose.
 *
 * @tparam T Type of scalar element.
 * @param a Specified vector.
 * @return Transposed vector.
 */
template <typename T, require_t<is_eigen_row_vector<T>>* = nullptr>
inline auto as_column_vector_or_scalar(T&& a) {
  return a.transpose();
}

/** \ingroup type_trait
 * Converts input argument to a column vector or a scalar. std::vector will be
 * converted to a column vector.
 *
 * @tparam T Type of scalar element.
 * @param a Specified vector.
 * @return input converted to a column vector.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto as_column_vector_or_scalar(T&& a) {
  using plain_vector = Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, 1>;
  using optionally_const_vector
      = std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                           const plain_vector, plain_vector>;
  return Eigen::Map<optionally_const_vector>(a.data(), a.size());
}

}  // namespace math
}  // namespace stan

#endif
