#ifndef STAN_MATH_PRIM_MAT_FUN_AS_COLUMN_VECTOR_OR_SCALAR_HPP
#define STAN_MATH_PRIM_MAT_FUN_AS_COLUMN_VECTOR_OR_SCALAR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup type_trait
 * Converts input argument to a column vector or a scalar. For column vector
 * inputs this is an identity function.
 *
 * @tparam T Type of scalar element.
 * @param a Specified vector.
 * @return Same vector.
 */
template <typename T, typename = require_t<is_eigen_col_vector<T>>>
inline const auto& as_column_vector_or_scalar(const T& a) {
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
template <typename T, typename = require_t<is_eigen_row_vector<T>>>
inline auto as_column_vector_or_scalar(const T& a) {
  return a.transpose();
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
