#ifndef STAN_MATH_PRIM_META_TO_VALUE_COLUMN_ARRAY_OR_SCALAR_HPP
#define STAN_MATH_PRIM_META_TO_VALUE_COLUMN_ARRAY_OR_SCALAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/holder.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <vector>

namespace stan {
namespace math {


/** \ingroup type_trait
 * Converts input argument to a column vector or a scalar. For scalar inputs
 * that is an identity function.
 *
 * @tparam T Type of scalar element.
 * @param a Specified scalar.
 * @return the scalar.
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline auto to_value_column_array_or_scalar(const T& a) {
  return value_of(a);
}

/** \ingroup type_trait
 * Converts input argument to a column vector or a scalar. For column vector
 * inputs this is an identity function.
 *
 * @tparam T Type of scalar element.
 * @param a Specified vector.
 * @return Same vector.
 */
template <typename T, require_eigen_col_vector_t<T>* = nullptr>
inline auto to_value_column_array_or_scalar(T&& a) {
  return value_of(std::forward<T>(a).array()).eval();
}

/** \ingroup type_trait
 * Converts input argument to a column vector or a scalar. For a row vector
 * input this is transpose.
 *
 * @tparam T Type of scalar element.
 * @param a Specified vector.
 * @return Transposed vector.
 */
template <typename T, require_eigen_row_vector_t<T>* = nullptr,
          require_not_eigen_col_vector_t<T>* = nullptr>
inline auto to_value_column_array_or_scalar(T&& a) {
  return value_of(std::forward<T>(a)).transpose().eval().array();
}

/** \ingroup type_trait
 * Accesses the values of a `var_value` with inner matrix type to
 * and casts them to an array.
 *
 * @tparam T Type of \c Eigen \c Matrix or expression
 * @param v Specified \c Eigen \c Matrix or expression.
 * @return Matrix converted to an array.
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto to_value_column_array_or_scalar(T&& v) {
  return v.val().array();
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
inline auto to_value_column_array_or_scalar(T&& a) {
  using plain_vector = Eigen::Array<value_type_t<T>, Eigen::Dynamic, 1>;
  using optionally_const_vector
      = std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                           const plain_vector, plain_vector>;
  using T_map = Eigen::Map<optionally_const_vector>;
  return value_of(T_map(a.data(), a.size())).eval();
}

}  // namespace math
}  // namespace stan

#endif
