#ifndef STAN_MATH_PRIM_META_AS_COLUMN_VECTOR_OR_SCALAR_HPP
#define STAN_MATH_PRIM_META_AS_COLUMN_VECTOR_OR_SCALAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/holder.hpp>
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
 * @return the scalar.
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline T as_column_vector_or_scalar(const T& a) {
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
template <typename T, require_eigen_col_vector_t<T>* = nullptr>
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
template <typename T, require_eigen_row_vector_t<T>* = nullptr,
          require_not_eigen_col_vector_t<T>* = nullptr>
inline auto as_column_vector_or_scalar(T&& a) {
  return make_holder([](auto& x) { return x.transpose(); }, std::forward<T>(a));
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
  using T_map = Eigen::Map<optionally_const_vector>;
  return make_holder([](auto& x) { return T_map(x.data(), x.size()); },
                     std::forward<T>(a));
}

}  // namespace math
}  // namespace stan

#endif
