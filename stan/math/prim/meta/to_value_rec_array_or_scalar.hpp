#ifndef STAN_MATH_PRIM_META_TO_VALUE_REC_ARRAY_OR_SCALAR_HPP
#define STAN_MATH_PRIM_META_TO_VALUE_REC_ARRAY_OR_SCALAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/holder.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup type_trait
 * Accesses the value of a scalar.
 *
 * @tparam T Type of element.
 * @param v Specified value.
 * @return Same value.
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline auto to_value_rec_array_or_scalar(T v) {
  return value_of_rec(v);
}

/** \ingroup type_trait
 * Access the values of a matrix type and casts them to an array.
 *
 * @tparam T Type of \c Eigen \c Matrix or expression
 * @param v Specified \c Eigen \c Matrix or expression.
 * @return Matrix converted to an array.
 */
template <typename T, require_eigen_t<T>* = nullptr,
          require_not_plain_type_t<T>* = nullptr>
inline auto to_value_rec_array_or_scalar(T&& v) {
  return value_of_rec(v.eval()).eval().array();
}

/** \ingroup type_trait
 * Access the values of a matrix type and casts them to an array.
 *
 * @tparam T Type of \c Eigen \c Matrix or expression
 * @param v Specified \c Eigen \c Matrix or expression.
 * @return Matrix converted to an array.
 */
template <typename T, require_eigen_t<T>* = nullptr,
          require_plain_type_t<T>* = nullptr>
inline auto to_value_rec_array_or_scalar(T&& v) {
  return value_of_rec(std::forward<T>(v)).eval().array();
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
inline auto to_value_rec_array_or_scalar(T&& v) {
  return v.val().array();
}

/** \ingroup type_trait
 * Access the values of an std::vector type and casts them to an array.
 *
 * @tparam T Type of scalar element.
 * @param v Specified vector.
 * @return Matrix converted to an array.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto to_value_rec_array_or_scalar(T&& v) {
  using T_map
      = Eigen::Map<const Eigen::Array<value_type_t<T>, Eigen::Dynamic, 1>>;
  return value_of_rec(T_map(v.data(), v.size())).eval();
}

}  // namespace math
}  // namespace stan

#endif
