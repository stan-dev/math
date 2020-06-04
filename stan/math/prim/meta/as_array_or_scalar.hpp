#ifndef STAN_MATH_PRIM_META_AS_ARRAY_OR_SCALAR_HPP
#define STAN_MATH_PRIM_META_AS_ARRAY_OR_SCALAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/holder.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
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
inline auto as_array_or_scalar(T&& v) {
  return make_holder([](auto& x) { return x.array(); }, std::forward<T>(v));
}

/** \ingroup type_trait
 * Converts a std::vector type to an array.
 *
 * @tparam T Type of scalar element.
 * @param v Specified vector.
 * @return Matrix converted to an array.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto as_array_or_scalar(T&& v) {
  using T_map
      = Eigen::Map<const Eigen::Array<value_type_t<T>, Eigen::Dynamic, 1>>;
  return make_holder([](auto& x) { return T_map(x.data(), x.size()); },
                     std::forward<T>(v));
}

}  // namespace math
}  // namespace stan

#endif
