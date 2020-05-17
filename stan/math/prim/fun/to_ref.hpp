#ifndef STAN_MATH_PRIM_FUN_TO_REF_HPP
#define STAN_MATH_PRIM_FUN_TO_REF_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * No-op that should be optimized away.
 * @tparam T non-Eigen argument type
 * @param a argument
 * @return argument
 */
template <typename T, require_not_eigen_t<T>* = nullptr>
inline T to_ref(T&& a) {
  return std::forward<T>(a);
}

/**
 * Converts Eigen argument into `Eigen::Ref`. This evaluate expensive
 * expressions.
 * @tparam T argument type (Eigen expression)
 * @param a argument
 * @return argument converted to `Eigen::Ref`
 */
template <typename T, require_eigen_t<T>* = nullptr>
inline Eigen::Ref<const plain_type_t<T>> to_ref(T&& a) {
  return std::forward<T>(a);
}

}  // namespace math
}  // namespace stan
#endif
