#ifndef STAN_MATH_PRIM_SCAL_META_ASSUME_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_META_ASSUME_TYPE_HPP

#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Assume which type we get. If actual type is convertible to assumed type or in
 * case of eigen types compile time rows and columns also match this is a no-op.
 * This is intended to be used in compile time branches that would otherwise
 * trigger compile error even though they are optimized away.
 * @tparam T_desired type of output we need to avoid compile time errors
 * @tparam T_actual actual type of the argument
 * @param a input value
 * @return the input value a
 */
template <typename T_desired, typename T_actual,
          typename
          = std::enable_if_t<std::is_convertible<T_actual, T_desired>::value
                             && !is_eigen<T_desired>::value>>
inline T_actual&& assume_type(T_actual&& a) {  // NOLINT
  return std::forward<T_actual>(a);
}

/**
 * Assume which type we get. If actual type is not convertible to assumed type
 * or in case of eigen types compile time rows and columns are not the same this
 * has return type of \c T_desired, but it only throws. This version should only
 * be used where it is optimized away so the throw should never happen. This is
 * intended to be used in compile time branches that would otherwise trigger
 * compile error even though they are optimized away.
 * @tparam T_desired type of output we need to avoid compile time errors
 * @tparam T_actual actual type of the argument
 * @param a input value
 * @return nothing, this always throws
 * @throw always throws std::runtime_error
 */
template <typename T_desired, typename T_actual,
          typename
          = std::enable_if_t<!std::is_convertible<T_actual, T_desired>::value>>
inline T_desired assume_type(const T_actual& a) {
  throw std::runtime_error("Wrong type assumed! Please file a bug report.");
}

}  // namespace math
}  // namespace stan

#endif
