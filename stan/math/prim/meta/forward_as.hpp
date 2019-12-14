#ifndef STAN_MATH_PRIM_META_FORWARD_AS_HPP
#define STAN_MATH_PRIM_META_FORWARD_AS_HPP

#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
namespace math {

/** \ingroup type_trait
 * Assume which type we get. If actual type is convertible to assumed type or in
 * case of eigen types compile time rows and columns also match this is a no-op.
 * Otherwise it throws std::runtime_error, which should never happen if used as
 * intended.
 *
 * This is intended to be used in compile time branches that would otherwise
 * trigger compile error even though they are never executed.
 *
 * @tparam T_desired type of output we need to avoid compile time errors
 * @tparam T_actual actual type of the argument
 * @param a input value
 * @return the input value a
 */
template <typename T_desired, typename T_actual,
          typename
          = std::enable_if_t<std::is_convertible<T_actual, T_desired>::value
                             && !is_eigen<T_desired>::value>>
inline T_actual&& forward_as(T_actual&& a) {  // NOLINT
  return std::forward<T_actual>(a);
}

/** \ingroup type_trait
 * Assume which type we get. If actual type is not convertible to assumed type
 * or in case of eigen types compile time rows and columns are not the same this
 * has return type of \c T_desired, but it only throws. This version should only
 * be used where it is optimized away so the throw should never happen.
 *
 * This is intended to be used in compile time branches that would otherwise
 * trigger compile error even though they are never executed.
 *
 * @tparam T_desired type of output we need to avoid compile time errors
 * @tparam T_actual actual type of the argument
 * @param a input value
 * @return nothing, this always throws
 * @throw always throws std::runtime_error
 */
template <typename T_desired, typename T_actual,
          typename
          = std::enable_if_t<!std::is_convertible<T_actual, T_desired>::value>>
inline T_desired forward_as(const T_actual& a) {
  throw std::runtime_error("Wrong type assumed! Please file a bug report.");
}

/** \ingroup type_trait
 * Assume which type we get. If actual type is convertible to assumed type or in
 * case of eigen types compile time rows and columns also match this is a no-op.
 * Otherwise it throws std::runtime_error, which should never happen if used as
 * intended.
 *
 * This is intended to be used in compile time branches that would otherwise
 * trigger compile error even though they are never executed.
 *
 * @tparam T_desired type of output we need to avoid compile time errors
 * @tparam T_actual actual type of the argument
 * @param a input value
 * @return the input value a
 */
// clang-format off
template <typename T_desired, typename T_actual,
          typename = std::enable_if_t<
              std::is_convertible<T_actual, T_desired>::value && static_cast<
                  int>(T_desired::RowsAtCompileTime)
                  == static_cast<int>(T_actual::RowsAtCompileTime)
              && static_cast<int>(T_desired::ColsAtCompileTime)
                     == static_cast<int>(T_actual::ColsAtCompileTime)>,
          typename = void>
// clang-format on
inline T_actual&& forward_as(T_actual&& a) {  // NOLINT
  return std::forward<T_actual>(a);
}

/**
 * Assume which type we get. If actual type is not convertible to assumed type
 * or in case of eigen types compile time rows and columns are not the same this
 * has return type of \c T_desired, but it only throws. This version should only
 * be used where it is optimized away so the throw should never happen.
 *
 * This is intended to be used in compile time branches that would otherwise
 * trigger compile error even though they are never executed.
 *
 * @tparam T_desired type of output we need to avoid compile time errors
 * @tparam T_actual actual type of the argument
 * @param a input value
 * @return nothing, this always throws
 * @throw always throws std::runtime_error
 */
template <typename T_desired, typename T_actual,
          typename = std::enable_if_t<
              !std::is_convertible<T_actual, T_desired>::value
              || static_cast<int>(T_desired::RowsAtCompileTime)
                     != static_cast<int>(T_actual::RowsAtCompileTime)
              || static_cast<int>(T_desired::ColsAtCompileTime)
                     != static_cast<int>(T_actual::ColsAtCompileTime)>,
          typename = void>
inline T_desired forward_as(const T_actual& a) {
  throw std::runtime_error("Wrong type assumed! Please file a bug report.");
}

}  // namespace math
}  // namespace stan
#endif
