#ifndef STAN_MATH_PRIM_META_FORWARD_AS_HPP
#define STAN_MATH_PRIM_META_FORWARD_AS_HPP

#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
namespace math {

namespace internal {
template <typename T1, typename T2>
constexpr bool eigen_static_size_match(T1 desired, T2 actual) {
  return static_cast<int>(desired) == Eigen::Dynamic
         || static_cast<int>(desired) == static_cast<int>(actual);
}
}  // namespace internal

/** \ingroup type_trait
 * Assume which type we get. If actual type is convertible to assumed type or in
 * case of eigen types compile time rows and columns also match or desired sizes
 * are dynamic this is a no-op. Otherwise it throws std::runtime_error, which
 * should never happen if used as intended.
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
          = std::enable_if_t<std::is_same<std::decay_t<T_actual>,
                                          std::decay_t<T_desired>>::value
                             && !is_eigen<T_desired>::value>>
inline T_actual&& forward_as(T_actual&& a) {  // NOLINT
  return std::forward<T_actual>(a);
}

/** \ingroup type_trait
 * Assume which type we get. If actual type is not convertible to assumed type
 * or in case of eigen types compile time rows and columns are not the same and
 * desired sizes are not dynamic this has return type of \c T_desired, but it
 * only throws. This version should only be used where it is optimized away so
 * the throw should never happen.
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
              !std::is_same<std::decay<T_actual>, std::decay<T_desired>>::value
              && (!is_eigen<T_desired>::value || !is_eigen<T_actual>::value)>>
inline T_desired forward_as(const T_actual& a) {
  throw std::runtime_error("Wrong type assumed! Please file a bug report.");
}

/** \ingroup type_trait
 * Assume which type we get. If actual type is convertible to assumed type or in
 * case of eigen types compile time rows and columns also match or desired sizes
 * are dynamic this is a no-op. Otherwise it throws std::runtime_error,
 * which should never happen if used as intended.
 *
 * This is intended to be used in compile time branches that would otherwise
 * trigger compile error even though they are never executed.
 *
 * @tparam T_desired type of output we need to avoid compile time errors
 * @tparam T_actual actual type of the argument
 * @param a input value
 * @return the input value a
 */
template <
    typename T_desired, typename T_actual,
    require_eigen_t<T_desired>* = nullptr,
    std::enable_if_t<
        std::is_same<value_type_t<T_actual>, value_type_t<T_desired>>::value
        && is_eigen<T_desired>::value && is_eigen<T_actual>::value
        && internal::eigen_static_size_match(
               T_desired::RowsAtCompileTime,
               std::decay_t<T_actual>::RowsAtCompileTime)
        && internal::eigen_static_size_match(
               T_desired::ColsAtCompileTime,
               std::decay_t<T_actual>::ColsAtCompileTime)>* = nullptr>
inline T_actual&& forward_as(T_actual&& a) {  // NOLINT
  return std::forward<T_actual>(a);
}

/**
 * Assume which type we get. If actual type is not convertible to assumed type
 * or in case of eigen types compile time rows and columns are not the same and
 * desired sizes are not dynamic this has return type of \c T_desired, but it
 * only throws. This version should only be used where it is optimized away so
 * the throw should never happen.
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
template <
    typename T_desired, typename T_actual,
    require_all_eigen_t<T_desired, T_actual>* = nullptr,
    std::enable_if_t<
        !std::is_same<value_type_t<T_actual>, value_type_t<T_desired>>::value
        || !internal::eigen_static_size_match(
               T_desired::RowsAtCompileTime,
               std::decay_t<T_actual>::RowsAtCompileTime)
        || !internal::eigen_static_size_match(
               T_desired::ColsAtCompileTime,
               std::decay_t<T_actual>::ColsAtCompileTime)>* = nullptr>
inline T_desired forward_as(const T_actual& a) {
  throw std::runtime_error("Wrong type assumed! Please file a bug report.");
}

}  // namespace math
}  // namespace stan
#endif
