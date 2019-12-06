#ifndef STAN_MATH_PRIM_SCAL_META_PLAIN_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_META_PLAIN_TYPE_HPP

#include <type_traits>

namespace stan {
namespace math {

/**
 * Determines plain (non expression) type associated with \c T. For non \c Eigen
 * types it is the input type.
 * @tparam T type to determine plain type of
 */
template <typename T, typename Enable = void>
struct plain_type {
  using type = std::decay_t<T>;
};

template <typename T>
using plain_type_t = typename plain_type<T>::type;

/**
 * Determines plain (non expression) type associated with \c T. For non \c Eigen
 * types it is the input type. Result will have same reference type as input.
 * @tparam T type to determine plain type of
 */
template <typename T>
struct plain_type_keep_constness_and_ref {
  using T1 = plain_type_t<T>;
  using T2
      = std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                           const T1, T1>;
  using T3 = std::conditional_t<std::is_lvalue_reference<T>::value, T2&, T2>;
  using type = std::conditional_t<std::is_rvalue_reference<T>::value, T3&&, T3>;
};

template <typename T>
using plain_type_keep_constness_and_ref_t =
    typename plain_type_keep_constness_and_ref<T>::type;

}  // namespace math
}  // namespace stan

#endif  // STAN_MATH_PRIM_SCAL_META_PLAIN_TYPE_HPP
