#ifndef STAN_MATH_PRIM_SCAL_META_PLAIN_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_META_PLAIN_TYPE_HPP

#include <type_traits>

namespace stan {

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
 * Determines return type of calling .eval() on eigen expression. This is the
 * same as \c plain_type, except if input is already a plain type. In such case
 * eval return type is a const reference to plain type.
 * @tparam T type to determine eval return type of
 */
template <typename T>
struct eval_return_type {
  using T1 = plain_type_t<T>;
  using type = std::conditional_t<std::is_same<std::decay_t<T>, T1>::value,
                                  const T1&, T1>;
};

template <typename T>
using eval_return_type_t = typename eval_return_type<T>::type;

}  // namespace stan

#endif  // STAN_MATH_PRIM_SCAL_META_PLAIN_TYPE_HPP
