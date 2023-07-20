#ifndef STAN_MATH_PRIM_META_PLAIN_TYPE_HPP
#define STAN_MATH_PRIM_META_PLAIN_TYPE_HPP

#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_detected.hpp>
#include <stan/math/prim/meta/is_var_matrix.hpp>
#include <type_traits>

namespace stan {

/**
 * Determines plain (non expression) type associated with \c T. For non \c Eigen
 * types it is the decayed input type.
 * @tparam T type to determine plain type of
 */
template <typename T, typename Enable = void>
struct plain_type {
  using type = std::decay_t<T>;
};

template <typename T>
using plain_type_t = typename plain_type<T>::type;

/**
 * Determines return type of calling \c .eval() on Eigen expression.
 *
 * If input type \c T is a plain type (\c plain_type_t<T> equals \c
 * std::decay<T>), than member \c type is defined as <code> const
 * plain_type_t<T>& </code>. Otherwise member \c type is defined as \c
 * plain_type_t<T>.
 *
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

namespace internal {
// primary template handles types that have no nested ::type member:
template <class, class = void>
struct has_plain_object : std::false_type {};

// specialization recognizes types that do have a nested ::type member:
template <class T>
struct has_plain_object<T, void_t<typename std::decay_t<T>::PlainObject>>
    : std::true_type {};

// primary template handles types that have no nested ::type member:
template <class, class = void>
struct has_eval : std::false_type {};

// specialization recognizes types that do have a nested ::type member:
template <class T>
struct has_eval<T, void_t<decltype(std::declval<std::decay_t<T>&>().eval())>>
    : std::true_type {};

}  // namespace internal

/**
 * Determines plain (non expression) type associated with \c T. For \c Eigen
 * expression it is a type the expression can be evaluated into.
 * @tparam T type to determine plain type of
 */
template <typename T>
struct plain_type<T, require_t<bool_constant<internal::has_eval<T>::value
                                             && is_eigen<T>::value>>> {
  using type = std::decay_t<decltype(std::declval<T&>().eval())>;
};

template <typename T>
struct plain_type<
    T, require_t<bool_constant<!internal::has_eval<T>::value
                               && internal::has_plain_object<T>::value
                               && is_eigen<T>::value>>> {
  using type = typename std::decay_t<T>::PlainObject;
};

}  // namespace stan

#endif  // STAN_MATH_PRIM_META_PLAIN_TYPE_HPP
