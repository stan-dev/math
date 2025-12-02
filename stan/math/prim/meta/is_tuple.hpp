#ifndef STAN_MATH_PRIM_META_IS_TUPLE_HPP
#define STAN_MATH_PRIM_META_IS_TUPLE_HPP

#include <cstddef>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <tuple>
#include <type_traits>

namespace stan {
namespace math {
namespace internal {

template <typename... Types>
struct is_tuple_impl : std::false_type {};

template <typename... Types>
struct is_tuple_impl<std::tuple<Types...>> : std::true_type {};
}  // namespace internal

template <typename T>
struct is_tuple : internal::is_tuple_impl<std::decay_t<T>> {};

template <typename T>
constexpr bool is_tuple_v = is_tuple<T>::value;
/*! \ingroup require_std */
/*! \defgroup tuple_types tuple  */
/*! \addtogroup tuple_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_tuple */
/*! @tparam T the type to check */
template <typename T>
using require_tuple_t = require_t<is_tuple<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_tuple */
/*! @tparam T the type to check */
template <typename T>
using require_not_tuple_t = require_not_t<is_tuple<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_tuple */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_tuple_t = require_all_t<is_tuple<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_tuple */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_tuple_t
    = require_all_not_t<is_tuple<std::decay_t<Types>>...>;
/*! @} */
}  // namespace math

template <typename T>
using is_tuple = math::is_tuple<T>;

template <typename T>
inline constexpr bool is_tuple_v = math::is_tuple_v<T>;

/**
 * Checks both that T is a tuple and that its size is N
 * @tparam T type to retrieve the element from
 * @tparam N expected size
 * @ingroup type_trait
 */
template <typename T, std::size_t N, typename = void>
struct is_tuple_of_size : std::false_type {};

template <typename T, std::size_t N>
struct is_tuple_of_size<T, N, std::enable_if_t<stan::is_tuple_v<T>>>
    : std::bool_constant<std::tuple_size_v<std::decay_t<T>> == N> {};

template <typename T, std::size_t N>
inline constexpr bool is_tuple_of_size_v = stan::is_tuple_of_size<T, N>::value;

}  // namespace stan

#endif
