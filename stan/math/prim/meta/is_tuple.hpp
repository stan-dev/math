#ifndef STAN_MATH_PRIM_META_IS_TUPLE_HPP
#define STAN_MATH_PRIM_META_IS_TUPLE_HPP

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
}  // namespace stan

#endif
