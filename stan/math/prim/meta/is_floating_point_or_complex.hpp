#ifndef STAN_MATH_PRIM_META_IS_FLOATING_POINT_OR_COMPLEX_HPP
#define STAN_MATH_PRIM_META_IS_FLOATING_POINT_OR_COMPLEX_HPP

#include <stan/math/prim/meta/is_complex.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

template <typename T>
struct is_floating_point_or_complex :
    std::bool_constant<std::is_floating_point_v<std::decay_t<T>> || stan::is_complex_v<std::decay_t<T>>> {};

template <typename T>
constexpr bool is_floating_point_or_complex_v = is_floating_point_or_complex<T>::value;

/*! \ingroup require_stan_scalar_floating_point_or_complex */
/*! \defgroup floating_point_or_complex_types complex  */
/*! \addtogroup floating_point_or_complex_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_floating_point_or_complex */
/*! @tparam T the type to check */
template <typename T>
using require_floating_point_or_complex_t = require_t<is_floating_point_or_complex<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_floating_point_or_complex */
/*! @tparam T the type to check */
template <typename T>
using require_not_floating_point_or_complex_t = require_not_t<is_floating_point_or_complex<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_floating_point_or_complex */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_floating_point_or_complex_t = require_all_t<is_floating_point_or_complex<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_floating_point_or_complex */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_floating_point_or_complex_t = require_any_t<is_floating_point_or_complex<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_floating_point_or_complex */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_floating_point_or_complex_t
    = require_all_not_t<is_floating_point_or_complex<std::decay_t<Types>>...>;
/*! @} */

}

#endif
