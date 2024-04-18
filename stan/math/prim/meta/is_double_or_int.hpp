#ifndef STAN_MATH_PRIM_META_IS_DOUBLE_OR_INT_HPP
#define STAN_MATH_PRIM_META_IS_DOUBLE_OR_INT_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <type_traits>

namespace stan {
/**
 * Checks if decayed type is a double or integer
 * @tparam The type to check
 * @ingroup type_trait
 */
template <typename T>
struct is_double_or_int
    : bool_constant<
          math::disjunction<std::is_same<double, std::decay_t<T>>,
                            std::is_same<int, std::decay_t<T>>>::value> {};

/*! \ingroup require_stan_scalar_real */
/*! \defgroup double_or_int_types double_or_int  */
/*! \addtogroup double_or_int_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_double_or_int */
/*! @tparam T the type to check */
template <typename T>
using require_double_or_int_t = require_t<is_double_or_int<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_double_or_int */
/*! @tparam T the type to check */
template <typename T>
using require_not_double_or_int_t
    = require_not_t<is_double_or_int<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_double_or_int */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_double_or_int_t
    = require_all_t<is_double_or_int<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_double_or_int */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_double_or_int_t
    = require_any_t<is_double_or_int<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_double_or_int */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_double_or_int_t
    = require_all_not_t<is_double_or_int<std::decay_t<Types>>...>;

/*! \brief Require at least one of the types do not satisfy @ref
 * is_double_or_int */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_not_double_or_int_t
    = require_any_not_t<is_double_or_int<std::decay_t<Types>>...>;
/*! @} */

}  // namespace stan
#endif
