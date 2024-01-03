#ifndef STAN_MATH_PRIM_META_IS_PLAIN_TYPE_HPP
#define STAN_MATH_PRIM_META_IS_PLAIN_TYPE_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/prim/meta/plain_type.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Checks whether the template type `T` is an assignable type. This is used
 * to detect whether a type is an Eigen matrix expression.
 */
template <typename S>
using is_plain_type = std::is_same<std::decay_t<S>, plain_type_t<S>>;

/*! \ingroup require_eigens_types */
/*! \defgroup plain_type_types plain_type  */
/*! \addtogroup plain_type_types */
/*! @{ */

/*! \brief Require type satisfies is_plain_type */
template <typename T>
using require_plain_type_t = require_t<is_plain_type<std::decay_t<T>>>;

/*! \brief Require type does not satisfy is_plain_type */
template <typename T>
using require_not_plain_type_t = require_not_t<is_plain_type<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy is_plain_type */
template <typename... Types>
using require_all_plain_type_t
    = require_all_t<is_plain_type<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy is_plain_type */
template <typename... Types>
using require_any_not_plain_type_t
    = require_any_not_t<is_plain_type<std::decay_t<Types>>...>;
/*! @} */

// STAN_ADD_REQUIRE_UNARY_INNER(plain_type, is_plain_type,
// require_eigens_types);
}  // namespace stan
#endif
