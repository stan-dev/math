#ifndef STAN_MATH_PRIM_META_IS_FVAR_HPP
#define STAN_MATH_PRIM_META_IS_FVAR_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Defines a static member function type which is defined to be false
 * as the primitive scalar types cannot be a stan::math::fvar type.
 */
template <typename T, typename = void>
struct is_fvar : std::false_type {};

/*! \ingroup require_stan_scalar_real */
/*! \defgroup fvar_types fvar  */
/*! \addtogroup fvar_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_fvar */
/*! @tparam T the type to check */
template <typename T>
using require_fvar_t = require_t<is_fvar<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_fvar */
/*! @tparam T the type to check */
template <typename T>
using require_not_fvar_t = require_not_t<is_fvar<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_fvar */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_fvar_t = require_all_t<is_fvar<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_fvar */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_fvar_t = require_any_t<is_fvar<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_fvar */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_fvar_t
    = require_all_not_t<is_fvar<std::decay_t<Types>>...>;

/*! \brief Require at least one of the types do not satisfy @ref is_fvar */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_not_fvar_t
    = require_any_not_t<is_fvar<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup require_stan_scalar_real */
/*! \addtogroup fvar_types */
/*! @{ */

/*! \brief Require any of the value types satisfy @ref is_fvar */
/*! @tparam Types The types with a valid overload of @ref value_type available
 */
template <typename... Types>
using require_any_vt_fvar
    = require_any_t<is_fvar<value_type_t<std::decay_t<Types>>>...>;

/*! \brief Require none of the value types satisfy @ref is_fvar */
/*! @tparam Types The types with a valid overload of @ref value_type available
 */
template <typename... Types>
using require_all_not_vt_fvar
    = require_all_not_t<is_fvar<value_type_t<std::decay_t<Types>>>...>;

/*! \brief Require scalar type satisfies @ref is_fvar */
/*! @tparam T A type with a valid overload of @ref scalar_type available */
template <typename T>
using require_st_fvar = require_t<is_fvar<scalar_type_t<std::decay_t<T>>>>;

/*! \brief Require scalar type does not satisfy @ref is_fvar */
/*! @tparam T A type with a valid overload of @ref scalar_type available */
template <typename T>
using require_not_st_fvar
    = require_not_t<is_fvar<scalar_type_t<std::decay_t<T>>>>;

/*! \brief Require any of the scalar types satisfy @ref is_fvar */
/*! @tparam Types The types with a valid overload of @ref scalar_type available
 */
template <typename... Types>
using require_any_st_fvar
    = require_any_t<is_fvar<scalar_type_t<std::decay_t<Types>>>...>;

/*! \brief Require none of the scalar types satisfy @ref is_fvar */
/*! @tparam Types The types with a valid overload of @ref scalar_type available
 */
template <typename... Types>
using require_all_not_st_fvar
    = require_all_not_t<is_fvar<scalar_type_t<std::decay_t<Types>>>...>;
/*! @} */

}  // namespace stan
#endif
