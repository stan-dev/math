#ifndef STAN_MATH_PRIM_META_IS_VAR_HPP
#define STAN_MATH_PRIM_META_IS_VAR_HPP

#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Defines a static member named value which is defined to be false
 * as the primitive scalar types cannot be a stan::math::var type.
 */
template <typename T, typename = void>
struct is_var : std::false_type {};

/*! \ingroup require_stan_scalar_real */
/*! \defgroup var_types var  */
/*! \addtogroup var_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_var */
/*! @tparam T the type to check */
template <typename T>
using require_var_t = require_t<is_var<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_var */
/*! @tparam T the type to check */
template <typename T>
using require_not_var_t = require_not_t<is_var<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_var */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_var_t = require_all_t<is_var<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_var */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_var_t = require_any_t<is_var<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_var */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_var_t = require_all_not_t<is_var<std::decay_t<Types>>...>;

/*! \brief Require at least one of the types do not satisfy @ref is_var */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_not_var_t = require_any_not_t<is_var<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup require_stan_scalar_real */
/*! \defgroup var_types var  */
/*! \addtogroup var_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_var */
/*! and `value_type` satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_var for and whose @ref value_type is
 * checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_var_vt = require_t<
    container_type_check_base<is_var, value_type_t, TypeCheck, Check...>>;

/*! \brief Require any of the types satisfy is_var */
/*! and any of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_var for and whose @ref value_type is
 * checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_var_vt = require_any_t<
    container_type_check_base<is_var, value_type_t, TypeCheck, Check>...>;

/*! \brief Require all of the types satisfy is_var */
/*! and all of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_var for and whose @ref value_type is
 * checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_var_vt = require_all_t<
    container_type_check_base<is_var, value_type_t, TypeCheck, Check>...>;
/*! @} */

/*! \ingroup require_stan_scalar_real */
/*! \addtogroup var_types */
/*! @{ */

/*! \brief Require value type does not satisfy @ref is_var */
/*! @tparam T A type with a valid overload of @ref value_type available */
template <typename T>
using require_not_vt_var = require_not_t<is_var<value_type_t<std::decay_t<T>>>>;

/*! \brief Require any of the @ref value_types satisfy @ref is_var */
/*! @tparam Types The types with a valid overload of @ref value_type available
 */
template <typename... Types>
using require_any_vt_var
    = require_any_t<is_var<value_type_t<std::decay_t<Types>>>...>;

/*! \brief Require none of the @ref value_types satisfy @ref is_var */
/*! @tparam Types The types with a valid overload of @ref value_type available
 */
template <typename... Types>
using require_all_not_vt_var
    = require_all_not_t<is_var<value_type_t<std::decay_t<Types>>>...>;

/*! \brief Require @ref scalar_type satisfies @ref is_var */
/*! @tparam T A type with a valid overload of @ref scalar_type available */
template <typename T>
using require_st_var = require_t<is_var<scalar_type_t<std::decay_t<T>>>>;

/*! \brief Require scalar type does not satisfy @ref is_var */
/*! @tparam T A type with a valid overload of @ref scalar_type available */
template <typename T>
using require_not_st_var
    = require_not_t<is_var<scalar_type_t<std::decay_t<T>>>>;

/*! \brief Require all of the scalar types satisfy @ref is_var */
/*! @tparam Types The types with a valid overload of @ref scalar_type available
 */
template <typename... Types>
using require_all_st_var
    = require_all_t<is_var<scalar_type_t<std::decay_t<Types>>>...>;

/*! \brief Require any of the scalar types satisfy @ref is_var */
/*! @tparam Types The types with a valid overload of @ref scalar_type available
 */
template <typename... Types>
using require_any_st_var
    = require_any_t<is_var<scalar_type_t<std::decay_t<Types>>>...>;

/*! \brief Require none of the scalar types satisfy @ref is_var */
/*! @tparam Types The types with a valid overload of @ref scalar_type available
 */
template <typename... Types>
using require_all_not_st_var
    = require_all_not_t<is_var<scalar_type_t<std::decay_t<Types>>>...>;
/*! @} */

template <typename T>
struct value_type<T, std::enable_if_t<is_var<T>::value>> {
  using type = typename std::decay_t<T>::value_type;
};

}  // namespace stan
#endif
