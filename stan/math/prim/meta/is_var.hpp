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

// STAN_ADD_REQUIRE_UNARY(var, is_var, require_stan_scalar_real);
template <typename T>
using require_var_t = require_t<is_var<std::decay_t<T>>>;

template <typename T>
using require_not_var_t = require_not_t<is_var<std::decay_t<T>>>;

template <typename... Types>
using require_all_var_t = require_all_t<is_var<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_var_t = require_any_t<is_var<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_var_t = require_all_not_t<is_var<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_var_t = require_any_not_t<is_var<std::decay_t<Types>>...>;

// STAN_ADD_REQUIRE_CONTAINER(var, is_var, require_stan_scalar_real);
template <template <class...> class TypeCheck, class... Check>
using require_var_vt = require_t<
    container_type_check_base<is_var, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_var_vt = require_any_t<
    container_type_check_base<is_var, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_var_vt = require_all_t<
    container_type_check_base<is_var, value_type_t, TypeCheck, Check>...>;

// STAN_ADD_REQUIRE_UNARY_INNER(var, is_var, require_stan_scalar_real);
template <typename T>
using require_not_vt_var = require_not_t<is_var<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_any_vt_var
    = require_any_t<is_var<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_var
    = require_all_not_t<is_var<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_var = require_t<is_var<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_var
    = require_not_t<is_var<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_var
    = require_all_t<is_var<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_var
    = require_any_t<is_var<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_var
    = require_all_not_t<is_var<scalar_type_t<std::decay_t<Types>>>...>;

template <typename T>
struct value_type<T, std::enable_if_t<is_var<T>::value>> {
  using type = typename std::decay_t<T>::value_type;
};

}  // namespace stan
#endif
