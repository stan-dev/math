#ifndef STAN_MATH_PRIM_META_IS_VARI_HPP
#define STAN_MATH_PRIM_META_IS_VARI_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Specialization for checking if value of T minus cv qualifier and pointer is a
 * vari.
 */
template <typename T, typename = void>
struct is_vari : std::false_type {};

//STAN_ADD_REQUIRE_UNARY(vari, is_vari, require_stan_scalar_real);
template <typename T>
using require_vari_t = require_t<is_vari<std::decay_t<T>>>;

template <typename T>
using require_not_vari_t
    = require_not_t<is_vari<std::decay_t<T>>>;

template <typename... Types>
using require_all_vari_t
    = require_all_t<is_vari<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_vari_t
    = require_any_t<is_vari<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_vari_t
    = require_all_not_t<is_vari<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_vari_t
    = require_any_not_t<is_vari<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_CONTAINER(vari, is_vari, require_stan_scalar_real);
template <template <class...> class TypeCheck, class... Check>
using require_vari_vt = require_t<
    container_type_check_base<is_vari, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_vari_vt = require_not_t<
    container_type_check_base<is_vari, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_vari_vt = require_any_t<
    container_type_check_base<is_vari, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_vari_vt = require_any_not_t<
    container_type_check_base<is_vari, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vari_vt = require_all_t<
    container_type_check_base<is_vari, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_vari_vt = require_all_not_t<
    container_type_check_base<is_vari, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_vari_st = require_t<
    container_type_check_base<is_vari, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_vari_st = require_not_t<
    container_type_check_base<is_vari, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_vari_st = require_any_t<
    container_type_check_base<is_vari, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_vari_st = require_any_not_t<
    container_type_check_base<is_vari, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vari_st = require_all_t<
    container_type_check_base<is_vari, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_vari_st = require_all_not_t<
    container_type_check_base<is_vari, scalar_type_t, TypeCheck, Check>...>;


template <typename T>
struct value_type<T, require_t<is_vari<T>>> {
  using type = typename std::decay_t<T>::value_type;
};

}  // namespace stan
#endif
