#ifndef STAN_MATH_PRIM_META_IS_ARENA_MATRIX_HPP
#define STAN_MATH_PRIM_META_IS_ARENA_MATRIX_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Defines a static member named value which is defined to be true
 * if the type is `arena_matrix<T>`
 */
template <typename T, typename = void>
struct is_arena_matrix : std::false_type {};

// STAN_ADD_REQUIRE_UNARY(arena_matrix, is_arena_matrix, require_eigens_types);
template <typename T>
using require_arena_matrix_t = require_t<is_arena_matrix<std::decay_t<T>>>;

template <typename T>
using require_not_arena_matrix_t
    = require_not_t<is_arena_matrix<std::decay_t<T>>>;

template <typename... Types>
using require_all_arena_matrix_t
    = require_all_t<is_arena_matrix<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_arena_matrix_t
    = require_any_t<is_arena_matrix<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_arena_matrix_t
    = require_all_not_t<is_arena_matrix<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_arena_matrix_t
    = require_any_not_t<is_arena_matrix<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_CONTAINER(arena_matrix, is_arena_matrix, require_eigens_types);
template <template <class...> class TypeCheck, class... Check>
using require_arena_matrix_vt = require_t<
    container_type_check_base<is_arena_matrix, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_arena_matrix_vt = require_not_t<
    container_type_check_base<is_arena_matrix, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_arena_matrix_vt = require_any_t<
    container_type_check_base<is_arena_matrix, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_arena_matrix_vt = require_any_not_t<
    container_type_check_base<is_arena_matrix, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_arena_matrix_vt = require_all_t<
    container_type_check_base<is_arena_matrix, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_arena_matrix_vt = require_all_not_t<
    container_type_check_base<is_arena_matrix, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_arena_matrix_st = require_t<
    container_type_check_base<is_arena_matrix, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_arena_matrix_st = require_not_t<
    container_type_check_base<is_arena_matrix, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_arena_matrix_st = require_any_t<
    container_type_check_base<is_arena_matrix, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_arena_matrix_st = require_any_not_t<
    container_type_check_base<is_arena_matrix, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_arena_matrix_st = require_all_t<
    container_type_check_base<is_arena_matrix, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_arena_matrix_st = require_all_not_t<
    container_type_check_base<is_arena_matrix, scalar_type_t, TypeCheck, Check>...>;

//STAN_ADD_REQUIRE_UNARY_INNER(arena_matrix, is_arena_matrix, require_eigens_types);
template <typename T>
using require_vt_arena_matrix
    = require_t<is_arena_matrix<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_arena_matrix
    = require_not_t<is_arena_matrix<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_arena_matrix
    = require_all_t<is_arena_matrix<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_arena_matrix
    = require_any_t<is_arena_matrix<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_arena_matrix
    = require_all_not_t<is_arena_matrix<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_arena_matrix
    = require_any_not_t<is_arena_matrix<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_arena_matrix
    = require_t<is_arena_matrix<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_arena_matrix
    = require_not_t<is_arena_matrix<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_arena_matrix
    = require_all_t<is_arena_matrix<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_arena_matrix
    = require_any_t<is_arena_matrix<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_arena_matrix
    = require_all_not_t<is_arena_matrix<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_arena_matrix
    = require_any_not_t<is_arena_matrix<scalar_type_t<std::decay_t<Types>>>...>;


}  // namespace stan
#endif
