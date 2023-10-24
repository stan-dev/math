#ifndef STAN_MATH_PRIM_META_IS_MATRIX_CL_HPP
#define STAN_MATH_PRIM_META_IS_MATRIX_CL_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {
namespace math {

template <typename T>
class arena_matrix_cl;

/**
 * Non-templated base class for `matrix_cl` simplifies checking if something is
 * matrix_cl.
 */
class matrix_cl_base {};

}  // namespace math

/** \ingroup type_traits
 * Checks if the decayed type of T is a matrix_cl.
 */
template <typename T>
struct is_matrix_cl
    : public std::is_base_of<math::matrix_cl_base, std::decay_t<T>> {};

// STAN_ADD_REQUIRE_UNARY(matrix_cl, is_matrix_cl, matrix_cl_group);
template <typename T>
using require_matrix_cl_t = require_t<is_matrix_cl<std::decay_t<T>>>;

template <typename T>
using require_not_matrix_cl_t = require_not_t<is_matrix_cl<std::decay_t<T>>>;

template <typename... Types>
using require_all_matrix_cl_t
    = require_all_t<is_matrix_cl<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_matrix_cl_t
    = require_any_t<is_matrix_cl<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_matrix_cl_t
    = require_all_not_t<is_matrix_cl<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_matrix_cl_t
    = require_any_not_t<is_matrix_cl<std::decay_t<Types>>...>;

// STAN_ADD_REQUIRE_CONTAINER(matrix_cl, is_matrix_cl, matrix_cl_group);
template <template <class...> class TypeCheck, class... Check>
using require_matrix_cl_vt = require_t<
    container_type_check_base<is_matrix_cl, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_matrix_cl_vt = require_not_t<
    container_type_check_base<is_matrix_cl, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_matrix_cl_vt = require_any_t<
    container_type_check_base<is_matrix_cl, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_matrix_cl_vt = require_any_not_t<
    container_type_check_base<is_matrix_cl, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_matrix_cl_vt = require_all_t<
    container_type_check_base<is_matrix_cl, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_matrix_cl_vt = require_all_not_t<
    container_type_check_base<is_matrix_cl, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_matrix_cl_st
    = require_t<container_type_check_base<is_matrix_cl, scalar_type_t,
                                          TypeCheck, Check...>>;

namespace internal {

template <typename T>
struct is_arena_matrix_cl_impl : public std::false_type {};

template <typename T>
struct is_arena_matrix_cl_impl<math::arena_matrix_cl<T>>
    : public std::true_type {};

}  // namespace internal

template <typename T>
struct is_arena_matrix_cl
    : public internal::is_arena_matrix_cl_impl<std::decay_t<T>> {};

// STAN_ADD_REQUIRE_UNARY(arena_matrix_cl, is_arena_matrix_cl, matrix_cl_group);
// STAN_ADD_REQUIRE_CONTAINER(arena_matrix_cl, is_arena_matrix_cl,
// matrix_cl_group);
}  // namespace stan
#endif
