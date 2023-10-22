#ifndef STAN_MATH_PRIM_META_IS_VECTOR_LIKE_HPP
#define STAN_MATH_PRIM_META_IS_VECTOR_LIKE_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_detected.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {

namespace internal {
/**
 * @brief Used to detect if object has operator[](int) defined
 */
template <typename T>
using operator_bracket_t = decltype(std::declval<T>()[int{}]);
}  // namespace internal

/** \ingroup type_trait
 * Template metaprogram indicates whether a type is vector_like.
 *
 * A type is vector_like if an instance can be accessed like a
 * vector, i.e. square brackets.
 *
 * Access is_vector_like::value for the result.
 *
 * Default behavior is to use the is_vector template metaprogram.
 *
 * @tparam T Type to test
 */
template <typename T>
struct is_vector_like
    : bool_constant<is_detected<T, internal::operator_bracket_t>::value> {};

//STAN_ADD_REQUIRE_UNARY(vector_like, is_vector_like, require_std);
template <typename T>
using require_vector_like_t = require_t<is_vector_like<std::decay_t<T>>>;

template <typename T>
using require_not_vector_like_t
    = require_not_t<is_vector_like<std::decay_t<T>>>;

template <typename... Types>
using require_all_vector_like_t
    = require_all_t<is_vector_like<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_vector_like_t
    = require_any_t<is_vector_like<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_vector_like_t
    = require_all_not_t<is_vector_like<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_vector_like_t
    = require_any_not_t<is_vector_like<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_CONTAINER(vector_like, is_vector_like, require_std);
template <template <class...> class TypeCheck, class... Check>
using require_vector_like_vt = require_t<
    container_type_check_base<is_vector_like, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_vector_like_vt = require_not_t<
    container_type_check_base<is_vector_like, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_vector_like_vt = require_any_t<
    container_type_check_base<is_vector_like, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_vector_like_vt = require_any_not_t<
    container_type_check_base<is_vector_like, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vector_like_vt = require_all_t<
    container_type_check_base<is_vector_like, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_vector_like_vt = require_all_not_t<
    container_type_check_base<is_vector_like, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_vector_like_st = require_t<
    container_type_check_base<is_vector_like, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_vector_like_st = require_not_t<
    container_type_check_base<is_vector_like, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_vector_like_st = require_any_t<
    container_type_check_base<is_vector_like, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_vector_like_st = require_any_not_t<
    container_type_check_base<is_vector_like, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vector_like_st = require_all_t<
    container_type_check_base<is_vector_like, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_vector_like_st = require_all_not_t<
    container_type_check_base<is_vector_like, scalar_type_t, TypeCheck, Check>...>;


}  // namespace stan
#endif
