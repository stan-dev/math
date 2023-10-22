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

//STAN_ADD_REQUIRE_UNARY(plain_type, is_plain_type, require_eigens_types);
template <typename T>
using require_plain_type_t = require_t<is_plain_type<std::decay_t<T>>>;

template <typename T>
using require_not_plain_type_t
    = require_not_t<is_plain_type<std::decay_t<T>>>;

template <typename... Types>
using require_all_plain_type_t
    = require_all_t<is_plain_type<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_plain_type_t
    = require_any_t<is_plain_type<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_plain_type_t
    = require_all_not_t<is_plain_type<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_plain_type_t
    = require_any_not_t<is_plain_type<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_UNARY_INNER(plain_type, is_plain_type, require_eigens_types);
template <typename T>
using require_vt_plain_type
    = require_t<is_plain_type<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_plain_type
    = require_not_t<is_plain_type<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_plain_type
    = require_all_t<is_plain_type<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_plain_type
    = require_any_t<is_plain_type<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_plain_type
    = require_all_not_t<is_plain_type<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_plain_type
    = require_any_not_t<is_plain_type<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_plain_type
    = require_t<is_plain_type<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_plain_type
    = require_not_t<is_plain_type<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_plain_type
    = require_all_t<is_plain_type<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_plain_type
    = require_any_t<is_plain_type<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_plain_type
    = require_all_not_t<is_plain_type<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_plain_type
    = require_any_not_t<is_plain_type<scalar_type_t<std::decay_t<Types>>>...>;


}  // namespace stan
#endif
