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

//STAN_ADD_REQUIRE_UNARY(double_or_int, is_double_or_int, require_stan_scalar_real);
template <typename T>
using require_double_or_int_t = require_t<is_double_or_int<std::decay_t<T>>>;

template <typename T>
using require_not_double_or_int_t
    = require_not_t<is_double_or_int<std::decay_t<T>>>;

template <typename... Types>
using require_all_double_or_int_t
    = require_all_t<is_double_or_int<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_double_or_int_t
    = require_any_t<is_double_or_int<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_double_or_int_t
    = require_all_not_t<is_double_or_int<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_double_or_int_t
    = require_any_not_t<is_double_or_int<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_UNARY_INNER(double_or_int, is_double_or_int, require_stan_scalar_real);
template <typename T>
using require_vt_double_or_int
    = require_t<is_double_or_int<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_double_or_int
    = require_not_t<is_double_or_int<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_double_or_int
    = require_all_t<is_double_or_int<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_double_or_int
    = require_any_t<is_double_or_int<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_double_or_int
    = require_all_not_t<is_double_or_int<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_double_or_int
    = require_any_not_t<is_double_or_int<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_double_or_int
    = require_t<is_double_or_int<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_double_or_int
    = require_not_t<is_double_or_int<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_double_or_int
    = require_all_t<is_double_or_int<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_double_or_int
    = require_any_t<is_double_or_int<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_double_or_int
    = require_all_not_t<is_double_or_int<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_double_or_int
    = require_any_not_t<is_double_or_int<scalar_type_t<std::decay_t<Types>>>...>;


}  // namespace stan
#endif
