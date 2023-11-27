#ifndef STAN_MATH_PRIM_META_REQUIRE_GENERICS_HPP
#define STAN_MATH_PRIM_META_REQUIRE_GENERICS_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

// STAN_ADD_REQUIRE_BINARY(same, std::is_same, require_std);
template <typename T, typename S>
using require_same_t
    = require_t<std::is_same<std::decay_t<T>, std::decay_t<S>>>;

template <typename T, typename S>
using require_not_same_t
    = require_not_t<std::is_same<std::decay_t<T>, std::decay_t<S>>>;

template <typename T, typename... Types>
using require_all_same_t
    = require_all_t<std::is_same<std::decay_t<T>, std::decay_t<Types>>...>;

template <typename T, typename... Types>
using require_any_same_t
    = require_any_t<std::is_same<std::decay_t<T>, std::decay_t<Types>>...>;

template <typename T, typename... Types>
using require_all_not_same_t
    = require_all_not_t<std::is_same<std::decay_t<T>, std::decay_t<Types>>...>;

// STAN_ADD_REQUIRE_BINARY_INNER(same, std::is_same, require_std);
template <typename T, typename S>
using require_st_same = require_t<std::is_same<scalar_type_t<std::decay_t<T>>,
                                               scalar_type_t<std::decay_t<S>>>>;

template <typename T, typename S>
using require_not_st_same
    = require_not_t<std::is_same<scalar_type_t<std::decay_t<T>>,
                                 scalar_type_t<std::decay_t<S>>>>;

template <typename T, typename... Types>
using require_all_st_same
    = require_all_t<std::is_same<scalar_type_t<std::decay_t<T>>,
                                 scalar_type_t<std::decay_t<Types>>>...>;

template <typename T, typename S>
using require_vt_same = require_t<
  std::is_same<value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<S>>>>;

template <typename T, typename... Types>
using require_all_vt_same
    = require_all_t<std::is_same<value_type_t<std::decay_t<T>>,
                                 value_type_t<std::decay_t<Types>>>...>;

// STAN_ADD_REQUIRE_BINARY(convertible, std::is_convertible, require_std);
template <typename T, typename S>
using require_convertible_t
    = require_t<std::is_convertible<std::decay_t<T>, std::decay_t<S>>>;

// STAN_ADD_REQUIRE_BINARY_INNER(convertible, std::is_convertible, require_std);
// STAN_ADD_REQUIRE_BINARY(assignable, std::is_assignable, require_std);
template <typename T, typename S>
using require_assignable_t
    = require_t<std::is_assignable<std::decay_t<T>, std::decay_t<S>>>;

// STAN_ADD_REQUIRE_BINARY_INNER(assignable, std::is_assignable, require_std);
// STAN_ADD_REQUIRE_BINARY(constructible, std::is_constructible, require_std);
template <typename T, typename S>
using require_constructible_t
    = require_t<std::is_constructible<std::decay_t<T>, std::decay_t<S>>>;

// STAN_ADD_REQUIRE_BINARY_INNER(constructible, std::is_constructible,
// require_std);
// STAN_ADD_REQUIRE_UNARY(arithmetic, std::is_arithmetic,
// require_stan_scalar_real);
template <typename T>
using require_arithmetic_t = require_t<std::is_arithmetic<std::decay_t<T>>>;

template <typename T>
using require_not_arithmetic_t
    = require_not_t<std::is_arithmetic<std::decay_t<T>>>;

template <typename... Types>
using require_all_arithmetic_t
    = require_all_t<std::is_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_arithmetic_t
    = require_any_t<std::is_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_arithmetic_t
    = require_any_not_t<std::is_arithmetic<std::decay_t<Types>>...>;

// STAN_ADD_REQUIRE_UNARY_INNER(arithmetic, std::is_arithmetic,
// require_stan_scalar_real);
template <typename... Types>
using require_all_vt_arithmetic
    = require_all_t<std::is_arithmetic<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_arithmetic = require_any_not_t<
    std::is_arithmetic<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_arithmetic
    = require_t<std::is_arithmetic<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_arithmetic
    = require_not_t<std::is_arithmetic<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_arithmetic
    = require_all_t<std::is_arithmetic<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_arithmetic
    = require_any_t<std::is_arithmetic<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_arithmetic = require_any_not_t<
    std::is_arithmetic<scalar_type_t<std::decay_t<Types>>>...>;

// STAN_ADD_REQUIRE_UNARY(floating_point, std::is_floating_point,
// require_stan_scalar_real);
template <typename T>
using require_floating_point_t
    = require_t<std::is_floating_point<std::decay_t<T>>>;

// STAN_ADD_REQUIRE_UNARY_INNER(floating_point, std::is_floating_point,
// require_stan_scalar_real);
// STAN_ADD_REQUIRE_UNARY(integral, std::is_integral, require_stan_scalar_real);
template <typename T>
using require_integral_t = require_t<std::is_integral<std::decay_t<T>>>;

// STAN_ADD_REQUIRE_UNARY_INNER(integral, std::is_integral,
// require_stan_scalar_real);
template <typename T>
using require_vt_integral
    = require_t<std::is_integral<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_st_integral
    = require_t<std::is_integral<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_integral
    = require_not_t<std::is_integral<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_not_st_integral = require_all_not_t<
    std::is_integral<scalar_type_t<std::decay_t<Types>>>...>;

}  // namespace stan
#endif
