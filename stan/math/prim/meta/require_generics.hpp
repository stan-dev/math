#ifndef STAN_MATH_PRIM_META_REQUIRE_GENERICS_HPP
#define STAN_MATH_PRIM_META_REQUIRE_GENERICS_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

STAN_ADD_REQUIRE_BINARY(same, std::is_same, require_std);
STAN_ADD_REQUIRE_BINARY_INNER(same, std::is_same, require_std);

STAN_ADD_REQUIRE_BINARY(convertible, std::is_convertible, require_std);
STAN_ADD_REQUIRE_BINARY_INNER(convertible, std::is_convertible, require_std);
STAN_ADD_REQUIRE_BINARY(assignable, std::is_assignable, require_std);
STAN_ADD_REQUIRE_BINARY_INNER(assignable, std::is_assignable, require_std);

STAN_ADD_REQUIRE_BINARY(constructible, std::is_constructible, require_std);
STAN_ADD_REQUIRE_BINARY_INNER(constructible, std::is_constructible,
                              require_std);

//STAN_ADD_REQUIRE_UNARY(arithmetic, std::is_arithmetic, require_stan_scalar_real);
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
using require_all_not_arithmetic_t
    = require_all_not_t<std::is_arithmetic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_arithmetic_t
    = require_any_not_t<std::is_arithmetic<std::decay_t<Types>>...>;

  
STAN_ADD_REQUIRE_UNARY_INNER(arithmetic, std::is_arithmetic,
                             require_stan_scalar_real);
//STAN_ADD_REQUIRE_UNARY(floating_point, std::is_floating_point, require_stan_scalar_real);
template <typename T>
using require_floating_point_t = require_t<std::is_floating_point<std::decay_t<T>>>;

template <typename T>
using require_not_floating_point_t
    = require_not_t<std::is_floating_point<std::decay_t<T>>>;

template <typename... Types>
using require_all_floating_point_t
    = require_all_t<std::is_floating_point<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_floating_point_t
    = require_any_t<std::is_floating_point<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_floating_point_t
    = require_all_not_t<std::is_floating_point<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_floating_point_t
    = require_any_not_t<std::is_floating_point<std::decay_t<Types>>...>;

  
STAN_ADD_REQUIRE_UNARY_INNER(floating_point, std::is_floating_point,
                             require_stan_scalar_real);
//STAN_ADD_REQUIRE_UNARY(integral, std::is_integral, require_stan_scalar_real);
template <typename T>
using require_integral_t = require_t<std::is_integral<std::decay_t<T>>>;

template <typename T>
using require_not_integral_t
    = require_not_t<std::is_integral<std::decay_t<T>>>;

template <typename... Types>
using require_all_integral_t
    = require_all_t<std::is_integral<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_integral_t
    = require_any_t<std::is_integral<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_integral_t
    = require_all_not_t<std::is_integral<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_integral_t
    = require_any_not_t<std::is_integral<std::decay_t<Types>>...>;

  
STAN_ADD_REQUIRE_UNARY_INNER(integral, std::is_integral,
                             require_stan_scalar_real);

}  // namespace stan
#endif
