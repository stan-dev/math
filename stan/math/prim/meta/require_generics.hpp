#ifndef STAN_MATH_PRIM_META_REQUIRE_GENERICS_HPP
#define STAN_MATH_PRIM_META_REQUIRE_GENERICS_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

STAN_ADD_REQUIRE_BINARY(same, std::is_same, require_std);
STAN_ADD_REQUIRE_BINARY_INNER(same, std::is_same, require_std);

STAN_ADD_REQUIRE_BINARY(convertible, std::is_convertible, require_std);
STAN_ADD_REQUIRE_BINARY_INNER(convertible, std::is_convertible, require_std);

/*! ingroup require_std */
/*! defgroup assignable_types check_type  */
/*! addtogroup assignable_types */
/*! @{ */
/*! brief Require types `T` and `S` satisfies `std::is_assignable` */
template <typename T, typename S>
using require_assignable_t = require_t<std::is_assignable<T, S>>;

/*! brief Require types `T` and `S` does not satisfy `std::is_assignable` */
template <typename T, typename S>
using require_not_assignable_t = require_not_t<std::is_assignable<T, S>>;

/*! brief Require `T` and all of the `Types` satisfy `std::is_assignable` */
template <typename T, typename... Types>
using require_all_assignable_t = require_all_t<std::is_assignable<T, Types>...>;

/*! brief Require any of the `Types` and `T` satisfy `std::is_assignable` */
template <typename T, typename... Types>
using require_any_assignable_t = require_any_t<std::is_assignable<T, Types>...>;

/*! brief Require none of the `Types` and `T` satisfy `std::is_assignable` */
template <typename T, typename... Types>
using require_all_not_assignable_t
    = require_all_not_t<std::is_assignable<T, Types>...>;

/*! brief Any one of the `Types` and `T` do not satisfy `std::is_assignable` */
template <typename T, typename... Types>
using require_any_not_assignable_t
    = require_any_not_t<std::is_assignable<T, Types>...>;
/*! @} */

STAN_ADD_REQUIRE_UNARY(arithmetic, std::is_arithmetic,
                       require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_INNER(arithmetic, std::is_arithmetic,
                             require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY(floating_point, std::is_floating_point,
                       require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_INNER(floating_point, std::is_floating_point,
                             require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY(integral, std::is_integral, require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_INNER(integral, std::is_integral,
                             require_stan_scalar_real);

}  // namespace stan
#endif
