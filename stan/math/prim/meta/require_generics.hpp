#ifndef STAN_MATH_PRIM_META_REQUIRE_GENERICS_HPP
#define STAN_MATH_PRIM_META_REQUIRE_GENERICS_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_container.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_vector_like.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <type_traits>
#include <string>

namespace stan {

STAN_ADD_REQUIRE_BINARY(same, std::is_same, require_std);
STAN_ADD_REQUIRE_BINARY_SCALAR(same, std::is_same, require_std);
STAN_ADD_REQUIRE_BINARY_VALUE(same, std::is_same, require_std);

STAN_ADD_REQUIRE_BINARY(convertible, std::is_convertible, require_std);
STAN_ADD_REQUIRE_BINARY_SCALAR(convertible, std::is_convertible, require_std);
STAN_ADD_REQUIRE_BINARY_VALUE(convertible, std::is_convertible, require_std);

STAN_ADD_REQUIRE_UNARY(arithmetic, std::is_arithmetic,
                       require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_SCALAR(arithmetic, std::is_arithmetic,
                              require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_VALUE(arithmetic, std::is_arithmetic,
                             require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY(floating_point, std::is_floating_point,
                       require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_SCALAR(floating_point, std::is_floating_point,
                              require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_VALUE(floating_point, std::is_floating_point,
                             require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY(index, std::is_integral, require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_SCALAR(index, std::is_integral,
                              require_stan_scalar_real);
STAN_ADD_REQUIRE_UNARY_VALUE(index, std::is_integral, require_stan_scalar_real);

template <typename Row, typename Col>
using require_eigen_row_and_col_t = require_t<
    math::conjunction<is_eigen_row_vector<Row>, is_eigen_col_vector<Col>>>;

template <typename Row, typename Col>
using require_not_eigen_row_and_col_t = require_not_t<
    math::conjunction<is_eigen_row_vector<Row>, is_eigen_col_vector<Col>>>;

/** @}*/
}  // namespace stan
#endif
