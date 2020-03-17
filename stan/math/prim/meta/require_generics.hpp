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
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>
#include <string>

namespace stan {

/** \addtogroup type_trait
 *  @{
 */



/**
 * Used as the base for checking whether a type is a container with
 * an underlying type
 *
 * @tparam ContainerCheck Templated struct or alias that wraps a static constant
 * value called type. Used to check the container satisfies a particular type
 * check.
 * @tparam CheckType Templated struct or alias that wraps a static constant
 * value called type. Used to check the container's underlying type satisfies a
 * particular type check.
 *
 */
template <template <class...> class ContainerCheck,
          template <class...> class TypeCheck, class... Check>
using container_value_type_check_base = bool_constant<
    math::conjunction<ContainerCheck<std::decay_t<Check>>...,
                      TypeCheck<value_type_t<Check>>...>::value>;

/**
 * Meta type trait to check if type is a standard vector and value type passes
 * additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_std_vector_value_check
    : container_value_type_check_base<is_std_vector, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type is a standard vector or eigen vector and
 * value type passes additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_vector_value_check
    : container_value_type_check_base<is_vector, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type is a standard vector, eigen vector, or array
 * and value type passes additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_vector_like_value_check
    : container_value_type_check_base<is_vector_like, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type inherits from EigenBase.
 */
template <template <class...> class TypeCheck, class... Check>
struct is_eigen_value_check
    : container_value_type_check_base<is_eigen, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type inherits from EigenBase and has row or col
 * values at compile time of 1.
 */
template <template <class...> class TypeCheck, class... Check>
struct is_eigen_vector_value_check
    : container_value_type_check_base<is_eigen_vector, TypeCheck, Check...> {};

/**
 * Used as the base for checking whether a type is a container with
 * an underlying scalar type
 *
 * @tparam ContainerCheck Templated struct or alias that wraps a static constant
 * scalar called type. Used to check the container satisfies a particular type
 * check.
 * @tparam CheckType Templated struct or alias that wraps a static constant
 * scalar called type. Used to check the container's underlying type satisfies a
 * particular type check.
 *
 */
template <template <class...> class ContainerCheck,
          template <class...> class TypeCheck, class... Check>
using container_scalar_type_check_base = bool_constant<
    math::conjunction<ContainerCheck<std::decay_t<Check>>...,
                      TypeCheck<scalar_type_t<Check>>...>::value>;

/**
 * Meta type trait to check if type is a standard vector and scalar type passes
 * additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_std_vector_scalar_check
    : container_scalar_type_check_base<is_std_vector, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type is a standard vector or eigen vector and
 * scalar type passes additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_vector_scalar_check
    : container_scalar_type_check_base<is_vector, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type is a standard vector, eigen vector, or array
 * and scalar type passes additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_vector_like_scalar_check
    : container_scalar_type_check_base<is_vector_like, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type inherits from EigenBase and scalar type pass
 * additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_eigen_scalar_check
    : container_scalar_type_check_base<is_eigen, TypeCheck, Check...> {};

/**
 * Meta type trait to check if type inherits from EigenBase, has 1 row or column
 * compile time size, and scalar type passes additional type trait check
 */
template <template <class...> class TypeCheck, class... Check>
struct is_eigen_vector_scalar_check
    : container_scalar_type_check_base<is_eigen_vector, TypeCheck, Check...> {};
/** @}*/

/** \addtogroup require_base_types
 *  @{
 */

STAN_ADD_REQUIRE_BINARY(same, std::is_same);
STAN_ADD_REQUIRE_BINARY_SCALAR(same, std::is_same);
STAN_ADD_REQUIRE_BINARY_VALUE(same, std::is_same);

STAN_ADD_REQUIRE_BINARY(convertible, std::is_convertible);
STAN_ADD_REQUIRE_BINARY_SCALAR(convertible, std::is_convertible);
STAN_ADD_REQUIRE_BINARY_VALUE(convertible, std::is_convertible);

STAN_ADD_REQUIRE_UNARY(arithmetic, std::is_arithmetic);
STAN_ADD_REQUIRE_UNARY_SCALAR(arithmetic, std::is_arithmetic);
STAN_ADD_REQUIRE_UNARY_VALUE(arithmetic, std::is_arithmetic);
STAN_ADD_REQUIRE_UNARY(floating_point, std::is_floating_point);
STAN_ADD_REQUIRE_UNARY_SCALAR(floating_point, std::is_floating_point);
STAN_ADD_REQUIRE_UNARY_VALUE(floating_point, std::is_floating_point);
STAN_ADD_REQUIRE_UNARY(index, std::is_integral);
STAN_ADD_REQUIRE_UNARY_SCALAR(index, std::is_integral);
STAN_ADD_REQUIRE_UNARY_VALUE(index, std::is_integral);




/** @}*/

/**
 * Requires for containers
 */

/** \addtogroup require_container_types
 *  @{
 */



/**
 * Check a templated type to see if it and its inner type pass a conditional
 * test.
 * @tparam ContainerCheck Templated struct or alias that wraps a static constant
 * value called type. Used to check the container satisfies a particular type
 * check. used like template <typename T, require_container_vt<is_std_vector,
 * is_var, T>...>
 */
template <template <class...> class ContainerCheck,
          template <class...> class TypeCheck, class... Check>
using require_container_vt = require_t<
    container_value_type_check_base<ContainerCheck, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_std_vector_vt
    = require_t<is_std_vector_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_std_vector_vt
    = require_not_t<is_std_vector_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_std_vector_vt
    = require_any_t<is_std_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_std_vector_vt
    = require_any_not_t<is_std_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_std_vector_vt
    = require_all_t<is_std_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_std_vector_vt
    = require_all_not_t<is_std_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_vector_vt = require_t<is_vector_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_vector_vt
    = require_not_t<is_vector_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_vector_vt
    = require_any_t<is_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_vector_vt
    = require_any_not_t<is_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vector_vt
    = require_all_t<is_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_vector_vt
    = require_all_not_t<is_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_vector_like_vt
    = require_t<is_vector_like_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_vector_like_vt
    = require_not_t<is_vector_like_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_vector_like_vt
    = require_any_t<is_vector_like_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_vector_like_vt
    = require_any_not_t<is_vector_like_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vector_like_vt
    = require_all_t<is_vector_like_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_vector_like_vt
    = require_all_not_t<is_vector_like_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_eigen_vt = require_t<is_eigen_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_vt
    = require_not_t<is_eigen_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_vt
    = require_any_t<is_eigen_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_eigen_vt
    = require_any_not_t<is_eigen_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_vt
    = require_all_t<is_eigen_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_eigen_vt
    = require_all_not_t<is_eigen_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_eigen_vector_vt
    = require_t<is_eigen_vector_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_vector_vt
    = require_not_t<is_eigen_vector_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_vector_vt
    = require_any_t<is_eigen_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_eigen_vector_vt
    = require_any_not_t<is_eigen_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_vector_vt
    = require_all_t<is_eigen_vector_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_eigen_vector_vt
    = require_all_not_t<is_eigen_vector_value_check<TypeCheck, Check>...>;

/**
 * Container Scalar check
 */

/**
 * Check a templated type to see if it and its inner type pass a conditional
 * test.
 * @tparam ContainerCheck Templated struct or alias that wraps a static constant
 * scalar called type. Used to check the container satisfies a particular type
 * check. Used like template <typename T, require_container_st<is_std_vector,
 * is_var, T>...>
 */
template <template <class...> class ContainerCheck,
          template <class...> class TypeCheck, class... Check>
using require_container_st = require_t<
    container_scalar_type_check_base<ContainerCheck, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_std_vector_st
    = require_t<is_std_vector_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_std_vector_st
    = require_not_t<is_std_vector_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_std_vector_st
    = require_any_t<is_std_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_std_vector_st
    = require_any_not_t<is_std_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_std_vector_st
    = require_all_t<is_std_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_std_vector_st
    = require_all_not_t<is_std_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_vector_st
    = require_t<is_vector_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_vector_st
    = require_not_t<is_vector_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_vector_st
    = require_any_t<is_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_vector_st
    = require_any_not_t<is_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vector_st
    = require_all_t<is_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_vector_st
    = require_all_not_t<is_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_vector_like_st
    = require_t<is_vector_like_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_vector_like_st
    = require_not_t<is_vector_like_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_vector_like_st
    = require_any_t<is_vector_like_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_vector_like_st
    = require_any_not_t<is_vector_like_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vector_like_st
    = require_all_t<is_vector_like_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_vector_like_st
    = require_all_not_t<is_vector_like_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_eigen_st = require_t<is_eigen_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_st
    = require_not_t<is_eigen_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_st
    = require_any_t<is_eigen_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_eigen_st
    = require_any_not_t<is_eigen_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_st
    = require_all_t<is_eigen_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_eigen_st
    = require_all_not_t<is_eigen_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_eigen_vector_st
    = require_t<is_eigen_vector_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_vector_st
    = require_not_t<is_eigen_vector_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_vector_st
    = require_any_t<is_eigen_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_eigen_vector_st
    = require_any_not_t<is_eigen_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_vector_st
    = require_all_t<is_eigen_vector_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_eigen_vector_st
    = require_all_not_t<is_eigen_vector_scalar_check<TypeCheck, Check>...>;

template <typename Row, typename Col>
using require_eigen_row_and_col_t = require_t<
    math::conjunction<is_eigen_row_vector<Row>, is_eigen_col_vector<Col>>>;

template <typename Row, typename Col>
using require_not_eigen_row_and_col_t = require_not_t<
    math::conjunction<is_eigen_row_vector<Row>, is_eigen_col_vector<Col>>>;

/** @}*/
}  // namespace stan
#endif
