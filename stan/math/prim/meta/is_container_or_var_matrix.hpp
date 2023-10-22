#ifndef STAN_MATH_PRIM_META_IS_CONTAINER_OR_VAR_MATRIX_HPP
#define STAN_MATH_PRIM_META_IS_CONTAINER_OR_VAR_MATRIX_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_container.hpp>
#include <stan/math/prim/meta/is_var_matrix.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {

/**
 * Deduces whether type is eigen matrix, standard vector, or var<Matrix>.
 * @tparam Container type to check
 */
template <typename Container>
using is_container_or_var_matrix
    = bool_constant<math::disjunction<is_container<Container>,
                                      is_var_matrix<Container>>::value>;

// STAN_ADD_REQUIRE_UNARY(container_or_var_matrix, is_container_or_var_matrix,
// general_types);
template <typename T>
using require_container_or_var_matrix_t = require_t<is_container_or_var_matrix<std::decay_t<T>>>;

template <typename T>
using require_not_container_or_var_matrix_t
    = require_not_t<is_container_or_var_matrix<std::decay_t<T>>>;

template <typename... Types>
using require_all_container_or_var_matrix_t
    = require_all_t<is_container_or_var_matrix<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_container_or_var_matrix_t
    = require_any_t<is_container_or_var_matrix<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_container_or_var_matrix_t
    = require_all_not_t<is_container_or_var_matrix<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_container_or_var_matrix_t
    = require_any_not_t<is_container_or_var_matrix<std::decay_t<Types>>...>;

  
  
//STAN_ADD_REQUIRE_CONTAINER(container_or_var_matrix, is_container_or_var_matrix, general_types);
template <template <class...> class TypeCheck, class... Check>
using require_container_or_var_matrix_vt = require_t<
    container_type_check_base<is_container_or_var_matrix, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_container_or_var_matrix_vt = require_not_t<
    container_type_check_base<is_container_or_var_matrix, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_container_or_var_matrix_vt = require_any_t<
    container_type_check_base<is_container_or_var_matrix, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_container_or_var_matrix_vt = require_any_not_t<
    container_type_check_base<is_container_or_var_matrix, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_container_or_var_matrix_vt = require_all_t<
    container_type_check_base<is_container_or_var_matrix, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_container_or_var_matrix_vt = require_all_not_t<
    container_type_check_base<is_container_or_var_matrix, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_container_or_var_matrix_st = require_t<
    container_type_check_base<is_container_or_var_matrix, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_container_or_var_matrix_st = require_not_t<
    container_type_check_base<is_container_or_var_matrix, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_container_or_var_matrix_st = require_any_t<
    container_type_check_base<is_container_or_var_matrix, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_container_or_var_matrix_st = require_any_not_t<
    container_type_check_base<is_container_or_var_matrix, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_container_or_var_matrix_st = require_all_t<
    container_type_check_base<is_container_or_var_matrix, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_container_or_var_matrix_st = require_all_not_t<
    container_type_check_base<is_container_or_var_matrix, scalar_type_t, TypeCheck, Check>...>;


}  // namespace stan

#endif
