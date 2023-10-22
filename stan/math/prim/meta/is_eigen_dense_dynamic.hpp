#ifndef STAN_MATH_PRIM_META_IS_EIGEN_DENSE_DYNAMIC_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_DENSE_DYNAMIC_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_eigen_matrix.hpp>
#include <stan/math/prim/meta/is_eigen_dense_base.hpp>
#include <stan/math/prim/meta/is_base_pointer_convertible.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

/**
 * Checks whether type T is derived from Eigen::DenseBase and has dynamic rows
 * and columns. If true this will have a static member function named value with
 * a type of true, else value is false.
 * @tparam T Type to check
 * @ingroup type_trait
 */
template <typename T>
using is_eigen_dense_dynamic = stan::internal::is_eigen_matrix_dynamic_impl<
    std::decay_t<T>, stan::is_eigen_dense_base<std::decay_t<T>>::value>;

//STAN_ADD_REQUIRE_UNARY(eigen_dense_dynamic, is_eigen_dense_dynamic, require_eigens_types);
template <typename T>
using require_eigen_dense_dynamic_t = require_t<is_eigen_dense_dynamic<std::decay_t<T>>>;

template <typename T>
using require_not_eigen_dense_dynamic_t
    = require_not_t<is_eigen_dense_dynamic<std::decay_t<T>>>;

template <typename... Types>
using require_all_eigen_dense_dynamic_t
    = require_all_t<is_eigen_dense_dynamic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_eigen_dense_dynamic_t
    = require_any_t<is_eigen_dense_dynamic<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_eigen_dense_dynamic_t
    = require_all_not_t<is_eigen_dense_dynamic<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_eigen_dense_dynamic_t
    = require_any_not_t<is_eigen_dense_dynamic<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_CONTAINER(eigen_dense_dynamic, is_eigen_dense_dynamic, require_eigens_types);
template <template <class...> class TypeCheck, class... Check>
using require_eigen_dense_dynamic_vt = require_t<
    container_type_check_base<is_eigen_dense_dynamic, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_dense_dynamic_vt = require_not_t<
    container_type_check_base<is_eigen_dense_dynamic, value_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_dense_dynamic_vt = require_any_t<
    container_type_check_base<is_eigen_dense_dynamic, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_eigen_dense_dynamic_vt = require_any_not_t<
    container_type_check_base<is_eigen_dense_dynamic, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_dense_dynamic_vt = require_all_t<
    container_type_check_base<is_eigen_dense_dynamic, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_eigen_dense_dynamic_vt = require_all_not_t<
    container_type_check_base<is_eigen_dense_dynamic, value_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_eigen_dense_dynamic_st = require_t<
    container_type_check_base<is_eigen_dense_dynamic, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_dense_dynamic_st = require_not_t<
    container_type_check_base<is_eigen_dense_dynamic, scalar_type_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_dense_dynamic_st = require_any_t<
    container_type_check_base<is_eigen_dense_dynamic, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_eigen_dense_dynamic_st = require_any_not_t<
    container_type_check_base<is_eigen_dense_dynamic, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_dense_dynamic_st = require_all_t<
    container_type_check_base<is_eigen_dense_dynamic, scalar_type_t, TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_eigen_dense_dynamic_st = require_all_not_t<
    container_type_check_base<is_eigen_dense_dynamic, scalar_type_t, TypeCheck, Check>...>;


}  // namespace stan

#endif
