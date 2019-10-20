#ifndef STAN_MATH_OPENCL_IS_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_IS_MATRIX_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {

namespace math {
// Dummy class to instantiate matrix_cl to enable for specific types.
template <typename T, typename = void>
class matrix_cl {
 public:
  using Scalar = T;
  using type = T;
};
}  // namespace math

namespace internal {

/**
 * This underlying implimentation is used when the type is not an std vector.
 */
template <typename T>
struct is_matrix_cl_impl : std::false_type {};

/**
 * This specialization implimentation has a static member named value when the
 * template type is an std vector.
 */
template <typename... Args>
struct is_matrix_cl_impl<stan::math::matrix_cl<Args...>> : std::true_type {};

}  // namespace internal

template <typename T, typename = void>
struct is_matrix_cl : std::false_type {};

/**
 * Checks if the decayed type of T is a matrix_cl.
 */
template <typename T>
struct is_matrix_cl<
    T, std::enable_if_t<internal::is_matrix_cl_impl<std::decay_t<T>>::value>>
    : std::true_type {};

template <typename T>
using require_matrix_cl_t = require_t<is_matrix_cl<T>>;

template <typename T>
using require_not_matrix_cl_t = require_not_t<is_matrix_cl<T>>;

template <typename... Types>
using require_all_matrix_cl_t = require_all_t<is_matrix_cl<Types>...>;

template <typename... Types>
using require_any_matrix_cl_t = require_any_t<is_matrix_cl<Types>...>;

template <typename... Types>
using require_all_not_matrix_cl_t = require_all_not_t<is_matrix_cl<Types>...>;

template <typename... Types>
using require_any_not_matrix_cl_t = require_any_not_t<is_matrix_cl<Types>...>;

/**
 * matrix_cl
 */
template <template <class...> class TypeCheck, class... Check>
struct is_matrix_cl_value_check
    : container_value_type_check_base<is_matrix_cl, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_matrix_cl_vt
    = require_t<is_matrix_cl_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_matrix_cl_vt
    = require_not_t<is_matrix_cl_value_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_matrix_cl_vt
    = require_any_t<is_matrix_cl_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_matrix_cl_vt
    = require_any_not_t<is_matrix_cl_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_matrix_cl_vt
    = require_all_t<is_matrix_cl_value_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_matrix_cl_vt
    = require_all_not_t<is_matrix_cl_value_check<TypeCheck, Check>...>;

// Scalar Check

/**
 * matrix_cl
 */
template <template <class...> class TypeCheck, class... Check>
struct is_matrix_cl_scalar_check
    : container_scalar_type_check_base<is_matrix_cl, TypeCheck, Check...> {};

template <template <class...> class TypeCheck, class... Check>
using require_matrix_cl_st
    = require_t<is_matrix_cl_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_not_matrix_cl_st
    = require_not_t<is_matrix_cl_scalar_check<TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_any_matrix_cl_st
    = require_any_t<is_matrix_cl_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_any_not_matrix_cl_st
    = require_any_not_t<is_matrix_cl_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_matrix_cl_st
    = require_all_t<is_matrix_cl_scalar_check<TypeCheck, Check>...>;

template <template <class...> class TypeCheck, class... Check>
using require_all_not_matrix_cl_st
    = require_all_not_t<is_matrix_cl_scalar_check<TypeCheck, Check>...>;

}  // namespace stan
#endif
#endif
