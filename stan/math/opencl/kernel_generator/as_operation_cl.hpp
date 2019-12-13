#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_AS_OPERATION_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_AS_OPERATION_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/load.hpp>
#include <stan/math/opencl/kernel_generator/scalar.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/is_matrix_cl.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Converts any valid kernel generator expression into an operation. This is an
 * overload for operations - a no-op
 * @tparam T_operation type of the input operation
 * @param a an operation
 * @return operation
 */
template <typename T_operation,
          typename = std::enable_if_t<std::is_base_of<
              operation_cl_base, std::remove_reference_t<T_operation>>::value>>
inline T_operation&& as_operation_cl(T_operation&& a) {
  return std::forward<T_operation>(a);
}

/**
 * Converts any valid kernel generator expression into an operation. This is an
 * overload for scalars (arithmetic types). It wraps them into \c scalar_.
 * @tparam T_scalar type of the input scalar
 * @param a scalar
 * @return \c scalar_ wrapping the input
 */
template <typename T_scalar, typename = require_arithmetic_t<T_scalar>>
inline scalar_<T_scalar> as_operation_cl(const T_scalar a) {
  return scalar_<T_scalar>(a);
}

/**
 * Converts any valid kernel generator expression into an operation. This is an
 * overload for \c matrix_cl. It wraps them into into \c load_.
 * @tparam T_matrix_cl \c matrix_cl
 * @param a \c matrix_cl
 * @return \c load_ wrapping the input
 */
template <typename T_matrix_cl, typename = require_matrix_cl_t<T_matrix_cl>>
inline load_<T_matrix_cl> as_operation_cl(T_matrix_cl&& a) {
  return load_<T_matrix_cl>(std::forward<T_matrix_cl>(a));
}

/**
 * Type that results when converting any valid kernel generator expression into
 * operation. If a function accepts a forwarding reference T&& a, the result of
 * as_operation_cl(a) should be stored in a variable of type
 * as_operation_cl_t<T>. If the return value of \c as_operation_cl() would be a
 * rvalue reference, the reference is removed, so that a variable of this type
 * actually stores the value.
 */
template <typename T>
using as_operation_cl_t = std::conditional_t<
    std::is_lvalue_reference<T>::value,
    decltype(as_operation_cl(std::declval<T>())),
    std::remove_reference_t<decltype(as_operation_cl(std::declval<T>()))>>;

}  // namespace math
}  // namespace stan

#endif
#endif
