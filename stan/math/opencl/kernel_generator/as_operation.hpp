#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_AS_OPERATION_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_AS_OPERATION_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/operation.hpp>
#include <stan/math/opencl/kernel_generator/load.hpp>
#include <stan/math/opencl/kernel_generator/scalar.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Converts any valid kernel generator expression into an operation. This is an
 * overload for operations - a no-op
 * @tparam T_operation type of the input operation
 * @param an operation
 * @return operation
 */
template <typename T_operation, typename = std::enable_if_t<std::is_base_of<
                          operation_base, std::remove_reference_t<T_operation>>::value>>
inline T_operation&& as_operation(T_operation&& a) {
  return std::forward<T_operation>(a);
}

/**
 * Converts any valid kernel generator expression into an operation. This is an
 * overload for scalars (atihmetic types). It wraps them into \c scalar__.
 * @tparam T_scalar type of the input scalar
 * @param a scalar
 * @return \c scalar__ wrapping the input
 */
template <typename T_scalar, typename = enable_if_arithmetic<T_scalar>>
inline scalar__<T_scalar> as_operation(const T_scalar a) {
  return scalar__<T_scalar>(a);
}

/**
 * Converts any valid kernel generator expression into an operation. This is an
 * overload for \c matrix_cl. It wraps them into into \c load__.
 * @tparam T_matrix_cl \c matrix_cl
 * @param a \c matrix_cl
 * @return \c load__ wrapping the input
 */
template <typename T_matrix_cl, typename = std::enable_if_t<std::is_base_of<
                          matrix_cl<typename std::remove_reference_t<T_matrix_cl>::type>,
                          std::remove_reference_t<T_matrix_cl>>::value>>
inline load__<T_matrix_cl> as_operation(T_matrix_cl&& a) {
  return load__<T_matrix_cl>(std::forward<T_matrix_cl>(a));
}

/**
 * Type that results when converting any valid kernel generator expression into
 * operation. If a function accepts a forwarding reference T&& a, the result of
 * as_operation(a) should be stored in a variable of type as_operation_t<T>. If
 * the return value of \c as_operation() would be a rvalue reference, the
 * reference is removed, so that a variable of this type actually stores the value.
 */
template <typename T>
using as_operation_t = std::conditional_t<
    std::is_lvalue_reference<T>::value,
    decltype(as_operation(std::declval<T>())),
    std::remove_reference_t<decltype(as_operation(std::declval<T>()))>>;

}  // namespace math
}  // namespace stan

#endif
#endif
