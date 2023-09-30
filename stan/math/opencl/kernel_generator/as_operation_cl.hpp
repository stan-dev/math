#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_AS_OPERATION_CL_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_AS_OPERATION_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/assignment_ops.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/load.hpp>
#include <stan/math/opencl/kernel_generator/scalar.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */

/**
 * Converts any valid kernel generator expression into an operation. This is an
 * overload for operations - a no-op
 * @tparam AssignOp ignored
 * @tparam T_operation type of the input operation
 * @param a an operation
 * @return operation
 */
template <assign_op_cl AssignOp = assign_op_cl::equals, typename T_operation,
          typename = std::enable_if_t<std::is_base_of<
              operation_cl_base, std::remove_reference_t<T_operation>>::value>>
inline T_operation&& as_operation_cl(T_operation&& a) {
  return std::forward<T_operation>(a);
}

/**
 * Converts any valid kernel generator expression into an operation. This is an
 * overload for scalars (arithmetic types). It wraps them into \c scalar_.
 * @tparam AssignOp ignored
 * @tparam T_scalar type of the input scalar
 * @param a scalar
 * @return \c scalar_ wrapping the input
 */
template <assign_op_cl AssignOp = assign_op_cl::equals, typename T_scalar,
          typename = require_arithmetic_t<T_scalar>,
          require_not_same_t<T_scalar, bool>* = nullptr>
inline scalar_<T_scalar> as_operation_cl(const T_scalar a) {
  return scalar_<T_scalar>(a);
}

/**
 * Converts any valid kernel generator expression into an operation. This is an
 * overload for bool scalars. It wraps them into \c scalar_<char> as \c bool can
 * not be used as a type of a kernel argument.
 * @tparam AssignOp ignored
 * @param a scalar
 * @return \c scalar_<char> wrapping the input
 */
template <assign_op_cl AssignOp = assign_op_cl::equals>
inline scalar_<char> as_operation_cl(const bool a) {
  return scalar_<char>(a);
}

/**
 * Converts any valid kernel generator expression into an operation. This is an
 * overload for \c matrix_cl. It wraps them into into \c load_.
 * @tparam AssignOp an optional `assign_op_cl` that dictates whether the object
 *  is assigned using standard or compound assign.
 * @tparam T_matrix_cl \c matrix_cl
 * @param a \c matrix_cl
 * @return \c load_ wrapping the input
 */
template <assign_op_cl AssignOp = assign_op_cl::equals, typename T_matrix_cl,
          typename = require_any_t<is_matrix_cl<T_matrix_cl>,
                                   is_arena_matrix_cl<T_matrix_cl>>>
inline load_<T_matrix_cl, AssignOp> as_operation_cl(T_matrix_cl&& a) {
  return load_<T_matrix_cl, AssignOp>(std::forward<T_matrix_cl>(a));
}

/**
 * Type that results when converting any valid kernel generator expression into
 * operation. If a function accepts a forwarding reference T&& a, the result of
 * as_operation_cl(a) should be stored in a variable of type
 * as_operation_cl_t<T>. If the return value of \c as_operation_cl() would be a
 * rvalue reference, the reference is removed, so that a variable of this type
 * actually stores the value.
 * @tparam T a `matrix_cl` or `Scalar` type
 * @tparam AssignOp an optional `assign_op_cl` that dictates whether the object
 *  is assigned using standard or compound assign.
 */
template <typename T, assign_op_cl AssignOp = assign_op_cl::equals>
using as_operation_cl_t
    = std::conditional_t<std::is_lvalue_reference<T>::value,
                         decltype(as_operation_cl<AssignOp>(std::declval<T>())),
                         std::remove_reference_t<decltype(
                             as_operation_cl<AssignOp>(std::declval<T>()))>>;

/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
