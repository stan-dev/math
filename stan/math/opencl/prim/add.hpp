#ifndef STAN_MATH_OPENCL_PRIM_ADD_HPP
#define STAN_MATH_OPENCL_PRIM_ADD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {
/** \ingroup opencl
 * Computes the sum of the input matrices.
 *
 * @param A first matrix
 * @param B second matrix
 * @return the product of the first and second matrix
 *
 * @throw <code>std::invalid_argument</code> if the
 *   number of columns in A and rows in B do not match
 */
template <typename T_a, typename T_b,
          typename = require_all_valid_expressions_and_none_scalar_t<T_a, T_b>>
inline addition_<as_operation_cl_t<T_a>, as_operation_cl_t<T_b>> add(
    T_a&& a, T_b&& b) {  // NOLINT
  return {as_operation_cl(std::forward<T_a>(a)),
          as_operation_cl(std::forward<T_b>(b))};
}

template <typename T_a, typename T_b,
          typename = require_all_valid_expressions_and_none_scalar_t<T_a>,
          typename = require_arithmetic_t<T_b>>
inline addition_<as_operation_cl_t<T_a>, scalar_<T_b>> add(
    T_a&& a, T_b&& b) {  // NOLINT
  return {as_operation_cl(std::forward<T_a>(a)),
          as_operation_cl(b)};
}

template <typename T_a, typename T_b,
          typename = require_arithmetic_t<T_a>,
          typename = require_all_valid_expressions_and_none_scalar_t<T_b>>
inline addition_<scalar_<T_a>, as_operation_cl_t<T_b>> add(
    T_a&& a, T_b&& b) {  // NOLINT
  return {as_operation_cl(a),
          as_operation_cl(std::forward<T_b>(b))};
}
}  // namespace math
}  // namespace stan
#endif
#endif
