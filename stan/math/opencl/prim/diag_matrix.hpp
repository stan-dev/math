#ifndef STAN_MATH_OPENCL_PRIM_DIAG_MATRIX_HPP
#define STAN_MATH_OPENCL_PRIM_DIAG_MATRIX_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/add_diag.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Return a square diagonal matrix with the specified vector of
 * coefficients as the diagonal values.
 *
 * @tparam T_x type of input kernel generator expression for the
 * diagonal
 *
 * @param x input kernel generator expression for the diagonal
 *
 * @return a kernel generator expression
 */
template <typename T_x,
          require_all_kernel_expressions_and_none_scalar_t<T_x>* = nullptr>
inline auto diag_matrix(T_x&& x) {  // NOLINT
  return add_diag(constant(value_type_t<T_x>(0), x.size(), x.size()).eval(), x);
}
}  // namespace math
}  // namespace stan

#endif
#endif
