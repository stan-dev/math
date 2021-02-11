#ifndef STAN_MATH_OPENCL_REV_AS_COLUMN_VECTOR_OR_SCALAR_HPP
#define STAN_MATH_OPENCL_REV_AS_COLUMN_VECTOR_OR_SCALAR_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Converts kernel generator expression row or column vector to a column vector.
 *
 * @tparam T kernel generator expression.
 * @param m Specified input.
 * @return input converted to a column vector.
 */
template <typename T,
          require_all_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
          require_any_var_t<T>* = nullptr>
inline auto as_column_vector_or_scalar(const T& m) {
  return m.as_column_vector_or_scalar();
}

}  // namespace math
}  // namespace stan

#endif
#endif
