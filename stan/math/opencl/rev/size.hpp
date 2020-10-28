#ifndef STAN_MATH_OPENCL_REV_SIZE_HPP
#define STAN_MATH_OPENCL_REV_SIZE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/is_kernel_expression.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
namespace stan {
namespace math {

/**
 * Returns the size (number of the elements) of a `var_vlaue<matrix_cl>`.
 * @param a `var_value` to determine size of
 * @return number of elements in a
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
int size(const var_value<T>& a) {
  return a.vi_->rows() * a.vi_->cols();
}

}  // namespace math
}  // namespace stan

#endif
#endif
