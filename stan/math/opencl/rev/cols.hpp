#ifndef STAN_MATH_OPENCL_REV_COLS_HPP
#define STAN_MATH_OPENCL_REV_COLS_HPP
#ifdef STAN_OPENCL

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta.hpp>
namespace stan {
namespace math {

/**
 * Returns the number of columns of a `var_vlaue<matrix_cl>`.
 * @param a `var_value` to determine the number of columns
 * @return number of columns
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
size_t cols(const var_value<T>& a) {
  return a.vi_->cols();
}

}  // namespace math
}  // namespace stan

#endif
#endif
