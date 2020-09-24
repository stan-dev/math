#ifndef STAN_MATH_OPENCL_PRIM_SIZE_HPP
#define STAN_MATH_OPENCL_PRIM_SIZE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/is_kernel_expression.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Returns the size (number of the elements) of a `matrix_cl`.
 * @param a `matric_cl` to determine size of
 * @return number of elements in a
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
int size(const T& a) {
  return a.rows() * a.cols();
}

}  // namespace math
}  // namespace stan

#endif
#endif
