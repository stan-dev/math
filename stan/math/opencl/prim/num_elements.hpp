#ifndef STAN_MATH_OPENCL_PRIM_NUM_ELEMENTS_HPP
#define STAN_MATH_OPENCL_PRIM_NUM_ELEMENTS_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Returns the number of the elements of a `matrix_cl` or
 * `var_value<matrix_cl<T>>`.
 * @param m input to determine size of
 * @return number of elements in m
 */
template <typename T,
          require_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
size_t num_elements(const T& m) {
  return size(m);
}

}  // namespace math
}  // namespace stan

#endif
#endif
