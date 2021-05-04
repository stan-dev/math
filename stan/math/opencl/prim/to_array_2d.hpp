#ifndef STAN_MATH_OPENCL_PRIM_TO_ARRAY_2D_HPP
#define STAN_MATH_OPENCL_PRIM_TO_ARRAY_2D_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns input matrix converted into a nested std vector. With matrix_cl that
 * is the same type.
 *
 * @tparam T_x type of the matrix
 *
 * @param x matrix
 * @return the vector representation of the input
 */
template <typename T_x,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_x>* = nullptr>
inline T_x to_array_2d(T_x&& x) {
  return std::forward<T_x>(x);
}

}  // namespace math
}  // namespace stan
#endif
#endif
