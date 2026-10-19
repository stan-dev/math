#ifndef STAN_MATH_OPENCL_PRIM_APPEND_ARRAY_HPP
#define STAN_MATH_OPENCL_PRIM_APPEND_ARRAY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/**
 * Return the concatenation of two specified vectors in the order of
 *   the arguments.
 *
 * @tparam T1 type of elements in first vector
 * @tparam T2 type of elements in second vector
 * @param x First vector
 * @param y Second vector
 * @return A vector of x and y concatenated together (in that order)
 */
template <
    typename T_x, typename T_y,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_x, T_y>* = nullptr>
inline auto append_array(T_x&& x, T_y&& y) {
  return append_row(std::forward<T_x>(x), std::forward<T_y>(y));
}
}  // namespace math
}  // namespace stan
#endif
#endif
