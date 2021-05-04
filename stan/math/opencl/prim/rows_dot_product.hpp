#ifndef STAN_MATH_OPENCL_PRIM_ROWS_DOT_PRODUCT_HPP
#define STAN_MATH_OPENCL_PRIM_ROWS_DOT_PRODUCT_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_vector.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of rows of the specified matrices.
 *
 * @tparam T_a type of the first matrix
 * @tparam T_b type of the second matrix
 *
 * @param a Matrix of first vectors.
 * @param b Matrix of second vectors.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the matrices are not the same
 * size
 */
template <typename T_a, typename T_b,
          require_all_kernel_expressions_and_none_scalar_t<T_a, T_b>* = nullptr>
inline auto rows_dot_product(T_a&& a, T_b&& b) {
  check_matching_sizes("rows_dot_product", "a", a, "b", b);
  return rowwise_sum(elt_multiply(std::forward<T_a>(a), std::forward<T_b>(b)));
}

}  // namespace math
}  // namespace stan

#endif
#endif
