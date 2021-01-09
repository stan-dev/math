#ifndef STAN_MATH_OPENCL_PRIM_COLUMNS_DOT_PRODUCT_HPP
#define STAN_MATH_OPENCL_PRIM_COLUMNS_DOT_PRODUCT_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_vector.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of columns of the specified matrices.
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
inline auto columns_dot_product(const T_a& a, const T_b& b) {
  using res_scal = std::common_type_t<value_type_t<T_a>, value_type_t<T_b>>;
  check_matching_sizes("columns_dot_product", "a", a, "b", b);
  matrix_cl<res_scal> res;

  if (size_zero(a, b)) {
    res = constant(res_scal(0), 1, a.cols());
    return res;
  }

  res = colwise_sum(elt_multiply(a, b));
  while (res.rows() > 1) {
    res = colwise_sum(res).eval();
  }
  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
