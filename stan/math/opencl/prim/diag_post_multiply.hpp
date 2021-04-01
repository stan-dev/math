#ifndef STAN_MATH_OPENCL_PRIM_DIAG_POST_MULTIPLY_HPP
#define STAN_MATH_OPENCL_PRIM_DIAG_POST_MULTIPLY_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Return the product of a matrix and the diagonal matrix formed from the vector
 * or row_vector.
 *
 * @tparam T1 type of the matrix
 * @tparam T2 type of the vector/row_vector
 * @param m1 input matrix
 * @param m2 input vector/row_vector
 *
 * @return product of a matrix and the diagonal matrix formed from the
 * vector or row_vector.
 */
template <typename T1, typename T2,
          require_all_kernel_expressions_and_none_scalar_t<T1, T2>* = nullptr>
inline auto diag_post_multiply(const T1& m1, const T2& m2) {
  check_size_match("diag_post_multiply", "m1.cols()", m1.cols(), "m2.size()",
                   m2.size());
  check_vector("diag_post_multiply", "m2", m2);
  // we need to eval - the branches would otherwise return different types
  if (m2.rows() == 1) {
    return eval(elt_multiply(m1, colwise_broadcast(m2)));
  } else {
    return eval(elt_multiply(m1, colwise_broadcast(transpose(m2))));
  }
}
}  // namespace math
}  // namespace stan

#endif
#endif
