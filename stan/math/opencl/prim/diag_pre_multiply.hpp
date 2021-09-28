#ifndef STAN_MATH_OPENCL_PRIM_DIAG_PRE_MULTIPLY_HPP
#define STAN_MATH_OPENCL_PRIM_DIAG_PRE_MULTIPLY_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Return the product of the diagonal matrix formed from the vector
 * or row_vector and a matrix.
 *
 * @tparam T1 type of the vector/row_vector
 * @tparam T2 type of the matrix
 * @param m1 input vector/row_vector
 * @param m2 input matrix
 *
 * @return product of the diagonal matrix formed from the
 * vector or row_vector and a matrix.
 */
template <typename T1, typename T2,
          require_all_kernel_expressions_and_none_scalar_t<T1, T2>* = nullptr>
inline auto diag_pre_multiply(const T1& m1, const T2& m2) {
  check_size_match("diag_pre_multiply", "m1.size()", m1.size(), "m2.rows()",
                   m2.rows());
  check_vector("diag_pre_multiply", "m1", m1);
  // we need to eval - the branches would otherwise return different types
  if (m1.cols() == 1) {
    return eval(elt_multiply(rowwise_broadcast(m1), m2));
  } else {
    return eval(elt_multiply(rowwise_broadcast(transpose(m1)), m2));
  }
}
}  // namespace math
}  // namespace stan

#endif
#endif
