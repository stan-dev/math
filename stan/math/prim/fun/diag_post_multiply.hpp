#ifndef STAN_MATH_PRIM_FUN_DIAG_POST_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_DIAG_POST_MULTIPLY_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the product of the matrix and the diagonal matrix
 * formed from the vector or row_vector.
 *
 * @tparam T1 type of the matrix
 * @tparam T2 type of the vector/row_vector
 * @param m1 input matrix
 * @param m2 input vector/row_vector
 *
 * @return product of the matrix and the diagonal matrix formed from the
 * vector or row_vector.
 */
template <typename T1, typename T2, require_eigen_t<T1>* = nullptr,
          require_eigen_vector_t<T2>* = nullptr,
          require_all_not_st_var<T1, T2>* = nullptr>
auto diag_post_multiply(const T1& m1, const T2& m2) {
  check_size_match("diag_post_multiply", "m2.size()", m2.size(), "m1.cols()",
                   m1.cols());
  return m1 * m2.asDiagonal();
}

}  // namespace math
}  // namespace stan
#endif
