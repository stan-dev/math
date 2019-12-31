#ifndef STAN_MATH_PRIM_MAT_FUN_ELT_MULTIPLY_HPP
#define STAN_MATH_PRIM_MAT_FUN_ELT_MULTIPLY_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise multiplication of the specified
 * matrices.
 *
 * @tparam T1 type of elements in first matrix
 * @tparam T2 type of elements in second matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m1 First matrix
 * @param m2 Second matrix
 * @return Elementwise product of matrices.
 */
template <typename T1, typename T2, int R, int C>
Eigen::Matrix<return_type_t<T1, T2>, R, C> elt_multiply(
    const Eigen::Matrix<T1, R, C>& m1, const Eigen::Matrix<T2, R, C>& m2) {
  check_matching_dims("elt_multiply", "m1", m1, "m2", m2);
  return m1.cwiseProduct(m2);
}

}  // namespace math
}  // namespace stan

#endif
