#ifndef STAN_MATH_PRIM_MAT_FUN_TCROSSPROD_HPP
#define STAN_MATH_PRIM_MAT_FUN_TCROSSPROD_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of post-multiplying a matrix by its
 * own transpose.
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param M Matrix to multiply.
 * @return M times its transpose.
 */
template <int R, int C>
inline Eigen::MatrixXd tcrossprod(const Eigen::Matrix<double, R, C>& M) {
  if (M.rows() == 0) {
    return matrix_d(0, 0);
  }
  if (M.rows() == 1) {
    return M * M.transpose();
  }
  matrix_d result(M.rows(), M.rows());
  return result.setZero().selfadjointView<Eigen::Upper>().rankUpdate(M);
}

}  // namespace math
}  // namespace stan

#endif
