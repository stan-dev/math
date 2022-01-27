#ifndef STAN_MATH_LAPLACE_BLOCK_MATRIX_SQRT_HPP
#define STAN_MATH_LAPLACE_BLOCK_MATRIX_SQRT_HPP

#include <stan/math/mix.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the matrix square-root for a block diagonal matrix.
 */
inline Eigen::SparseMatrix<double> block_matrix_sqrt(
    const Eigen::SparseMatrix<double>& W, const Eigen::Index block_size) {
  const Eigen::Index n_block = W.cols() / block_size;
  Eigen::MatrixXd local_block(block_size, block_size);
  Eigen::SparseMatrix<double> W_root(W.rows(), W.cols());
  return W.sqrt();
}

}  // namespace math
}  // namespace stan

#endif
