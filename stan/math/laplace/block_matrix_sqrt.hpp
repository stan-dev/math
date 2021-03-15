#ifndef STAN_MATH_LAPLACE_BLOCK_MATRIX_SQRT_HPP
#define STAN_MATH_LAPLACE_BLOCK_MATRIX_SQRT_HPP

#include <stan/math/prim/fun/cholesky_decompose.hpp>

#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>

namespace stan {
namespace math {

/**
 * Return the matrix square-root for a block diagonal matrix.
 */
 Eigen::SparseMatrix<double>
 block_matrix_sqrt(Eigen::SparseMatrix<double> W,
                   int block_size) {
  int n_block = W.cols() / block_size;
  Eigen::MatrixXd local_block(block_size, block_size);
  Eigen::MatrixXd local_block_sqrt(block_size, block_size);
  Eigen::SparseMatrix<double> W_root(W.rows(), W.cols());
  W_root.reserve(Eigen::VectorXi::Constant(W_root.cols(), block_size));

  // No block operation available for sparse matrices, so we have to loop.
  // See https://eigen.tuxfamily.org/dox/group__TutorialSparse.html#title7
  for (int i = 0; i < n_block; i++) {
    std::cout << "block number: " << i << std::endl;
    for (int j = 0; j < block_size; j++) {
      for (int k = 0; k < block_size; k++) {
        local_block(j, k) = W.coeffRef(i * block_size + j, i * block_size + k);
      }
    }

    local_block_sqrt = local_block.sqrt();
    // local_block_sqrt = cholesky_decompose(local_block);

    for (int j = 0; j < block_size; j++) {
      for (int k = 0; k < block_size; k++) {
        W_root.insert(i * block_size + j, i * block_size + k)
          = local_block_sqrt(j, k);
      }
    }
  }

  W_root.makeCompressed();

  return W_root;
}


}  // namespace math
}  // namespace stan

#endif
