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
  W_root.reserve(Eigen::VectorXi::Constant(W_root.cols(), block_size));

  // No block operation available for sparse matrices, so we have to loop.
  // See https://eigen.tuxfamily.org/dox/group__TutorialSparse.html#title7
  for (Eigen::Index k = 0; k < W.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(W, k); it; ++it) {
      W_root.insert(it.row(), it.col()) = std::sqrt(it.value());
    }
  }
  W_root.makeCompressed();

  return W_root;
}

}  // namespace math
}  // namespace stan

#endif
