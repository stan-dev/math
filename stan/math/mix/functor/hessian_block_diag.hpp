#ifndef STAN_MATH_MIX_FUNCTOR_HESSIAN_BLOCK_DIAG_HPP
#define STAN_MATH_MIX_FUNCTOR_HESSIAN_BLOCK_DIAG_HPP

#include <stan/math/mix/functor/hessian_times_vector.hpp>
#include <Eigen/Sparse>

namespace stan {
namespace math {


/**
 * Returns a block diagonal Hessian by computing the relevant directional
 * derivatives and storing them in a matrix.
 * For m the size of each block, the operations const m calls to
 * hessian_times_vector, that is m forward sweeps and m reverse sweeps.
 * @tparam F
 * @tparam Eta
 * @tparam Args
 * @param f
 * @param x
 * @param eta
 * @param hessian_block_size
 * @param args
 */
template <typename F, typename Eta, typename... Args,
          require_eigen_t<Eta>* = nullptr>
inline Eigen::SparseMatrix<double> hessian_block_diag(
    F&& f, const Eigen::VectorXd& x, const Eta& eta,
    const Eigen::Index hessian_block_size, Args&&... args) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  const Eigen::Index x_size = x.size();
  Eigen::SparseMatrix<double> H(x_size, x_size);
  H.reserve(Eigen::VectorXi::Constant(x_size, hessian_block_size));
  VectorXd v(x_size);
  Eigen::Index n_blocks = x_size / hessian_block_size;
  for (Eigen::Index i = 0; i < hessian_block_size; ++i) {
    v.setZero();
    for (Eigen::Index j = i; j < x_size; j += hessian_block_size) {
      v(j) = 1;
    }
    VectorXd Hv = hessian_times_vector(f, x, eta, v, args...);
    for (int j = 0; j < n_blocks; ++j) {
      for (int k = 0; k < hessian_block_size; ++k) {
        H.insert(k + j * hessian_block_size, i + j * hessian_block_size)
            = Hv(k + j * hessian_block_size);
      }
    }
  }
  return H;
}

}  // namespace math
}  // namespace stan

#endif
