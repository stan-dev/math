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
 */
/*
template <typename F>
inline Eigen::MatrixXd hessian_block_diag(F&& f, const Eigen::VectorXd& x,
                              const Eigen::VectorXd& eta,
                              const Eigen::VectorXd& delta,
                              const std::vector<int>& delta_int,
                              const Eigen::Index hessian_block_size,
                              std::ostream* pstream = 0) {
 using Eigen::MatrixXd;
 using Eigen::VectorXd;

 const Eigen::Index x_size = x.size();
 VectorXd v;
 Eigen::MatrixXd H = Eigen::MatrixXd::Zero(x_size, x_size);
 const Eigen::Index n_blocks = x_size / hessian_block_size;
 for (Eigen::Index i = 0; i < hessian_block_size; ++i) {
   v.setZero();
   for (Eigen::Index j = i; j < x_size; j += hessian_block_size)
     v(j) = 1;
   VectorXd Hv = hessian_times_vector(f, x, eta, delta, delta_int, v, pstream);
   for (Eigen::Index j = 0; j < n_blocks; ++j) {
     for (Eigen::Index k = 0; k < hessian_block_size; ++k)
       H(k + j * hessian_block_size, i + j * hessian_block_size)
           = Hv(k + j * hessian_block_size);
   }
 }
 return H;
}
*/
/**
 * Overload for case where hessian is stored as a sparse matrix.
 */
template <typename F, typename Eta, require_eigen_t<Eta>* = nullptr>
inline Eigen::SparseMatrix<double> hessian_block_diag(
    F&& f, const Eigen::VectorXd& x, const Eta& eta,
    const Eigen::VectorXd& delta, const std::vector<int>& delta_int,
    const Eigen::Index hessian_block_size, std::ostream* pstream = 0) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  const Eigen::Index x_size = x.size();
  // H = MatrixXd::Zero(x_size, x_size);
  Eigen::SparseMatrix<double> H(x_size, x_size);
  // H.reserve(Eigen::VectorXi::Constant(x_size, hessian_block_size));

  VectorXd v(x_size);
  Eigen::Index n_blocks = x_size / hessian_block_size;
  for (Eigen::Index i = 0; i < hessian_block_size; ++i) {
    v.setZero();
    for (Eigen::Index j = i; j < x_size; j += hessian_block_size) {
      v(j) = 1;
    }
    VectorXd Hv = hessian_times_vector(f, x, eta, delta, delta_int, v, pstream);
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
