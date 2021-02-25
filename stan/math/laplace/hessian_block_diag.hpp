#ifndef STAN_MATH_LAPLACE_HESSIAN_BLOCK_DIAG_HPP
#define STAN_MATH_LAPLACE_HESSIAN_BLOCK_DIAG_HPP

// TODO: refine include.
#include <stan/math/mix.hpp>
#include <stan/math/laplace/hessian_times_vector.hpp>

namespace stan {
namespace math {

  /**
   * Returns a block diagonal Hessian by computing the relevant directional
   * derivatives and storing them in a matrix.
   * For m the size of each block, the operations const m calls to
   * hessian_times_vector, that is m forward sweeps and m reverse sweeps.
   */
   template <typename F>
   void hessian_block_diag(const F& f,
                           const Eigen::VectorXd& x,
                           const Eigen::VectorXd& eta,
                           const Eigen::VectorXd& delta,
                           const std::vector<int>& delta_int,
                           int m,
                           double& fx,
                           Eigen::MatrixXd& H,
                           std::ostream* pstream = 0) {
    using Eigen::VectorXd;
    using Eigen::MatrixXd;

    int x_size = x.size();
    VectorXd v;
    H = MatrixXd::Zero(x_size, x_size);
    int n_blocks = x_size / m;
    for (int i = 0; i < m; ++i) {
      v = VectorXd::Zero(x_size);
      for (int j = i; j < x_size; j += m) v(j) = 1;
      VectorXd Hv;
      hessian_times_vector(f, x, eta, delta, delta_int, v, fx, Hv, pstream);
      std::cout << "Hv: " << Hv << std::endl;
      for (int j = 0; j < n_blocks; ++j) {
        for (int k = 0; k < m; ++k) H(k + j * m, i + j * m) = Hv(k + j * m);
      }
    }
  }

}  // namespace math
}  // namespace stan

#endif
