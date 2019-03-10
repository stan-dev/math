#ifndef STAN_MATH_REV_MAT_FUN_INV_SQRT_SPD_HPP
#define STAN_MATH_REV_MAT_FUN_INV_SQRT_SPD_HPP

#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/mat/functor/adj_jac_apply.hpp>
#include <unsupported/Eigen/KroneckerProduct>
#include <cstddef>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace {
class inv_sqrt_spd_op {
  int K_;
  double *J_;  // Holds the Jacobian, which has K_^2 rows & columns

 public:
  inv_sqrt_spd_op() : K_(0), J_(nullptr) {}

  template <std::size_t size>
  Eigen::MatrixXd operator()(const std::array<bool, size> &needs_adj,
                             const Eigen::MatrixXd &m) {
    using Eigen::kroneckerProduct;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    K_ = m.rows();
    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(m);
    MatrixXd V = solver.eigenvectors();
    VectorXd sqrt_ev = solver.eigenvalues().cwiseSqrt();
    MatrixXd inv_sqrt_m
        = V * sqrt_ev.cwiseInverse().asDiagonal() * V.transpose();
    MatrixXd sqrt_m = V * sqrt_ev.asDiagonal() * V.transpose();
    int K2 = K_ * K_;
    VectorXd d(K2);
    //  https://math.stackexchange.com/questions/540361/
    //  and then simplifying the inverse of a Kronecker sum
    //  and multiplying by the Jacobian of the inverse of a SPD matrix
    int pos = 0;
    for (int i = 0; i < K_; i++)
      for (int j = 0; j < K_; j++)
        d.coeffRef(pos++) = sqrt_ev.coeff(i) * sqrt_ev.coeff(j)
                            * (sqrt_ev.coeff(i) + sqrt_ev.coeff(j));
    MatrixXd U = kroneckerProduct(V, V);
    MatrixXd J = -U * d.cwiseInverse().asDiagonal() * U.transpose();
    int K4 = K2 * K2;
    J_ = ChainableStack::instance().memalloc_.alloc_array<double>(K4);
    for (int k = 0; k < K4; ++k)
      J_[k] = J.coeff(k);
    return inv_sqrt_m;
  }

  template <std::size_t size>
  std::tuple<Eigen::MatrixXd> multiply_adjoint_jacobian(
      const std::array<bool, size> &needs_adj,
      const Eigen::MatrixXd &adj) const {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    MatrixXd output(K_, K_);
    Map<VectorXd> map_output(output.data(), output.size());
    Map<const VectorXd> map_input(adj.data(), adj.size());
    Map<const MatrixXd> J(J_, adj.size(), adj.size());
    // J is symmetric so do not need to transpose it
    map_output = J.selfadjointView<Eigen::Lower>() * map_input;
    return std::make_tuple(output);
  }
};
}  // namespace

inline Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> inv_sqrt_spd(
    const Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> &m) {
  check_nonzero_size("inv_sqrt_spd", "m", m);
  check_symmetric("inv_sqrt_spd", "m", m);
  return adj_jac_apply<inv_sqrt_spd_op>(m);
}

}  // namespace math
}  // namespace stan
#endif
