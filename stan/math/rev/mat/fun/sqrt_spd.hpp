#ifndef STAN_MATH_REV_MAT_FUN_SQRT_SPD_HPP
#define STAN_MATH_REV_MAT_FUN_SQRT_SPD_HPP

#include <cstddef>
#include <stan/math/prim/mat/fun/sqrt_spd.hpp>
#include <stan/math/rev/mat/functor/adj_jac_apply.hpp>
#include <tuple>
#include <unsupported/Eigen/KroneckerProduct>
#include <vector>

namespace stan {
namespace math {

namespace {
class sqrt_spd_op {
  int K_;
  double *y_;  // Holds the output, which has K_^2 elements

 public:
  sqrt_spd_op() : K_(0), y_(nullptr) {}

  template <std::size_t size>
  Eigen::MatrixXd operator()(const std::array<bool, size> &needs_adj,
                             const Eigen::MatrixXd &m) {
    K_ = m.rows();
    Eigen::MatrixXd sqrt_m = sqrt_spd(m);
    int KK = K_ * K_;
    y_ = ChainableStack::instance().memalloc_.alloc_array<double>(KK);
    for (int k = 0; k < KK; ++k)
      y_[k] = sqrt_m.coeff(k);
    return sqrt_m;
  }

  //  https://math.stackexchange.com/questions/540361/
  //  derivative-or-differential-of-symmetric-square-root-of-a-matrix

  template <std::size_t size>
  std::tuple<Eigen::MatrixXd> multiply_adjoint_jacobian(
      const std::array<bool, size> &needs_adj,
      const Eigen::MatrixXd &adj) const {
    using Eigen::MatrixXd;
    using Eigen::kroneckerProduct;
    Eigen::MatrixXd output(K_, K_);
    Eigen::Map<const Eigen::VectorXd> map_input(adj.data(), adj.size());
    Eigen::Map<Eigen::VectorXd> map_output(output.data(), output.size());
    Eigen::Map<MatrixXd> sqrt_m(y_, K_, K_);

    map_output = (kroneckerProduct(sqrt_m, MatrixXd::Identity(K_, K_))
                  + kroneckerProduct(MatrixXd::Identity(K_, K_), sqrt_m))
                     // the above is symmetric so can skip the transpose
                     .ldlt()
                     .solve(map_input);

    return std::make_tuple(output);
  }
};
}  // namespace

inline Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> sqrt_spd(
    const Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> &m) {
  return adj_jac_apply<sqrt_spd_op>(m);
}

}  // namespace math
}  // namespace stan
#endif
