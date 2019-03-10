#ifndef STAN_MATH_REV_MAT_FUN_INVERSE_SPD_HPP
#define STAN_MATH_REV_MAT_FUN_INVERSE_SPD_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/rev/mat/functor/adj_jac_apply.hpp>
#include <cstddef>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace {
class inverse_spd_op {
  int K_;
  double *y_;  // Holds the inverse, which has K_ rows & columns

 public:
  inverse_spd_op() : K_(0), y_(nullptr) {}

  template <std::size_t size>
  Eigen::MatrixXd operator()(const std::array<bool, size> &needs_adj,
                             const Eigen::MatrixXd &m) {
    using Eigen::MatrixXd;
    K_ = m.rows();
    Eigen::LDLT<MatrixXd> ldlt(m);
    if (ldlt.info() != Eigen::Success)
      domain_error("invese_spd", "LDLT factor failed", "", "");
    if (!ldlt.isPositive() || ldlt.vectorD().coeff(K_ - 1) <= 0)
      domain_error("invese_spd", "matrix not positive definite", "", "");
    MatrixXd inv_m = ldlt.solve(MatrixXd::Identity(m.rows(), m.cols()));

    int K2 = K_ * K_;
    y_ = ChainableStack::instance().memalloc_.alloc_array<double>(K2);
    for (int k = 0; k < K2; ++k)
      y_[k] = inv_m.coeff(k);
    return inv_m;
  }

  template <std::size_t size>
  std::tuple<Eigen::MatrixXd> multiply_adjoint_jacobian(
      const std::array<bool, size> &needs_adj,
      const Eigen::MatrixXd &adj) const {
    using Eigen::MatrixXd;
    Eigen::Map<const MatrixXd> y(y_, K_, K_);
    MatrixXd output = y.selfadjointView<Eigen::Lower>() * -adj
                      * y.selfadjointView<Eigen::Lower>();
    return std::make_tuple(output);
  }
};
}  // namespace

inline Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> inverse_spd(
    const Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> &m) {
  check_symmetric("inverse_spd", "m", m);
  return adj_jac_apply<inverse_spd_op>(m);
}

}  // namespace math
}  // namespace stan
#endif
