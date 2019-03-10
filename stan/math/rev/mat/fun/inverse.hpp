#ifndef STAN_MATH_REV_MAT_FUN_INVERSE_HPP
#define STAN_MATH_REV_MAT_FUN_INVERSE_HPP

#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/rev/mat/functor/adj_jac_apply.hpp>
#include <cstddef>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace {
class inverse_op {
  int K_;
  double *y_;  // Holds the inverse, which has K_ rows & columns

 public:
  inverse_op() : K_(0), y_(nullptr) {}

  template <std::size_t size>
  Eigen::MatrixXd operator()(const std::array<bool, size> &needs_adj,
                             const Eigen::MatrixXd &m) {
    using Eigen::MatrixXd;
    K_ = m.rows();
    MatrixXd inv_m = m.inverse();
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
    Eigen::Map<const MatrixXd> y(y_, adj.rows(), adj.cols());
    MatrixXd output = -y.transpose() * adj * y.transpose();
    return std::make_tuple(output);
  }
};
}  // namespace

inline Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> inverse(
    const Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> &m) {
  check_square("inverse", "m", m);
  return adj_jac_apply<inverse_op>(m);
}

}  // namespace math
}  // namespace stan
#endif
