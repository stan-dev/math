#ifndef STAN_MATH_REV_FUN_SIMPLEX_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_SIMPLEX_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/adj_jac_apply.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <cmath>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace internal {
template <typename T>
class simplex_constrain_op {
public:
  adj_op<T> diag_;
  adj_op<T> z_;
  explicit simplex_constrain_op(const T& x) : diag_(x.size()), z_(x.size()) {}
  /**
   * Return the simplex corresponding to the specified free vector.
   * A simplex is a vector containing values greater than or equal
   * to 0 that sum to 1.  A vector with (K-1) unconstrained values
   * will produce a simplex of size K.
   *
   * The transform is based on a centered stick-breaking process.
   *
   * @tparam size Number of adjoints to return
   * @param needs_adj Boolean indicators of if adjoints of arguments will be
   * needed
   * @param y Free vector input of dimensionality K - 1
   * @return Simplex of dimensionality K
   */
  Eigen::VectorXd operator()(const Eigen::VectorXd& y) {
    Eigen::Matrix<double, Eigen::Dynamic, 1> x(y.size() + 1);
    double stick_len(1.0);
    for (int k = 0; k < y.size(); ++k) {
      double log_N_minus_k = std::log(y.size() - k);
      z_(k) = inv_logit(y(k) - log_N_minus_k);
      diag_(k) = stick_len * z_(k) * inv_logit(log_N_minus_k - y(k));
      x(k) = stick_len * z_(k);
      stick_len -= x(k);
    }
    x(y.size()) = stick_len;
    return x;
  }

  /**
   * Compute the result of multiply the transpose of the adjoint vector times
   * the Jacobian of the simplex_constrain operator.
   *
   * @tparam size Number of adjoints to return
   * @param needs_adj Boolean indicators of if adjoints of arguments will be
   * needed
   * @param adj Eigen::VectorXd of adjoints at the output of the softmax
   * @return Eigen::VectorXd of adjoints propagated through softmax operation
   */
  auto multiply_adjoint_jacobian(const Eigen::VectorXd& adj) const {
    const auto z_size = z_.size();
    Eigen::VectorXd adj_times_jac(z_size);
    double acc = adj(z_size);

    if (z_size > 0) {
      adj_times_jac(z_size - 1) = diag_(z_size - 1) * (adj(z_size - 1) - acc);
      for (int n = z_size - 1; --n >= 0;) {
        acc = adj(n + 1) * z_(n + 1) + (1 - z_(n + 1)) * acc;
        adj_times_jac(n) = diag_(n) * (adj(n) - acc);
      }
    }
    return std::make_tuple(adj_times_jac);
  }
};
}  // namespace internal

/**
 * Return the simplex corresponding to the specified free vector.
 * A simplex is a vector containing values greater than or equal
 * to 0 that sum to 1.  A vector with (K-1) unconstrained values
 * will produce a simplex of size K.
 *
 * The transform is based on a centered stick-breaking process.
 *
 * @param y Free vector input of dimensionality K - 1
 * @return Simplex of dimensionality K
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline auto simplex_constrain(const T& y) {
  return adj_jac_apply<internal::simplex_constrain_op<T>>(y);
}

}  // namespace math
}  // namespace stan
#endif
