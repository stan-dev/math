#ifndef STAN_MATH_REV_FUN_POSITIVE_ORDERED_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_POSITIVE_ORDERED_CONSTRAIN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/adj_jac_apply.hpp>
#include <cmath>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace internal {
template <typename T>
class positive_ordered_constrain_op {
public:
  adj_op<T> exp_x_;
  explicit positive_ordered_constrain_op(const T& x) : exp_x_(exp(x.val())) {}
  /**
   * Return an increasing positive ordered vector derived from the specified
   * free vector.  The returned constrained vector will have the
   * same dimensionality as the specified free vector.
   *
   * @tparam size Number of arguments
   * @param needs_adj Boolean indicators of if adjoints of arguments will be
   * needed
   * @param x Free vector of scalars
   * @return Positive, increasing ordered vector
   */
  template <std::size_t size>
  Eigen::VectorXd operator()(const std::array<bool, size>& needs_adj,
                             const Eigen::VectorXd& x) {
    using std::exp;
    Eigen::Matrix<double, Eigen::Dynamic, 1> y(x.size());
    if (x.size() == 0) {
      return y;
    }
    y[0] = exp_x_(0);
    for (int n = 1; n < x.size(); ++n) {
      y[n] = y[n - 1] + exp_x_(n);
    }
    return y;
  }

  /**
   * Compute the result of multiply the transpose of the adjoint vector times
   * the Jacobian of the positive_ordered_constrain operator.
   *
   * @tparam size Number of adjoints to return
   * @param needs_adj Boolean indicators of if adjoints of arguments will be
   * needed
   * @param adj Eigen::VectorXd of adjoints at the output of the softmax
   * @return Eigen::VectorXd of adjoints propagated through softmax operation
   */
  template <std::size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::VectorXd& adj) const {
    const auto x_size = exp_x_.size();
    Eigen::VectorXd adj_times_jac(x_size);
    double rolling_adjoint_sum = 0.0;

    for (int n = x_size; --n >= 0;) {
      rolling_adjoint_sum += adj(n);
      adj_times_jac(n) = exp_x_(n) * rolling_adjoint_sum;
    }

    return std::make_tuple(adj_times_jac);
  }
};
}  // namespace internal

/**
 * Return an increasing positive ordered vector derived from the specified
 * free vector.  The returned constrained vector will have the
 * same dimensionality as the specified free vector.
 *
 * @param x Free vector of scalars
 * @return Positive, increasing ordered vector
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline auto positive_ordered_constrain(
    const T& x) {
  return adj_jac_apply<internal::positive_ordered_constrain_op<T>>(x);
}

}  // namespace math
}  // namespace stan
#endif
