#ifndef STAN_MATH_LAPLACE_THIRD_DIFF_DIRECTIONAL_HPP
#define STAN_MATH_LAPLACE_THIRD_DIFF_DIRECTIONAL_HPP

// TODO: refine include.
#include <stan/math/mix.hpp>

namespace stan {
namespace math {

  /**
   * Return the third-order directional derivative of a function
   * which maps to a scalar. The derivative is taken with respect
   * to do two directions: v and w.
   */
  template <typename F>
  void third_diff_directional(
    const F& f, const Eigen::VectorXd& x,
    const Eigen::VectorXd& eta,
    const Eigen::VectorXd& delta,
    const std::vector<int>& delta_int,
    double& fx,
    Eigen::VectorXd& third_diff,
    Eigen::VectorXd& v,
    Eigen::VectorXd& w,
    std::ostream* pstream = 0) {
    using Eigen::Matrix;
    using Eigen::Dynamic;
    nested_rev_autodiff nested;

    int x_size = x.size();
    Matrix<var, Dynamic, 1> x_var = x;
    Matrix<fvar<var>, Dynamic, 1> x_fvar(x_size);
    for (int i = 0; i < x_size; ++i) {
      x_fvar(i) = fvar<var>(x_var(i), v(i));
    }
    fvar<var> fx_fvar = f(x_fvar, eta, delta, delta_int, pstream);

    Matrix<fvar<fvar<var>>, Dynamic, 1> x_ffvar(x_size);
    for (int i = 0; i < x_size; ++i) {
      x_ffvar(i) = fvar<fvar<var>>(x_fvar(i), w(i));
    }
    fvar<fvar<var>> fx_ffvar = f(x_ffvar, eta, delta, delta_int, pstream);

    grad(fx_ffvar.d_.d_.vi_);

    third_diff.resize(x_size);
    for (int i = 0; i < x_size; ++i) {
      third_diff(i) = x_var(i).adj();
    }
  }

}  // namespace math
}  // namespace stan

#endif
