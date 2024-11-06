#ifndef STAN_MATH_MIX_FUNCTOR_FINITE_DIFF_GRAD_HESSIAN_AUTO_HPP
#define STAN_MATH_MIX_FUNCTOR_FINITE_DIFF_GRAD_HESSIAN_AUTO_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/mix/functor/hessian.hpp>
#include <stan/math/prim/functor/finite_diff_gradient_auto.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Calculate the value, Hessian, and the gradient of the Hessian of
 * the specified function at the specified argument using second-order
 * autodiff and first-order finite difference.
 *
 * <p>The functor must implement
 *
 * <code>
 * double operator()(const Eigen::VectorXd&) const;
 * </code>
 *
 * <p>Reference for finite difference to compute gradient:
 *
 * <br />De Levie: An improved numerical approximation
 * for the first derivative, page 3
 *
 * <p>Step size for dimension `i` is set automatically using
 * `stan::math::finite_diff_stepsize(H(i, j))`; the nested
 * finite differences are over entries in the Hessian.
 *
 * <p>Evaluating this function involves 6 calls to the Hessian
 * autodiff function for each entry in the Hessian.
 *
 * @tparam F Type of function
 * @param[in] f Function
 * @param[in] x Argument to function
 * @param[out] fx Function applied to argument
 * @param[out] hess Hessian matrix
 * @param[out] grad_hess_fx gradient of Hessian of function at argument
 */
template <typename F>
void finite_diff_grad_hessian_auto(const F& f, const Eigen::VectorXd& x,
                                   double& fx, Eigen::MatrixXd& hess,
                                   std::vector<Eigen::MatrixXd>& grad_hess_fx) {
  int d = x.size();

  grad_hess_fx.clear();
  grad_hess_fx.reserve(d);

  Eigen::VectorXd x_temp(x);
  Eigen::VectorXd grad_auto(d);
  Eigen::MatrixXd hess_auto(d, d);
  Eigen::MatrixXd hess_diff(d, d);

  hessian(f, x, fx, grad_auto, hess);

  for (int i = 0; i < d; ++i) {
    double dummy_fx_eval;
    double epsilon = finite_diff_stepsize(x(i));
    hess_diff.setZero();

    x_temp(i) = x(i) + 2 * epsilon;
    hessian(f, x_temp, dummy_fx_eval, grad_auto, hess_auto);
    hess_diff = -hess_auto;

    x_temp(i) = x(i) + -2 * epsilon;
    hessian(f, x_temp, dummy_fx_eval, grad_auto, hess_auto);
    hess_diff += hess_auto;

    x_temp(i) = x(i) + epsilon;
    hessian(f, x_temp, dummy_fx_eval, grad_auto, hess_auto);
    hess_diff += 8 * hess_auto;

    x_temp(i) = x(i) + -epsilon;
    hessian(f, x_temp, dummy_fx_eval, grad_auto, hess_auto);
    hess_diff -= 8 * hess_auto;

    x_temp(i) = x(i);
    hess_diff /= 12 * epsilon;

    grad_hess_fx.push_back(hess_diff);
  }
  fx = f(x);
}

}  // namespace math
}  // namespace stan
#endif
