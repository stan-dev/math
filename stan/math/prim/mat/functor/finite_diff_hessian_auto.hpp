#ifndef STAN_MATH_PRIM_MAT_FUNCTOR_FINITE_DIFF_HESSIAN_AUTO_HPP
#define STAN_MATH_PRIM_MAT_FUNCTOR_FINITE_DIFF_HESSIAN_AUTO_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/functor/finite_diff_gradient_auto.hpp>
#include <stan/math/prim/mat/functor/finite_diff_hessian_helper.hpp>
#include <stan/math/prim/scal/fun/finite_diff_stepsize.hpp>

namespace stan {
namespace math {

/**
 * Calculate the value and the Hessian of the specified function at
 * the specified argument using second-order finite difference,
 * automatically setting the stepsize between the function evaluations
 * along a dimension.
 *
 * <p>The functor must implement
 *
 * <code>
 * double operator()(const Eigen::VectorXd&)
 * </code>
 *
 * <p>For details of the algorithm, see
 * <br />Eberly, D., 2008. Derivative approximation by finite
 * differences. Magic Software, Inc., Page 6.
 *
 * <p>Step size for dimension `i` is set automatically using
 * `stan::math::finite_diff_stepsize(x(i))`.
 *
 * <p>For each non-diagonal entry in the Hessian, the function is
 * evaluated 16 times; the diagonal entries require 4 function evaluations.
 *
 * @tparam F Type of function
 * @param[in] f Function
 * @param[in] x Argument to function
 * @param[out] fx Function applied to argument
 * @param[out] grad_fx Gradient of function at argument
 * @param[out] hess_fx Hessian of function at argument
 */
template <typename F>
void finite_diff_hessian_auto(const F& f, const Eigen::VectorXd& x, double& fx,
                              Eigen::VectorXd& grad_fx,
                              Eigen::MatrixXd& hess_fx) {
  int d = x.size();

  Eigen::VectorXd x_temp(x);
  hess_fx.resize(d, d);

  finite_diff_gradient_auto(f, x, fx, grad_fx);
  double f_diff = 0;
  for (int i = 0; i < d; ++i) {
    for (int j = i; j < d; ++j) {
      double epsilon = finite_diff_stepsize(x(i));
      x_temp(i) += 2 * epsilon;
      if (i != j) {
        f_diff = -finite_diff_hessian_helper(f, x_temp, j, epsilon);
        x_temp(i) = x(i) + -2 * epsilon;
        f_diff += finite_diff_hessian_helper(f, x_temp, j, epsilon);
        x_temp(i) = x(i) + epsilon;
        f_diff += 8 * finite_diff_hessian_helper(f, x_temp, j, epsilon);
        x_temp(i) = x(i) + -epsilon;
        f_diff -= 8 * finite_diff_hessian_helper(f, x_temp, j, epsilon);
        f_diff /= 12 * epsilon * 12 * epsilon;
      } else {
        f_diff = -f(x_temp);
        f_diff -= 30 * fx;
        x_temp(i) = x(i) + -2 * epsilon;
        f_diff -= f(x_temp);
        x_temp(i) = x(i) + epsilon;
        f_diff += 16 * f(x_temp);
        x_temp(i) = x(i) - epsilon;
        f_diff += 16 * f(x_temp);
        f_diff /= 12 * epsilon * epsilon;
      }
      x_temp(i) = x(i);
      hess_fx(j, i) = f_diff;
      hess_fx(i, j) = hess_fx(j, i);
    }
  }
}
}  // namespace math
}  // namespace stan
#endif
