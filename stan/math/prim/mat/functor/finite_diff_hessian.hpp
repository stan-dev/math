#ifndef STAN_MATH_PRIM_MAT_FUNCTOR_FINITE_DIFF_HESSIAN_HPP
#define STAN_MATH_PRIM_MAT_FUNCTOR_FINITE_DIFF_HESSIAN_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/functor/finite_diff_gradient.hpp>

namespace stan {
namespace math {
namespace internal {

/**
 * Return the finite difference approximation to the gradient of the
 * specified function at the specified index of the specified argument
 * using the specified step size.
 *
 * <p>The functor must implement
 * <code>
 * double operator()(const Eigen::VectorXd&) const;
 * </code>
 *
 * @tparam F type of function
 * @param f function to differentiate
 * @param x argument to function
 * @param i dimensio of argument for derivative
 * @param epsilon step size for finite differences
 * @return derivative of f(x) with respect to x(i)
 */
template <typename F>
double finite_diff_grad4(const F& f, const Eigen::VectorXd& x, int i,
                         double epsilon = 1e-03) {
  Eigen::VectorXd x_temp(x);

  x_temp(i) = x(i) + 2 * epsilon;
  double grad = -f(x_temp);

  x_temp(i) = x(i) + -2 * epsilon;
  grad += f(x_temp);

  x_temp(i) = x(i) + epsilon;
  grad += 8 * f(x_temp);

  x_temp(i) = x(i) + -epsilon;
  grad -= 8 * f(x_temp);

  return grad;
}
}  // namespace internal

/**
 * Calculate the value and the Hessian of the specified function at
 * the specified argument using second-order finite difference with
 * the specified perturbation step size.
 *
 * <p>The functor must implement
 * <code>
 * double
 * operator()(const Eigen::VectorXd&) const;
 * </code>
 *
 * <p>For details of the algorithm, see
 * <br />Eberly, D., 2008. Derivative approximation by finite
 * differences. Magic Software, Inc., Page 6.
 *
 * @tparam F Type of function
 * @param[in] f Function
 * @param[in] x Argument to function
 * @param[out] fx Function applied to argument
 * @param[out] grad_fx Gradient of function at argument
 * @param[out] hess_fx Hessian of function at argument
 * @param[in] epsilon perturbation step size
 */
template <typename F>
void finite_diff_hessian(const F& f, const Eigen::VectorXd& x, double& fx,
                         Eigen::VectorXd& grad_fx, Eigen::MatrixXd& hess_fx,
                         double epsilon = 1e-03) {
  int d = x.size();
  Eigen::VectorXd x_temp(x);
  hess_fx.resize(d, d);
  finite_diff_gradient(f, x, fx, grad_fx);
  for (int i = 0; i < d; ++i) {
    for (int j = i; j < d; ++j) {
      double f_diff = 0;
      x_temp(i) += 2 * epsilon;
      if (i != j) {
        f_diff = -internal::finite_diff_grad4(f, x_temp, j);
        x_temp(i) = x(i) + -2 * epsilon;
        f_diff += internal::finite_diff_grad4(f, x_temp, j);
        x_temp(i) = x(i) + epsilon;
        f_diff += 8 * internal::finite_diff_grad4(f, x_temp, j);
        x_temp(i) = x(i) + -epsilon;
        f_diff -= 8 * internal::finite_diff_grad4(f, x_temp, j);
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
