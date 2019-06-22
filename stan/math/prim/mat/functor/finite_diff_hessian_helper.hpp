#ifndef STAN_MATH_PRIM_MAT_FUNCTOR_FINITE_DIFF_HESSIAN_HELPER_HPP
#define STAN_MATH_PRIM_MAT_FUNCTOR_FINITE_DIFF_HESSIAN_HELPER_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the subcalculation required by `finite_diff_hessian` and
 * `finite_diff_hessian_auto`.  The calculation is like a partial
 * derivative calculation, but returns different values.  The only
 * utility of this function is as a subroutine for the two Hessians.
 *
 * <p>The functor must implement
 * <code>
 * double operator()(const Eigen::VectorXd&) const;
 * </code>
 *
 * <p>This function evaluations the functor four times.
 *
 * @tparam F type of function
 * @param f function to differentiate
 * @param x argument to function
 * @param i dimension of argument for derivative
 * @param epsilon step size for finite differences
 * @return derivative of f(x) with respect to x(i)
 */
template <typename F>
double finite_diff_hessian_helper(const F& f, const Eigen::VectorXd& x, int i,
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
}  // namespace math
}  // namespace stan
#endif
