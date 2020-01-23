#ifndef STAN_MATH_PRIM_FUNCTOR_FINITE_DIFF_GRADIENT_AUTO_HPP
#define STAN_MATH_PRIM_FUNCTOR_FINITE_DIFF_GRADIENT_AUTO_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/finite_diff_stepsize.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Calculate the value and the gradient of the specified function
 * at the specified argument using finite difference.
 *
 * <p>The functor must implement
 *
 * <code>
 * double operator()(const Eigen::Matrix<double, -, 1>&) const;
 * </code>
 *
 * <p>Error of derivative in dimension `i` should be on the should be on
 * order of `epsilon(i)^6`, where `epsilon(i) = sqrt(delta) * abs(x(i))`
 * for input `x` at dimension `i`.
 *
 * The reference for this algorithm is:
 *
 * <br />Robert de Levie. 2009. An improved numerical approximation
 * for the first derivative. Journal of Chemical Sciences 121(5), page
 * 3.
 *
 * <p>The reference for automatically setting the difference is this
 * section of the Wikipedia,
 *
 * <br /><a
 * href="https://en.wikipedia.org/wiki/Numerical_differentiation#Practical_considerations_using_floating-point_arithmetic">Numerical
 * differentiation: practical considerations using floating point
 * arithmetic</a>.
 *
 * <p>Evaluating this function involves 6 calls to the function being
 * differentiated for each dimension in the input, plus one global
 * evaluation.  All evaluations will be for double-precision inputs.
 *
 * @tparam F Type of function
 * @param[in] f function
 * @param[in] x argument to function
 * @param[out] fx function applied to argument
 * @param[out] grad_fx gradient of function at argument
 */
template <typename F>
void finite_diff_gradient_auto(const F& f, const Eigen::VectorXd& x, double& fx,
                               Eigen::VectorXd& grad_fx) {
  Eigen::VectorXd x_temp(x);
  fx = f(x);
  grad_fx.resize(x.size());
  for (int i = 0; i < x.size(); ++i) {
    double h = finite_diff_stepsize(x(i));

    double delta_f = 0;

    x_temp(i) = x(i) + 3 * h;
    delta_f += f(x_temp);

    x_temp(i) = x(i) + 2 * h;
    delta_f -= 9 * f(x_temp);

    x_temp(i) = x(i) + h;
    delta_f += 45 * f(x_temp);

    x_temp(i) = x(i) + -3 * h;
    delta_f -= f(x_temp);

    x_temp(i) = x(i) + -2 * h;
    delta_f += 9 * f(x_temp);

    x_temp(i) = x(i) - h;
    delta_f -= 45 * f(x_temp);

    delta_f /= 60 * h;

    x_temp(i) = x(i);
    grad_fx(i) = delta_f;
  }
}

}  // namespace math
}  // namespace stan
#endif
