#ifndef STAN_MATH_PRIM_SCAL_FUNCTOR_COMPLEX_STEP_DERIVATIVE_HPP
#define STAN_MATH_PRIM_SCAL_FUNCTOR_COMPLEX_STEP_DERIVATIVE_HPP

#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>

#include <vector>
#include <iostream>

namespace stan {
namespace math {

/**
 * Return a double that has value of given functor F with signature
 * (complex, std::vector<double>, std::vector<int>, stream*) : complex
 *
 * @tparam F type of functor F
 * @param[in] f functor for the complex number evaluation,
 * must support @c std::complex<double> as arg.
 * @param[in] theta parameter where f and df/d(theta) is requested.
 * @param[in] x_r continuous data vector for the ODE.
 * @param[in] x_i integer data vector for the ODE.
 * @param[out] msgs the print stream for warning messages.
 * @return a var with value f(theta.val()) and derivative at theta.
 */
template <typename F>
double complex_step_derivative(const F& f, const double& theta,
                               const std::vector<double>& x_r,
                               const std::vector<int>& x_i,
                               std::ostream* msgs) {
  return f(theta, x_r, x_i, msgs);
}

}  // namespace math
}  // namespace stan

#endif
