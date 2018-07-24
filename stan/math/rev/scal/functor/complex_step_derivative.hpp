#ifndef STAN_MATH_REV_SCAL_FUNCTOR_COMPLEX_STEP_DERIVATIVE_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_COMPLEX_STEP_DERIVATIVE_HPP

#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/rev/core/precomp_v_vari.hpp>

#include <complex>
#include <iostream>

namespace stan {
namespace math {

template <typename F>
stan::math::var complex_step_derivative(const F& f,
                                        const stan::math::var& theta,
                                        const std::vector<double> &x_r,
                                        const std::vector<int> &x_i,
                                        std::ostream* msgs) {
  using namespace std::complex_literals;
  using std::complex;
  using stan::math::var;
  static double h = 1.e-32;
  const double theta_d = theta.val();
  const double res = f(theta_d, x_r, x_i, msgs);
  const double g = std::imag(f(theta_d + h * 1i, x_r, x_i, msgs)) / h;
  return var(new stan::math::precomp_v_vari(res, theta.vi_, g));
}

}
}

#endif
