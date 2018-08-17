#ifndef STAN_MATH_PRIM_SCAL_FUN_BINORMAL_INTEGRAL_OWENS_HPP
#define STAN_MATH_PRIM_SCAL_FUN_BINORMAL_INTEGRAL_OWENS_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/fun/Phi.hpp>
#include <stan/math/prim/scal/fun/owens_t.hpp>
#include <stan/math/prim/scal/fun/fmin.hpp>
#include <stan/math/prim/scal/fun/exp.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <iostream>

namespace stan {
namespace math {

template <typename T_z1, typename T_z2, typename T_rho>
typename return_type<T_z1, T_z2, T_rho>::type 
binormal_integral_owens(const T_z1& z1, const T_z2& z2, const T_rho& rho) {
  static const char* function = "binormal_integral_owens";
  typedef typename return_type<T_z1, T_z2, T_rho>::type T_return;
  using std::sqrt;
  using std::fabs;
  using std::min;
  using std::exp;
  using std::asin;

  check_not_nan(function, "Random variable 1", z1);
  check_not_nan(function, "Random variable 2", z2);
  check_bounded(function, "Correlation coefficient", rho, -1,1);
  
  if (z1 > 10 && z2 > 10)
    return 1;
  if (z1 < -37.5 || z2 < -37.5)
    return 0;
  if (rho == 0) 
    return Phi(z1) * Phi(z2);
  if (z1 == rho * z2) 
    return 0.5 / pi() * exp(-0.5 * z2 * z2) * asin(rho) + Phi(z1) * Phi(z2);
  if (z1 * rho == z2) 
    return 0.5 / pi() * exp(-0.5 * z1 * z1) * asin(rho) + Phi(z1) * Phi(z2);
  if (z2 > 40)
    return Phi(z1);
  if (z1 > 40)
    return Phi(z2);
  if (fabs(rho) < 0.7 && fabs(z1) < 1.0 && fabs(z2) < 1.0) {
    T_rho denom = sqrt((1 + rho) * (1 - rho));
    T_return a1 = (z2 / z1 - rho) / denom;
    T_return a2 = (z1 / z2 - rho) / denom;
    T_return product = z1 * z2;
    T_return delta = product < 0 || (product == 0 && (z1 + z2) < 0);
    return 0.5 * (Phi(z1) + Phi(z2) - delta) - owens_t(z1, a1) - owens_t(z2, a2); 
  } 
  if (fabs(rho) < 1) {
    int s = rho > 0 ? 1 : -1;
    const T_z1 h = -z1;
    const T_z2 k = -z2;

    auto f = [&h, &k](const T_return& x) {
      return exp(0.5 * (-h * h - k * k + h * k * 2 * sin(x))/(cos(x) * cos(x)));
    };
    boost::math::quadrature::tanh_sinh<T_return> integrator;
    T_return accum = s == 1 ? integrator.integrate(f, 0, asin(rho))
                            : integrator.integrate(f, static_cast<T_return>(asin(rho)),
                                                   static_cast<T_return>(0));

    accum *= s * 0.5 / pi();
    accum += Phi(z1) * Phi(z2); 
    return accum;
  }
  if (rho == 1) 
    return fmin(Phi(z1),Phi(z2));
  return z2 > -z1 > 0 ? Phi(z1) + Phi(z2) - 1 : 0;
}

}  // namespace math
}  // namespace stan
#endif
