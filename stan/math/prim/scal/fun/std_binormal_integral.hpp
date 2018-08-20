#ifndef STAN_MATH_PRIM_SCAL_FUN_STD_BINORMAL_INTEGRAL_HPP
#define STAN_MATH_PRIM_SCAL_FUN_STD_BINORMAL_INTEGRAL_HPP

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
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/prob/std_normal_lpdf.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <cmath>
#include <iostream>

namespace stan {
namespace math {

double std_binormal_integral(const double z1, 
                             const double z2, 
                             const double rho) {
  static const char* function = "std_binormal_integral";
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
  if (fabs(rho) < 0.7 && fabs(z1) <= 1 && fabs(z2) <= 1) {
    double denom = sqrt((1 + rho) * (1 - rho));
    double a1 = (z2 / z1 - rho) / denom;
    double a2 = (z1 / z2 - rho) / denom;
    double product = z1 * z2;
    double delta = product < 0 || (product == 0 && (z1 + z2) < 0);
    return 0.5 * (Phi(z1) + Phi(z2) - delta) - owens_t(z1, a1) - owens_t(z2, a2); 
  } 
  if (fabs(rho) < 1) {
    int s = rho > 0 ? 1 : -1;
    double h = -z1;
    double k = -z2;

    auto f = [&h, &k](const double x) {
      return exp(0.5 * (-h * h - k * k + h * k * 2 * sin(x))/(cos(x) * cos(x)));
    };
    double L1, error;
    boost::math::quadrature::tanh_sinh<double> integrator;
    double accum = s == 1 ? integrator.integrate(f, 0, asin(rho), 
                                                   sqrt(stan::math::EPSILON),
                                                   &error, &L1)
                            : integrator.integrate(f, static_cast<double>(asin(rho)),
                                                   static_cast<double>(0), 
                                                   sqrt(stan::math::EPSILON),
                                                   &error, &L1);
    if (exp(log(L1) - log(accum)) > 1.5) {
      std::cout << "L1: " << L1 << " integral: " << accum << std::endl;
      domain_error(function, "the numeric integration of bivariate normal density condition number", exp(log(L1) - log(accum)),
                   "",
                   " indicates the integral in poorly conditioned");
    }


    accum *= s * 0.5 / pi();
    accum += Phi(z1) * Phi(z2); 
    return accum;
  }
  if (rho == 1) 
    return fmin(Phi(z1),Phi(z2));
  return z2 > -z1 > 0 ? Phi(z1) + Phi(z2) - 1 : 0;
}

template <typename T_z1, typename T_z2, typename T_rho>
typename return_type<T_z1, T_z2, T_rho>::type
                                     std_binormal_integral(
    const T_z1& z1, const T_z2& z2, const T_rho& rho) {
  typedef typename partials_return_type<T_z1, T_z2, T_rho>::type partials_type;
  using std::exp;
  using std::sqrt;
  using std::log;
  using std::asin;
  using stan::math::std_binormal_integral;
  using stan::math::std_normal_lpdf;

  const partials_type z1_dbl = value_of(z1);
  const partials_type z2_dbl = value_of(z2);
  const partials_type rho_dbl = value_of(rho);

  partials_type cdf = std_binormal_integral(z1_dbl, z2_dbl, rho_dbl);
  operands_and_partials<T_z1, T_z2, T_rho> ops_partials(z1, z2, rho);

  if (contains_nonconstant_struct<T_z1, T_z2, T_rho>::value) {
    const partials_type one_minus_rho_sq = (1 + rho_dbl) * (1 - rho_dbl);
    const partials_type sqrt_one_minus_rho_sq = sqrt(one_minus_rho_sq);
    const partials_type rho_times_z2 = rho_dbl * z2_dbl;
    const partials_type z1_minus_rho_times_z2 = z1_dbl - rho_times_z2;
    if (!is_constant_struct<T_z1>::value)
      ops_partials.edge1_.partials_[0] 
        += cdf > 0 && cdf < 1 ? exp(std_normal_lpdf(z1_dbl))
            * Phi((z2_dbl - rho_dbl * z1_dbl) / sqrt_one_minus_rho_sq)
            : 0;
    if (!is_constant_struct<T_z2>::value)
      ops_partials.edge2_.partials_[0] 
        += cdf > 0 && cdf < 1 ? exp(std_normal_lpdf(z2_dbl))
            * Phi(z1_minus_rho_times_z2 / sqrt_one_minus_rho_sq)
            : 0;
    if (!is_constant_struct<T_rho>::value)
      ops_partials.edge3_.partials_[0]
          += cdf > 0 && cdf < 1 ? 0.5 / (stan::math::pi() * 
                                         sqrt_one_minus_rho_sq)
          * exp(-0.5 / one_minus_rho_sq 
                * z1_minus_rho_times_z2 *  z1_minus_rho_times_z2
                -0.5 * z2_dbl * z2_dbl)
            : 0;
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
