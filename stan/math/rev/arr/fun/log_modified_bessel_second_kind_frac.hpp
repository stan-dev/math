#ifndef STAN_MATH_REV_SCAL_FUN_MODIFIED_BESSEL_SECOND_KIND_FRAC_HPP
#define STAN_MATH_REV_SCAL_FUN_MODIFIED_BESSEL_SECOND_KIND_FRAC_HPP

#include <stan/math/rev/arr/functor/integrate_1d.hpp>
#include <stan/math/rev/core.hpp>
#include <vector>

namespace {
// Note: This HAS to be outside the stan::math namespace, otherwise resolution
// between std::pow and stan::math::pow fails
struct integrate_func {
  template <typename T1, typename T2, typename T3>
  inline typename stan::return_type<T1, T2, T3>::type operator()(
      const T1 &u, const T2 &, const std::vector<T3> &theta,
      const std::vector<double> &, const std::vector<int> &,
      std::ostream *) const {
    typedef typename stan::return_type<T1, T2, T3>::type T_Ret;

    auto v = theta[0];
    auto z = theta[1];

    auto v_ = fabs(v);
    auto v_mhalf = v_ - 0.5;
    auto neg2v_m1 = -2 * v_ - 1;
    auto beta = 16.0 / (2 * v_ + 1);

    T_Ret uB = pow(u, beta);
    T_Ret first
        = beta * exp(-uB) * pow(2 * z + uB, v_mhalf) * boost::math::pow<7>(u);
    T_Ret second = exp(-1.0 / u);
    if (second > 0)
      second = second * pow(u, neg2v_m1) * pow(2 * z * u + 1, v_mhalf);
    return first + second;
  }
};
}  // namespace

namespace stan {
namespace math {

template <typename T0__, typename T1__>
typename boost::math::tools::promote_args<T0__, T1__>::type
log_modified_bessel_second_kind_frac(const T0__ &v, const T1__ &z,
                                     std::ostream &msgs) {
  typedef typename boost::math::tools::promote_args<T0__, T1__>::type T_Ret;

  // Based on
  // https://github.com/stan-dev/stan/wiki/Stan-Development-Meeting-Agenda/0ca4e1be9f7fc800658bfbd97331e800a4f50011
  // Which is in turn based on Equation 26 of Rothwell: Computation of the
  // logarithm of Bessel functions of complex argument and fractional order
  // https://scholar.google.com/scholar?cluster=2908870453394922596&hl=en&as_sdt=5,33&sciodt=0,33

  using std::exp;
  using std::fabs;
  using std::log;
  using std::pow;

  if (v == 0.5)
    return 0.5 * log(M_PI / (2 * z)) - z;
  const T0__ v_ = fabs(v);
  const T_Ret lead = 0.5 * log(M_PI) - lgamma(v_ + 0.5) - v_ * log(2 * z) - z;
  if (is_inf(lead))
    return -z + 0.5 * log(0.5 * M_PI / z);
  const T0__ beta = 16.0 / (2 * v_ + 1);

  long double error;
  long double L1;
  size_t levels;
  long double condition_number;

  std::vector<T_Ret> theta = {v, z};

  // Note: we increase the relative tolerance, as the integrator fails with low
  // tolerance, although the results of the full function still passes the tests
  T_Ret Q = integrate_1d(integrate_func{}, 0.0, 1.0, theta,
                         std::vector<double>(), std::vector<int>(), msgs, 1e-4);

  return lead + log(Q);
}

}  // namespace math
}  // namespace stan
#endif
