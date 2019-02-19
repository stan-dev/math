#ifndef STAN_MATH_REV_SCAL_FUN_MODIFIED_BESSEL_SECOND_KIND_FRAC_HPP
#define STAN_MATH_REV_SCAL_FUN_MODIFIED_BESSEL_SECOND_KIND_FRAC_HPP

#include <stan/math/rev/arr/functor/integrate_1d.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
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


double
compute_inner_integral_with_gradient(const double &v, const double &z) {
  double relative_tolerance = std::sqrt(std::numeric_limits<double>::epsilon());

  double error;
  double L1;
  size_t levels;

  boost::math::quadrature::tanh_sinh<double> integrator;
  
  std::vector<double> theta = {v, z};

  auto f_wrap = [&](double x, double xc) { 
    integrate_func f;  
    return f(x, xc, theta, std::vector<double>(), std::vector<int>(), 0); 
  };
  double integral = integrator.integrate(f_wrap, 0.0, 1.0, relative_tolerance, 
                                &error, &L1, &levels);

  //TODO: check error

  return integral;
}


template <typename T_v, typename T_z>
typename boost::enable_if_c<boost::is_same<T_v, stan::math::var>::value
                                       || boost::is_same<T_z, stan::math::var>::value,
                                   stan::math::var>::type
compute_inner_integral_with_gradient(const T_v &v, const T_z &z) {
  typedef typename boost::math::tools::promote_args<T_v, T_z>::type T_Ret;

  double error;
  double L1;
  size_t levels;
  double condition_number;
  double relative_tolerance = std::sqrt(std::numeric_limits<double>::epsilon());

  std::vector<T_Ret> theta = {v, z};

  boost::math::quadrature::tanh_sinh<double> integrator;
  
  double integral = compute_inner_integral_with_gradient(stan::math::value_of(v), stan::math::value_of(z));

  std::vector<stan::math::var> theta_concat;
  std::vector<double> dintegral_dtheta;

  std::vector<double> theta_vals = stan::math::value_of(theta);

  std::ostringstream msgs;

  if(stan::is_var<T_v>::value) {
    theta_concat.push_back(v);
    auto gradient1_wrap = [&](double x, double xc) { return 
        stan::math::gradient_of_f(integrate_func{}, 
            x, xc, theta_vals, 
                              std::vector<double>(), std::vector<int>(), 0,
                              std::ref(msgs));
    };
    dintegral_dtheta.push_back(
      integrator.integrate(
            gradient1_wrap,
            0.0, 1.0, relative_tolerance,
            &error, &L1, &levels
            )
    );
    //TODO check error
  }  

  if(stan::is_var<T_z>::value) {
    theta_concat.push_back(z);
    auto gradient2_wrap = [&](double x, double xc) { return 
        stan::math::gradient_of_f(integrate_func{}, 
            x, xc, theta_vals, 
                              std::vector<double>(), std::vector<int>(), 1,
                              std::ref(msgs));
    };
    dintegral_dtheta.push_back(
      integrator.integrate(
            gradient2_wrap,
            0.0, 1.0, relative_tolerance,
            &error, &L1, &levels
            )
    );
    //TODO check error
  }    

  return stan::math::precomputed_gradients(integral, theta_concat, dintegral_dtheta);
}

}  // namespace

namespace stan {
namespace math {

template <typename T_v, typename T_z>
typename boost::math::tools::promote_args<T_v, T_z>::type
log_modified_bessel_second_kind_frac(const T_v &v, const T_z &z,
                                     std::ostream &msgs) {
  typedef typename boost::math::tools::promote_args<T_v, T_z>::type T_Ret;

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
  const T_v v_ = fabs(v);
  const T_Ret lead = 0.5 * log(M_PI) - lgamma(v_ + 0.5) - v_ * log(2 * z) - z;
  if (is_inf(lead))
    return -z + 0.5 * log(0.5 * M_PI / z);
  const T_v beta = 16.0 / (2 * v_ + 1);

  T_Ret Q = compute_inner_integral_with_gradient(v, z);
  
  return lead + log(Q);
}

}  // namespace math
}  // namespace stan
#endif
