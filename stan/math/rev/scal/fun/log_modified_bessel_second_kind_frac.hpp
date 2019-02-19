#ifndef STAN_MATH_REV_SCAL_FUN_LOG_MODIFIED_BESSEL_SECOND_KIND_FRAC_HPP
#define STAN_MATH_REV_SCAL_FUN_LOG_MODIFIED_BESSEL_SECOND_KIND_FRAC_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <vector>
#include <iostream>

// The formulas and code are based on
// https://github.com/stan-dev/stan/wiki/Stan-Development-Meeting-Agenda/0ca4e1be9f7fc800658bfbd97331e800a4f50011
// Which is in turn based on Equation 26 of Rothwell: Computation of the
// logarithm of Bessel functions of complex argument and fractional order
// https://scholar.google.com/scholar?cluster=2908870453394922596&hl=en&as_sdt=5,33&sciodt=0,33

namespace stan {
namespace math {

namespace {

template <typename T_v, typename T_z, typename T_u>
class inner_integral {
 private:
  T_v v;
  T_z z;

 public:
  typedef typename boost::math::tools::promote_args<T_v, T_z, T_u>::type T_Ret;

  inner_integral(const T_v &v, const T_z &z) : v(v), z(z) {}

  inline T_Ret operator()(const T_u &u, const T_u &) const {
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

// Uses nested autodiff to get gradient with respect to v
// Gradient with respect to z can be computed analytically
// Code modified from gradient_of_f
template <typename T>
class inner_integral_grad_v {
 private:
  T v;
  T z;

 public:
  inner_integral_grad_v(const T &v, const T &z) : v(v), z(z) {}

  inline T operator()(const T &u, const T &uc) const {
    double gradient = 0.0;
    start_nested();
    var v_var(v);
    try {
      auto f = inner_integral<var, T, T>(v_var, z);
      var fx = f(u, uc);
      fx.grad();
      gradient = v_var.adj();
      if (is_nan(gradient)) {
        if (fx.val() == 0) {
          gradient = 0;
        } else {
          domain_error("inner_integral_grad_v",
                       "The gradient of inner_integral is nan", 0, "", "");
        }
      }
    } catch (const std::exception &e) {
      recover_memory_nested();
      throw;
    }
    recover_memory_nested();
    return gradient;
  }
};

double compute_inner_integral_with_gradient(const double &v, const double &z) {
  double relative_tolerance = std::sqrt(std::numeric_limits<double>::epsilon());

  double error;
  double L1;
  size_t levels;

  boost::math::quadrature::tanh_sinh<double> integrator;

  auto f = inner_integral<double, double, double>(v, z);
  double integral = integrator.integrate(f, 0.0, 1.0, relative_tolerance,
                                         &error, &L1, &levels);

  if (error > 1e-6 * L1) {
    domain_error("compute_inner_integral_with_gradient",
                 "error estimate of integral / L1 ", error / L1, "",
                 "is larger than 1e-6");
  }
  return integral;
}

var compute_inner_integral_with_gradient(const var &v, const double &z) {
  double integral = compute_inner_integral_with_gradient(
      stan::math::value_of(v), stan::math::value_of(z));

  std::vector<stan::math::var> theta_concat;
  std::vector<double> dintegral_dtheta;

  theta_concat.push_back(v);
  auto f = inner_integral_grad_v<double>(v.val(), z);

  double error;
  double L1;
  size_t levels;
  double condition_number;
  double relative_tolerance = std::sqrt(std::numeric_limits<double>::epsilon());
  boost::math::quadrature::tanh_sinh<double> integrator;

  dintegral_dtheta.push_back(integrator.integrate(
      f, 0.0, 1.0, relative_tolerance, &error, &L1, &levels));

  if (error > 1e-6 * L1) {
    domain_error("compute_inner_integral_with_gradient",
                 "error estimate of integral / L1 ", error / L1, "",
                 "is larger than 1e-6");
  }

  return stan::math::precomputed_gradients(integral, theta_concat,
                                           dintegral_dtheta);
}

template <typename T_v, typename T_z>
typename boost::math::tools::promote_args<T_v, T_z>::type compute_lead(
    const T_v &v, const T_z &z) {
  typedef typename boost::math::tools::promote_args<T_v, T_z>::type T_Ret;

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

  return lead;
}
}  // namespace

template <typename T_v>
T_v log_modified_bessel_second_kind_frac(const T_v &v, const double &z) {
  T_v lead = compute_lead(v, z);
  T_v Q = compute_inner_integral_with_gradient(v, z);
  return lead + log(Q);
}

template <typename T_v>
var log_modified_bessel_second_kind_frac(const T_v &v, const var &z) {
  T_v lead = compute_lead(v, z.val());

  T_v Q = compute_inner_integral_with_gradient(v, z.val());

  T_v value = lead + log(Q);

  double value_vm1
      = log_modified_bessel_second_kind_frac(value_of(v) - 1, z.val());
  double gradient_dz
      = -exp(value_vm1 - value_of(value)) - value_of(v) / z.val();

  std::vector<var> operands;
  std::vector<double> gradients;
  if (is_var<T_v>::value) {
    // A trick to combine the autodiff gradient with precomputed_gradients
    operands.push_back(value);
    gradients.push_back(1);
  }
  operands.push_back(z);
  gradients.push_back(gradient_dz);

  return precomputed_gradients(value_of(value), operands, gradients);
}

}  // namespace math
}  // namespace stan
#endif
