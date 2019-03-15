#ifndef STAN_MATH_REV_SCAL_FUN_LOG_MODIFIED_BESSEL_SECOND_KIND_FRAC_HPP
#define STAN_MATH_REV_SCAL_FUN_LOG_MODIFIED_BESSEL_SECOND_KIND_FRAC_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>

#include <stan/math/rev/scal/fun/pow.hpp>
#include <stan/math/rev/scal/fun/fabs.hpp>

#include <stan/math/rev/scal/fun/cos.hpp>

#include <stan/math/rev/arr/fun/log_sum_exp.hpp>

#include <stan/math/rev/scal/fun/exp.hpp>
#include <stan/math/prim/scal/fun/exp.hpp>

#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>

#include <boost/math/special_functions/bessel.hpp>


#include <stan/math/prim/scal/err/domain_error.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <vector>
#include <limits>

// The formulas and code for small z are based on
// https://github.com/stan-dev/stan/wiki/Stan-Development-Meeting-Agenda/0ca4e1be9f7fc800658bfbd97331e800a4f50011
// Which is in turn based on Equation 26 of Rothwell: Computation of the
// logarithm of Bessel functions of complex argument and fractional order
// https://scholar.google.com/scholar?cluster=2908870453394922596&hl=en&as_sdt=5,33&sciodt=0,33
// 
// http://mathworld.wolfram.com/ModifiedBesselFunctionoftheSecondKind.html

// Assuming all v's positive (flipped at the top-level call)


namespace stan {
namespace math {

namespace besselk_internal {

template <typename T_v, typename T_z, typename T_u>
class inner_integral_rothwell {
 private:
  T_v v;
  T_z z;

 public:
  typedef typename boost::math::tools::promote_args<T_v, T_z, T_u>::type T_Ret;

  inner_integral_rothwell(const T_v &v, const T_z &z) : v(v), z(z) {}

  inline T_Ret operator()(const T_u &u) const {
    using std::pow;
    using std::exp;

    auto v_mhalf = v - 0.5;
    auto neg2v_m1 = -2 * v - 1;
    auto beta = 16.0 / (2 * v + 1);

    T_Ret value;
    T_Ret uB = pow(u, beta);
    T_Ret first
        = beta * exp(-uB) * pow(2 * z + uB, v_mhalf) * boost::math::pow<7>(u);
    T_Ret second = exp(-1.0 / u);
    if (second > 0) {
      second = second * pow(u, neg2v_m1);
      if(is_inf(second)) {
        second = exp(-1.0 / u + neg2v_m1 * log(u));
      }
      second = second * pow(2 * z * u + 1, v_mhalf);
    }
    value = first + second;

    //std::cout << std::setprecision(22) << "U: " << u  << " val: " << value << std::endl;
    
    return value;
  }

  template<class F> 
  static double integrate(const F f, double tolerance, double* error, double* L1, std::size_t* levels) {
    boost::math::quadrature::tanh_sinh<double> integrator;
    return integrator.integrate(f, 0.0, 1.0, tolerance, error, L1, levels);
  };  
};



template <typename T_v, typename T_z, typename T_u>
class inner_integral_mathematica {
 private:
  T_v v;
  T_z z;

 public:
  typedef typename boost::math::tools::promote_args<T_v, T_z, T_u>::type T_Ret;

  inner_integral_mathematica(const T_v &v, const T_z &z) : v(v), z(z) {}

  inline T_Ret operator()(const T_u &u) const {
    using std::pow;
    using std::cos;

    return cos(u) / pow(u * u + z * z, v + 0.5);
  }

  template<class F>
  static double integrate(const F f, double tolerance, double* error, double* L1, std::size_t* levels) {
    boost::math::quadrature::exp_sinh<double> integrator;
    return integrator.integrate(f, tolerance, error, L1, levels);
  };

};

template <typename T_v, typename T_z, typename T_u>
class inner_integral_mathematica_large {
 private:
  T_v v;
  T_z z;

 public:
  typedef typename boost::math::tools::promote_args<T_v, T_z, T_u>::type T_Ret;

  inner_integral_mathematica_large(const T_v &v, const T_z &z) : v(v), z(z) {}

  inline static T_Ret log_factor(const T_u &u, const T_v &v, const T_z &z) {
    //return v * (log(2) + log(z)) - (v + 0.5) * (log(u*u + z *z));
    return v * (log(2) + log(z)) + lgamma(v + 0.5) - (v + 0.5) * (log(u*u + z *z));
  }

  inline T_Ret operator()(const T_u &u) const {
    using std::log;
    using std::exp;
    using std::cos;

    return cos(u) * exp(log_factor(u, v, z));
  }

  template<class F>
  static double integrate(const F f, double tolerance, double* error, double* L1, std::size_t* levels) {
    boost::math::quadrature::exp_sinh<double> integrator;
    return integrator.integrate(f, tolerance, error, L1, levels);
  };

};

// Uses nested autodiff to get gradient with respect to v
// Gradient with respect to z can be computed analytically
// Code modified from gradient_of_f
template <typename F, typename T>
class inner_integral_grad_v {
 private:
  T v;
  T z;

 public:
  inner_integral_grad_v(const T &v, const T &z) : v(v), z(z) {}

  inline T operator()(const T &u) const {
    double gradient = 0.0;
    start_nested();
    var v_var(stan::math::value_of(v));
    try {
      auto f = F(v_var, z);
      var fx = f(u);
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

template <template<typename,typename,typename> class INTEGRAL>
double compute_inner_integral_with_gradient(const double &v, const double &z) {
  double relative_tolerance = std::sqrt(std::numeric_limits<double>::epsilon());

  double error;
  double L1;
  size_t levels;


  auto f = INTEGRAL<double,double,double>(v, z);
  double integral = INTEGRAL<double,double,double>::integrate(f, relative_tolerance,
                                         &error, &L1, &levels);

  if (error > 1e-6 * L1) {
    domain_error("compute_inner_integral_with_gradient(double, double)",
                 "error estimate of integral / L1 ", error / L1, "",
                 "is larger than 1e-6");
  }
  return integral;
}

template <template<typename,typename,typename> class INTEGRAL>
var compute_inner_integral_with_gradient(const var &v, const double &z) {
  double integral = compute_inner_integral_with_gradient<INTEGRAL>(
      stan::math::value_of(v), stan::math::value_of(z));

  std::vector<stan::math::var> theta_concat;
  std::vector<double> dintegral_dtheta;

  theta_concat.push_back(v);
  auto f = inner_integral_grad_v<INTEGRAL<var, double, double>, double>(v.val(), z);

  double error;
  double L1;
  size_t levels;
  double condition_number;
  double relative_tolerance = std::sqrt(std::numeric_limits<double>::epsilon());
  

  dintegral_dtheta.push_back(INTEGRAL<var, double, double>::integrate(
      f, relative_tolerance, &error, &L1, &levels));

  if (error > 1e-6 * L1) {
    domain_error("compute_inner_integral_with_gradient(var, double)",
                 "error estimate of integral / L1 ", error / L1, "",
                 "is larger than 1e-6");
  }

  return stan::math::precomputed_gradients(integral, theta_concat,
                                           dintegral_dtheta);
}

template <typename T_v, typename T_z>
typename boost::math::tools::promote_args<T_v, T_z>::type compute_lead_rothwell(
    const T_v &v, const T_z &z) {
  typedef typename boost::math::tools::promote_args<T_v, T_z>::type T_Ret;

  using std::exp;
  using std::log;
  using std::pow;

  const T_Ret lead = 0.5 * log(pi()) - lgamma(v + 0.5) - v * log(2 * z) - z;
  if (is_inf(lead))
    return -z + 0.5 * log(0.5 * pi() / z);

  return lead;
}

template <typename T_v, typename T_z>
typename boost::math::tools::promote_args<T_v, T_z>::type compute_lead_mathematica(
    const T_v &v, const T_z &z) {
  typedef typename boost::math::tools::promote_args<T_v, T_z>::type T_Ret;

  using std::log;
  
  return lgamma(v + 0.5) + v * (log(2) + log(z)) - 0.5 * log(pi());
}

template <typename T_v, typename T_z>
typename boost::math::tools::promote_args<T_v, T_z>::type compute_lead_mathematica_large(
    const T_v &v, const T_z &z) {
  typedef typename boost::math::tools::promote_args<T_v, T_z>::type T_Ret;

  using std::log;  
  //return lgamma(v + 0.5) - 0.5 * log(pi());
  //return v * (log(2) + log(z)) - 0.5 * log(pi());
  return - 0.5 * log(pi());
}

// Using the first two terms from
// Sidi & Hoggan 2011 "ASYMPTOTICS OF MODIFIED BESSELFUNCTIONS OF HIGH ORDER"
// https://pdfs.semanticscholar.org/d659/ecd3e66fefb2c2fafe16bf3086fcf41c573d.pdf
// Without the second term, the approximation is slightly worse and third term
// did actually do some harm.
// With only the first term the formula reduces to 1.10 of
// Temme, Journal of Computational Physics, vol 19, 324 (1975)
// https://doi.org/10.1016/0021-9991(75)90082-0
template <typename T_v>
T_v
asymptotic_large_v(const T_v &v, const double &z) {
  using std::log;

  //return 0.5 * (log(stan::math::pi()) - log(2) - log(v)) - v * (log(z) - log(2) - log(v));
  return stan::math::LOG_2 - v * (log(z) - stan::math::LOG_2) + lgamma(v) 
      //+ log(1 + (0.25 * boost::math::pow<2>(z) / v) ) 
      //Third term, currently removed
      //+ (- 0.25 * boost::math::pow<2>(z) +  
      //    0.5 * 0.125 * boost::math::pow<4>(z)) / boost::math::pow<2>(v)             
      ;
}

// https://dlmf.nist.gov/10.40
// does not really work
template <typename T_v>
T_v
asymptotic_large_z(const T_v &v, const double &z) {
  using std::pow;
  using std::log;

  const int max_terms = 50;
  int n_terms = std::min(max_terms, static_cast<int>(value_of(v) + 0.5));

  T_v log_series_sum;
  if(n_terms > 1) {
    std::vector<T_v> log_terms;
    log_terms.reserve(max_terms - 1);

    T_v log_a_k(0);
    double log_z = log(z);
    T_v v_squared_4 = v * v * 4;
    double log_8 = log(8);
    
    for(int k = 1; k < n_terms; k++) {
      log_a_k = log_a_k + log(v_squared_4 - (2 * k - 1)) - log(k) - k * log_8;
      log_terms.push_back(log_a_k - k * log_z);
      if(log_terms.back() < -20) {
        break;
      }
    }

    log_series_sum = log_sum_exp(log_terms);
  } else {
    log_series_sum = 0;
  }
  return 0.5 * (log(pi()) - log(2) - log(z)) -z + log_series_sum;
}

//The code to choose computation method is separate, because it is 
//referenced from the test code.
enum class ComputationType { 
  Rothwell, Mathematica, Mathematica_Large, Asymp_v, Asymp_z
  };

const double rothwell_max_v = 75;
const double rothwell_max_log_z_over_v = 600;
const double mathematica_min_log_factor = -200;
const double mathematica_max_log_factor = 100;

inline ComputationType 
choose_computation_type(const double& v, const double& z) {
  using std::fabs;
  using std::pow;
  const double v_ = fabs(v);
  const double rothwell_log_z_boundary = 
    rothwell_max_log_z_over_v * (v_ - 0.5) - log(2);
  const double mathematica_log_factor = -(v_ + 0.5) * 2 * log(z);
  const double mathematica_large_log_factor = 
    inner_integral_mathematica_large<double,double,double>::log_factor(
      0, v_, z);

  if(v_ < rothwell_max_v && 
    (v_ <= 0.5 || log(z) < rothwell_log_z_boundary))
  {
    return ComputationType::Rothwell;
  } else if(mathematica_log_factor < mathematica_max_log_factor && 
    mathematica_log_factor > mathematica_min_log_factor) {
      return ComputationType::Mathematica;
  } else if(mathematica_large_log_factor < mathematica_max_log_factor && 
    mathematica_large_log_factor > mathematica_min_log_factor) {
      return ComputationType::Mathematica_Large;
  } else if(v_ > z){
    return ComputationType::Asymp_v;
  }
  else {
    return ComputationType::Asymp_z;
  }
}

void check_params(const double& v, const double& z) {
  const char *function = "log_modified_bessel_second_kind_frac";
  if(!std::isfinite(v)) {
    stan::math::domain_error(function, 
      "v must be finite", v, "");
  }
  if(!std::isfinite(z)) {
    stan::math::domain_error(function, 
      "z must be finite", z, "");
  }
  if(z < 0) {
    stan::math::domain_error(function, 
      "z is negative", z, "");
  }
}

}  // namespace besselk_internal


template <typename T_v>
T_v log_modified_bessel_second_kind_frac(const T_v &v, const double &z) {
  using std::pow;
  using std::fabs;
  using namespace besselk_internal;
  check_params(value_of(v), value_of(z));

  if(z == 0) {
    return std::numeric_limits<double>::infinity();
  }

  T_v v_ = fabs(v);
  switch(choose_computation_type(value_of(v_), value_of(z))) {
    case ComputationType::Rothwell : {
      T_v lead = compute_lead_rothwell(v_, z);
      T_v Q = compute_inner_integral_with_gradient<inner_integral_rothwell>(
        v_, z);
      return lead + log(Q);
    }
    case ComputationType::Mathematica_Large : {
      T_v lead = compute_lead_mathematica_large(v_, z);
      T_v Q = compute_inner_integral_with_gradient<inner_integral_mathematica_large>(v_, z);
      return lead + log(Q);    
    }
    case ComputationType::Mathematica : {
      T_v lead = compute_lead_mathematica(v_, z);
      T_v Q = compute_inner_integral_with_gradient<inner_integral_mathematica>(
        v_, z);
      return lead + log(Q);    
    }
    case ComputationType::Asymp_v : {
      return asymptotic_large_v(v_, z);
    }
    case ComputationType::Asymp_z : {
      return asymptotic_large_z(v_, z);
    } 
    default : {
          stan::math::domain_error("log_modified_bessel_second_kind_frac", 
      "Invalid computation type ", 0, "");
      return asymptotic_large_v(v_, z);
    }
  }
}

template <typename T_v>
var log_modified_bessel_second_kind_frac(const T_v &v, const var &z) {  
  T_v value = log_modified_bessel_second_kind_frac(v, z.val());

  double value_vm1
      = log_modified_bessel_second_kind_frac(value_of(v) - 1, z.val());
  double gradient_dz
      = -std::exp(value_vm1 - value_of(value)) - value_of(v) / z.val();
  // Compute using boost, seems to be less stable
  // double gradient_dz = 
  //   -boost::math::cyl_bessel_k(value_of(v) - 1, value_of(z)) / 
  //    boost::math::cyl_bessel_k(value_of(v), value_of(z))
  //    - value_of(v) / value_of(z);


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
