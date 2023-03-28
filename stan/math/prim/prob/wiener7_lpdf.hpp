#ifndef STAN_MATH_PRIM_FUN_WIENER7_LPDF_HPP
#define STAN_MATH_PRIM_FUN_WIENER7_LPDF_HPP

#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/functor/hcubature.hpp>
#include <stan/math/prim/prob/wiener5_lpdf.hpp>

namespace stan {
namespace math {
namespace internal {

// tools
/*struct my_params {
  double y, a, v, w, t0, sv, sw, st0, lerr;
};

/*struct my_params2 {
  double y, a, v, w, w_lower, w_upper, t0, sv, sw, sw_mean, st0, lerr;
};

struct my_params3 {
  double y, a, v, w, t0, t0_mean, sv, sw, st0, st0_mean, lerr;
};*/

// calculate derivative of density of wiener5 with respect to y (version for
// wiener7)
double dtdwiener5_for_7(const double& y, const double& a, const double& v,
                        const double& w, const double& sv, const double& err) {
  double kll, kss, ans;

  // prepare some variables
  double y_asq = y / square(a);
  double la = 2.0 * log(a);
  double ans0, lg1;
  if (sv != 0) {
    double sv_sqr = square(sv);
    double one_plus_svsqr_y = (1 + sv_sqr * y);
    ans0 = -0.5
           * (square(sv_sqr) * (y + square(a * w))
              + sv_sqr * (1 - 2 * a * v * w) + square(v))
           / square(one_plus_svsqr_y);
    lg1 = (sv_sqr * square(a * w) - 2 * a * v * w - square(v) * y) / 2.0
              / one_plus_svsqr_y
          - la - 0.5 * log(one_plus_svsqr_y);
  } else {
    ans0 = -0.5 * square(v);
    lg1 = (-2 * a * v * w - square(v) * y) / 2.0 - la;
  }
  double factor = lg1 - la;
  double ld = dwiener5(y, a, -v, 1 - w, sv,
                       err - log(max(fabs(ans0), fabs(ans0 - 1.5 / y))));

  // calculate the number of terms kss needed for small y
  double es = err - lg1;
  es = es + la;
  double K1 = (sqrt(3.0 * y_asq) + w) / 2.0;
  double u_eps = fmin(
      -1.0, (log(8.0 / 27.0) + LOG_PI + 4.0 * log(y_asq) + 2.0 * es) / 3.0);
  double arg = -3.0 * y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  double K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate number of terms kll needed for large y
  double el = err - lg1;
  el = el + la;
  K1 = sqrt(3.0 / y_asq) / pi();
  u_eps = fmin(-1.0, el + log(0.6) + LOG_PI + 2.0 * log(y_asq));
  static const double PISQ = square(pi());  // pi*pi
  arg = -2.0 / PISQ / y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  double kl = (arg > 0) ? sqrt(arg) : K1;
  kll = ceil(fmax(kl, K1));

  double erg;
  double newsign = 1;
  double fplus = NEGATIVE_INFTY;
  double fminus = NEGATIVE_INFTY;

  // if small y is better
  if (2 * kss < kll) {
    // calculate terms of the sum for small y
    double twot = 2.0 * y_asq;
    if (static_cast<size_t>(kss) > 0) {
      for (size_t k = static_cast<size_t>(kss); k >= 1; k--) {
        double w_plus_2k = w + 2.0 * k;
        double w_minus_2k = w - 2.0 * k;
        fplus = log_sum_exp(3.0 * log(w_plus_2k) - w_plus_2k * w_plus_2k / twot,
                            fplus);
        fminus = log_sum_exp(
            3.0 * log(-w_minus_2k) - w_minus_2k * w_minus_2k / twot, fminus);
      }
    }
    fplus = log_sum_exp(3.0 * log(w) - w * w / twot, fplus);
    if (fplus < fminus) {
      newsign = -1;
      erg = log_diff_exp(fminus, fplus);
    } else {
      erg = log_diff_exp(fplus, fminus);
    }
    ans = ans0 - 1.5 / y
          + newsign
                * exp(factor - 1.5 * LOG_TWO - LOG_SQRT_PI - 3.5 * log(y_asq)
                      + erg - ld);
    // if large y is better
  } else {
    // calculate terms of the sum for large y
    double halfy = y_asq / 2.0;
    for (size_t k = static_cast<size_t>(kll); k >= 1; k--) {
      double pi_k = pi() * k;
      double zwi = sin(pi_k * w);
      if (zwi > 0) {
        fplus
            = log_sum_exp(3.0 * log(k) - pi_k * pi_k * halfy + log(zwi), fplus);
      }
      if (zwi < 0) {
        fminus = log_sum_exp(3.0 * log(k) - pi_k * pi_k * halfy + log(-zwi),
                             fminus);
      }
    }
    if (fplus < fminus) {
      erg = log_diff_exp(fminus, fplus);
      newsign = -1;
    } else {
      erg = log_diff_exp(fplus, fminus);
    }
    ans = ans0 - newsign * exp(factor + 3.0 * LOG_PI - LOG_TWO + erg - ld);
  }
  return ans * exp(ld);
}
//-----------------------------------------------

// integrand density (on normal scale)
template <typename T_x, typename T_p>
return_type_t<T_x, T_p> int_ddiff(const T_x& x, const T_p& p) {
	struct my_params {
  double y, a, v, w, t0, sv, sw, st0, lerr;
};
  my_params* params = static_cast<my_params*>(p);
  double y = (params->y);
  double a = (params->a);
  double v = (params->v);
  double t0 = (params->t0);
  double w = (params->w);
  double sw = (params->sw);
  double sv = (params->sv);
  double st0 = (params->st0);
  double lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);
  double retval;
  double omega = sw ? w + sw * (x_vec[0] - 0.5) : w;
  double t0_ = sw ? (st0 ? t0 + st0 * x_vec[1] : t0)
                  : (st0 ? t0 + st0 * x_vec[0] : t0);
  if (y - t0_ <= 0) {
    retval = 0.0;
  } else {
    double ldW = dwiener5(y - t0_, a, v, omega, sv, lerr);
    retval = exp(ldW);
  }
  return retval;
}
//-----------------------------------------------

// integrand d/dt (on normal scale)
template <typename T_x, typename T_p>
return_type_t<T_x, T_p> int_dtddiff(const T_x& x, const T_p& p) {
	struct my_params {
  double y, a, v, w, t0, sv, sw, st0, lerr;
};
  my_params* params = static_cast<my_params*>(p);
  double y = (params->y);
  double a = (params->a);
  double v = (params->v);
  double t0 = (params->t0);
  double w = (params->w);
  double sw = (params->sw);
  double sv = (params->sv);
  double st0 = (params->st0);
  double lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);

  double retval;
  double omega = sw ? w + sw * (x_vec[0] - 0.5) : w;
  double t0_ = sw ? (st0 ? t0 + st0 * x_vec[1] : t0)
                  : (st0 ? t0 + st0 * x_vec[0] : t0);
  if (y - t0_ <= 0) {
    retval = 0.0;
  } else {
    retval = dtdwiener5_for_7(y - t0_, a, -v, 1 - omega, sv, lerr);
  }
  return retval;
}
//-----------------------------------------------

// integrand d/da (on normal scale)
template <typename T_x, typename T_p>
return_type_t<T_x, T_p> int_daddiff(const T_x& x, const T_p& p) {
	struct my_params {
  double y, a, v, w, t0, sv, sw, st0, lerr;
};
  my_params* params = static_cast<my_params*>(p);
  double y = (params->y);
  double a = (params->a);
  double v = (params->v);
  double t0 = (params->t0);
  double w = (params->w);
  double sw = (params->sw);
  double sv = (params->sv);
  double st0 = (params->st0);
  double lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);

  double retval;
  double omega = sw ? w + sw * (x_vec[0] - 0.5) : w;
  double t0_ = sw ? (st0 ? t0 + st0 * x_vec[1] : t0)
                  : (st0 ? t0 + st0 * x_vec[0] : t0);
  if (y - t0_ <= 0) {
    retval = 0.0;
  } else {
    // prepare some variables for error estimation in density
    double sv_sqr = square(sv);
    double ans0 = (v * (1 - omega) + sv_sqr * square(1 - omega) * a)
                  / (1 + sv_sqr * (y - t0_));
    double factor
        = max(ans0 + 1.0 / a, ans0 - 2.0 / a);  // factor from small and large
                                                // representation, for same
                                                // computation as in dadwiener5
    // compute partial derivative
    retval = dadwiener5(y - t0_, a, v, omega, sv, lerr, 1);
  }
  return retval;
}
//-----------------------------------------------

// integrand d/dv (on normal scale)
template <typename T_x, typename T_p>
return_type_t<T_x, T_p> int_dvddiff(const T_x& x, const T_p& p) {
	struct my_params {
  double y, a, v, w, t0, sv, sw, st0, lerr;
};
  my_params* params = static_cast<my_params*>(p);
  double y = (params->y);
  double a = (params->a);
  double v = (params->v);
  double t0 = (params->t0);
  double w = (params->w);
  double sw = (params->sw);
  double sv = (params->sv);
  double st0 = (params->st0);
  double lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);

  double retval;
  double omega = sw ? w + sw * (x_vec[0] - 0.5) : w;
  double t0_ = sw ? (st0 ? t0 + st0 * x_vec[1] : t0)
                  : (st0 ? t0 + st0 * x_vec[0] : t0);
  if (y - t0_ <= 0) {
    retval = 0.0;
  } else {
    retval = dvdwiener5(y - t0_, a, v, omega, sv)
             * exp(dwiener5(y - t0_, a, v, omega, sv, lerr));
  }
  return retval;
}
//-----------------------------------------------

// integrand d/dw (on normal scale)
template <typename T_x, typename T_p>
return_type_t<T_x, T_p> int_dwddiff(const T_x& x, const T_p& p) {
	struct my_params {
  double y, a, v, w, t0, sv, sw, st0, lerr;
};
  my_params* params = static_cast<my_params*>(p);
  double y = (params->y);
  double a = (params->a);
  double v = (params->v);
  double t0 = (params->t0);
  double w = (params->w);
  double sw = (params->sw);
  double sv = (params->sv);
  double st0 = (params->st0);
  double lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);

  double retval;
  double omega = sw ? w + sw * (x_vec[0] - 0.5) : w;
  double t0_ = sw ? (st0 ? t0 + st0 * x_vec[1] : t0)
                  : (st0 ? t0 + st0 * x_vec[0] : t0);
  if (y - t0_ <= 0) {
    retval = 0.0;
  } else {
    // prepare some variables for error estimation in density
    double sv_sqr = square(sv);
    double ans0
        = (v * a + sv_sqr * square(a) * (1 - omega)) / (1 + sv_sqr * (y - t0_));
    // compute partial derivative
    retval = dwdwiener5(y - t0_, a, v, omega, sv, lerr, 1);
  }
  return retval;
}
//-----------------------------------------------

// integrand d/dsv (on normal scale)
template <typename T_x, typename T_p>
return_type_t<T_x, T_p> int_dsvddiff(const T_x& x, const T_p& p) {
	struct my_params {
  double y, a, v, w, t0, sv, sw, st0, lerr;
};
  my_params* params = static_cast<my_params*>(p);
  double y = (params->y);
  double a = (params->a);
  double v = (params->v);
  double t0 = (params->t0);
  double w = (params->w);
  double sw = (params->sw);
  double sv = (params->sv);
  double st0 = (params->st0);
  double lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);

  double retval;
  double omega = sw ? w + sw * (x_vec[0] - 0.5) : w;
  double t0_ = sw ? (st0 ? t0 + st0 * x_vec[1] : t0)
                  : (st0 ? t0 + st0 * x_vec[0] : t0);
  if (y - t0_ <= 0) {
    retval = 0.0;
  } else {
    retval = dsvdwiener5(y - t0_, a, v, omega, sv)
             * exp(dwiener5(y - t0_, a, v, omega, sv, lerr));
  }
  return retval;
}
//-----------------------------------------------

// integrand d/dsw (on normal scale)
template <typename T_x, typename T_p>
return_type_t<T_x, T_p> int_dswddiff(const T_x& x, const T_p& p) {
struct my_params2 {
  double y, a, v, w, w_lower, w_upper, t0, sv, sw, sw_mean, st0, lerr;
};
  my_params2* params = static_cast<my_params2*>(p);
  double y = (params->y);
  double a = (params->a);
  double v = (params->v);
  double t0 = (params->t0);
  double w = (params->w);
  double w_lower = (params->w_lower);
  double w_upper = (params->w_upper);
  double sw = (params->sw);
  double sw_mean = (params->sw_mean);
  double sv = (params->sv);
  double st0 = (params->st0);
  double lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);
  double t0_ = sw ? (st0 ? t0 + st0 * x_vec[1] : t0)
                  : (st0 ? t0 + st0 * x_vec[0] : t0);
  double f, fl, fu;
  if (y - t0_ <= 0) {
    return 0.0;
  } else {
    fl = exp(internal::dwiener5(y - t0_, a, v, w_lower, sv, lerr));
    fu = exp(internal::dwiener5(y - t0_, a, v, w_upper, sv, lerr));
  }
  return 1 / sw_mean * 0.5 * (fl + fu);
}
//-----------------------------------------------

// integrand d/dst0 (on normal scale)
template <typename T_x, typename T_p>
return_type_t<T_x, T_p> int_dst0ddiff(const T_x& x, const T_p& p) {
  struct my_params3 {
  double y, a, v, w, t0, t0_mean, sv, sw, st0, st0_mean, lerr;
};
  my_params3* params = static_cast<my_params3*>(p);
  double y = (params->y);
  double a = (params->a);
  double v = (params->v);
  double t0 = (params->t0);
  double t0_mean = (params->t0_mean);
  double w = (params->w);
  double sw = (params->sw);
  double sv = (params->sv);
  double st0 = (params->st0);
  double st0_mean = (params->st0_mean);
  double lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);

  double t0_ = sw ? (st0_mean ? t0_mean + st0_mean * x_vec[1] : t0_mean)
                  : (st0_mean ? t0_mean + st0_mean * x_vec[0] : t0_mean);
  double omega = sw ? w + sw * (x_vec[0] - 0.5) : w;
  double f;
  if (y - t0_ <= 0) {
    return 0.0;
  } else {
    f = exp(internal::dwiener5(y - t0_, a, v, omega, sv, lerr));
  }
  return 1 / st0 * f;
}
//-----------------------------------------------

}  // namespace internal

template <bool propto, typename T_y, typename T_a, typename T_t0, typename T_w,
          typename T_v, typename T_sv, typename T_sw, typename T_st0>
inline return_type_t<T_y, T_a, T_t0, T_w, T_v, T_sv, T_sw, T_st0> wiener7_lpdf(
    const T_y& y, const T_a& a, const T_t0& t0, const T_w& w, const T_v& v,
    const T_sv& sv, const T_sw& sw, const T_st0& st0, const double& prec) {
  using T_y_ref = ref_type_t<T_y>;
  using T_a_ref = ref_type_t<T_a>;
  using T_v_ref = ref_type_t<T_v>;
  using T_w_ref = ref_type_t<T_w>;
  using T_t0_ref = ref_type_t<T_t0>;
  using T_sv_ref = ref_type_t<T_sv>;
  using T_sw_ref = ref_type_t<T_sw>;
  using T_st0_ref = ref_type_t<T_st0>;

  const char* function_name = "wiener7_lpdf";
  check_consistent_sizes(function_name, "Random variable", y,
                         "Boundary separation", a, "Drift rate", v,
                         "A-priori bias", w, "Nondecision time", t0,
                         "Inter-trial variability in drift rate", sv,
                         "Inter-trial variability in A-priori bias", sw,
                         "Inter-trial variability in Nondecision time", st0);
  check_consistent_size(function_name, "Random variable", y, 1);
  check_consistent_size(function_name, "Boundary separation", a, 1);
  check_consistent_size(function_name, "Drift rate", v, 1);
  check_consistent_size(function_name, "A-priori bias", w, 1);
  check_consistent_size(function_name, "Nondecision time", t0, 1);
  check_consistent_size(function_name, "Inter-trial variability in drift rate",
                        sv, 1);
  check_consistent_size(function_name,
                        "Inter-trial variability in A-priori bias", sw, 1);
  check_consistent_size(function_name,
                        "Inter-trial variability in Nondecision time", st0, 1);

  T_y_ref y_ref = y;
  T_a_ref a_ref = a;
  T_v_ref v_ref = v;
  T_w_ref w_ref = w;
  T_t0_ref t0_ref = t0;
  T_sv_ref sv_ref = sv;
  T_sw_ref sw_ref = sw;
  T_st0_ref st0_ref = st0;

  check_positive_finite(function_name, "Random variable", value_of(y_ref));
  check_positive_finite(function_name, "Boundary separation", value_of(a_ref));
  check_finite(function_name, "Drift rate", value_of(v_ref));
  check_less(function_name, "A-priori bias", value_of(w_ref), 1);
  check_greater(function_name, "A-priori bias", value_of(w_ref), 0);
  check_nonnegative(function_name, "Nondecision time", value_of(t0_ref));
  check_finite(function_name, "Nondecision time", value_of(t0_ref));
  check_nonnegative(function_name, "Inter-trial variability in drift rate",
                    value_of(sv_ref));
  check_finite(function_name, "Inter-trial variability in drift rate",
               value_of(sv_ref));
  check_bounded(function_name, "Inter-trial variability in A-priori bias",
                value_of(sw_ref), 0, 1);
  check_nonnegative(function_name,
                    "Inter-trial variability in Nondecision time",
                    value_of(st0_ref));
  check_finite(function_name, "Inter-trial variability in Nondecision time",
               value_of(st0_ref));

  if (size_zero(y, a, v, w, t0) || size_zero(sv, sw, st0)) {
    return 0;
  }
  size_t N = max_size(y, a, v, w, t0, sv, sw, st0);
  if (!N) {
    return 0;
  }
  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_a_ref> a_vec(a_ref);
  scalar_seq_view<T_v_ref> v_vec(v_ref);
  scalar_seq_view<T_w_ref> w_vec(w_ref);
  scalar_seq_view<T_t0_ref> t0_vec(t0_ref);
  scalar_seq_view<T_sv_ref> sv_vec(sv_ref);
  scalar_seq_view<T_sw_ref> sw_vec(sw_ref);
  scalar_seq_view<T_st0_ref> st0_vec(st0_ref);
  size_t N_y_t0 = max_size(y, t0, st0);

  for (size_t i = 0; i < N_y_t0; ++i) {
    if (y_vec[i] <= t0_vec[i]) {
      std::stringstream msg;
      msg << ", but must be greater than nondecision time = " << t0_vec[i];
      std::string msg_str(msg.str());
      throw_domain_error(function_name, "Random variable", y_vec[i], " = ",
                         msg_str.c_str());
    }
  }
  size_t N_beta_sw = max_size(w, sw);
  for (size_t i = 0; i < N_beta_sw; ++i) {
    if (w_vec[i] - .5 * sw_vec[i] <= 0) {
      std::stringstream msg;
      msg << ", but must be smaller than 2*(A-priori bias) = " << 2 * w_vec[i];
      std::string msg_str(msg.str());
      throw_domain_error(function_name,
                         "Inter-trial variability in A-priori bias", sw_vec[i],
                         " = ", msg_str.c_str());
    }
    if (w_vec[i] + .5 * sw_vec[i] >= 1) {
      std::stringstream msg;
      msg << ", but must be smaller than 2*(1-A-priori bias) = "
          << 2 * (1 - w_vec[i]);
      std::string msg_str(msg.str());
      throw_domain_error(function_name,
                         "Inter-trial variability in A-priori bias", sw_vec[i],
                         " = ", msg_str.c_str());
    }
  }
  if (!include_summand<propto, T_y, T_a, T_v, T_w, T_t0, T_sv, T_sw,
                       T_st0>::value) {
    return 0;
  }

  double error_bound_dens = 1e-6;  // precision for density
  double lerror_bound_dens = log(error_bound_dens);
  double error_bound = prec;  // precision for
  // derivatives (controllable by user)
  double lerror_bound = log(error_bound);  // log(alpha)
  double abstol = 0.0;
  double reltol = .9 * error_bound;  // eps_rel(Integration)
  double abstol_wiener5 = 1e-12;     // eps_abs(wiener5)
  double labstol_wiener5 = log(abstol_wiener5);
  // log(eps_abs(wiener5)
  int Meval = 6000;
  double dens = 0.0;
  double ld = 0.0;
  operands_and_partials<T_y_ref, T_a_ref, T_t0_ref, T_w_ref, T_v_ref, T_sv_ref,
                        T_sw_ref, T_st0_ref>
      ops_partials(y_ref, a_ref, t0_ref, w_ref, v_ref, sv_ref, sw_ref, st0_ref);
  static constexpr double LOG_FOUR = LOG_TWO + LOG_TWO;
  static constexpr double LOG_POINT1 = -1;
  struct my_params {
  double y, a, v, w, t0, sv, sw, st0, lerr;
};
  struct my_params2 {
  double y, a, v, w, w_lower, w_upper, t0, sv, sw, sw_mean, st0, lerr;
};
  struct my_params3 {
  double y, a, v, w, t0, t0_mean, sv, sw, st0, st0_mean, lerr;
};

  // calculate density and partials
  for (size_t i = 0; i < N; i++) {
    const double y_val = y_vec.val(i);
    const double a_val = a_vec.val(i);
    const double v_val = v_vec.val(i);
    const double w_val = w_vec.val(i);
    const double t0_val = t0_vec.val(i);
    const double sv_val = sv_vec.val(i);
    const double sw_val = sw_vec.val(i);
    const double st0_val = st0_vec.val(i);
    my_params params = {y_val,  a_val,   v_val,
                                  w_val,  t0_val,  sv_val,
                                  sw_val, st0_val, labstol_wiener5 - LOG_TWO};
    int dim = (sw_val != 0) + (st0_val != 0);
    check_positive(function_name,
                   "(Inter-trial variability in A-priori bias) + "
                   "(Inter-trial variability in nondecision time)",
                   dim);

    std::vector<double> xmin(dim, 0);
    std::vector<double> xmax(dim, 1);
    if (st0_val) {
      xmax[dim - 1] = fmin(1.0, (y_val - t0_val) / st0_val);
    }

    dens = hcubature(internal::int_ddiff<std::vector<double>, void*>, &params,
                     dim, xmin, xmax, Meval, abstol, reltol / 2);
    if (labstol_wiener5
        > fabs(log(dens)) + LOG_POINT1 + lerror_bound_dens - LOG_TWO) {
      double new_error
          = LOG_POINT1 + lerror_bound_dens - LOG_TWO + log(fabs(dens));
      my_params params_new_error
          = {y_val,  a_val,  v_val,   w_val,    t0_val,
             sv_val, sw_val, st0_val, new_error};
      dens = hcubature(internal::int_ddiff<std::vector<double>, void*>,
                       &params_new_error, dim, xmin, xmax, Meval, abstol,
                       reltol);
    }
    ld += log(dens);

    // computation of derivative for t and precision check in order to give
    // the value as deriv_t to edge1 and as -deriv_t to edge5
    double deriv_t_7;
    deriv_t_7
        = 1 / dens
          * hcubature(internal::int_dtddiff<std::vector<double>, void*>,
                      &params, dim, xmin, xmax, Meval, abstol, reltol / 2);
    if (labstol_wiener5
        > log(fabs(deriv_t_7)) + LOG_POINT1 + lerror_bound - LOG_TWO) {
      double new_error
          = LOG_POINT1 + lerror_bound - LOG_FOUR + log(fabs(deriv_t_7));
      my_params params_new_error
          = {y_val,  a_val,  v_val,   w_val,    t0_val,
             sv_val, sw_val, st0_val, new_error};
      deriv_t_7 = 1 / dens
                  * hcubature(internal::int_dtddiff<std::vector<double>, void*>,
                              &params_new_error, dim, xmin, xmax, Meval, abstol,
                              reltol / 2);
    }

    // computation of derivatives and precision checks
    double deriv;
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[i] = deriv_t_7;
    }
    if (!is_constant_all<T_a>::value) {
      deriv = 1 / dens
              * hcubature(internal::int_daddiff<std::vector<double>, void*>,
                          &params, dim, xmin, xmax, Meval, abstol, reltol / 2);
      if (labstol_wiener5
          > log(fabs(deriv)) + LOG_POINT1 + lerror_bound - LOG_TWO) {
        double new_error
            = LOG_POINT1 + lerror_bound - LOG_FOUR + log(fabs(deriv));
        my_params params_new_error
            = {y_val,  a_val,  v_val,   w_val,    t0_val,
               sv_val, sw_val, st0_val, new_error};
        deriv = 1 / dens
                * hcubature(internal::int_daddiff<std::vector<double>, void*>,
                            &params_new_error, dim, xmin, xmax, Meval, abstol,
                            reltol / 2);
      }
      ops_partials.edge2_.partials_[i] = deriv;
    }
    if (!is_constant_all<T_t0>::value) {
      ops_partials.edge3_.partials_[i] = -deriv_t_7;
    }
    if (!is_constant_all<T_w>::value) {
      deriv = 1 / dens
              * hcubature(internal::int_dwddiff<std::vector<double>, void*>,
                          &params, dim, xmin, xmax, Meval, abstol, reltol / 2);
      if (labstol_wiener5
          > log(fabs(deriv)) + LOG_POINT1 + lerror_bound - LOG_TWO) {
        double new_error
            = LOG_POINT1 + lerror_bound - LOG_FOUR + log(fabs(deriv));
        my_params params_new_error
            = {y_val,  a_val,  v_val,   w_val,    t0_val,
               sv_val, sw_val, st0_val, new_error};
        deriv = 1 / dens
                * hcubature(internal::int_dwddiff<std::vector<double>, void*>,
                            &params_new_error, dim, xmin, xmax, Meval, abstol,
                            reltol / 2);
      }
      ops_partials.edge4_.partials_[i] = deriv;
    }
    if (!is_constant_all<T_v>::value) {
      deriv = 1 / dens
              * hcubature(internal::int_dvddiff<std::vector<double>, void*>,
                          &params, dim, xmin, xmax, Meval, abstol, reltol / 2);
      if (labstol_wiener5
          > log(fabs(deriv)) + LOG_POINT1 + lerror_bound - LOG_TWO) {
        double new_error
            = LOG_POINT1 + lerror_bound - LOG_TWO + log(fabs(deriv));
        my_params params_new_error
            = {y_val,  a_val,  v_val,   w_val,    t0_val,
               sv_val, sw_val, st0_val, new_error};
        deriv = 1 / dens
                * hcubature(internal::int_dvddiff<std::vector<double>, void*>,
                            &params_new_error, dim, xmin, xmax, Meval, abstol,
                            reltol / 2);
      }
      ops_partials.edge5_.partials_[i] = deriv;
    }
    if (!is_constant_all<T_sv>::value) {
      deriv = 1 / dens
              * hcubature(internal::int_dsvddiff<std::vector<double>, void*>,
                          &params, dim, xmin, xmax, Meval, abstol, reltol / 2);
      if (labstol_wiener5
          > log(fabs(deriv)) + LOG_POINT1 + lerror_bound - LOG_TWO) {
        double new_error
            = LOG_POINT1 + lerror_bound - LOG_TWO + log(fabs(deriv));
        my_params params_new_error
            = {y_val,  a_val,  v_val,   w_val,    t0_val,
               sv_val, sw_val, st0_val, new_error};
        deriv = 1 / dens
                * hcubature(internal::int_dsvddiff<std::vector<double>, void*>,
                            &params_new_error, dim, xmin, xmax, Meval, abstol,
                            reltol / 2);
      }
      ops_partials.edge6_.partials_[i] = deriv;
    }
    if (!is_constant_all<T_sw>::value) {
      if (sw_val == 0) {
        ops_partials.edge7_.partials_[i] = 0;
      } else {
        double lower, upper, width, fl, fu;
        lower = w_val - sw_val / 2;
        lower = (0 > lower) ? 0 : lower;
        upper = w_val + sw_val / 2;
        upper = (1 < upper) ? 1 : upper;
        width = upper - lower;

        int dim_ = (st0_val != 0);
        if (dim_ == 0) {
          fl = exp(internal::dwiener5(y_val - t0_val, a_val, v_val, lower,
                                      sv_val, labstol_wiener5));
          fu = exp(internal::dwiener5(y_val - t0_val, a_val, v_val, upper,
                                      sv_val, labstol_wiener5));
          if (labstol_wiener5 > log(fabs(fl + fu) + log(sw_val) - LOG_TWO
                                    + lerror_bound - LOG_TWO)) {
            fl = exp(internal::dwiener5(y_val - t0_val, a_val, v_val, lower,
                                        sv_val,
                                        lerror_bound - LOG_TWO + log(sw_val)
                                            - LOG_TWO + log(fabs(fl + fu))));
            fu = exp(internal::dwiener5(y_val - t0_val, a_val, v_val, upper,
                                        sv_val,
                                        lerror_bound - LOG_TWO + log(sw_val)
                                            - LOG_TWO + log(fabs(fl + fu))));
          }
          deriv = 1 / width * 0.5 * (fl + fu);
        } else {
          my_params2 params_sw
              = {y_val, a_val,  v_val,   w_val,
                 lower, upper,  t0_val,  sv_val,
                 0,     sw_val, st0_val, labstol_wiener5 - LOG_TWO};
          deriv = hcubature(internal::int_dswddiff<std::vector<double>, void*>,
                            &params_sw, dim_, xmin, xmax, Meval, abstol,
                            reltol / 2);
          if (labstol_wiener5
              > log(fabs(deriv) + LOG_POINT1 + lerror_bound - LOG_TWO)) {
            double new_error
                = log(fabs(deriv)) + LOG_POINT1 + lerror_bound - LOG_TWO;
            my_params2 params_new_error_sw
                = {y_val,  a_val,  v_val, w_val,  lower,   upper,
                   t0_val, sv_val, 0,     sw_val, st0_val, new_error};
            deriv
                = hcubature(internal::int_dswddiff<std::vector<double>, void*>,
                            &params_new_error_sw, dim_, xmin, xmax, Meval,
                            abstol, reltol / 2);
          }
        }
        ops_partials.edge7_.partials_[i] = deriv / dens - 1 / sw_val;
      }
    }
    if (!is_constant_all<T_st0>::value) {
      double f;
      if (st0_val == 0) {
        ops_partials.edge8_.partials_[i] = 0;
      } else if (y_val - (t0_val + st0_val) <= 0) {
        ops_partials.edge8_.partials_[i] = -1 / st0_val;
      } else {
        int dim_ = (sw_val != 0);
        if (dim_ == 0) {
          f = exp(internal::dwiener5(y_val - (t0_val + st0_val), a_val, v_val,
                                     w_val, sv_val, labstol_wiener5));
          if (labstol_wiener5 > log(fabs(f) + log(st0_val) - LOG_TWO
                                    + lerror_bound - LOG_TWO)) {
            f = exp(internal::dwiener5(y_val - (t0_val + st0_val), a_val, v_val,
                                       w_val, sv_val,
                                       lerror_bound - LOG_TWO + log(st0_val)
                                           - LOG_TWO + log(fabs(f))));
          }
          deriv = 1 / st0_val * f;
        } else {
          double new_error = labstol_wiener5 - LOG_TWO;
          my_params3 params_st
              = {y_val,  a_val,  v_val,   w_val, t0_val,   t0_val + st0_val,
                 sv_val, sw_val, st0_val, 0,     new_error};
          deriv = hcubature(internal::int_dst0ddiff<std::vector<double>, void*>,
                            &params_st, dim_, xmin, xmax, Meval, abstol,
                            reltol / 2);
          if (labstol_wiener5
              > log(fabs(deriv) + LOG_POINT1 + lerror_bound - LOG_TWO)) {
            double new_error
                = log(fabs(deriv)) + LOG_POINT1 + lerror_bound - LOG_TWO;
            my_params3 params_new_error_st
                = {y_val,  a_val,  v_val,   w_val, t0_val,   t0_val + st0_val,
                   sv_val, sw_val, st0_val, 0,     new_error};
            deriv
                = hcubature(internal::int_dst0ddiff<std::vector<double>, void*>,
                            &params_new_error_st, dim_, xmin, xmax, Meval,
                            abstol, reltol / 2);
          }
        }
        ops_partials.edge8_.partials_[i] = -1 / st0_val + deriv / dens;
      }
    }
    std::vector<double>().swap(xmin);
    std::vector<double>().swap(xmax);
  }
  return ops_partials.build(ld);
}
}  // namespace math
}  // namespace stan
#endif
