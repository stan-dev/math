#ifndef STAN_MATH_PRIM_FUN_WIENER5_LPDF_HPP
#define STAN_MATH_PRIM_FUN_WIENER5_LPDF_HPP

#include <stan/math/prim/fun.hpp>

namespace stan {
namespace math {
namespace internal {

// calculate density in log
double dwiener5(const double& y, const double& a, const double& vn,
                const double& wn, const double& sv, const double& err) {
  double kll, kss, ans, v, w;
  w = 1.0 - wn;
  v = -vn;
  double y_asq = y / square(a);
  ans = 0.0;

  // calculate the number of terms needed for short t
  double lg1;
  if (sv != 0) {
    double sv_sqr = square(sv);
    double one_plus_svsqr_y = 1 + sv_sqr * y;
    lg1 = (sv_sqr * square(a * w) - 2 * a * v * w - square(v) * y) / 2.0
              / one_plus_svsqr_y
          - 2 * log(a) - 0.5 * log(one_plus_svsqr_y);
  } else {
    lg1 = (-2 * a * v * w - square(v) * y) / 2.0 - 2 * log(a);
  }
  double es = (err - lg1);
  double K1 = (sqrt(2.0 * y_asq) + w) / 2.0;
  double u_eps = fmin(-1.0, LOG_TWO + LOG_PI + 2.0 * log(y_asq) + 2.0 * (es));
  double arg = -y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  double K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate the number of terms needed for large t
  double el = es;
  K1 = 1.0 / (pi() * sqrt(y_asq));
  K2 = 0.0;
  double two_log_piy = -2.0 * (log(pi() * y_asq) + el);
  static const double PISQ = square(pi());  // pi*pi
  if (two_log_piy >= 0) {
    K2 = sqrt(two_log_piy / (PISQ * y_asq));
  }
  kll = ceil(fmax(K1, K2));

  // if small t is better
  if (2 * kss <= kll) {
    double fplus = NEGATIVE_INFTY;
    double fminus = NEGATIVE_INFTY;
    double twoy = 2.0 * y_asq;
    if (static_cast<size_t>(kss) > 0) {
      for (size_t k = static_cast<size_t>(kss); k >= 1; k--) {
        double w_plus_2k = w + 2.0 * k;
        double w_minus_2k = w - 2.0 * k;

        fplus = log_sum_exp(log(w_plus_2k) - square(w_plus_2k) / twoy, fplus);
        fminus
            = log_sum_exp(log(-w_minus_2k) - square(w_minus_2k) / twoy, fminus);
      }
    }
    fplus = log_sum_exp(log(w) - square(w) / twoy, fplus);
    ans = lg1
          + (-0.5 * LOG_TWO - LOG_SQRT_PI - 1.5 * log(y_asq)
             + log_diff_exp(fplus, fminus));
    // if large t is better
  } else {
    double fplus = NEGATIVE_INFTY;
    double fminus = NEGATIVE_INFTY;
    double halfy = y_asq / 2.0;
    for (size_t k = static_cast<size_t>(kll); k >= 1; k--) {
      double pi_k = k * pi();
      double check = sin(pi_k * w);
      if (check > 0) {
        fplus = log_sum_exp(log(k) - square(pi_k) * halfy + log(check), fplus);
      } else {
        fminus
            = log_sum_exp(log(k) - square(pi_k) * halfy + log(-check), fminus);
      }
    }
    if (fplus < fminus) {
      ans = NEGATIVE_INFTY;
    } else {
      ans = lg1 + log_diff_exp(fplus, fminus) + LOG_PI;
    }
  }
  return ans;
}
//-----------------------------------------------

// d/dt DENSITY
// calculate derivative of density with respect to t (in log, ans =
// d/dt(log(f))=d/dt f'/f; ans*exp(ld)=f' on normal scale)
double dtdwiener5(const double& y, const double& a, const double& vn,
                  const double& wn, const double& sv, const double& err) {
  double kll, kss, ans, v, w;
  w = 1.0 - wn;
  v = -vn;

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
  double ld = dwiener5(y, a, vn, wn, sv,
                       err - log(max(fabs(ans0 - 1.5 / y), fabs(ans0))));

  // calculate the number of terms kss needed for small t
  double es = err - lg1;
  es = es + la;
  double K1 = (sqrt(3.0 * y_asq) + w) / 2.0;
  double u_eps = fmin(
      -1.0, (log(8.0 / 27.0) + LOG_PI + 4.0 * log(y_asq) + 2.0 * es) / 3.0);
  double arg = -3.0 * y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  double K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate number of terms kll needed for large t
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

  // if small t is better
  if (2 * kss < kll) {
    // calculate terms of the sum for small t
    double twoy = 2.0 * y_asq;
    if (static_cast<size_t>(kss) > 0) {
      for (size_t k = static_cast<size_t>(kss); k >= 1; k--) {
        double w_plus_2k = w + 2.0 * k;
        double w_minus_2k = w - 2.0 * k;
        fplus = log_sum_exp(3.0 * log(w_plus_2k) - w_plus_2k * w_plus_2k / twoy,
                            fplus);
        fminus = log_sum_exp(
            3.0 * log(-w_minus_2k) - w_minus_2k * w_minus_2k / twoy, fminus);
      }
    }
    fplus = log_sum_exp(3.0 * log(w) - w * w / twoy, fplus);
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
    // if large t is better
  } else {
    // calculate terms of the sum for large t
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
  return ans;
}
//-----------------------------------------------

// d/da DENSITY
// calculate derivative of density with respect to a (in log, ans =
// d/da(log(f))=d/da f'/f; ans*exp(ld)=f' on normal scale)
double dadwiener5(const double& y, const double& a, const double& vn,
                  const double& wn, const double& sv, const double& err,
                  const int& normal_or_log) {
  double kll, kss, ans, v, w;

  double la = log(a);
  double ly = log(y);
  w = 1.0 - wn;
  v = -vn;

  // prepare some variables
  double y_asq = y / square(a);
  double ans0, lg1;
  if (sv != 0) {
    double sv_sqr = square(sv);
    double one_plus_svsqr_y = (1 + sv_sqr * y);
    ans0 = (-v * w + sv_sqr * square(w) * a) / one_plus_svsqr_y;
    lg1 = (sv_sqr * square(a * w) - 2 * a * v * w - square(v) * y) / 2.0
              / one_plus_svsqr_y
          - 2 * la - 0.5 * log(one_plus_svsqr_y);
  } else {
    ans0 = -v * w;
    lg1 = (-2 * a * v * w - square(v) * y) / 2.0 - 2 * la;
  }
  double factor = lg1 - 3 * la;
  double ld
      = dwiener5(y, a, vn, wn, sv,
                 err - log(max(fabs(ans0 + 1.0 / a), fabs(ans0 - 2.0 / a))));

  // calculate the number of terms kss needed for small t
  double es = err - lg1;
  es = es + la;
  es = es - LOG_TWO + 2.0 * la - ly;
  double K1 = (sqrt(3.0 * y_asq) + w) / 2.0;
  double u_eps = fmin(
      -1.0, (log(8.0 / 27.0) + LOG_PI + 4.0 * log(y_asq) + 2.0 * es) / 3.0);
  double arg = -3.0 * y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  double K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate number of terms kll needed for large t
  double el = err - lg1;
  el = el + la;
  el = el - LOG_TWO + 2.0 * la - ly;
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

  // if small t is better
  if (2 * kss < kll) {
    // calculate terms of the sum for short t
    double twoy = 2.0 * y_asq;
    if (static_cast<int>(kss) > 0) {
      for (size_t k = static_cast<size_t>(kss); k >= 1; k--) {
        double w_plus_2k = w + 2.0 * k;
        double w_minus_2k = w - 2.0 * k;
        fplus = log_sum_exp(3.0 * log(w_plus_2k) - w_plus_2k * w_plus_2k / twoy,
                            fplus);
        fminus = log_sum_exp(
            3.0 * log(-w_minus_2k) - w_minus_2k * w_minus_2k / twoy, fminus);
      }
    }
    fplus = log_sum_exp(3.0 * log(w) - w * w / twoy, fplus);
    if (fplus < fminus) {
      newsign = -1;
      erg = log_diff_exp(fminus, fplus);
    } else {
      erg = log_diff_exp(fplus, fminus);
    }
    ans = ans0 + 1.0 / a
          - newsign
                * exp(-0.5 * LOG_TWO - LOG_SQRT_PI - 2.5 * ly + 4.0 * la + lg1
                      + erg - ld);
    // if large t is better
  } else {
    // calculate terms of the sum for large t
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
    if (fplus > fminus) {
      erg = log_diff_exp(fplus, fminus);
    } else {
      erg = log_diff_exp(fminus, fplus);
      newsign = -1;
    }
    ans = ans0 - 2.0 / a + newsign * exp(ly + factor + 3.0 * LOG_PI + erg - ld);
  }
  if (normal_or_log == 1) {
    return ans * exp(ld);  // derivative of f for hcubature
  } else {
    return ans;  // derivative of log(f)
  }
}
//-----------------------------------------------

// d/dv DENSITY
// calculate derivative of density with respect to v (in log, ans =
// d/dv(log(f))=d/dv f'/f; ans*exp(ld)=f' on normal scale)
double dvdwiener5(const double& y, const double& a, const double& vn,
                  const double& wn, const double& sv) {
  double ans;
  if (sv != 0) {
    ans = 1 + square(sv) * y;
    ans = (a * (1 - wn) - vn * y) / ans;
  } else {
    ans = (a * (1 - wn) - vn * y);
  }
  return ans;
}
//-----------------------------------------------

// d/dw DENSITY
// calculate derivative of density with respect to w (in log, ans =
// d/dw(log(f))=d/dw f'/f; ans*exp(ld)=f' on normal scale)
double dwdwiener5(const double& y, const double& a, const double& vn,
                  const double& wn, const double& sv, const double& err,
                  const int& normal_or_log) {
  double kll, kss, ans, v, w;
  double sign = -1;
  w = 1.0 - wn;
  v = -vn;

  // prepare some variables
  double y_asq = y / square(a);
  double ans0, lg1;
  if (sv != 0) {
    double sv_sqr = square(sv);
    double one_plus_svsqr_y = (1 + sv_sqr * y);
    ans0 = (-v * a + sv_sqr * square(a) * w) / one_plus_svsqr_y;
    lg1 = (sv_sqr * square(a * w) - 2 * a * v * w - square(v) * y) / 2.0
              / one_plus_svsqr_y
          - 2 * log(a) - 0.5 * log(one_plus_svsqr_y);
  } else {
    ans0 = -v * a;
    lg1 = (-2 * a * v * w - square(v) * y) / 2.0 - 2 * log(a);
  }
  double ld = dwiener5(y, a, vn, wn, sv, err - log(fabs(ans0)));
  double ls = -lg1 + ld;
  double ll = -lg1 + ld;

  // calculate the number of terms kss needed for small t
  double K1 = (sqrt(3.0 * y_asq) + w) / 2.0;
  double u_eps
      = fmin(-1.0, 2.0 * (err - lg1) + LOG_TWO + LOG_PI + 2.0 * log(y_asq));
  double arg = -y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  double K2 = (arg > 0) ? 0.5 * (sqrt(arg) + w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate number of terms kll needed for large t
  K1 = sqrt(2.0 / y_asq) / pi();
  static const double TWO_LOG_PI = 2.0 * LOG_PI;
  u_eps = fmin(
      -1.0, log(4.0 / 9.0) + TWO_LOG_PI + 3.0 * log(y_asq) + 2.0 * (err - lg1));
  arg = -(u_eps - sqrt(-2.0 * u_eps - 2.0));
  K2 = (arg > 0) ? 1.0 / pi() * sqrt(arg / y_asq) : K1;
  kll = ceil(fmax(K1, K2));

  double erg;
  double newsign = 1;
  double fplus = NEGATIVE_INFTY;
  double fminus = NEGATIVE_INFTY;

  // if small t is better
  if (2 * kss < kll) {
    // calculate terms of the sum for short t
    double twoy = 2.0 * y_asq;
    for (size_t k = static_cast<size_t>(kss); k >= 1; k--) {
      double sqrt_w_plus_2k = square(w + 2 * k);
      double sqrt_w_minus_2k = square(w - 2 * k);
      double wp2k_minusy = sqrt_w_plus_2k - y_asq;
      double wm2k_minusy = sqrt_w_minus_2k - y_asq;
      if (wp2k_minusy > 0) {
        fplus = log_sum_exp(log(wp2k_minusy) - sqrt_w_plus_2k / twoy, fplus);
      } else if (wp2k_minusy < 0) {
        fminus
            = log_sum_exp(log(-(wp2k_minusy)) - sqrt_w_plus_2k / twoy, fminus);
      }
      if (wm2k_minusy > 0) {
        fplus = log_sum_exp(log(wm2k_minusy) - sqrt_w_minus_2k / twoy, fplus);
      } else if (wm2k_minusy < 0) {
        fminus
            = log_sum_exp(log(-(wm2k_minusy)) - sqrt_w_minus_2k / twoy, fminus);
      }
    }
    double sqr_w = square(w);
    double sqrt_w_plus_2k = sqr_w - y_asq;
    if (sqrt_w_plus_2k > 0) {
      fplus = log_sum_exp(log(sqrt_w_plus_2k) - sqr_w / twoy, fplus);
    } else if (sqrt_w_plus_2k < 0) {
      fminus = log_sum_exp(log(-(sqrt_w_plus_2k)) - sqr_w / twoy, fminus);
    }
    if (fplus < fminus) {
      newsign = -1;
      erg = log_diff_exp(fminus, fplus);
    } else {
      erg = log_diff_exp(fplus, fminus);
    }
    ans = ans0
          - newsign
                * exp(erg - ls - 2.5 * log(y_asq) - 0.5 * LOG_TWO
                      - 0.5 * LOG_PI);
    // if large t is better
  } else {
    // calculate terms of the sum for large t
    double halfy = y_asq / 2.0;
    for (size_t k = static_cast<size_t>(kll); k >= 1; k--) {
      double pi_k = pi() * k;
      double x = cos(pi_k * w);
      if (x > 0) {
        fplus
            = log_sum_exp(2.0 * log(k) - square(pi_k) * halfy + log(x), fplus);
      } else if (x < 0) {
        fminus = log_sum_exp(2.0 * log(k) - square(pi_k) * halfy + log(-x),
                             fminus);
      }
    }
    if (fplus < fminus) {
      erg = log_diff_exp(fminus, fplus);
      newsign = -1;
    } else {
      erg = log_diff_exp(fplus, fminus);
    }
    ans = ans0 + newsign * exp(erg - ll + TWO_LOG_PI);
  }
  if (normal_or_log == 1) {
    return ans * sign * exp(ld);  // derivative of f for hcubature
  } else {
    return ans * sign;  // derivative of log(f)
  }
}
//-----------------------------------------------

// d/dsv DENSITY
// calculate derivative of density with respect to sv (in log, ans =
// d/dsv(log(f))=d/dsv f'/f; ans*exp(ld)=f' on normal scale)
double dsvdwiener5(const double& y, const double& a, const double& vn,
                   const double& wn, const double& sv) {
  double v, w;
  v = -vn;
  w = 1 - wn;
  double one_sqrsv_y = 1 + square(sv) * y;
  double t1 = -y / one_sqrsv_y;
  double t2 = (square(a * w) + 2 * a * v * w * y + square(v * y))
              / square(one_sqrsv_y);
  return sv * (t1 + t2);
}
//-----------------------------------------------

}  // namespace internal

template <bool propto, typename T_y, typename T_a, typename T_t0, typename T_w,
          typename T_v, typename T_sv>
inline return_type_t<T_y, T_a, T_t0, T_w, T_v, T_sv> wiener5_lpdf(
    const T_y& y, const T_a& a, const T_t0& t0, const T_w& w, const T_v& v,
    const T_sv& sv, const double& prec) {
  using T_y_ref = ref_type_t<T_y>;
  using T_a_ref = ref_type_t<T_a>;
  using T_t0_ref = ref_type_t<T_t0>;
  using T_w_ref = ref_type_t<T_w>;
  using T_v_ref = ref_type_t<T_v>;
  using T_sv_ref = ref_type_t<T_sv>;

  const char* function_name = "wiener5_lpdf";
  check_consistent_sizes(function_name, "Random variable", y,
                         "Boundary separation", a, "Drift rate", v,
                         "A-priori bias", w, "Nondecision time", t0,
                         "Inter-trial variability in drift rate", sv);
  check_consistent_size(function_name, "Random variable", y, 1);
  check_consistent_size(function_name, "Boundary separation", a, 1);
  check_consistent_size(function_name, "Nondecision time", t0, 1);
  check_consistent_size(function_name, "A-priori bias", w, 1);
  check_consistent_size(function_name, "Drift rate", v, 1);
  check_consistent_size(function_name, "Inter-trial variability in drift rate",
                        sv, 1);

  T_y_ref y_ref = y;
  T_a_ref a_ref = a;
  T_t0_ref t0_ref = t0;
  T_w_ref w_ref = w;
  T_v_ref v_ref = v;
  T_sv_ref sv_ref = sv;

  check_positive_finite(function_name, "Random variable", value_of(y_ref));
  check_positive_finite(function_name, "Boundary separation", value_of(a_ref));
  check_nonnegative(function_name, "Nondecision time", value_of(t0_ref));
  check_finite(function_name, "Nondecision time", value_of(t0_ref));
  check_less(function_name, "A-priori bias", value_of(w_ref), 1);
  check_greater(function_name, "A-priori bias", value_of(w_ref), 0);
  check_finite(function_name, "Drift rate", value_of(v_ref));
  check_nonnegative(function_name, "Inter-trial variability in drift rate",
                    value_of(sv_ref));
  check_finite(function_name, "Inter-trial variability in drift rate",
               value_of(sv_ref));

  if (size_zero(y, a, t0, w, v) || size_zero(sv)) {
    return 0;
  }
  size_t N = max_size(y, a, t0, w, v, sv);
  if (!N) {
    return 0;
  }
  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_a_ref> a_vec(a_ref);
  scalar_seq_view<T_t0_ref> t0_vec(t0_ref);
  scalar_seq_view<T_w_ref> w_vec(w_ref);
  scalar_seq_view<T_v_ref> v_vec(v_ref);
  scalar_seq_view<T_sv_ref> sv_vec(sv_ref);
  size_t N_y_t0 = max_size(y, t0);

  for (size_t i = 0; i < N_y_t0; ++i) {
    if (y_vec[i] <= t0_vec[i]) {
      std::stringstream msg;
      msg << ", but must be greater than nondecision time = " << t0_vec[i];
      std::string msg_str(msg.str());
      throw_domain_error(function_name, "Random variable", y_vec[i], " = ",
                         msg_str.c_str());
    }
  }

  if (!include_summand<propto, T_y, T_a, T_t0, T_w, T_v, T_sv>::value) {
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
  operands_and_partials<T_y_ref, T_a_ref, T_t0_ref, T_w_ref, T_v_ref, T_sv_ref>
      ops_partials(y_ref, a_ref, t0_ref, w_ref, v_ref, sv_ref);

  static constexpr double LOG_FOUR = LOG_TWO + LOG_TWO;
  static constexpr double LOG_POINT1 = -1;

  // calculate density and partials
  for (size_t i = 0; i < N; i++) {
    // Calculate 4-parameter model without inter-trial variabilities (if
    // sv_vec[i] == 0) or 5-parameter model with inter-trial variability in
    // drift rate (if sv_vec[i] != 0)
    const double y_val = y_vec.val(i);
    const double a_val = a_vec.val(i);
    const double t0_val = t0_vec.val(i);
    const double w_val = w_vec.val(i);
    const double v_val = v_vec.val(i);
    const double sv_val = sv_vec.val(i);

    dens = internal::dwiener5(y_val - t0_val, a_val, v_val, w_val, sv_val,
                              labstol_wiener5);
    if (labstol_wiener5 > fabs(dens) + lerror_bound_dens - LOG_TWO) {
      dens = internal::dwiener5(y_val - t0_val, a_val, v_val, w_val, sv_val,
                                fabs(dens) + lerror_bound_dens - LOG_TWO);
    }
    ld += dens;

    // computation of derivative for t and precision check in order to give
    // the value as deriv_y to edge1 and as -deriv_y to edge5
    double deriv_y = internal::dtdwiener5(y_val - t0_val, a_val, v_val, w_val,
                                          sv_val, labstol_wiener5);
    if (labstol_wiener5 > log(fabs(deriv_y)) + dens + lerror_bound - LOG_TWO) {
      deriv_y = internal::dtdwiener5(
          y_val - t0_val, a_val, v_val, w_val, sv_val,
          log(fabs(deriv_y)) + dens + lerror_bound - LOG_FOUR);
    }

    // computation of derivatives and precision checks
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[i] = deriv_y;
    }
    if (!is_constant_all<T_a>::value) {
      double deriv_a = internal::dadwiener5(y_val - t0_val, a_val, v_val, w_val,
                                            sv_val, labstol_wiener5, 0);
      if (labstol_wiener5
          > log(fabs(deriv_a)) + dens + lerror_bound - LOG_TWO) {
        deriv_a = internal::dadwiener5(
            y_val - t0_val, a_val, v_val, w_val, sv_val,
            log(fabs(deriv_a)) + dens + lerror_bound - LOG_FOUR, 0);
      }
      ops_partials.edge2_.partials_[i] = deriv_a;
    }
    if (!is_constant_all<T_t0>::value) {
      ops_partials.edge3_.partials_[i] = -deriv_y;
    }
    if (!is_constant_all<T_w>::value) {
      double deriv_w = internal::dwdwiener5(y_val - t0_val, a_val, v_val, w_val,
                                            sv_val, labstol_wiener5, 0);
      if (labstol_wiener5
          > log(fabs(deriv_w)) + dens + lerror_bound - LOG_TWO) {
        deriv_w = internal::dwdwiener5(
            y_val - t0_val, a_val, v_val, w_val, sv_val,
            log(fabs(deriv_w)) + dens + lerror_bound - LOG_FOUR, 0);
      }
      ops_partials.edge4_.partials_[i] = deriv_w;
    }
    if (!is_constant_all<T_v>::value) {
      ops_partials.edge5_.partials_[i]
          = internal::dvdwiener5(y_val - t0_val, a_val, v_val, w_val, sv_val);
    }
    if (!is_constant_all<T_sv>::value) {
      ops_partials.edge6_.partials_[i]
          = internal::dsvdwiener5(y_val - t0_val, a_val, v_val, w_val, sv_val);
    }
  }  // end for loop
  return ops_partials.build(ld);
}  // end wiener5_lpdf
}  // namespace math
}  // namespace stan
#endif
