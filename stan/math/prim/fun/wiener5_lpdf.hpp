#ifndef STAN_MATH_PRIM_FUN_WIENER5_LPDF_HPP
#define STAN_MATH_PRIM_FUN_WIENER5_LPDF_HPP

#include <stan/math/prim/fun/ceil.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>

namespace stan {
namespace math {
namespace internal {

// calculate density in log 
template <typename T_y, typename T_a, typename T_v, typename T_w,
          typename T_sv>
return_type_t<T_y, T_a, T_v, T_w, T_sv> dwiener5(
    const T_y& y, const T_a& a, const T_v& vn, const T_w& wn,
    const T_sv& sv, const double& err) {
  using T_return_type = return_type_t<T_y, T_a, T_v, T_w, T_sv>;

  T_return_type kll, kss, ans, v, w;

  w = 1.0 - wn;
  v = -vn;

  T_return_type y_asq = y / square(a);
  ans = 0.0;

  // calculate the number of terms needed for short t
  T_return_type lg1;
  if (sv != 0) {
    T_return_type sv_sqr = square(sv);
    T_return_type one_plus_svsqr_y = 1 + sv_sqr * y;
    lg1 = (sv_sqr * square(a * w) - 2 * a * v * w - square(v) * y) / 2.0 / one_plus_svsqr_y
          - 2 * log(a) - 0.5 * log(one_plus_svsqr_y);
  } else {
    lg1 = (-2 * a * v * w - square(v) * y) / 2.0 - 2 * log(a);
  }
  T_return_type es = (err - lg1);
  T_return_type K1 = (sqrt(2.0 * y_asq) + w) / 2.0;
  T_return_type u_eps
      = fmin(-1.0, LOG_TWO + LOG_PI + 2.0 * log(y_asq) + 2.0 * (es));
  T_return_type arg = -y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  T_return_type K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate the number of terms needed for large t
  T_return_type el = es;
  K1 = 1.0 / (pi() * sqrt(y_asq));
  K2 = 0.0;
  T_return_type two_log_piy = -2.0 * (log(pi() * y_asq) + el);
  static const double PISQ = square(pi());  // pi*pi
  if (two_log_piy >= 0) {
    K2 = sqrt(two_log_piy / (PISQ * y_asq));
  }
  kll = ceil(fmax(K1, K2));

  // if small t is better
  if (2 * kss <= kll) {
    T_return_type fplus = NEGATIVE_INFTY;
    T_return_type fminus = NEGATIVE_INFTY;
    T_return_type twoy = 2.0 * y_asq;
    if (static_cast<size_t>(kss) > 0) {
      for (size_t k = static_cast<size_t>(kss); k >= 1; k--) {
        T_return_type w_plus_2k = w + 2.0 * k;
        T_return_type w_minus_2k = w - 2.0 * k;

        fplus = log_sum_exp(log(w_plus_2k) - square(w_plus_2k) / twoy, fplus);
        fminus = log_sum_exp(log(-w_minus_2k) - square(w_minus_2k) / twoy, fminus);
      }
    }
    fplus = log_sum_exp(log(w) - square(w) / twoy, fplus);
    ans = lg1
          + (-0.5 * LOG_TWO - LOG_SQRT_PI - 1.5 * log(y_asq)
             + log_diff_exp(fplus, fminus));
    // if large t is better
  } else {
    T_return_type fplus = NEGATIVE_INFTY;
    T_return_type fminus = NEGATIVE_INFTY;
    T_return_type halfy = y_asq / 2.0;
    for (size_t k = static_cast<size_t>(kll); k >= 1; k--) {
      T_return_type pi_k = k * pi();
      T_return_type check = sin(pi_k * w);
      if (check > 0) {
        fplus = log_sum_exp(log(k)
                                - square(pi_k) * halfy + log(check),
                            fplus);
	  }
      else {
        fminus = log_sum_exp(log(k)
                                 - square(pi_k) * halfy + log(-check),
                             fminus);
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
template <typename T_y, typename T_a, typename T_v, typename T_w,
          typename T_sv>
return_type_t<T_y, T_a, T_v, T_sv, T_w> dtdwiener5(
    const T_y& y, const T_a& a, const T_v& vn, const T_w& wn,
    const T_sv& sv, const double& err) {
  using T_return_type = return_type_t<T_y, T_a, T_v, T_w, T_sv>;

  T_return_type kll, kss, ans, v, w;

  w = 1.0 - wn;
  v = -vn;

  // prepare some variables
  T_return_type y_asq = y / square(a);
  T_return_type la = 2.0 * log(a);
  T_return_type ans0, lg1;
  if (sv != 0) {
    T_return_type sv_sqr = square(sv);
    T_return_type one_plus_svsqr_y = (1 + sv_sqr * y);
    ans0 = -0.5
           * (square(sv_sqr) * (y + square(a * w))
              + sv_sqr * (1 - 2 * a * v * w) + square(v))
           / square(one_plus_svsqr_y);
    lg1 = (sv_sqr * square(a * w) - 2 * a * v * w - square(v) * y) / 2.0 / one_plus_svsqr_y
          - la - 0.5 * log(one_plus_svsqr_y);
  } else {
    ans0 = -0.5 * square(v);
    lg1 = (-2 * a * v * w - square(v) * y) / 2.0 - la;
  }
  T_return_type factor = lg1 - la;

  T_return_type ld = dwiener5(y, a, vn, wn, sv,
                              err - log(max(fabs(ans0 - 1.5 / y), fabs(ans0))));

  // calculate the number of terms kss needed for small t
  T_return_type es = err - lg1;
  es = es + la;
  T_return_type K1 = (sqrt(3.0 * y_asq) + w) / 2.0;
  T_return_type u_eps = fmin(
      -1.0, (log(8.0 / 27.0) + LOG_PI + 4.0 * log(y_asq) + 2.0 * es) / 3.0);
  T_return_type arg = -3.0 * y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  T_return_type K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate number of terms kll needed for large t
  T_return_type el = err - lg1;
  el = el + la;
  K1 = sqrt(3.0 / y_asq) / pi();
  u_eps = fmin(-1.0, el + log(0.6) + LOG_PI + 2.0 * log(y_asq));
  static const double PISQ = square(pi());  // pi*pi
  arg = -2.0 / PISQ / y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  T_return_type kl = (arg > 0) ? sqrt(arg) : K1;
  kll = ceil(fmax(kl, K1));

  T_return_type erg;
  T_return_type newsign = 1;
  T_return_type fplus = NEGATIVE_INFTY;
  T_return_type fminus = NEGATIVE_INFTY;

  // if small t is better
  if (2 * kss < kll) {
    // calculate terms of the sum for small t
    T_return_type twoy = 2.0 * y_asq;
    if (static_cast<size_t>(kss) > 0) {
      for (size_t k = static_cast<size_t>(kss); k >= 1; k--) {
        T_return_type w_plus_2k = w + 2.0 * k;
        T_return_type w_minus_2k = w - 2.0 * k;
        fplus = log_sum_exp(3.0 * log(w_plus_2k) - w_plus_2k * w_plus_2k / twoy, fplus);
        fminus = log_sum_exp(3.0 * log(-w_minus_2k) - w_minus_2k * w_minus_2k / twoy, fminus);
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
    T_return_type halfy = y_asq / 2.0;
    for (size_t k = static_cast<size_t>(kll); k >= 1; k--) {
      T_return_type pi_k = pi() * k;
      T_return_type zwi = sin(pi_k * w);
      if (zwi > 0) {
        fplus = log_sum_exp(3.0 * log(k)
                                - pi_k * pi_k * halfy + log(zwi),
                            fplus);
      }
      if (zwi < 0) {
        fminus = log_sum_exp(3.0 * log(k)
                                 - pi_k * pi_k * halfy + log(-zwi),
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
template <typename T_y, typename T_a, typename T_v, typename T_w,
          typename T_sv>
return_type_t<T_y, T_a, T_v, T_w, T_sv> dadwiener5(
    const T_y& y, const T_a& a, const T_v& vn, const T_w& wn,
    const T_sv& sv, const double& err, const int& normal_or_log) {
  using T_return_type = return_type_t<T_y, T_a, T_v, T_w>;

  T_return_type kll, kss, ans, v, w;

  T_return_type la = log(a);
  T_return_type ly = log(y);

  w = 1.0 - wn;
  v = -vn;

  // prepare some variables
  T_return_type y_asq = y / square(a);
  T_return_type ans0, lg1;
  if (sv != 0) {
    T_return_type sv_sqr = square(sv);
    T_return_type one_plus_svsqr_y = (1 + sv_sqr * y);
    ans0 = (-v * w + sv_sqr * square(w) * a) / one_plus_svsqr_y;
    lg1 = (sv_sqr * square(a * w) - 2 * a * v * w - square(v) * y) / 2.0 / one_plus_svsqr_y
          - 2 * la - 0.5 * log(one_plus_svsqr_y);
  } else {
    ans0 = -v * w;
    lg1 = (-2 * a * v * w - square(v) * y) / 2.0 - 2 * la;
  }
  T_return_type factor = lg1 - 3 * la;

  T_return_type ld
      = dwiener5(y, a, vn, wn, sv,
                 err - log(max(fabs(ans0 + 1.0 / a), fabs(ans0 - 2.0 / a))));

  // calculate the number of terms kss needed for small t
  T_return_type es = err - lg1;
  es = es + la;
  es = es - LOG_TWO + 2.0 * la - ly;
  T_return_type K1 = (sqrt(3.0 * y_asq) + w) / 2.0;
  T_return_type u_eps = fmin(
      -1.0, (log(8.0 / 27.0) + LOG_PI + 4.0 * log(y_asq) + 2.0 * es) / 3.0);
  T_return_type arg = -3.0 * y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  T_return_type K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate number of terms kll needed for large t
  T_return_type el = err - lg1;
  el = el + la;
  el = el - LOG_TWO + 2.0 * la - ly;
  K1 = sqrt(3.0 / y_asq) / pi();
  u_eps = fmin(-1.0, el + log(0.6) + LOG_PI + 2.0 * log(y_asq));
  static const double PISQ = square(pi());  // pi*pi
  arg = -2.0 / PISQ / y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  T_return_type kl = (arg > 0) ? sqrt(arg) : K1;
  kll = ceil(fmax(kl, K1));

  T_return_type erg;
  T_return_type newsign = 1;
  T_return_type fplus = NEGATIVE_INFTY;
  T_return_type fminus = NEGATIVE_INFTY;

  // if small t is better
  if (2 * kss < kll) {
    // calculate terms of the sum for short t
    T_return_type twoy = 2.0 * y_asq;
    if (static_cast<int>(kss) > 0) {
      for (size_t k = static_cast<size_t>(kss); k >= 1; k--) {
        T_return_type w_plus_2k = w + 2.0 * k;
        T_return_type w_minus_2k = w - 2.0 * k;
        fplus = log_sum_exp(3.0 * log(w_plus_2k) - w_plus_2k * w_plus_2k / twoy, fplus);
        fminus = log_sum_exp(3.0 * log(-w_minus_2k) - w_minus_2k * w_minus_2k / twoy, fminus);
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
    T_return_type halfy = y_asq / 2.0;
    for (size_t k = static_cast<size_t>(kll); k >= 1; k--) {
      T_return_type pi_k = pi() * k;
      T_return_type zwi = sin(pi_k * w);
      if (zwi > 0) {
        fplus = log_sum_exp(3.0 * log(k)
                                - pi_k * pi_k * halfy + log(zwi),
                            fplus);
      }
      if (zwi < 0) {
        fminus = log_sum_exp(3.0 * log(k)
                                 - pi_k * pi_k * halfy + log(-zwi),
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
  }
  else {
    return ans;  // derivative of log(f)
  }
}
//-----------------------------------------------

// d/dv DENSITY
// calculate derivative of density with respect to v (in log, ans =
// d/dv(log(f))=d/dv f'/f; ans*exp(ld)=f' on normal scale)
template <typename T_y, typename T_a, typename T_v, typename T_w,
          typename T_sv>
return_type_t<T_y, T_a, T_v, T_w, T_sv> dvdwiener5(const T_y& y,
                                                              const T_a& a,
                                                              const T_v& vn,
                                                              const T_w& wn,
                                                              const T_sv& sv) {
  using T_return_type = return_type_t<T_y, T_a, T_v, T_w, T_sv>;

  T_return_type ans;
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
template <typename T_y, typename T_a, typename T_v, typename T_w,
          typename T_sv>
return_type_t<T_y, T_a, T_v, T_w, T_sv> dwdwiener5(
    const T_y& y, const T_a& a, const T_v& vn, const T_w& wn,
    const T_sv& sv, const double& err, const int& normal_or_log) {
  using T_return_type = return_type_t<T_y, T_a, T_v, T_w, T_sv>;

  T_return_type kll, kss, ans, v, w;
  T_return_type sign = -1;

  w = 1.0 - wn;
  v = -vn;

  // prepare some variables
  T_return_type y_asq = y / square(a);
  T_return_type ans0, lg1;
  if (sv != 0) {
    T_return_type sv_sqr = square(sv);
    T_return_type one_plus_svsqr_y = (1 + sv_sqr * y);
    ans0 = (-v * a + sv_sqr * square(a) * w) / one_plus_svsqr_y;
    lg1 = (sv_sqr * square(a * w) - 2 * a * v * w - square(v) * y) / 2.0 / one_plus_svsqr_y
          - 2 * log(a) - 0.5 * log(one_plus_svsqr_y);
  } else {
    ans0 = -v * a;
    lg1 = (-2 * a * v * w - square(v) * y) / 2.0 - 2 * log(a);
  }
  T_return_type ld = dwiener5(y, a, vn, wn, sv, err - log(fabs(ans0)));

  T_return_type ls = -lg1 + ld;
  T_return_type ll = -lg1 + ld;

  // calculate the number of terms kss needed for small t
  T_return_type K1 = (sqrt(3.0 * y_asq) + w) / 2.0;
  T_return_type u_eps
      = fmin(-1.0, 2.0 * (err - lg1) + LOG_TWO + LOG_PI + 2.0 * log(y_asq));
  T_return_type arg = -y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  T_return_type K2 = (arg > 0) ? 0.5 * (sqrt(arg) + w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate number of terms kll needed for large t
  K1 = sqrt(2.0 / y_asq) / pi();
  static const double TWO_LOG_PI = 2.0 * LOG_PI;
  u_eps = fmin(
      -1.0, log(4.0 / 9.0) + TWO_LOG_PI + 3.0 * log(y_asq) + 2.0 * (err - lg1));
  arg = -(u_eps - sqrt(-2.0 * u_eps - 2.0));
  K2 = (arg > 0) ? 1.0 / pi() * sqrt(arg / y_asq) : K1;
  kll = ceil(fmax(K1, K2));

  T_return_type erg;
  T_return_type newsign = 1;
  T_return_type fplus = NEGATIVE_INFTY;
  T_return_type fminus = NEGATIVE_INFTY;

  // if small t is better
  if (2 * kss < kll) {
    // calculate terms of the sum for short t
    T_return_type twoy = 2.0 * y_asq;
    for (size_t k = static_cast<size_t>(kss); k >= 1; k--) {
      T_return_type sqrt_w_plus_2k = square(w + 2 * k);
      T_return_type sqrt_w_minus_2k = square(w - 2 * k);
      T_return_type wp2k_minusy = sqrt_w_plus_2k - y_asq;
      T_return_type wm2k_minusy = sqrt_w_minus_2k - y_asq;
      if (wp2k_minusy > 0) {
        fplus = log_sum_exp(log(wp2k_minusy) - sqrt_w_plus_2k / twoy, fplus);
	  }
      else if (wp2k_minusy < 0) {
        fminus = log_sum_exp(log(-(wp2k_minusy)) - sqrt_w_plus_2k / twoy, fminus);
	  }
      if (wm2k_minusy > 0) {
        fplus = log_sum_exp(log(wm2k_minusy) - sqrt_w_minus_2k / twoy, fplus);
	  }
      else if (wm2k_minusy < 0) {
        fminus = log_sum_exp(log(-(wm2k_minusy)) - sqrt_w_minus_2k / twoy, fminus);
	  }
    }
    T_return_type sqr_w = square(w);
    T_return_type sqrt_w_plus_2k = sqr_w - y_asq;
    if (sqrt_w_plus_2k > 0) {
      fplus = log_sum_exp(log(sqrt_w_plus_2k) - sqr_w / twoy, fplus);
	}
    else if (sqrt_w_plus_2k < 0) {
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
    T_return_type halfy = y_asq / 2.0;
    for (size_t k = static_cast<size_t>(kll); k >= 1; k--) {
      T_return_type pi_k = pi() * k;
      T_return_type x = cos(pi_k * w);
      if (x > 0) {
        fplus = log_sum_exp(2.0 * log(k)
                                - square(pi_k) * halfy + log(x),
                            fplus);
	  }
      else if (x < 0) {
        fminus = log_sum_exp(2.0 * log(k)
                                 - square(pi_k) * halfy + log(-x),
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
  }
  else {
    return ans * sign;  // derivative of log(f)
  }
}
//-----------------------------------------------

// d/dsv DENSITY
// calculate derivative of density with respect to sv (in log, ans =
// d/dsv(log(f))=d/dsv f'/f; ans*exp(ld)=f' on normal scale)
template <typename T_y, typename T_a, typename T_v, typename T_w,
          typename T_sv>
return_type_t<T_y, T_a, T_v, T_w, T_sv> dsvdwiener5(
    const T_y& y, const T_a& a, const T_v& vn, const T_w& wn,
    const T_sv& sv) {
  using T_return_type = return_type_t<T_y, T_a, T_v, T_w, T_sv>;

  T_return_type v, w;

  v = -vn;
  w = 1 - wn;

  T_return_type one_sqrsv_y = 1 + square(sv) * y;
  T_return_type t1 = -y / one_sqrsv_y;
  T_return_type t2
      = (square(a * w) + 2 * a * v * w * y + square(v * y)) / square(one_sqrsv_y);
  return sv * (t1 + t2);
}
//-----------------------------------------------

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
