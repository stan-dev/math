#ifndef STAN_MATH_PRIM_FUN_WIENER7_LPDF_HPP
#define STAN_MATH_PRIM_FUN_WIENER7_LPDF_HPP

#include <stan/math/prim/fun/ceil.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/wiener5_lpdf.hpp>

namespace stan {
namespace math {
namespace internal {

// tools
template <typename T_y, typename T_a, typename T_v, typename T_w, typename T_t0,
          typename T_sv, typename T_sw, typename T_st0, typename T_err>
struct my_params {
  T_y y;
  T_a a;
  T_v v;
  T_w w;
  T_t0 t0;
  T_sv sv;
  T_sw sw;
  T_st0 st0;
  T_err lerr;
};

template <typename T_y, typename T_a, typename T_v, typename T_w,
          typename T_w_lower, typename T_w_upper, typename T_t0, typename T_sv,
          typename T_sw, typename T_sw_mean, typename T_st0, typename T_err>
struct my_params2 {
  T_y y;
  T_a a;
  T_v v;
  T_w w;
  T_w_lower w_lower;
  T_w_upper w_upper;
  T_t0 t0;
  T_sv sv;
  T_sw sw;
  T_sw_mean sw_mean;
  T_st0 st0;
  T_err lerr;
};

template <typename T_y, typename T_a, typename T_v, typename T_w, typename T_t0,
          typename T_t0_mean, typename T_sv, typename T_sw, typename T_st0,
          typename T_st0_mean, typename T_err>
struct my_params3 {
  T_y y;
  T_a a;
  T_v v;
  T_w w;
  T_t0 t0;
  T_t0_mean t0_mean;
  T_sv sv;
  T_sw sw;
  T_st0 st0;
  T_st0_mean st0_mean;
  T_err lerr;
};

// calculate derivative of density of wiener5 with respect to y (version for
// wiener7)
template <typename T_y, typename T_a, typename T_v, typename T_w, typename T_sv>
return_type_t<T_y, T_a, T_v, T_w, T_sv> dtdwiener5_for_7(
    const T_y& y, const T_a& a, const T_v& v, const T_w& w, const T_sv& sv,
    const double& err) {
  using T_return_type = return_type_t<T_y, T_a, T_v, T_w, T_sv>;

  T_return_type kll, kss, ans;

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
    lg1 = (sv_sqr * square(a * w) - 2 * a * v * w - square(v) * y) / 2.0
              / one_plus_svsqr_y
          - la - 0.5 * log(one_plus_svsqr_y);
  } else {
    ans0 = -0.5 * square(v);
    lg1 = (-2 * a * v * w - square(v) * y) / 2.0 - la;
  }
  T_return_type factor = lg1 - la;

  T_return_type ld = dwiener5(y, a, -v, 1 - w, sv,
                              err - log(max(fabs(ans0), fabs(ans0 - 1.5 / y))));

  // calculate the number of terms kss needed for small y
  T_return_type es = err - lg1;
  es = es + la;
  T_return_type K1 = (sqrt(3.0 * y_asq) + w) / 2.0;
  T_return_type u_eps = fmin(
      -1.0, (log(8.0 / 27.0) + LOG_PI + 4.0 * log(y_asq) + 2.0 * es) / 3.0);
  T_return_type arg = -3.0 * y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  T_return_type K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate number of terms kll needed for large y
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

  // if small y is better
  if (2 * kss < kll) {
    // calculate terms of the sum for small y
    T_return_type twot = 2.0 * y_asq;
    if (static_cast<size_t>(kss) > 0) {
      for (size_t k = static_cast<size_t>(kss); k >= 1; k--) {
        T_return_type w_plus_2k = w + 2.0 * k;
        T_return_type w_minus_2k = w - 2.0 * k;
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
    T_return_type halfy = y_asq / 2.0;
    for (size_t k = static_cast<size_t>(kll); k >= 1; k--) {
      T_return_type pi_k = pi() * k;
      T_return_type zwi = sin(pi_k * w);
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
  using T_return_type = return_type_t<T_x, T_p>;
  my_params<T_return_type, T_return_type, T_return_type, T_return_type,
            T_return_type, T_return_type, T_return_type, T_return_type,
            T_return_type>* params
      = static_cast<my_params<T_return_type, T_return_type, T_return_type,
                              T_return_type, T_return_type, T_return_type,
                              T_return_type, T_return_type, T_return_type>*>(p);
  T_return_type y = (params->y);
  T_return_type a = (params->a);
  T_return_type v = (params->v);
  T_return_type t0 = (params->t0);
  T_return_type w = (params->w);
  T_return_type sw = (params->sw);
  T_return_type sv = (params->sv);
  T_return_type st0 = (params->st0);
  T_return_type lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);

  T_return_type retval;

  T_return_type omega = sw ? w + sw * (x_vec[0] - 0.5) : w;
  T_return_type t0_ = sw ? (st0 ? t0 + st0 * x_vec[1] : t0)
                         : (st0 ? t0 + st0 * x_vec[0] : t0);
  if (y - t0_ <= 0) {
    retval = 0.0;
  } else {
    T_return_type ldW = dwiener5(y - t0_, a, v, omega, sv, lerr);
    retval = exp(ldW);
  }
  return retval;
}
//-----------------------------------------------

// integrand d/dt (on normal scale)
template <typename T_x, typename T_p>
return_type_t<T_x, T_p> int_dtddiff(const T_x& x, const T_p& p) {
  using T_return_type = return_type_t<T_x, T_p>;
  my_params<T_return_type, T_return_type, T_return_type, T_return_type,
            T_return_type, T_return_type, T_return_type, T_return_type,
            T_return_type>* params
      = static_cast<my_params<T_return_type, T_return_type, T_return_type,
                              T_return_type, T_return_type, T_return_type,
                              T_return_type, T_return_type, T_return_type>*>(p);
  T_return_type y = (params->y);
  T_return_type a = (params->a);
  T_return_type v = (params->v);
  T_return_type t0 = (params->t0);
  T_return_type w = (params->w);
  T_return_type sw = (params->sw);
  T_return_type sv = (params->sv);
  T_return_type st0 = (params->st0);
  T_return_type lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);

  T_return_type retval;

  T_return_type omega = sw ? w + sw * (x_vec[0] - 0.5) : w;
  T_return_type t0_ = sw ? (st0 ? t0 + st0 * x_vec[1] : t0)
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
  using T_return_type = return_type_t<T_x, T_p>;
  my_params<T_return_type, T_return_type, T_return_type, T_return_type,
            T_return_type, T_return_type, T_return_type, T_return_type,
            T_return_type>* params
      = static_cast<my_params<T_return_type, T_return_type, T_return_type,
                              T_return_type, T_return_type, T_return_type,
                              T_return_type, T_return_type, T_return_type>*>(p);
  T_return_type y = (params->y);
  T_return_type a = (params->a);
  T_return_type v = (params->v);
  T_return_type t0 = (params->t0);
  T_return_type w = (params->w);
  T_return_type sw = (params->sw);
  T_return_type sv = (params->sv);
  T_return_type st0 = (params->st0);
  T_return_type lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);

  T_return_type retval;

  T_return_type omega = sw ? w + sw * (x_vec[0] - 0.5) : w;
  T_return_type t0_ = sw ? (st0 ? t0 + st0 * x_vec[1] : t0)
                         : (st0 ? t0 + st0 * x_vec[0] : t0);

  if (y - t0_ <= 0) {
    retval = 0.0;
  } else {
    // prepare some variables for error estimation in density
    T_return_type sv_sqr = square(sv);
    T_return_type ans0 = (v * (1 - omega) + sv_sqr * square(1 - omega) * a)
                         / (1 + sv_sqr * (y - t0_));
    T_return_type factor
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
  using T_return_type = return_type_t<T_x, T_p>;
  my_params<T_return_type, T_return_type, T_return_type, T_return_type,
            T_return_type, T_return_type, T_return_type, T_return_type,
            T_return_type>* params
      = static_cast<my_params<T_return_type, T_return_type, T_return_type,
                              T_return_type, T_return_type, T_return_type,
                              T_return_type, T_return_type, T_return_type>*>(p);
  T_return_type y = (params->y);
  T_return_type a = (params->a);
  T_return_type v = (params->v);
  T_return_type t0 = (params->t0);
  T_return_type w = (params->w);
  T_return_type sw = (params->sw);
  T_return_type sv = (params->sv);
  T_return_type st0 = (params->st0);
  T_return_type lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);

  T_return_type retval;

  T_return_type omega = sw ? w + sw * (x_vec[0] - 0.5) : w;
  T_return_type t0_ = sw ? (st0 ? t0 + st0 * x_vec[1] : t0)
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
  using T_return_type = return_type_t<T_x, T_p>;
  my_params<T_return_type, T_return_type, T_return_type, T_return_type,
            T_return_type, T_return_type, T_return_type, T_return_type,
            T_return_type>* params
      = static_cast<my_params<T_return_type, T_return_type, T_return_type,
                              T_return_type, T_return_type, T_return_type,
                              T_return_type, T_return_type, T_return_type>*>(p);
  T_return_type y = (params->y);
  T_return_type a = (params->a);
  T_return_type v = (params->v);
  T_return_type t0 = (params->t0);
  T_return_type w = (params->w);
  T_return_type sw = (params->sw);
  T_return_type sv = (params->sv);
  T_return_type st0 = (params->st0);
  T_return_type lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);

  T_return_type retval;

  T_return_type omega = sw ? w + sw * (x_vec[0] - 0.5) : w;
  T_return_type t0_ = sw ? (st0 ? t0 + st0 * x_vec[1] : t0)
                         : (st0 ? t0 + st0 * x_vec[0] : t0);

  if (y - t0_ <= 0) {
    retval = 0.0;
  } else {
    // prepare some variables for error estimation in density
    T_return_type sv_sqr = square(sv);
    T_return_type ans0
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
  using T_return_type = return_type_t<T_x, T_p>;
  my_params<T_return_type, T_return_type, T_return_type, T_return_type,
            T_return_type, T_return_type, T_return_type, T_return_type,
            T_return_type>* params
      = static_cast<my_params<T_return_type, T_return_type, T_return_type,
                              T_return_type, T_return_type, T_return_type,
                              T_return_type, T_return_type, T_return_type>*>(p);
  T_return_type y = (params->y);
  T_return_type a = (params->a);
  T_return_type v = (params->v);
  T_return_type t0 = (params->t0);
  T_return_type w = (params->w);
  T_return_type sw = (params->sw);
  T_return_type sv = (params->sv);
  T_return_type st0 = (params->st0);
  T_return_type lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);

  T_return_type retval;

  T_return_type omega = sw ? w + sw * (x_vec[0] - 0.5) : w;
  T_return_type t0_ = sw ? (st0 ? t0 + st0 * x_vec[1] : t0)
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
  using T_return_type = return_type_t<T_x, T_p>;
  my_params2<T_return_type, T_return_type, T_return_type, T_return_type,
             T_return_type, T_return_type, T_return_type, T_return_type,
             T_return_type, T_return_type, T_return_type, T_return_type>* params
      = static_cast<my_params2<T_return_type, T_return_type, T_return_type,
                               T_return_type, T_return_type, T_return_type,
                               T_return_type, T_return_type, T_return_type,
                               T_return_type, T_return_type, T_return_type>*>(
          p);
  T_return_type y = (params->y);
  T_return_type a = (params->a);
  T_return_type v = (params->v);
  T_return_type t0 = (params->t0);
  T_return_type w = (params->w);
  T_return_type w_lower = (params->w_lower);
  T_return_type w_upper = (params->w_upper);
  T_return_type sw = (params->sw);
  T_return_type sw_mean = (params->sw_mean);
  T_return_type sv = (params->sv);
  T_return_type st0 = (params->st0);
  T_return_type lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);

  T_return_type t0_ = sw ? (st0 ? t0 + st0 * x_vec[1] : t0)
                         : (st0 ? t0 + st0 * x_vec[0] : t0);

  T_return_type f, fl, fu;
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
  using T_return_type = return_type_t<T_x, T_p>;
  my_params3<T_return_type, T_return_type, T_return_type, T_return_type,
             T_return_type, T_return_type, T_return_type, T_return_type,
             T_return_type, T_return_type, T_return_type>* params
      = static_cast<
          my_params3<T_return_type, T_return_type, T_return_type, T_return_type,
                     T_return_type, T_return_type, T_return_type, T_return_type,
                     T_return_type, T_return_type, T_return_type>*>(p);
  T_return_type y = (params->y);
  T_return_type a = (params->a);
  T_return_type v = (params->v);
  T_return_type t0 = (params->t0);
  T_return_type t0_mean = (params->t0_mean);
  T_return_type w = (params->w);
  T_return_type sw = (params->sw);
  T_return_type sv = (params->sv);
  T_return_type st0 = (params->st0);
  T_return_type st0_mean = (params->st0_mean);
  T_return_type lerr = (params->lerr);

  using T_x_ref = ref_type_t<T_x>;
  T_x_ref x_ref = x;
  scalar_seq_view<T_x_ref> x_vec(x_ref);

  T_return_type t0_ = sw ? (st0_mean ? t0_mean + st0_mean * x_vec[1] : t0_mean)
                         : (st0_mean ? t0_mean + st0_mean * x_vec[0] : t0_mean);

  T_return_type omega = sw ? w + sw * (x_vec[0] - 0.5) : w;

  T_return_type f;
  if (y - t0_ <= 0) {
    return 0.0;
  } else {
    f = exp(internal::dwiener5(y - t0_, a, v, omega, sv, lerr));
  }
  return 1 / st0 * f;
}
//-----------------------------------------------

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
