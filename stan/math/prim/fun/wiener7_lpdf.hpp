// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef STAN_MATH_PRIM_FUN_WIENER7_LPDF_HPP
#define STAN_MATH_PRIM_FUN_WIENER7_LPDF_HPP

#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/wiener5_lpdf.hpp>

namespace stan {
namespace math {

/*
 * Helper functions for the log of the first passage time density function for a
 * (Wiener) drift diffusion model with 7 parameters. -boundary separation
 * (alpha), -drift rate (delta), -relative starting point (beta), -non-decision
 * time (tau), -inter-trial variability in drift rate (sv). -inter-trial
 * variability in relative starting point (sw). -inter-trial variability in
 * non-decition time (st0). For details see wiener_full_lpdf.
 */

namespace internal {

// tools
template <typename T_y, typename T_alpha, typename T_delta, typename T_beta,
          typename T_t0, typename T_sv, typename T_sw, typename T_st0,
          typename T_err>
struct my_params {
  T_y y;
  T_alpha a;
  T_delta v;
  T_beta w;
  T_t0 t0;
  T_sv sv;
  T_sw sw;
  T_st0 st0;
  T_err lerr;
};

template <typename T_y, typename T_alpha, typename T_delta, typename T_beta,
          typename T_beta_lower, typename T_beta_upper, typename T_t0,
          typename T_sv, typename T_sw, typename T_sw_mean, typename T_st0,
          typename T_err>
struct my_params2 {
  T_y y;
  T_alpha a;
  T_delta v;
  T_beta w;
  T_beta_lower w_lower;
  T_beta_upper w_upper;
  T_t0 t0;
  T_sv sv;
  T_sw sw;
  T_sw_mean sw_mean;
  T_st0 st0;
  T_err lerr;
};

template <typename T_y, typename T_alpha, typename T_delta, typename T_beta,
          typename T_t0, typename T_t0_mean, typename T_sv, typename T_sw,
          typename T_st0, typename T_st0_mean, typename T_err>
struct my_params3 {
  T_y y;
  T_alpha a;
  T_delta v;
  T_beta w;
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
template <typename T_y, typename T_alpha, typename T_delta, typename T_beta,
          typename T_sv>
return_type_t<T_y, T_alpha, T_delta, T_beta, T_sv> dtdwiener5_for_7(
    const T_y& q, const T_alpha& a, const T_delta& v, const T_beta& w,
    const T_sv& sv, const double& err) {
  using T_return_type = return_type_t<T_y, T_alpha, T_delta, T_beta, T_sv>;
  using std::ceil;
  using std::exp;
  using std::log;
  using std::pow;
  using std::sqrt;

  static const double PISQ = square(pi());  // pi*pi

  T_return_type kll, kss, ans;

  // prepare some variables
  T_return_type q_asq = q / pow(a, 2);
  T_return_type la = 2.0 * log(a);
  T_return_type ans0, lg1;
  if (sv != 0) {
    T_return_type eta_sqr = pow(sv, 2);
    T_return_type temp = (1 + eta_sqr * q);
    ans0 = -0.5
           * (pow(eta_sqr, 2) * (q + pow(a * w, 2))
              + eta_sqr * (1 - 2 * a * v * w) + pow(v, 2))
           / pow(temp, 2);
    lg1 = (eta_sqr * pow(a * w, 2) - 2 * a * v * w - pow(v, 2) * q) / 2.0 / temp
          - la - 0.5 * log(temp);
  } else {
    ans0 = -0.5 * pow(v, 2);
    lg1 = (-2 * a * v * w - pow(v, 2) * q) / 2.0 - la;
  }
  T_return_type factor = lg1 - la;

  T_return_type ld = dwiener5(q, a, -v, 1 - w, sv,
                              err - log(max(fabs(ans0), fabs(ans0 - 1.5 / q))));

  // calculate the number of terms kss needed for small y
  T_return_type es = err - lg1;
  es = es + la;
  T_return_type K1 = (sqrt(3.0 * q_asq) + w) / 2.0;
  T_return_type u_eps = fmin(
      -1.0, (log(8.0 / 27.0) + LOG_PI + 4.0 * log(q_asq) + 2.0 * es) / 3.0);
  T_return_type arg = -3.0 * q_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  T_return_type K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate number of terms kll needed for large y
  T_return_type el = err - lg1;
  el = el + la;
  K1 = sqrt(3.0 / q_asq) / pi();
  u_eps = fmin(-1.0, el + log(0.6) + LOG_PI + 2.0 * log(q_asq));
  arg = -2.0 / PISQ / q_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  T_return_type kl = (arg > 0) ? sqrt(arg) : K1;
  kll = ceil(fmax(kl, K1));

  T_return_type erg;
  T_return_type newsign = 1;
  T_return_type fplus = NEGATIVE_INFTY;
  T_return_type fminus = NEGATIVE_INFTY;

  // if small y is better
  if (2 * kss < kll) {
    // calculate terms of the sum for small y
    T_return_type twot = 2.0 * q_asq;
    if (static_cast<size_t>(kss) > 0) {
      for (size_t k = static_cast<size_t>(kss); k >= 1; k--) {
        T_return_type temp1 = w + 2.0 * k;
        T_return_type temp2 = w - 2.0 * k;
        fplus = log_sum_exp(3.0 * log(temp1) - temp1 * temp1 / twot, fplus);
        fminus = log_sum_exp(3.0 * log(-temp2) - temp2 * temp2 / twot, fminus);
      }
    }
    fplus = log_sum_exp(3.0 * log(w) - w * w / twot, fplus);
    if (fplus < fminus) {
      newsign = -1;
      erg = log_diff_exp(fminus, fplus);
    } else {
      erg = log_diff_exp(fplus, fminus);
    }

    ans = ans0 - 1.5 / q
          + newsign
                * exp(factor - 1.5 * LOG_TWO - LOG_SQRT_PI - 3.5 * log(q_asq)
                      + erg - ld);
    // if large y is better
  } else {
    // calculate terms of the sum for large y
    T_return_type halfq = q_asq / 2.0;
    for (size_t k = static_cast<size_t>(kll); k >= 1; k--) {
      T_return_type temp = pi() * k;
      T_return_type zwi = sin(temp * w);
      if (zwi > 0) {
        fplus = log_sum_exp(3.0 * log(static_cast<T_return_type>(k))
                                - temp * temp * halfq + log(zwi),
                            fplus);
      }
      if (zwi < 0) {
        fminus = log_sum_exp(3.0 * log(static_cast<T_return_type>(k))
                                 - temp * temp * halfq + log(-zwi),
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
    T_return_type eta_sqr = pow(sv, 2);
    T_return_type temp = (1 + eta_sqr * (y - t0_));
    T_return_type ans0
        = (v * (1 - omega) + eta_sqr * pow((1 - omega), 2) * a) / temp;
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
    T_return_type eta_sqr = pow(sv, 2);
    T_return_type temp = (1 + eta_sqr * (y - t0_));
    T_return_type ans0 = (v * a + eta_sqr * pow(a, 2) * (1 - omega)) / temp;
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
