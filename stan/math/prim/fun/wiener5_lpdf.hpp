// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef STAN_MATH_PRIM_FUN_WIENER5_LPDF_HPP
#define STAN_MATH_PRIM_FUN_WIENER5_LPDF_HPP

#include <stan/math/prim/fun/pow.hpp>
#include <stan/math/prim/fun/square.hpp>

namespace stan {
namespace math {

/*
 * The log of the first passage time density function for a (Wiener)
 * drift diffusion model with 5 parameters.
 * -boundary separation (alpha),
 * -drift rate (delta),
 * -relative starting point (beta),
 * -non-decision time (tau),
 * -inter-trial variability in drift rate (sv).
 * For details see wiener_full_lpdf.
 */

//-----------------------------------------------
// DENSITY and derivatives
//-----------------------------------------------

namespace internal {

// calculate density in log
template <typename T_y, typename T_alpha, typename T_delta, typename T_beta,
          typename T_sv>
return_type_t<T_y, T_alpha, T_delta, T_beta, T_sv> dwiener5(
    const T_y& q, const T_alpha& a, const T_delta& vn, const T_beta& wn,
    const T_sv& sv, const double& err) {
  using T_return_type = return_type_t<T_y, T_alpha, T_delta, T_beta, T_sv>;
  using std::ceil;
  using std::log;
  using std::pow;
  using std::sqrt;

  static const double PISQ = square(pi());  // pi*pi

  T_return_type kll, kss, ans, v, w;

  w = 1.0 - wn;
  v = -vn;

  T_return_type q_asq = q / pow(a, 2);
  ans = 0.0;

  // calculate the number of terms needed for short t
  T_return_type lg1;
  if (sv != 0) {
    T_return_type eta_sqr = pow(sv, 2);
    T_return_type temp = 1 + eta_sqr * q;
    lg1 = (eta_sqr * pow(a * w, 2) - 2 * a * v * w - pow(v, 2) * q) / 2.0 / temp
          - 2 * log(a) - 0.5 * log(temp);
  } else {
    lg1 = (-2 * a * v * w - pow(v, 2) * q) / 2.0 - 2 * log(a);
  }
  T_return_type es = (err - lg1);
  T_return_type K1 = (sqrt(2.0 * q_asq) + w) / 2.0;
  T_return_type u_eps
      = fmin(-1.0, LOG_TWO + LOG_PI + 2.0 * log(q_asq) + 2.0 * (es));
  T_return_type arg = -q_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  T_return_type K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate the number of terms needed for large t
  T_return_type el = es;
  K1 = 1.0 / (pi() * sqrt(q_asq));
  K2 = 0.0;
  T_return_type temp_ = -2.0 * (log(pi() * q_asq) + el);
  if (temp_ >= 0)
    K2 = sqrt(temp_ / (pow(pi(), 2) * q_asq));
  kll = ceil(fmax(K1, K2));

  // if small t is better
  if (2 * kss <= kll) {
    T_return_type fplus = NEGATIVE_INFTY;
    T_return_type fminus = NEGATIVE_INFTY;
    T_return_type twot = 2.0 * q_asq;
    if (static_cast<size_t>(kss) > 0) {
      for (size_t k = static_cast<size_t>(kss); k >= 1; k--) {
        T_return_type temp1 = w + 2.0 * k;
        T_return_type temp2 = w - 2.0 * k;

        fplus = log_sum_exp(log(temp1) - pow(temp1, 2) / twot, fplus);
        fminus = log_sum_exp(log(-temp2) - pow(temp2, 2) / twot, fminus);
      }
    }
    fplus = log_sum_exp(log(w) - pow(w, 2) / twot, fplus);
    ans = lg1
          + (-0.5 * LOG_TWO - LOG_SQRT_PI - 1.5 * log(q_asq)
             + log_diff_exp(fplus, fminus));
  // if large t is better
  } else {
    T_return_type fplus = NEGATIVE_INFTY;
    T_return_type fminus = NEGATIVE_INFTY;
    T_return_type halfq = q_asq / 2.0;
    for (size_t k = static_cast<size_t>(kll); k >= 1; k--) {
      T_return_type temp = k * pi();
      T_return_type check = sin(temp * w);
      if (check > 0)
        fplus = log_sum_exp(log(static_cast<T_return_type>(k))
                                - pow(temp, 2) * halfq + log(check),
                            fplus);
      else
        fminus = log_sum_exp(log(static_cast<T_return_type>(k))
                                 - pow(temp, 2) * halfq + log(-check),
                             fminus);
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
template <typename T_y, typename T_alpha, typename T_delta, typename T_beta,
          typename T_sv>
return_type_t<T_y, T_alpha, T_delta, T_sv, T_beta> dtdwiener5(
    const T_y& q, const T_alpha& a, const T_delta& vn, const T_beta& wn,
    const T_sv& sv, const double& err) {
  using T_return_type = return_type_t<T_y, T_alpha, T_delta, T_beta, T_sv>;
  using std::ceil;
  using std::exp;
  using std::log;
  using std::pow;
  using std::sqrt;

  static const double PISQ = square(pi());  // pi*pi

  T_return_type kll, kss, ans, v, w;

  w = 1.0 - wn;
  v = -vn;

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

  T_return_type ld = dwiener5(q, a, vn, wn, sv,
                              err - log(max(fabs(ans0 - 1.5 / q), fabs(ans0))));

  // calculate the number of terms kss needed for small t
  T_return_type es = err - lg1;
  es = es + la;
  T_return_type K1 = (sqrt(3.0 * q_asq) + w) / 2.0;
  T_return_type u_eps = fmin(
      -1.0, (log(8.0 / 27.0) + LOG_PI + 4.0 * log(q_asq) + 2.0 * es) / 3.0);
  T_return_type arg = -3.0 * q_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  T_return_type K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate number of terms kll needed for large t
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

  // if small t is better
  if (2 * kss < kll) {
    // calculate terms of the sum for small t
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
  // if large t is better
  } else {
    // calculate terms of the sum for large t
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
  return ans;
}
//-----------------------------------------------

// d/da DENSITY
// calculate derivative of density with respect to a (in log, ans =
// d/da(log(f))=d/da f'/f; ans*exp(ld)=f' on normal scale)
template <typename T_y, typename T_alpha, typename T_delta, typename T_beta,
          typename T_sv>
return_type_t<T_y, T_alpha, T_delta, T_beta, T_sv> dadwiener5(
    const T_y& q, const T_alpha& a, const T_delta& vn, const T_beta& wn,
    const T_sv& sv, const double& err, const int& normal_or_log) {
  using T_return_type = return_type_t<T_y, T_alpha, T_delta, T_beta>;
  using std::ceil;
  using std::exp;
  using std::log;
  using std::pow;
  using std::sqrt;

  static const double PISQ = square(pi());  // pi*pi

  T_return_type kll, kss, ans, v, w;

  T_return_type la = log(a);
  T_return_type lq = log(q);

  w = 1.0 - wn;
  v = -vn;

  // prepare some variables
  T_return_type q_asq = q / pow(a, 2);
  T_return_type ans0, lg1;
  if (sv != 0) {
    T_return_type eta_sqr = pow(sv, 2);
    T_return_type temp = (1 + eta_sqr * q);
    ans0 = (-v * w + eta_sqr * pow(w, 2) * a) / temp;
    lg1 = (eta_sqr * pow(a * w, 2) - 2 * a * v * w - pow(v, 2) * q) / 2.0 / temp
          - 2 * la - 0.5 * log(temp);
  } else {
    ans0 = -v * w;
    lg1 = (-2 * a * v * w - pow(v, 2) * q) / 2.0 - 2 * la;
  }
  T_return_type factor = lg1 - 3 * la;

  T_return_type ld
      = dwiener5(q, a, vn, wn, sv,
                 err - log(max(fabs(ans0 + 1.0 / a), fabs(ans0 - 2.0 / a))));

  // calculate the number of terms kss needed for small t
  T_return_type es = err - lg1;
  es = es + la;
  es = es - LOG_TWO + 2.0 * la - lq;
  T_return_type K1 = (sqrt(3.0 * q_asq) + w) / 2.0;
  T_return_type u_eps = fmin(
      -1.0, (log(8.0 / 27.0) + LOG_PI + 4.0 * log(q_asq) + 2.0 * es) / 3.0);
  T_return_type arg = -3.0 * q_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  T_return_type K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate number of terms kll needed for large t
  T_return_type el = err - lg1;
  el = el + la;
  el = el - LOG_TWO + 2.0 * la - lq;
  K1 = sqrt(3.0 / q_asq) / pi();
  u_eps = fmin(-1.0, el + log(0.6) + LOG_PI + 2.0 * log(q_asq));
  arg = -2.0 / PISQ / q_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  T_return_type kl = (arg > 0) ? sqrt(arg) : K1;
  kll = ceil(fmax(kl, K1));

  T_return_type erg;
  T_return_type newsign = 1;
  T_return_type fplus = NEGATIVE_INFTY;
  T_return_type fminus = NEGATIVE_INFTY;

  // if small t is better
  if (2 * kss < kll) {
    // calculate terms of the sum for short t
    T_return_type twot = 2.0 * q_asq;
    if (static_cast<int>(kss) > 0) {
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

    ans = ans0 + 1.0 / a
          - newsign
                * exp(-0.5 * LOG_TWO - LOG_SQRT_PI - 2.5 * lq + 4.0 * la + lg1
                      + erg - ld);
  // if large t is better
  } else {
    // calculate terms of the sum for large t
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
    if (fplus > fminus) {
      erg = log_diff_exp(fplus, fminus);
    } else {
      erg = log_diff_exp(fminus, fplus);
      newsign = -1;
    }

    ans = ans0 - 2.0 / a + newsign * exp(lq + factor + 3.0 * LOG_PI + erg - ld);
  }
  if (normal_or_log == 1)
    return ans * exp(ld);  // derivative of f for hcubature
  else
    return ans;  // derivative of log(f)
}
//-----------------------------------------------

// d/dv DENSITY
// calculate derivative of density with respect to v (in log, ans =
// d/dv(log(f))=d/dv f'/f; ans*exp(ld)=f' on normal scale)
template <typename T_y, typename T_alpha, typename T_delta, typename T_beta,
          typename T_sv>
return_type_t<T_y, T_alpha, T_delta, T_beta, T_sv> dvdwiener5(const T_y& q,
                                                              const T_alpha& a,
                                                              const T_delta& vn,
                                                              const T_beta& wn,
                                                              const T_sv& sv) {
  using T_return_type = return_type_t<T_y, T_alpha, T_delta, T_beta, T_sv>;
  using std::pow;

  T_return_type temp;
  if (sv != 0) {
    temp = 1 + pow(sv, 2) * q;
    temp = (a * (1 - wn) - vn * q) / temp;
  } else {
    temp = (a * (1 - wn) - vn * q);
  }
  return temp;
}
//-----------------------------------------------

// d/dw DENSITY
// calculate derivative of density with respect to w (in log, ans =
// d/dw(log(f))=d/dw f'/f; ans*exp(ld)=f' on normal scale)
template <typename T_y, typename T_alpha, typename T_delta, typename T_beta,
          typename T_sv>
return_type_t<T_y, T_alpha, T_delta, T_beta, T_sv> dwdwiener5(
    const T_y& q, const T_alpha& a, const T_delta& vn, const T_beta& wn,
    const T_sv& sv, const double& err, const int& normal_or_log) {
  using T_return_type = return_type_t<T_y, T_alpha, T_delta, T_beta, T_sv>;
  using std::ceil;
  using std::exp;
  using std::log;
  using std::pow;
  using std::sqrt;

  T_return_type kll, kss, ans, v, w;
  T_return_type sign = -1;

  w = 1.0 - wn;
  v = -vn;

  // prepare some variables
  T_return_type q_asq = q / pow(a, 2);
  T_return_type ans0, lg1;
  if (sv != 0) {
    T_return_type eta_sqr = pow(sv, 2);
    T_return_type temp = (1 + eta_sqr * q);
    ans0 = (-v * a + eta_sqr * pow(a, 2) * w) / temp;
    lg1 = (eta_sqr * pow(a * w, 2) - 2 * a * v * w - pow(v, 2) * q) / 2.0 / temp
          - 2 * log(a) - 0.5 * log(temp);
  } else {
    ans0 = -v * a;
    lg1 = (-2 * a * v * w - pow(v, 2) * q) / 2.0 - 2 * log(a);
  }
  T_return_type ld = dwiener5(q, a, vn, wn, sv, err - log(fabs(ans0)));

  T_return_type ls = -lg1 + ld;
  T_return_type ll = -lg1 + ld;

  // calculate the number of terms kss needed for small t
  T_return_type K1 = (sqrt(3.0 * q_asq) + w) / 2.0;
  T_return_type u_eps
      = fmin(-1.0, 2.0 * (err - lg1) + LOG_TWO + LOG_PI + 2.0 * log(q_asq));
  T_return_type arg = -q_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  T_return_type K2 = (arg > 0) ? 0.5 * (sqrt(arg) + w) : K1;
  kss = ceil(fmax(K1, K2));

  // calculate number of terms kll needed for large t
  K1 = sqrt(2.0 / q_asq) / pi();
  u_eps = fmin(-1.0, log(4.0 / 9.0) + 2.0 * LOG_PI + 3.0 * log(q_asq)
                         + 2.0 * (err - lg1));
  arg = -(u_eps - sqrt(-2.0 * u_eps - 2.0));
  K2 = (arg > 0) ? 1.0 / pi() * sqrt(arg / q_asq) : K1;
  kll = ceil(fmax(K1, K2));

  T_return_type erg;
  T_return_type newsign = 1;
  T_return_type fplus = NEGATIVE_INFTY;
  T_return_type fminus = NEGATIVE_INFTY;

  // if small t is better
  if (2 * kss < kll) {
    // calculate terms of the sum for short t
    T_return_type twot = 2.0 * q_asq;
    for (size_t k = static_cast<size_t>(kss); k >= 1; k--) {
      T_return_type temp1 = pow((w + 2 * k), 2);
      T_return_type temp2 = pow((w - 2 * k), 2);
      T_return_type temp3 = temp1 - q_asq;
      T_return_type temp4 = temp2 - q_asq;
      if (temp3 > 0)
        fplus = log_sum_exp(log(temp3) - temp1 / twot, fplus);
      else if (temp3 < 0)
        fminus = log_sum_exp(log(-(temp3)) - temp1 / twot, fminus);
      if (temp4 > 0)
        fplus = log_sum_exp(log(temp4) - temp2 / twot, fplus);
      else if (temp4 < 0)
        fminus = log_sum_exp(log(-(temp4)) - temp2 / twot, fminus);
    }
    T_return_type temp = pow(w, 2);
    T_return_type temp1 = temp - q_asq;
    if (temp1 > 0)
      fplus = log_sum_exp(log(temp1) - temp / twot, fplus);
    else if (temp1 < 0)
      fminus = log_sum_exp(log(-(temp1)) - temp / twot, fminus);

    if (fplus < fminus) {
      newsign = -1;
      erg = log_diff_exp(fminus, fplus);
    } else {
      erg = log_diff_exp(fplus, fminus);
    }

    ans = ans0
          - newsign
                * exp(erg - ls - 2.5 * log(q_asq) - 0.5 * LOG_TWO
                      - 0.5 * LOG_PI);
  // if large t is better
  } else {
    // calculate terms of the sum for large t
    T_return_type halfq = q_asq / 2.0;
    for (size_t k = static_cast<size_t>(kll); k >= 1; k--) {
      T_return_type temp = pi() * k;
      T_return_type x = cos(temp * w);
      if (x > 0)
        fplus = log_sum_exp(2.0 * log(static_cast<T_return_type>(k))
                                - pow(temp, 2) * halfq + log(x),
                            fplus);
      else if (x < 0)
        fminus = log_sum_exp(2.0 * log(static_cast<T_return_type>(k))
                                 - pow(temp, 2) * halfq + log(-x),
                             fminus);
    }
    if (fplus < fminus) {
      erg = log_diff_exp(fminus, fplus);
      newsign = -1;
    } else {
      erg = log_diff_exp(fplus, fminus);
    }

    ans = ans0 + newsign * exp(erg - ll + 2.0 * LOG_PI);
  }
  if (normal_or_log == 1)
    return ans * sign * exp(ld);  // derivative of f for hcubature
  else
    return ans * sign;  // derivative of log(f)
}
//-----------------------------------------------

// d/dsv DENSITY
// calculate derivative of density with respect to sv (in log, ans =
// d/dsv(log(f))=d/dsv f'/f; ans*exp(ld)=f' on normal scale)
template <typename T_y, typename T_alpha, typename T_delta, typename T_beta,
          typename T_sv>
return_type_t<T_y, T_alpha, T_delta, T_beta, T_sv> dsvdwiener5(
    const T_y& q, const T_alpha& a, const T_delta& vn, const T_beta& wn,
    const T_sv& sv) {
  using T_return_type = return_type_t<T_y, T_alpha, T_delta, T_beta, T_sv>;
  using std::pow;

  T_return_type v, w;

  v = -vn;
  w = 1 - wn;

  T_return_type temp = 1 + pow(sv, 2) * q;
  T_return_type t1 = -q / temp;
  T_return_type t2
      = (pow(a * w, 2) + 2 * a * v * w * q + pow(v * q, 2)) / pow(temp, 2);
  return sv * (t1 + t2);
}
//-----------------------------------------------

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
