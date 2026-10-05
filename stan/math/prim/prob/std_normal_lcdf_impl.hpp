#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_LCDF_IMPL_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_LCDF_IMPL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/pow.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>
#include <limits>
#include <utility>

namespace stan {
namespace math {
namespace internal {

/**
 * The piecewise value and derivative of the standard normal log CDF. This is
 * the body of `std_normal_lcdf`, moved here so that other distributions can
 * share it; the inverse Gaussian cdf family is the first.
 *
 * Both are written in terms of the scaled variable `scaled`, related to the
 * argument of \f$\log\Phi\f$ by \f$\texttt{scaled} = z/\sqrt{2}\f$.
 * `std_normal_lcdf` passes `y * INV_SQRT_TWO`. `std_normal_lcdf_grad` returns
 * the derivative with respect to `scaled`, so callers apply the chain rule for
 * their own parameterization.
 *
 * The piecewise structure is the same as `normal_lcdf`, and the cutoffs,
 * their provenance and the measurements behind them are documented once, on
 * that function. See prim/prob/normal_lcdf.hpp.
 *
 * Callers reject NaN before calling. A NaN `scaled` fails every branch test
 * in `std_normal_lcdf_value` and returns `-inf`.
 */

/**
 * Return the Cody (1969) lower-tail approximation to
 * \f$\log\Phi(\sqrt{2}\,\texttt{scaled})\f$ with the Gaussian exponent
 * \f$-\texttt{scaled}^2\f$ left off.
 *
 * The full branch value is `std_normal_lcdf_cody_correction(scaled) -
 * square(scaled)`. The two are exposed separately because a caller that must
 * add a large positive quantity to the log CDF -- the inverse Gaussian CDF
 * family adds \f$2\lambda/\mu\f$ to \f$\log\Phi(-z_2)\f$ -- can cancel that
 * quantity against the exponent analytically rather than forming both and
 * subtracting two large opposed terms.
 *
 * Valid for `scaled < -4`; the caller is responsible for the branch test.
 *
 * @tparam T a floating point type
 * @param scaled the standardized argument divided by sqrt(2)
 * @return the Cody branch value excluding `-square(scaled)`
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline T std_normal_lcdf_cody_correction(const T& scaled) {
  using std::log;
  const T x2 = square(scaled);
  const T x4 = pow(scaled, 4);
  const T x6 = pow(scaled, 6);
  const T x8 = pow(scaled, 8);
  const T x10 = pow(scaled, 10);
  const T temp_p = 0.000658749161529837803157 + 0.0160837851487422766278 / x2
                   + 0.125781726111229246204 / x4 + 0.360344899949804439429 / x6
                   + 0.305326634961232344035 / x8
                   + 0.0163153871373020978498 / x10;
  const T temp_q = -0.00233520497626869185443 - 0.0605183413124413191178 / x2
                   - 0.527905102951428412248 / x4 - 1.87295284992346047209 / x6
                   - 2.56852019228982242072 / x8 - 1.0 / x10;
  return LOG_HALF + log(INV_SQRT_PI + (temp_p / temp_q) / x2) - log(-scaled);
}

/**
 * Return \f$\log\Phi(\sqrt{2}\,\texttt{scaled})\f$.
 *
 * Rigorous numerical approximations are applied here to deal with values of
 * `|scaled| >> 0`. This is needed to deal with rare base-rate logistic
 * regression problems where it is useful to use an alternative link function
 * instead. `erfc()` is used instead of `erf()` in order to retain precision,
 * since for `x > 0` `erfc() -> 0`.
 *
 * See the branch documentation above.
 *
 * @tparam T a floating point type
 * @param scaled the standardized argument divided by sqrt(2)
 * @return the log of the standard normal cdf
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline T std_normal_lcdf_value(const T& scaled) {
  using std::fabs;
  using std::log;
  const T x2 = square(scaled);
  if (scaled > 0.0) {
    // CDF(x) = 1/2 + 1/2erf(x) = 1 - 1/2erfc(x)
    return log1p(-0.5 * erfc(scaled));
  } else if (scaled > -4.0) {
    // CDF(x) = 1/2 - 1/2erf(-x) = 1/2erfc(-x); -4 is R pnorm's M_SQRT_32
    // Since we scale by sqrt(2), we use sqrt(32)/sqrt(2) = 4
    return log(erfc(-scaled)) + LOG_HALF;
  } else if (10.0 * log(fabs(scaled)) < log(std::numeric_limits<T>::max())) {
    // entering territory where erfc(-x)~0
    // need to use direct numerical approximation of the log cdf instead
    // the following based on W. J. Cody, Math. Comp. 23(107):631-638 (1969)
    // CDF(x) = 1/2erfc(-x)
    return std_normal_lcdf_cody_correction(scaled) - x2;
  }
  // scaled^10 term will overflow
  return NEGATIVE_INFTY;
}

struct std_normal_lcdf_value_fun {
  template <typename T>
  static inline auto fun(T&& scaled) {
    return std_normal_lcdf_value(std::forward<T>(scaled));
  }
};

/**
 * A vectorized version of std_normal_lcdf_value().
 *
 * @tparam T a container type
 * @param scaled the standardized arguments divided by sqrt(2)
 * @return elementwise log of the standard normal cdf
 */
template <typename T, require_container_t<T>* = nullptr>
inline auto std_normal_lcdf_value(T&& scaled) {
  return apply_scalar_unary<std_normal_lcdf_value_fun, T>::apply(
      std::forward<T>(scaled));
}

/**
 * Return the derivative of \f$\log\Phi(\sqrt{2}\,\texttt{scaled})\f$ with
 * respect to `scaled`.
 *
 * Based on the analytic form `dln(CDF)/dx = exp(-x^2)/(sqrt(pi)*(1/2 +
 * erf(x)/2))`, calculated using a piecewise function due to instability and
 * inaccuracy in the various approximations. See the branch documentation
 * above.
 *
 * @tparam T a floating point type
 * @param scaled the standardized argument divided by sqrt(2)
 * @return d/d(scaled) of the log of the standard normal cdf
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline T std_normal_lcdf_grad(const T& scaled) {
  using std::exp;
  using std::fabs;
  using std::log;
  using std::pow;
  const T x2 = square(scaled);
  T t = 0.0;
  T t2 = 0.0;
  T t4 = 0.0;

  if (scaled > 2.9) {
    // approximation derived from Abramowitz and Stegun (1964) 7.1.26
    t = 1.0 / (1.0 + 0.3275911 * scaled);
    t2 = square(t);
    t4 = pow(t, 4);
    // A&S 7.1.26 puts exp(-x2) in the numerator; keep it there so it
    // underflows to zero instead of overflowing inside a denominator
    const T exp_m_x2 = exp(-x2);
    return INV_SQRT_PI * exp_m_x2
           / (1.0
              - exp_m_x2
                    * (0.254829592 - 0.284496736 * t + 1.421413741 * t2
                       - 1.453152027 * t2 * t + 1.061405429 * t4));
  } else if (scaled > 2.5) {
    // in the trouble area where all of the standard numerical
    // approximations are unstable - bridge the gap using Taylor
    // expansions of the analytic function
    // use Taylor expansion centred around x=2.7
    t = scaled - 2.7;
    t2 = square(t);
    t4 = pow(t, 4);
    return 0.0003849882382 - 0.002079084702 * t + 0.005229340880 * t2
           - 0.008029540137 * t2 * t + 0.008232190507 * t4
           - 0.005692364250 * t4 * t + 0.002399496363 * pow(t, 6);
  } else if (scaled > 2.1) {
    // use Taylor expansion centred around x=2.3
    t = scaled - 2.3;
    t2 = square(t);
    t4 = pow(t, 4);
    return 0.002846135439 - 0.01310032351 * t + 0.02732189391 * t2
           - 0.03326906904 * t2 * t + 0.02482478940 * t4
           - 0.009883071924 * t4 * t - 0.0002771362254 * pow(t, 6);
  } else if (scaled > 1.5) {
    // use Taylor expansion centred around x=1.85
    t = scaled - 1.85;
    t2 = square(t);
    t4 = pow(t, 4);
    return 0.01849212058 - 0.06876280470 * t + 0.1099906382 * t2
           - 0.09274533184 * t2 * t + 0.03543327418 * t4
           + 0.005644855518 * t4 * t - 0.01111434424 * pow(t, 6);
  } else if (scaled > 0.8) {
    // use Taylor expansion centred around x=1.15
    t = scaled - 1.15;
    t2 = square(t);
    t4 = pow(t, 4);
    return 0.1585747034 - 0.3898677543 * t + 0.3515963775 * t2
           - 0.09748053605 * t2 * t - 0.04347986191 * t4
           + 0.02182506378 * t4 * t + 0.01074751427 * pow(t, 6);
  } else if (scaled > 0.1) {
    // use Taylor expansion centred around x=0.45
    t = scaled - 0.45;
    t2 = square(t);
    t4 = pow(t, 4);
    return 0.6245634904 - 0.9521866949 * t + 0.3986215682 * t2
           + 0.04700850676 * t2 * t - 0.03478651979 * t4
           - 0.01772675404 * t4 * t + 0.0006577254811 * pow(t, 6);
  } else if (scaled < -29.0) {
    // asymptotic Mills ratio, DLMF 7.12.1: the derivative grows linearly as
    // -2*scaled, so no quadratic residual fit can track it
    const T inv_x2 = 1.0 / x2;
    return -2.0 * scaled
           / (1.0 + inv_x2 * (-0.5 + inv_x2 * (0.75 + inv_x2 * -1.875)));
  } else if (10.0 * log(fabs(scaled)) < log(std::numeric_limits<T>::max())) {
    // approximation derived from Abramowitz and Stegun (1964) 7.1.26
    // use fact that erf(x)=-erf(-x)
    // Abramowitz and Stegun define this for -inf<x<0 but seems to be
    // accurate for -inf<x<0.1
    t = 1.0 / (1.0 - 0.3275911 * scaled);
    t2 = square(t);
    t4 = pow(t, 4);
    T grad = 2.0 * INV_SQRT_PI
             / (0.254829592 * t - 0.284496736 * t2 + 1.421413741 * t2 * t
                - 1.453152027 * t4 + 1.061405429 * t4 * t);
    // check if we need to add a correction term
    // (from cubic fit of residuals)
    if (scaled < -17.0) {
      grad += 0.0001263257217272 * x2 * scaled + 0.0123586859488623 * x2
              - 0.0860505264736028 * scaled - 1.252783383752970;
    } else if (scaled < -7.0) {
      grad += 0.000471585349920831 * x2 * scaled + 0.0296839305424034 * x2
              + 0.207402143352332 * scaled + 0.425316974683324;
    } else if (scaled < -3.9) {
      grad += -0.0006972280656443 * x2 * scaled + 0.0068218494628567 * x2
              + 0.0585761964460277 * scaled + 0.1034397670201370;
    } else if (scaled < -2.1) {
      grad += -0.0018742199480885 * x2 * scaled - 0.0097119598291202 * x2
              - 0.0170137970924080 * scaled - 0.0100428567412041;
    }
    return grad;
  }
  return INFTY;
}

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
