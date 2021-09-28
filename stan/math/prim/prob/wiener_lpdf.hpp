// Original code from which Stan's code is derived:
// Copyright (c) 2013, Joachim Vandekerckhove.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted
// provided that the following conditions are met:
//
//   * Redistributions of source code must retain the above copyright notice,
//   * this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright notice,
//   * this list of conditions and the following disclaimer in the
//   * documentation and/or other materials provided with the distribution.
//   * Neither the name of the University of California, Irvine nor the names
//   * of its contributors may be used to endorse or promote products derived
//   * from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
// THE POSSIBILITY OF SUCH DAMAGE.

#ifndef STAN_MATH_PRIM_PROB_WIENER_LPDF_HPP
#define STAN_MATH_PRIM_PROB_WIENER_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/ceil.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <algorithm>
#include <cmath>
#include <string>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the first passage time density function for a (Wiener)
 *  drift diffusion model for the given \f$y\f$,
 * boundary separation \f$\alpha\f$, nondecision time \f$\tau\f$,
 * relative bias \f$\beta\f$, and drift rate \f$\delta\f$.
 * \f$\alpha\f$ and \f$\tau\f$ must be greater than 0, and
 * \f$\beta\f$ must be between 0 and 1. \f$y\f$ should contain
 * reaction times in seconds (strictly positive) with
 * upper-boundary responses.
 *
 * @tparam T_y type of scalar
 * @tparam T_alpha type of alpha parameter
 * @tparam T_tau type of tau parameter
 * @tparam T_beta type of beta parameter
 * @tparam T_delta type of delta parameter
 *
 * @param y A scalar variate.
 * @param alpha The boundary separation.
 * @param tau The nondecision time.
 * @param beta The relative bias.
 * @param delta The drift rate.
 * @return The log of the Wiener first passage time density of
 *  the specified arguments.
 */
template <bool propto, typename T_y, typename T_alpha, typename T_tau,
          typename T_beta, typename T_delta>
return_type_t<T_y, T_alpha, T_tau, T_beta, T_delta> wiener_lpdf(
    const T_y& y, const T_alpha& alpha, const T_tau& tau, const T_beta& beta,
    const T_delta& delta) {
  using T_return_type = return_type_t<T_y, T_alpha, T_tau, T_beta, T_delta>;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_alpha>;
  using T_tau_ref = ref_type_t<T_tau>;
  using T_beta_ref = ref_type_t<T_beta>;
  using T_delta_ref = ref_type_t<T_delta>;
  using std::ceil;
  using std::exp;
  using std::floor;
  using std::log;
  using std::sin;
  using std::sqrt;
  static const char* function = "wiener_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Boundary separation",
                         alpha, "A-priori bias", beta, "Nondecision time", tau,
                         "Drift rate", delta);

  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_tau_ref tau_ref = tau;
  T_beta_ref beta_ref = beta;
  T_delta_ref delta_ref = delta;

  check_positive(function, "Random variable", value_of(y_ref));
  check_positive_finite(function, "Boundary separation", value_of(alpha_ref));
  check_positive_finite(function, "Nondecision time", value_of(tau_ref));
  check_bounded(function, "A-priori bias", value_of(beta_ref), 0, 1);
  check_finite(function, "Drift rate", value_of(delta_ref));

  if (size_zero(y, alpha, beta, tau, delta)) {
    return 0;
  }

  T_return_type lp(0.0);

  size_t N = max_size(y, alpha, beta, tau, delta);
  if (!N) {
    return 0.0;
  }

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  scalar_seq_view<T_tau_ref> tau_vec(tau_ref);
  scalar_seq_view<T_delta_ref> delta_vec(delta_ref);
  size_t N_y_tau = max_size(y, tau);

  for (size_t i = 0; i < N_y_tau; ++i) {
    if (y_vec[i] <= tau_vec[i]) {
      std::stringstream msg;
      msg << ", but must be greater than nondecision time = " << tau_vec[i];
      std::string msg_str(msg.str());
      throw_domain_error(function, "Random variable", y_vec[i], " = ",
                         msg_str.c_str());
    }
  }

  if (!include_summand<propto, T_y, T_alpha, T_tau, T_beta, T_delta>::value) {
    return 0;
  }

  static const double WIENER_ERR = 0.000001;
  static const double PI_TIMES_WIENER_ERR = pi() * WIENER_ERR;
  static const double LOG_PI_LOG_WIENER_ERR = LOG_PI + log(WIENER_ERR);
  static const double TWO_TIMES_SQRT_TWO_PI_TIMES_WIENER_ERR
      = 2.0 * SQRT_TWO_PI * WIENER_ERR;
  static const double LOG_TWO_OVER_TWO_PLUS_LOG_SQRT_PI
      = LOG_TWO / 2 + LOG_SQRT_PI;
  static const double SQUARE_PI_OVER_TWO = square(pi()) * 0.5;
  static const double TWO_TIMES_LOG_SQRT_PI = 2.0 * LOG_SQRT_PI;

  for (size_t i = 0; i < N; i++) {
    typename scalar_type<T_beta>::type one_minus_beta = 1.0 - beta_vec[i];
    typename scalar_type<T_alpha>::type alpha2 = square(alpha_vec[i]);
    T_return_type x = (y_vec[i] - tau_vec[i]) / alpha2;
    T_return_type kl, ks, tmp = 0;
    T_return_type k, K;
    T_return_type sqrt_x = sqrt(x);
    T_return_type log_x = log(x);
    T_return_type one_over_pi_times_sqrt_x = 1.0 / pi() * sqrt_x;

    // calculate number of terms needed for large t:
    // if error threshold is set low enough
    if (PI_TIMES_WIENER_ERR * x < 1) {
      // compute bound
      kl = sqrt(-2.0 * SQRT_PI * (LOG_PI_LOG_WIENER_ERR + log_x)) / sqrt_x;
      // ensure boundary conditions met
      kl = (kl > one_over_pi_times_sqrt_x) ? kl : one_over_pi_times_sqrt_x;
    } else {
      kl = one_over_pi_times_sqrt_x;  // set to boundary condition
    }
    // calculate number of terms needed for small t:
    // if error threshold is set low enough
    T_return_type tmp_expr0 = TWO_TIMES_SQRT_TWO_PI_TIMES_WIENER_ERR * sqrt_x;
    if (tmp_expr0 < 1) {
      // compute bound
      ks = 2.0 + sqrt_x * sqrt(-2 * log(tmp_expr0));
      // ensure boundary conditions are met
      T_return_type sqrt_x_plus_one = sqrt_x + 1.0;
      ks = (ks > sqrt_x_plus_one) ? ks : sqrt_x_plus_one;
    } else {     // if error threshold was set too high
      ks = 2.0;  // minimal kappa for that case
    }
    if (ks < kl) {   // small t
      K = ceil(ks);  // round to smallest integer meeting error
      T_return_type tmp_expr1 = (K - 1.0) / 2.0;
      T_return_type tmp_expr2 = ceil(tmp_expr1);
      for (k = -floor(tmp_expr1); k <= tmp_expr2; k++) {
        tmp += (one_minus_beta + 2.0 * k)
               * exp(-(square(one_minus_beta + 2.0 * k)) * 0.5 / x);
      }
      tmp = log(tmp) - LOG_TWO_OVER_TWO_PLUS_LOG_SQRT_PI - 1.5 * log_x;
    } else {         // if large t is better...
      K = ceil(kl);  // round to smallest integer meeting error
      for (k = 1; k <= K; ++k) {
        tmp += k * exp(-(square(k)) * (SQUARE_PI_OVER_TWO * x))
               * sin(k * pi() * one_minus_beta);
      }
      tmp = log(tmp) + TWO_TIMES_LOG_SQRT_PI;
    }

    // convert to f(t|v,a,w) and return result
    lp += delta_vec[i] * alpha_vec[i] * one_minus_beta
          - square(delta_vec[i]) * x * alpha2 / 2.0 - log(alpha2) + tmp;
  }
  return lp;
}

template <typename T_y, typename T_alpha, typename T_tau, typename T_beta,
          typename T_delta>
inline return_type_t<T_y, T_alpha, T_tau, T_beta, T_delta> wiener_lpdf(
    const T_y& y, const T_alpha& alpha, const T_tau& tau, const T_beta& beta,
    const T_delta& delta) {
  return wiener_lpdf<false>(y, alpha, tau, beta, delta);
}

}  // namespace math
}  // namespace stan
#endif
