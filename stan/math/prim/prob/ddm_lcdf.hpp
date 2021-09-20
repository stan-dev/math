#ifndef STAN_MATH_PRIM_PROB_DDM_LCDF_HPP
#define STAN_MATH_PRIM_PROB_DDM_LCDF_HPP

#include <cmath>
#include <string>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/meta/ref_type.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/invalid_argument.hpp>

// Open the Namespace
namespace stan {
namespace math {
using stan::return_type_t;

/**
 * The log of the first passage time distribution function for a
 * (Ratcliff, 1978) drift diffusion model with intrinsic trial-trial variability
 * for the given response time \f$rt\f$, response \f$response\f$, boundary
 * separation \f$a\f$, mean drift rate across trials \f$v\f$, non-decision
 * time \f$t0\f$, relative bias \f$w\f$, standard deviation of drift rate
 * across trials \f$sv\f$.
 *
 * @tparam T_rt type of parameter `rt`
 * @tparam T_response type of parameter `response`
 * @tparam T_a type of parameter `a`
 * @tparam T_v type of parameter `v`
 * @tparam T_t0 type of parameter `t0`
 * @tparam T_w type of parameter `w`
 * @tparam T_sv type of parameter `sv`
 *
 * @param rt The response time; rt >= 0.
 * @param response The response; response in {1, 2}.
 * @param a The threshold separation; a > 0.
 * @param v The mean drift rate across trials.
 * @param t0 The non-decision time; t0 >= 0.
 * @param w The relative a priori bias; 0 < w < 1.
 * @param sv The standard deviation of drift rate across trials; sv >= 0.
 * @return The log of the Wiener first passage time density of
 *  the specified arguments.
 */
template <bool propto, typename T_rt, typename T_response, typename T_a,
          typename T_v, typename T_t0, typename T_w, typename T_sv>
return_type_t<T_rt, T_response, T_a, T_v, T_t0, T_w, T_sv> ddm_lcdf(
    const T_rt& rt, const T_response& response, const T_a& a, const T_v& v,
    const T_t0& t0, const T_w& w, const T_sv& sv) {
  using T_return_type
      = return_type_t<T_rt, T_response, T_a, T_v, T_t0, T_w, T_sv>;
  using stan::ref_type_t;
  using stan::scalar_seq_view;
  using stan::math::include_summand;
  using stan::math::invalid_argument;
  using stan::math::throw_domain_error;
  using std::erf;
  using std::exp;
  using std::isfinite;
  using std::isnan;
  using std::log;
  using std::max;
  using std::sqrt;
  using std::vector;
  using T_rt_ref = ref_type_t<T_rt>;
  using T_response_ref = ref_type_t<T_response>;
  using T_a_ref = ref_type_t<T_a>;
  using T_v_ref = ref_type_t<T_v>;
  using T_t0_ref = ref_type_t<T_t0>;
  using T_w_ref = ref_type_t<T_w>;
  using T_sv_ref = ref_type_t<T_sv>;

  // Constants
  static const char* function = "ddm_lcdf";
  static const double ERR_TOL = 0.000001;  // error tolerance for PDF approx
  static const double PI_CONST = 3.14159265358979323846;  // define pi like C++
  static const double SQRT_2PI = sqrt(2 * PI_CONST);
  static const double SQRT_2PI_INV = 1 / SQRT_2PI;
  static const double SQRT_2_INV_NEG = -1 / sqrt(2);

  // Convert Inputs
  T_rt_ref rt_ref = rt;
  T_response_ref response_ref = response;
  T_a_ref a_ref = a;
  T_v_ref v_ref = v;
  T_t0_ref t0_ref = t0;
  T_w_ref w_ref = w;
  T_sv_ref sv_ref = sv;
  scalar_seq_view<T_rt_ref> rt_vec(rt_ref);
  scalar_seq_view<T_response_ref> response_vec(response_ref);
  scalar_seq_view<T_a_ref> a_vec(a_ref);
  scalar_seq_view<T_v_ref> v_vec(v_ref);
  scalar_seq_view<T_t0_ref> t0_vec(t0_ref);
  scalar_seq_view<T_w_ref> w_vec(w_ref);
  scalar_seq_view<T_sv_ref> sv_vec(sv_ref);

  // Parameter Checks
  size_t Nrt = rt_vec.size();
  size_t Nres = response_vec.size();
  size_t Na = a_vec.size();
  size_t Nv = v_vec.size();
  size_t Nt0 = t0_vec.size();
  size_t Nw = w_vec.size();
  size_t Nsv = sv_vec.size();
  size_t Nmax = max({Nrt, Nres, Na, Nv, Nt0, Nw, Nsv});
  vector<int> out(Nmax);  // initialize output-checking vector

  if (Nrt
      < 1) {  // rt, invalid inputs will be handled in calculation of the CDF
    return 0;
  }

  if (Nres < 1) {  // response
    return 0;
  } else {
    for (size_t i = 0; i < Nres; i++) {
      if (response_vec[i] == 1) {  // lower
        for (size_t j = i; j < Nmax; j += Nres) {
          out[j] = 1;
        }
      } else if (response_vec[i] == 2) {  // upper
        for (size_t j = i; j < Nmax; j += Nres) {
          out[j] = 2;
        }
      } else {  // {NaN, NA} evaluate to FALSE
        throw_domain_error(function, "response", response_vec[i], " = ",
                           ", but it must be either 1 (lower) or 2 (upper)");
      }
    }
  }

  if (Na < 1) {  // a
    return 0;
  } else {
    for (size_t i = 0; i < Na; i++) {
      if (a_vec[i] > 0) {
        if (isfinite(a_vec[i])) {
          continue;
        } else {  // a = Inf implies PDF = log(0) and CDF problems
          throw_domain_error(function, "a (threshold separation)", a_vec[i],
                             " = ", ", but it must be finite");
        }
      } else {  // {NaN, NA} evaluate to FALSE
        throw_domain_error(function, "a (threshold separation)", a_vec[i],
                           " = ", ", but it must be positive and finite");
      }
    }
  }

  if (Nv < 1) {  // v
    return 0;
  } else {
    for (size_t i = 0; i < Nv; i++) {
      if (isfinite(v_vec[i])) {
        continue;
      } else {  // NaN, NA, Inf, -Inf are not finite
        throw_domain_error(function, "v (drift rate)", v_vec[i], " = ",
                           ", but it must be finite");
      }
    }
  }

  if (Nt0 < 1) {  // t0
    return 0;
  } else {
    for (size_t i = 0; i < Nt0; i++) {
      if (t0_vec[i] >= 0) {
        if (isfinite(t0_vec[i])) {  // this could also be handled in calculation
                                    // of CDF
          continue;
        } else {  // t0 = Inf implies rt - t0 < 0 implies CDF = log(0)
          throw_domain_error(function, "t0 (non-decision time)", t0_vec[i],
                             " = ", ", but it must be finite");
        }
      } else {  // {NaN, NA} evaluate to FALSE
        throw_domain_error(function, "t0 (non-decision time)", t0_vec[i], " = ",
                           ", but it must be positive and finite");
      }
    }
  }

  if (Nw < 1) {  // w
    return 0;
  } else {
    for (size_t i = 0; i < Nw; i++) {
      if (w_vec[i] > 0 && w_vec[i] < 1) {
        continue;
      } else {  // {NaN, NA} evaluate to FALSE
        throw_domain_error(function, "w (relative a priori bias)", w_vec[i],
                           " = ", ", but it must be that 0 < w < 1");
      }
    }
  }

  if (Nsv < 1) {  // sv
    return 0;
  } else {
    for (size_t i = 0; i < Nsv; i++) {
      if (sv_vec[i] >= 0) {
        if (isfinite(sv_vec[i])) {
          continue;
        } else {  // sv = Inf implies PDF = log(0) and CDF problems
          throw_domain_error(
              function, "sv (standard deviation of drift rate across trials)",
              sv_vec[i], " = ", ", but it must be finite");
        }
      } else {  // {NaN, NA} evaluate to FALSE
        throw_domain_error(
            function, "sv (standard deviation of drift rate across trials)",
            sv_vec[i], " = ", ", but it must be positive and finite");
      }
    }
  }

  if (!include_summand<propto, T_rt, T_response, T_a, T_v, T_t0, T_w,
                       T_sv>::value) {
    return 0;
  }

  // Calculate log(CDF)
  T_return_type lp(0.0);
  double t, a_i, v_i, w_i, sv_i;
  for (size_t i = 0; i < Nmax; i++) {
    // Check Parameter Values
    t = rt_vec[i % Nrt]
        - t0_vec[i % Nt0];  // response time minus non-decision time
    if (t > 0) {            // sort response and calculate density
      a_i = a_vec[i % Na];
      sv_i = sv_vec[i % Nsv];
      if (out[i] == 1) {  // response is "lower" so use unchanged parameters
        v_i = v_vec[i % Nv];
        w_i = w_vec[i % Nw];
      } else {  // response is "upper" so use alternate parameters
        v_i = -v_vec[i % Nv];
        w_i = 1 - w_vec[i % Nw];
      }

      if (t > 32) {  // approximation for t = +Infinity
        t = 32;
      }

      // Calculate sum multiplier
      double mult = (sv_i * sv_i * a_i * a_i * w_i * w_i - 2 * v_i * a_i * w_i
                     - v_i * v_i * t)
                    / (2 + 2 * sv_i * sv_i * t);

      // Scale error so it is valid inside the sum (not logged)
      double exp_err = ERR_TOL * exp(-mult);

      // Calculate sum
      double sum = 0;
      double gamma = v_i - sv_i * sv_i * a_i * w_i;
      double lambda = 1 + sv_i * sv_i * t;
      double rho = sqrt(t * lambda);

      int j = 0;
      double rj = a_i * j + a_i * w_i;
      double m1 = (lambda * rj - gamma * t) / rho;
      double m2 = (lambda * rj + gamma * t) / rho;
      double mills_1, mills_2;
      if (m1 < 6.5) {
        mills_1 = SQRT_2PI * 0.5 * exp(0.5 * m1 * m1)
                  * (1 + erf(SQRT_2_INV_NEG * m1));
      } else {
        double m1sq = m1 * m1;
        mills_1 = (1 - 1 / (m1sq + 2) + 1 / ((m1sq + 2) * (m1sq + 4))
                   - 5 / ((m1sq + 2) * (m1sq + 4) * (m1sq + 6))
                   + 9 / ((m1sq + 2) * (m1sq + 4) * (m1sq + 6) * (m1sq + 8))
                   - 129
                         / ((m1sq + 2) * (m1sq + 4) * (m1sq + 6) * (m1sq + 8)
                            * (m1sq + 10)))
                  / m1;
      }
      if (m2 < 6.5) {
        mills_2 = SQRT_2PI * 0.5 * exp(0.5 * m2 * m2)
                  * (1 + erf(SQRT_2_INV_NEG * m2));
      } else {
        double m2sq = m2 * m2;
        mills_2 = (1 - 1 / (m2sq + 2) + 1 / ((m2sq + 2) * (m2sq + 4))
                   - 5 / ((m2sq + 2) * (m2sq + 4) * (m2sq + 6))
                   + 9 / ((m2sq + 2) * (m2sq + 4) * (m2sq + 6) * (m2sq + 8))
                   - 129
                         / ((m2sq + 2) * (m2sq + 4) * (m2sq + 6) * (m2sq + 8)
                            * (m2sq + 10)))
                  / m2;
      }
      double term
          = SQRT_2PI_INV * exp(-0.5 * rj * rj / t) * (mills_1 + mills_2);
      sum += term;

      while (term > exp_err) {
        if (j > 1000) {
          // maybe include a warning here?
          break;
        }
        j++;
        rj = a_i * j + a_i * (1 - w_i);
        m1 = (lambda * rj - gamma * t) / rho;
        m2 = (lambda * rj + gamma * t) / rho;
        if (m1 < 6.5) {
          mills_1 = SQRT_2PI * 0.5 * exp(0.5 * m1 * m1)
                    * (1 + erf(SQRT_2_INV_NEG * m1));
        } else {
          double m1sq = m1 * m1;
          mills_1 = (1 - 1 / (m1sq + 2) + 1 / ((m1sq + 2) * (m1sq + 4))
                     - 5 / ((m1sq + 2) * (m1sq + 4) * (m1sq + 6))
                     + 9 / ((m1sq + 2) * (m1sq + 4) * (m1sq + 6) * (m1sq + 8))
                     - 129
                           / ((m1sq + 2) * (m1sq + 4) * (m1sq + 6) * (m1sq + 8)
                              * (m1sq + 10)))
                    / m1;
        }
        if (m2 < 6.5) {
          mills_2 = SQRT_2PI * 0.5 * exp(0.5 * m2 * m2)
                    * (1 + erf(SQRT_2_INV_NEG * m2));
        } else {
          double m2sq = m2 * m2;
          mills_2 = (1 - 1 / (m2sq + 2) + 1 / ((m2sq + 2) * (m2sq + 4))
                     - 5 / ((m2sq + 2) * (m2sq + 4) * (m2sq + 6))
                     + 9 / ((m2sq + 2) * (m2sq + 4) * (m2sq + 6) * (m2sq + 8))
                     - 129
                           / ((m2sq + 2) * (m2sq + 4) * (m2sq + 6) * (m2sq + 8)
                              * (m2sq + 10)))
                    / m2;
        }
        term = SQRT_2PI_INV * exp(-0.5 * rj * rj / t) * (mills_1 + mills_2);
        sum -= term;

        if (term <= exp_err)
          break;

        j++;
        rj = a_i * j + a_i * w_i;
        m1 = (lambda * rj - gamma * t) / rho;
        m2 = (lambda * rj + gamma * t) / rho;
        if (m1 < 6.5) {
          mills_1 = SQRT_2PI * 0.5 * exp(0.5 * m1 * m1)
                    * (1 + erf(SQRT_2_INV_NEG * m1));
        } else {
          double m1sq = m1 * m1;
          mills_1 = (1 - 1 / (m1sq + 2) + 1 / ((m1sq + 2) * (m1sq + 4))
                     - 5 / ((m1sq + 2) * (m1sq + 4) * (m1sq + 6))
                     + 9 / ((m1sq + 2) * (m1sq + 4) * (m1sq + 6) * (m1sq + 8))
                     - 129
                           / ((m1sq + 2) * (m1sq + 4) * (m1sq + 6) * (m1sq + 8)
                              * (m1sq + 10)))
                    / m1;
        }
        if (m2 < 6.5) {
          mills_2 = SQRT_2PI * 0.5 * exp(0.5 * m2 * m2)
                    * (1 + erf(SQRT_2_INV_NEG * m2));
        } else {
          double m2sq = m2 * m2;
          mills_2 = (1 - 1 / (m2sq + 2) + 1 / ((m2sq + 2) * (m2sq + 4))
                     - 5 / ((m2sq + 2) * (m2sq + 4) * (m2sq + 6))
                     + 9 / ((m2sq + 2) * (m2sq + 4) * (m2sq + 6) * (m2sq + 8))
                     - 129
                           / ((m2sq + 2) * (m2sq + 4) * (m2sq + 6) * (m2sq + 8)
                              * (m2sq + 10)))
                    / m2;
        }
        term = SQRT_2PI_INV * exp(-0.5 * rj * rj / t) * (mills_1 + mills_2);
        sum += term;
      }

      // Add sum and multiplier to lp
      if (sum >= 0) {  // if result is negative, treat as 0 and don't add to lp
        lp += mult + log(sum);
      }
    } else {  // {NaN, NA} evaluate to FALSE
      if (isnan(t)) {
        throw_domain_error(function, "rt (response time)", rt_vec[i % Nrt],
                           "is a NaN and = ", ", but this value is invalid");
      } else {
        throw_domain_error(
            function, "rt (response time)", t0_vec[i % Nt0],
            "is not greater than t0 = ", ", but it must be that rt - t0 > 0");
      }
    }
  }

  return lp;
}

template <typename T_rt, typename T_response, typename T_a, typename T_v,
          typename T_t0, typename T_w, typename T_sv>
inline return_type_t<T_rt, T_response, T_a, T_v, T_t0, T_w, T_sv> ddm_lcdf(
    const T_rt& rt, const T_response& response, const T_a& a, const T_v& v,
    const T_t0& t0, const T_w& w, const T_sv& sv) {
  return ddm_lcdf<false>(rt, response, a, v, t0, w, sv);
}

}  // namespace math
}  // namespace stan
#endif
