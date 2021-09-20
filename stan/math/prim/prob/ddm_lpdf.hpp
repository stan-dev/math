#ifndef STAN_MATH_PRIM_PROB_DDM_LPDF_HPP
#define STAN_MATH_PRIM_PROB_DDM_LPDF_HPP

#define _USE_MATH_DEFINES
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
 * The log of the first passage time density function for a (Ratcliff, 1978)
 * drift diffusion model with intrinsic trial-trial variability
 * for the given response time \f$rt\f$, response \f$response\f$, boundary
 * separation \f$a\f$, mean drift rate across trials \f$v\f$, non-decision
 * time \f$t0\f$, relative bias \f$w\f$, and standard deviation of drift rate
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
return_type_t<T_rt, T_response, T_a, T_v, T_t0, T_w, T_sv> ddm_lpdf(
    const T_rt& rt, const T_response& response, const T_a& a, const T_v& v,
    const T_t0& t0, const T_w& w, const T_sv& sv) {
  using T_return_type
      = return_type_t<T_rt, T_response, T_a, T_v, T_t0, T_w, T_sv>;
  using stan::ref_type_t;
  using stan::scalar_seq_view;
  using stan::math::include_summand;
  using stan::math::invalid_argument;
  using stan::math::throw_domain_error;
  using std::ceil;
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
  static const char* function = "ddm_lpdf";
  static const int max_terms_large = 1;    // heuristic for switching mechanism
  static const double ERR_TOL = 0.000001;  // error tolerance for PDF approx
  static const double SV_THRESH = 0;  // threshold for using variable drift rate
  static const double LOG_PI = log(M_PI);
  static const double LOG_2PI_2 = 0.5 * log(2 * M_PI);

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
      < 1) {  // rt, invalid inputs will be handled in calculation of the pdf
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
        } else {  // a = Inf implies PDF = log(0)
          throw_domain_error(function, "a (threshold separation)", a_vec[i],
                             " = ", ", but it must be positive and finite");
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
        if (isfinite(
                t0_vec[i])) {  // this could also be handled in calculate_pdf()
          continue;
        } else {  // t0 = Inf implies rt - t0 < 0 implies PDF = log(0)
          throw_domain_error(function, "t0 (non-decision time)", t0_vec[i],
                             " = ", ", but it must be positive and finite");
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
        } else {  // sv = Inf implies PDF = log(0)
          throw_domain_error(
              function, "sv (standard deviation of drift rate across trials)",
              sv_vec[i], " = ", ", but it must be positive and finite");
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

  // Calculate log(PDF)
  T_return_type lp(0.0);
  double t, a_i, v_i, w_i, sv_i;
  for (size_t i = 0; i < Nmax; i++) {
    // Check Parameter Values
    t = rt_vec[i % Nrt]
        - t0_vec[i % Nt0];            // response time minus non-decision time
    if (t > 0 && std::isfinite(t)) {  // sort response and calculate density
      a_i = a_vec[i % Na];
      sv_i = sv_vec[i % Nsv];
      if (out[i] == 1) {  // response is "lower" so use unchanged parameters
        v_i = v_vec[i % Nv];
        w_i = w_vec[i % Nw];
      } else {  // response is "upper" so use alternate parameters
        v_i = -v_vec[i % Nv];
        w_i = 1 - w_vec[i % Nw];
      }

      // Approximate log(PDF)
      double mult;
      // Check large time
      if (sv_i <= SV_THRESH) {  // no sv
        mult = -v_i * a_i * w_i - v_i * v_i * t / 2 - 2 * log(a_i);
      } else {  // sv
        mult = (sv_i * sv_i * a_i * a_i * w_i * w_i - 2 * v_i * a_i * w_i
                - v_i * v_i * t)
                   / (2 + 2 * sv_i * sv_i * t)
               - 0.5 * log(1 + sv_i * sv_i * t) - 2 * log(a_i);
      }
      int kl;
      double exp_err = ERR_TOL * exp(-mult);
      double taa = t / (a_i * a_i);
      double bc = 1 / (M_PI * sqrt(taa));  // boundary conditions
      if (bc > INT_MAX)
        return INT_MAX;
      if (exp_err * M_PI * taa < 1) {  // error threshold is low enough
        double kl_tmp
            = sqrt(-2 * log(M_PI * taa * exp_err) / (M_PI * M_PI * taa));
        if (kl_tmp > INT_MAX) {
          kl = INT_MAX;
        } else {
          kl = ceil(max(kl_tmp, bc));  // ensure boundary conditions are met
        }
      } else {
        kl = ceil(bc);  // else set to boundary condition
      }

      // Compare kl (large time) to max_terms_large (small time)
      if (kl <= max_terms_large) {  // use large time
        double gamma = -0.5 * M_PI * M_PI * taa;
        double sum = 0.0;
        for (size_t j = 1; j <= kl; j++) {
          sum += j * sin(j * w_i * M_PI) * exp(gamma * j * j);
        }
        if (sum >= 0) {  // if result is negative, don't add to lp
          lp += LOG_PI + mult + log(sum);
        }
      } else {                    // use small time
        if (sv_i <= SV_THRESH) {  // no sv
          mult = log(a_i) - LOG_2PI_2 - 1.5 * log(t) - v_i * a_i * w_i
                 - v_i * v_i * t / 2;
        } else {  // sv
          mult = log(a_i) - 1.5 * log(t) - LOG_2PI_2
                 - 0.5 * log(1 + sv_i * sv_i * t)
                 + (sv_i * sv_i * a_i * a_i * w_i * w_i - 2 * v_i * a_i * w_i
                    - v_i * v_i * t)
                       / (2 + 2 * sv_i * sv_i * t);
        }
        exp_err = ERR_TOL / exp(mult);
        size_t minterms
            = sqrt(t) / a_i - w_i;  // min number of terms, truncates toward 0
        double gamma = -1 / (2 * taa);
        double sum = w_i * exp(gamma * w_i * w_i);  // initialize with j=0 term
        double term, rj;
        size_t j = 0;
        if (minterms % 2) {  // minterms is odd (and at least 1)
          j++;
          rj = j + 1 - w_i;
          term = rj * exp(gamma * rj * rj);
          sum -= term;
          while (j < minterms) {
            j++;
            rj = j + w_i;
            sum += rj * exp(gamma * rj * rj);
            j++;
            rj = j + 1 - w_i;
            term = rj * exp(gamma * rj * rj);
            sum -= term;
          }
          j++;
          rj = j + w_i;  // j is now even
          term = rj * exp(gamma * rj * rj);
          sum += term;
          while (term > exp_err) {
            j++;
            rj = j + 1 - w_i;
            term = rj * exp(gamma * rj * rj);
            sum -= term;
            if (term <= exp_err)
              break;
            j++;
            rj = j + w_i;
            term = rj * exp(gamma * rj * rj);
            sum += term;
          }
        } else {                  // minterms is even (and at least 0)
          while (j < minterms) {  // j is currently 0
            j++;
            rj = j + 1 - w_i;
            sum -= rj * exp(gamma * rj * rj);
            j++;
            rj = j + w_i;
            term = rj * exp(gamma * rj * rj);
            sum += term;
          }
          j++;
          rj = j + 1 - w_i;  // j is now odd
          term = rj * exp(gamma * rj * rj);
          sum -= term;
          while (term > exp_err) {
            j++;
            rj = j + w_i;
            term = rj * exp(gamma * rj * rj);
            sum += term;
            if (term <= exp_err)
              break;
            j++;
            rj = j + 1 - w_i;
            term = rj * exp(gamma * rj * rj);
            sum -= term;
          }
        }
        if (sum >= 0) {  // if result is negative, don't add to lp
          lp += mult + log(sum);
        }
      }
    } else {  // {NaN, NA} evaluate to FALSE
      if (isnan(t)) {
        throw_domain_error(
            function, "rt (response time)", rt_vec[i % Nrt],
            "is a NaN and = ", ", but rt must be positive and finite");
      } else {
        throw_domain_error(
            function, "rt (response time)", t0_vec[i % Nt0],
            "is not greater than t0 = ",
            ", but it must be that rt - t0 is positive and finite");
      }
    }
  }

  return lp;
}

template <typename T_rt, typename T_response, typename T_a, typename T_v,
          typename T_t0, typename T_w, typename T_sv>
inline return_type_t<T_rt, T_response, T_a, T_v, T_t0, T_w, T_sv> ddm_lpdf(
    const T_rt& rt, const T_response& response, const T_a& a, const T_v& v,
    const T_t0& t0, const T_w& w, const T_sv& sv) {
  return ddm_lpdf<false>(rt, response, a, v, t0, w, sv);
}

}  // namespace math
}  // namespace stan
#endif
