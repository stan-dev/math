#ifndef STAN_MATH_PRIM_PROB_WIENER_FULL_LCCDF_UNNORM_HPP
#define STAN_MATH_PRIM_PROB_WIENER_FULL_LCCDF_UNNORM_HPP

#include <stan/math/prim/prob/wiener4_lccdf_unnorm.hpp>
#include <stan/math/prim/prob/wiener_full_lcdf_unnorm.hpp>

namespace stan {
namespace math {
namespace internal {

/**
 * Calculate the derivative of the wiener7 density w.r.t. 'sw'
 *
 * @tparam T_y type of reaction time
 * @tparam T_a type of boundary separation
 * @tparam T_v type of drift rate
 * @tparam T_w type of relative starting point
 * @tparam T_sv type of inter-trial variability in v
 * @tparam T_err type of log error tolerance
 *
 * @param y The reaction time in seconds
 * @param a The boundary separation
 * @param v The drift rate
 * @param w The relative starting point
 * @param sw The inter-trial variability of the relative starting point
 * @param wildcard This parameter space is needed for a functor. Could be
 * deleted when another solution is found
 * @param log_error The log error tolerance
 * @return Gradient with respect to sw
 */
template <typename T_y, typename T_a, typename T_v, typename T_w, typename T_sw,
          typename T_err>
inline auto wiener7_ccdf_grad_sw(const T_y& y, const T_a& a, const T_v& v,
                                 const T_w& w, const T_sw& sw,
                                 T_err log_error) {
  auto low = fmax(0.0, w - sw / 2.0);
  auto high = fmin(1.0, w + sw / 2.0);

  const auto lower_value = wiener4_ccdf(y, a, v, low, log_error);
  const auto upper_value = wiener4_ccdf(y, a, v, high, log_error);
  return 0.5 * (lower_value + upper_value) / sw;
}

}  // namespace internal

/** \ingroup prob_dists
 * Returns the log CCDF of the Wiener distribution for a
 * (Wiener) drift diffusion model with up to 7 parameters. If containers
 * are supplied, returns the log sum of the probabilities.
 * See 'wiener_full_lpdf' for more comprehensive documentation
 * As the CDF goes to the probability to hit the upper bound
 * (instead of 1, as it is usually the case) when the reaction time
 * goes to infinity, the CCDF is defined as
 * ccdf = probability_to_hit_the_upper_bound - cdf.
 *
 * @tparam T_y type of reaction time
 * @tparam T_a type of boundary separation
 * @tparam T_t0 type of non-decision time
 * @tparam T_w type of relative starting point
 * @tparam T_v type of drift rate
 * @tparam T_sv type of inter-trial variability of drift rate
 * @tparam T_sw type of inter-trial variability of relative starting point
 * @tparam T_st0 type of inter-trial variability of non-decision time
 *
 * @param y The reaction time in seconds
 * @param a The boundary separation
 * @param t0 The non-decision time
 * @param w The relative starting point
 * @param v The drift rate
 * @param sv The inter-trial variability of the drift rate
 * @param sw The inter-trial variability of the relative starting point
 * @param st0 The inter-trial variability of the non-decision time
 * @param precision_derivatives Level of precision in estimation of partial
 * derivatives
 *
 * @return log probability or log sum of probabilities for upper
 * boundary responses
 * @throw std::domain_error if non-decision time \c t0 is greater than reaction
 * time \c y.
 * @throw std::domain_error if \c 1-sw/2 is smaller than or equal to \c w.
 * @throw std::domain_error if \c sw/2 is larger than or equal to \c w.
 */

template <bool propto = false, typename T_y, typename T_a, typename T_t0,
          typename T_w, typename T_v, typename T_sv, typename T_sw,
          typename T_st0>
inline auto wiener_lccdf_unnorm(const T_y& y, const T_a& a, const T_t0& t0,
                                const T_w& w, const T_v& v, const T_sv& sv,
                                const T_sw& sw, const T_st0& st0,
                                const double& precision_derivatives = 1e-8) {
  using ret_t = return_type_t<T_y, T_a, T_t0, T_w, T_v, T_sv, T_sw, T_st0>;
  using T_y_ref = ref_type_t<T_y>;
  using T_a_ref = ref_type_t<T_a>;
  using T_v_ref = ref_type_t<T_v>;
  using T_w_ref = ref_type_t<T_w>;
  using T_t0_ref = ref_type_t<T_t0>;
  using T_sv_ref = ref_type_t<T_sv>;
  using T_sw_ref = ref_type_t<T_sw>;
  using T_st0_ref = ref_type_t<T_st0>;
  using internal::GradientCalc;
  using T_partials_return
      = partials_return_t<T_y, T_a, T_t0, T_w, T_v, T_sv, T_sw, T_st0>;

  T_y_ref y_ref = y;
  T_a_ref a_ref = a;
  T_v_ref v_ref = v;
  T_w_ref w_ref = w;
  T_t0_ref t0_ref = t0;
  T_sv_ref sv_ref = sv;
  T_sw_ref sw_ref = sw;
  T_st0_ref st0_ref = st0;

  auto y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  auto a_val = to_ref(as_value_column_array_or_scalar(a_ref));
  auto v_val = to_ref(as_value_column_array_or_scalar(v_ref));
  auto w_val = to_ref(as_value_column_array_or_scalar(w_ref));
  auto t0_val = to_ref(as_value_column_array_or_scalar(t0_ref));
  auto sv_val = to_ref(as_value_column_array_or_scalar(sv_ref));
  auto sw_val = to_ref(as_value_column_array_or_scalar(sw_ref));
  auto st0_val = to_ref(as_value_column_array_or_scalar(st0_ref));

  if (!include_summand<propto, T_y, T_a, T_v, T_w, T_t0, T_sv, T_sw,
                       T_st0>::value) {
    return ret_t(0.0);
  }

  static constexpr const char* function_name = "wiener_lccdf_unnorm";
  check_consistent_sizes(function_name, "Random variable", y,
                         "Boundary separation", a, "Drift rate", v,
                         "A-priori bias", w, "Nondecision time", t0,
                         "Inter-trial variability in drift rate", sv,
                         "Inter-trial variability in A-priori bias", sw,
                         "Inter-trial variability in Nondecision time", st0);
  check_positive_finite(function_name, "Random variable", y_val);
  check_positive_finite(function_name, "Boundary separation", a_val);
  check_finite(function_name, "Drift rate", v_val);
  check_less(function_name, "A-priori bias", w_val, 1);
  check_greater(function_name, "A-priori bias", w_val, 0);
  check_nonnegative(function_name, "Nondecision time", t0_val);
  check_finite(function_name, "Nondecision time", t0_val);
  check_nonnegative(function_name, "Inter-trial variability in drift rate",
                    sv_val);
  check_finite(function_name, "Inter-trial variability in drift rate", sv_val);
  check_bounded(function_name, "Inter-trial variability in A-priori bias",
                sw_val, 0, 1);
  check_nonnegative(function_name,
                    "Inter-trial variability in Nondecision time", st0_val);
  check_finite(function_name, "Inter-trial variability in Nondecision time",
               st0_val);

  const size_t N = max_size(y, a, v, w, t0, sv, sw, st0);
  if (N == 0) {
    return ret_t(0.0);
  }
  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_a_ref> a_vec(a_ref);
  scalar_seq_view<T_v_ref> v_vec(v_ref);
  scalar_seq_view<T_w_ref> w_vec(w_ref);
  scalar_seq_view<T_t0_ref> t0_vec(t0_ref);
  scalar_seq_view<T_sv_ref> sv_vec(sv_ref);
  scalar_seq_view<T_sw_ref> sw_vec(sw_ref);
  scalar_seq_view<T_st0_ref> st0_vec(st0_ref);
  const size_t N_y_t0 = max_size(y, t0, st0);

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
    if (unlikely(w_vec[i] - .5 * sw_vec[i] <= 0)) {
      std::stringstream msg;
      msg << ", but must be smaller than 2*(A-priori bias) = "
          << 2.0 * w_vec[i];
      std::string msg_str(msg.str());
      throw_domain_error(function_name,
                         "Inter-trial variability in A-priori bias", sw_vec[i],
                         " = ", msg_str.c_str());
    }
    if (unlikely(w_vec[i] + .5 * sw_vec[i] >= 1)) {
      std::stringstream msg;
      msg << ", but must be smaller than 2*(1-A-priori bias) = "
          << 2.0 * (1.0 - w_vec[i]);
      std::string msg_str(msg.str());
      throw_domain_error(function_name,
                         "Inter-trial variability in A-priori bias", sw_vec[i],
                         " = ", msg_str.c_str());
    }
  }

  // for precs. 1e-6, 1e-12, see Hartmann et al. (2021), Henrich et al. (2023)
  const T_partials_return log_error_cdf = log(1e-6);  // precision for density
  const auto error_bound = precision_derivatives;     // precision for
  // derivatives (controllable by user)
  const auto log_error_derivatives = log(error_bound);
  const T_partials_return absolute_error_hcubature = 0.0;
  const T_partials_return relative_error_hcubature
      = .9 * error_bound;  // eps_rel(Integration)
  const T_partials_return log_error_absolute = log(1e-12);
  const int maximal_evaluations_hcubature = 6000;
  T_partials_return lccdf = 0.0;
  auto ops_partials = make_partials_propagator(y_ref, a_ref, t0_ref, w_ref,
                                               v_ref, sv_ref, sw_ref, st0_ref);
  ret_t result = 0.0;

  // calculate density and partials
  for (size_t i = 0; i < N; i++) {
    if (sv_vec[i] == 0 && sw_vec[i] == 0 && st0_vec[i] == 0) {
      result += wiener_lccdf_unnorm<propto>(y_vec[i], a_vec[i], t0_vec[i],
                                            w_vec[i], v_vec[i],
                                            precision_derivatives);
      continue;
    }
    const T_partials_return y_value = y_vec.val(i);
    const T_partials_return a_value = a_vec.val(i);
    const T_partials_return v_value = v_vec.val(i);
    const T_partials_return w_value = w_vec.val(i);
    const T_partials_return t0_value = t0_vec.val(i);
    const T_partials_return sv_value = sv_vec.val(i);
    const T_partials_return sw_value = sw_vec.val(i);
    const T_partials_return st0_value = st0_vec.val(i);
    const int dim = (sv_value != 0) + (sw_value != 0) + (st0_value != 0);
    check_positive(function_name,
                   "(Inter-trial variability in drift rate) + "
                   "(Inter-trial variability in A-priori bias) + "
                   "(Inter-trial variability in nondecision time)",
                   dim);

    Eigen::Matrix<T_partials_return, -1, 1> xmin = Eigen::VectorXd::Zero(dim);
    Eigen::Matrix<T_partials_return, -1, 1> xmax = Eigen::VectorXd::Ones(dim);
    for (int i = 0; i < dim; i++) {
      xmin[i] = 0;
      xmax[i] = 1;
    }
    if (sv_value != 0) {
      xmin[0] = -1;
      xmax[0] = 1;
    }
    if (st0_value != 0) {
      xmax[dim - 1] = fmin(1.0, (y_value - t0_value) / st0_value);
    }

    T_partials_return hcubature_err
        = log_error_absolute - log_error_cdf + LOG_TWO + 1;

    const auto params = std::make_tuple(y_value, a_value, v_value, w_value,
                                        t0_value, sv_value, sw_value, st0_value,
                                        log_error_absolute - LOG_TWO);

    const T_partials_return ccdf
        = internal::wiener7_integrate_cdf<GradientCalc::OFF, GradientCalc::OFF,
                                          GradientCalc::OFF, GradientCalc::OFF,
                                          GradientCalc::OFF, GradientCalc::OFF>(
            [&](auto&&... args) { return internal::wiener4_ccdf(args...); },
            hcubature_err, params, dim, xmin, xmax,
            maximal_evaluations_hcubature, absolute_error_hcubature,
            relative_error_hcubature / 2);
    lccdf += log(ccdf);

    hcubature_err = log_error_absolute - log_error_derivatives + log(fabs(ccdf))
                    + LOG_TWO + 1;

    // computation of derivative for t and precision check in order to give
    // the value as deriv_t to edge1 and as -deriv_t to edge5

    // computation of derivatives and precision checks
    if (!is_constant_all<T_y>::value || !is_constant_all<T_t0>::value) {
      const T_partials_return deriv_t_7
          = -internal::wiener7_integrate_cdf<
                GradientCalc::OFF, GradientCalc::OFF, GradientCalc::OFF,
                GradientCalc::OFF, GradientCalc::ON>(
                [&](auto&&... args) {
                  return internal::wiener5_density<GradientCalc::ON>(args...);
                },
                hcubature_err, params, dim, xmin, xmax,
                maximal_evaluations_hcubature, absolute_error_hcubature,
                relative_error_hcubature / 2)
            / ccdf;
      if (!is_constant_all<T_y>::value) {
        partials<0>(ops_partials)[i] = deriv_t_7;
      }
      if (!is_constant_all<T_t0>::value) {
        partials<2>(ops_partials)[i] = -deriv_t_7;
      }
    }
    T_partials_return deriv;
    if (!is_constant_all<T_a>::value) {
      partials<1>(ops_partials)[i]
          = internal::wiener7_integrate_cdf(
                [&](auto&&... args) {
                  return internal::wiener4_ccdf_grad_a(args...);
                },
                hcubature_err, params, dim, xmin, xmax,
                maximal_evaluations_hcubature, absolute_error_hcubature,
                relative_error_hcubature / 2)
            / ccdf;
    }
    if (!is_constant_all<T_w>::value) {
      partials<3>(ops_partials)[i]
          = internal::wiener7_integrate_cdf<GradientCalc::OFF,
                                            GradientCalc::ON>(
                [&](auto&&... args) {
                  return internal::wiener4_ccdf_grad_w(args...);
                },
                hcubature_err, params, dim, xmin, xmax,
                maximal_evaluations_hcubature, absolute_error_hcubature,
                relative_error_hcubature / 2)
            / ccdf;
    }
    if (!is_constant_all<T_v>::value) {
      partials<4>(ops_partials)[i]
          = internal::wiener7_integrate_cdf(
                [&](auto&&... args) {
                  return internal::wiener4_ccdf_grad_v(args...);
                },
                hcubature_err, params, dim, xmin, xmax,
                maximal_evaluations_hcubature, absolute_error_hcubature,
                relative_error_hcubature / 2)
            / ccdf;
    }
    if (!is_constant_all<T_sv>::value) {
      if (sv_value == 0) {
        partials<5>(ops_partials)[i] = 0.0;
      } else {
        partials<5>(ops_partials)[i]
            = internal::wiener7_integrate_cdf<
                  GradientCalc::OFF, GradientCalc::OFF, GradientCalc::ON>(
                  [&](auto&&... args) {
                    return internal::wiener4_ccdf_grad_v(args...);
                  },
                  hcubature_err, params, dim, xmin, xmax,
                  maximal_evaluations_hcubature, absolute_error_hcubature,
                  relative_error_hcubature / 2)
              / ccdf;
      }
    }
    if (!is_constant_all<T_sw>::value) {
      if (sw_value == 0) {
        partials<6>(ops_partials)[i] = 0.0;
      } else {
        if (st0_value == 0 && sv_value == 0) {
          deriv = internal::estimate_with_err_check<5, 0, GradientCalc::OFF,
                                                    GradientCalc::ON>(
              [](auto&&... args) {
                return internal::wiener7_ccdf_grad_sw(args...);
              },
              hcubature_err, y_value - t0_value, a_value, v_value, w_value,
              sw_value, log_error_absolute - LOG_TWO);
          deriv = deriv / ccdf - 1 / sw_value;
        } else {
          deriv = internal::wiener7_integrate_cdf<
                      GradientCalc::OFF, GradientCalc::OFF, GradientCalc::OFF,
                      GradientCalc::ON>(
                      [&](auto&&... args) {
                        return internal::wiener4_ccdf_grad_w(args...);
                      },
                      hcubature_err, params, dim, xmin, xmax,
                      maximal_evaluations_hcubature, absolute_error_hcubature,
                      relative_error_hcubature / 2)
                  / ccdf;
        }
        partials<6>(ops_partials)[i] = deriv;
      }
    }
    if (!is_constant_all<T_st0>::value) {
      if (st0_value == 0) {
        partials<7>(ops_partials)[i] = 0.0;
      } else if (y_value - (t0_value + st0_value) <= 0) {
        partials<7>(ops_partials)[i] = -1 / st0_value;
      } else {
        const auto t0_st0 = t0_value + st0_value;
        if (sw_value == 0 && sv_value == 0) {
          deriv = internal::estimate_with_err_check<4, 0>(
              [](auto&&... args) { return internal::wiener4_ccdf(args...); },
              log_error_derivatives + log(st0_value), y_value - t0_st0, a_value,
              v_value, w_value, log_error_absolute - LOG_TWO);
          deriv = deriv / st0_value / ccdf - 1 / st0_value;
        } else {
          const int dim_st = (sv_value != 0) + (sw_value != 0);
          const auto new_error = log_error_absolute - LOG_TWO;
          const auto& params_st
              = std::make_tuple(y_value, a_value, v_value, w_value, t0_st0,
                                sv_value, sw_value, 0.0, new_error);
          deriv = internal::wiener7_integrate_cdf<
              GradientCalc::OFF, GradientCalc::OFF, GradientCalc::OFF,
              GradientCalc::OFF, GradientCalc::OFF, GradientCalc::OFF>(
              [&](auto&&... args) { return internal::wiener4_ccdf(args...); },
              hcubature_err, params_st, dim_st, xmin, xmax,
              maximal_evaluations_hcubature, absolute_error_hcubature,
              relative_error_hcubature / 2);
          deriv = deriv / st0_value / ccdf - 1 / st0_value;
        }
        partials<7>(ops_partials)[i] = deriv;
      }
    }
  }
  return result + ops_partials.build(lccdf);
}
}  // namespace math
}  // namespace stan
#endif
