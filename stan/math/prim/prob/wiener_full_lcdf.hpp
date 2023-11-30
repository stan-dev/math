#ifndef STAN_MATH_PRIM_PROB_WIENER_FULL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_WIENER_FULL_LCDF_HPP

#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/prob.hpp>
#include <stan/math/prim/functor/hcubature.hpp>

namespace stan {
namespace math {
namespace internal {

/**
 * Calculate the derivative of the wiener7 density w.r.t. 'sw'
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v The drift rate
 * @param w The relative starting point
 * @param sw The inter-trial variability of the relative starting point
 * @param log_error The log error tolerance
 * @return Gradient w.r.t. sw
 */
template <typename Scalar, typename ReturnT = return_type_t<Scalar>>
inline ReturnT wiener7_cdf_grad_sw(Scalar& y, Scalar& a, Scalar& v, Scalar& w,
                                   Scalar& sw, Scalar wildcard,
                                   Scalar log_error) {
  Scalar low = w - sw / 2.0;
  low = (0 > low) ? 0 : low;
  Scalar high = w + sw / 2.0;
  high = (1 < high) ? 1 : high;

  const Scalar lower_value
      = wiener4_distribution<true, Scalar>(y, a, v, low, 0, log_error);
  const Scalar upper_value
      = wiener4_distribution<true, Scalar>(y, a, v, high, 0, log_error);
  return 0.5 * (lower_value + upper_value) / sw;
}

/**
 * Helper function for agnostically calling wiener5 functions
 * (to be integrated over) or directly calling wiener7 functions,
 * accounting for the different number of arguments.
 *
 * Specialisation for wiener5 functions
 *
 * @tparam GradSW Whether the wiener7 gradient function is passed
 * @tparam F Type of Gradient/density functor
 *
 * @param functor Gradient/density functor to apply
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v The drift rate
 * @param w The relative starting point
 * @param t0 The non-decision time
 * @param sv The inter-trial variability of the drift rate
 * @param sw The inter-trial variability of the relative starting point
 * @param cdf The value of the distribution
 * @param log_error The log error tolerance
 * @return Functor applied to arguments
 */
template <typename Scalar, bool GradSW, typename F,
          std::enable_if_t<!GradSW>* = nullptr>
inline Scalar conditionally_grad_sw_cdf(const F& functor, Scalar y_diff,
                                        Scalar a, Scalar v, Scalar w, Scalar sv,
                                        Scalar sw, Scalar cdf,
                                        Scalar log_error) {
  return functor(y_diff, a, v, w, cdf, log_error);
}

/**
 * Helper function for agnostically calling wiener5 functions
 * (to be integrated over) or directly calling wiener7 functions,
 * accounting for the different number of arguments.
 *
 * Specialisation for wiener7 functions
 *
 * @tparam GradSW Whether the wiener7 gradient function is passed
 * @tparam F Type of Gradient/density functor
 *
 * @param functor Gradient/density functor to apply
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v The drift rate
 * @param w The relative starting point
 * @param t0 The non-decision time
 * @param sv The inter-trial variability of the drift rate
 * @param sw The inter-trial variability of the relative starting point
 * @param cdf The value of the distribution
 * @param log_error The log error tolerance
 * @return Functor applied to arguments
 */
template <typename Scalar, bool GradSW, typename F,
          std::enable_if_t<GradSW>* = nullptr>
inline Scalar conditionally_grad_sw_cdf(const F& functor, Scalar y_diff,
                                        Scalar a, Scalar v, Scalar w, Scalar sv,
                                        Scalar sw, Scalar cdf,
                                        Scalar log_error) {
  return functor(y_diff, a, v, w, sw, cdf, log_error);
}

/**
 * Implementation function for preparing arguments and functor to be passed
 * to the hcubature() function for calculating wiener7 parameters via
 * integration
 *
 * @tparam GradSW Whether the wiener7 gradient function is passed
 * @tparam Wiener7FunctorT Type of functor
 * @tparam Targs... Types of arguments in parameter pack
 *
 * @param wiener7_functor Gradient/density functor to apply
 * @param hcubature_err Error tolerance for calculation
 * @param args Additional arguments to be passed to the hcubature function
 * @return Wiener7 density or gradient calculated by integration
 */
template <typename Scalar, bool Distribution = false, bool GradW7 = false,
          bool GradSV = false, bool GradSW = false, bool GradSW_ord = false,
          bool GradT = false, typename Wiener7FunctorT, typename... TArgs,
          typename ReturnT = return_type_t<Scalar>>
ReturnT wiener7_integrate_cdf(const Wiener7FunctorT& wiener7_functor,
                              Scalar hcubature_err, TArgs&&... args) {
  const auto wiener7_integrand_impl = [&](auto&& x, Scalar y, Scalar a,
                                          Scalar v, Scalar w, Scalar t0,
                                          Scalar sv, Scalar sw, Scalar st0,
                                          Scalar lerr) {
    scalar_seq_view<decltype(x)> x_vec(x);
    const Scalar temp = (sv != 0) ? square(x[0]) : 0;
    const Scalar factor = (sv != 0) ? x[0] / (1 - temp) : 0;
    const Scalar new_v = (sv != 0) ? v + sv * factor : v;
    const Scalar sv_val = (GradSW) ? 0 : sv;
    const Scalar sw_val = (GradSW) ? 0 : sw;
    const Scalar new_w = (sv_val != 0)
                             ? ((sw != 0) ? w + sw * (x[1] - 0.5) : w)
                             : ((sw_val != 0) ? w + sw * (x[0] - 0.5) : w);
    const Scalar new_t0
        = (sv_val != 0) ? ((sw != 0) ? ((st0 != 0) ? t0 + st0 * x[2] : t0)
                                     : ((st0 != 0) ? t0 + st0 * x[1] : t0))
                        : ((sw_val != 0) ? ((st0 != 0) ? t0 + st0 * x[1] : t0)
                                         : ((st0 != 0) ? t0 + st0 * x[0] : t0));
    if (y - new_t0 <= 0) {
      return static_cast<Scalar>(0.0);
    } else {
      const auto dist = GradT ? 0
                              : wiener4_distribution<true, Scalar>(
                                  y - new_t0, a, new_v, new_w, 0, lerr);
      const Scalar temp2 = (sv != 0) ? -0.5 * square(factor) - LOG_SQRT_PI
                                           - 0.5 * LOG_TWO + log1p(temp)
                                           - 2 * log1p(-temp)
                                     : 0;
      const Scalar factor_sv = GradSV ? factor : 1;
      const Scalar factor_sw
          = GradSW_ord ? ((sv_val != 0) ? (x[1] - 0.5) : (x[0] - 0.5)) : 1;
      const auto integrand
          = Distribution
                ? dist
                : GradT ? internal::conditionally_grad_sw<Scalar, GradSW>(
                      wiener7_functor, y - new_t0, a, v, new_w, sv, sw, lerr)
                        : factor_sv * factor_sw
                              * conditionally_grad_sw_cdf<Scalar, GradSW>(
                                  wiener7_functor, y - new_t0, a, new_v, new_w,
                                  sv, sw, dist, lerr);
      return integrand * exp(temp2);
    }
  };
  const auto functor = [&](auto&&... int_args) {
    return hcubature(wiener7_integrand_impl, int_args...);
  };
  return estimate_with_err_check<Scalar, 0, GradW7, 8>(functor, hcubature_err,
                                                       args...);
}
}  // namespace internal

/** \ingroup prob_dists
 * Returns the log CDF of the Wiener distribution for a
 * (Wiener) drift diffusion model with up to 7 parameters. If containers
 * are supplied, returns the log sum of the probabilities.
 * See 'wiener_full_lpdf' for more comprehensive documentation
 * If the reaction time goes to infinity, the CDF goes to the probability to
 * hit the upper bound (instead of 1, as it is usually the case)
 *
 * @tparam T_y type of scalar
 * @tparam T_a type of boundary separation
 * @tparam T_t0 type of non-decision time
 * @tparam T_w type of relative starting point
 * @tparam T_v type of drift rate
 * @tparam T_sv type of inter-trial variability of drift rate
 * @tparam T_sw type of inter-trial variability of relative starting point
 * @tparam T_st0 type of inter-trial variability of non-decision time
 *
 * @param y A scalar variable; the reaction time in seconds
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
          typename T_st0,
          typename ReturnT
          = return_type_t<T_y, T_a, T_t0, T_w, T_v, T_sv, T_sw, T_st0>>
inline ReturnT wiener_full_lcdf(const T_y& y, const T_a& a, const T_t0& t0,
                                const T_w& w, const T_v& v, const T_sv& sv,
                                const T_sw& sw, const T_st0& st0,
                                const double& precision_derivatives = 1e-8) {
  if (!include_summand<propto, T_y, T_a, T_v, T_w, T_t0, T_sv, T_sw,
                       T_st0>::value) {
    return 0;
  }
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_a_ref = ref_type_if_t<!is_constant<T_a>::value, T_a>;
  using T_v_ref = ref_type_if_t<!is_constant<T_v>::value, T_v>;
  using T_w_ref = ref_type_if_t<!is_constant<T_w>::value, T_w>;
  using T_t0_ref = ref_type_if_t<!is_constant<T_t0>::value, T_t0>;
  using T_sv_ref = ref_type_if_t<!is_constant<T_sv>::value, T_sv>;
  using T_sw_ref = ref_type_if_t<!is_constant<T_sw>::value, T_sw>;
  using T_st0_ref = ref_type_if_t<!is_constant<T_st0>::value, T_st0>;

  using T_partials_return
      = partials_return_t<T_y, T_a, T_t0, T_w, T_v, T_sv, T_sw, T_st0>;

  static constexpr const char* function_name = "wiener_full_lcdf";
  check_consistent_sizes(function_name, "Random variable", y,
                         "Boundary separation", a, "Drift rate", v,
                         "A-priori bias", w, "Nondecision time", t0,
                         "Inter-trial variability in drift rate", sv,
                         "Inter-trial variability in A-priori bias", sw,
                         "Inter-trial variability in Nondecision time", st0);

  T_y_ref y_ref = y;
  T_a_ref a_ref = a;
  T_v_ref v_ref = v;
  T_w_ref w_ref = w;
  T_t0_ref t0_ref = t0;
  T_sv_ref sv_ref = sv;
  T_sw_ref sw_ref = sw;
  T_st0_ref st0_ref = st0;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) a_val = to_ref(as_value_column_array_or_scalar(a_ref));
  decltype(auto) v_val = to_ref(as_value_column_array_or_scalar(v_ref));
  decltype(auto) w_val = to_ref(as_value_column_array_or_scalar(w_ref));
  decltype(auto) t0_val = to_ref(as_value_column_array_or_scalar(t0_ref));
  decltype(auto) sv_val = to_ref(as_value_column_array_or_scalar(sv_ref));
  decltype(auto) sw_val = to_ref(as_value_column_array_or_scalar(sw_ref));
  decltype(auto) st0_val = to_ref(as_value_column_array_or_scalar(st0_ref));
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

  if (size_zero(y, a, v, w, t0) || size_zero(sv, sw, st0)) {
    return 0;
  }
  const size_t N = max_size(y, a, v, w, t0, sv, sw, st0);
  if (!N) {
    return 0;
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
      msg << ", but must be smaller than 2*(A-priori bias) = " << 2 * w_vec[i];
      std::string msg_str(msg.str());
      throw_domain_error(function_name,
                         "Inter-trial variability in A-priori bias", sw_vec[i],
                         " = ", msg_str.c_str());
    }
    if (unlikely(w_vec[i] + .5 * sw_vec[i] >= 1)) {
      std::stringstream msg;
      msg << ", but must be smaller than 2*(1-A-priori bias) = "
          << 2 * (1 - w_vec[i]);
      std::string msg_str(msg.str());
      throw_domain_error(function_name,
                         "Inter-trial variability in A-priori bias", sw_vec[i],
                         " = ", msg_str.c_str());
    }
  }

  const T_partials_return log_error_cdf = log(1e-6);  // precision for density
  const T_partials_return error_bound = precision_derivatives;  // precision for
  // derivatives (controllable by user)
  const T_partials_return lerror_bound = log(error_bound);
  const T_partials_return absolute_error_hcubature = 0.0;
  const T_partials_return relative_error_hcubature
      = .9 * error_bound;  // eps_rel(Integration)
  const T_partials_return log_error_absolute = log(1e-12);
  const int maximal_evaluations_hcubature = 6000;
  T_partials_return lcdf = 0.0;
  auto ops_partials = make_partials_propagator(y_ref, a_ref, t0_ref, w_ref,
                                               v_ref, sv_ref, sw_ref, st0_ref);
  ReturnT result = 0;

  // calculate density and partials
  for (size_t i = 0; i < N; i++) {
    if (sv_vec[i] == 0 && sw_vec[i] == 0 && st0_vec[i] == 0) {
      result += wiener4_lcdf<propto>(y_vec[i], a_vec[i], t0_vec[i], w_vec[i],
                                     v_vec[i], precision_derivatives);
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

    const T_partials_return cdf
        = internal::wiener7_integrate_cdf<T_partials_return, true>(
            [&](auto&&... args) {
              return internal::wiener4_distribution<true, T_partials_return>(
                  args...);
            },
            hcubature_err, params, dim, xmin, xmax,
            maximal_evaluations_hcubature, absolute_error_hcubature,
            relative_error_hcubature / 2);
    lcdf += log(cdf);

    hcubature_err
        = log_error_absolute - lerror_bound + log(fabs(cdf)) + LOG_TWO + 1;

    // computation of derivative for t and precision check in order to give
    // the value as deriv_t to edge1 and as -deriv_t to edge5
    const auto params_dt7 = std::make_tuple(
        y_value, a_value, v_value, w_value, t0_value, sv_value, sw_value,
        st0_value, log_error_absolute - LOG_TWO - 9 * LOG_TWO);
    T_partials_return deriv_t_7
        = internal::wiener7_integrate_cdf<T_partials_return, false, false,
                                          false, false, false, true>(
              [&](auto&&... args) {
                return internal::wiener5_density<true>(args...);
              },
              hcubature_err, params, dim, xmin, xmax,
              maximal_evaluations_hcubature, absolute_error_hcubature,
              relative_error_hcubature / 2)
          / cdf;

    // computation of derivatives and precision checks
    T_partials_return deriv;
    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials)[i] = deriv_t_7;
    }
    if (!is_constant_all<T_a>::value) {
      partials<1>(ops_partials)[i]
          = internal::wiener7_integrate_cdf<T_partials_return>(
                [&](auto&&... args) {
                  return internal::wiener4_cdf_grad_a<T_partials_return>(
                      args...);
                },
                hcubature_err, params, dim, xmin, xmax,
                maximal_evaluations_hcubature, absolute_error_hcubature,
                relative_error_hcubature / 2)
            / cdf;
    }
    if (!is_constant_all<T_t0>::value) {
      partials<2>(ops_partials)[i] = -deriv_t_7;
    }
    if (!is_constant_all<T_w>::value) {
      partials<3>(ops_partials)[i]
          = internal::wiener7_integrate_cdf<T_partials_return, false, true>(
                [&](auto&&... args) {
                  return internal::wiener4_cdf_grad_w<T_partials_return>(
                      args...);
                },
                hcubature_err, params, dim, xmin, xmax,
                maximal_evaluations_hcubature, absolute_error_hcubature,
                relative_error_hcubature / 2)
            / cdf;
    }
    if (!is_constant_all<T_v>::value) {
      partials<4>(ops_partials)[i]
          = internal::wiener7_integrate_cdf<T_partials_return>(
                [&](auto&&... args) {
                  return internal::wiener4_cdf_grad_v<T_partials_return>(
                      args...);
                },
                hcubature_err, params, dim, xmin, xmax,
                maximal_evaluations_hcubature, absolute_error_hcubature,
                relative_error_hcubature / 2)
            / cdf;
    }
    if (!is_constant_all<T_sv>::value) {
      if (sv_value == 0) {
        partials<5>(ops_partials)[i] = 0;
      } else {
        partials<5>(ops_partials)[i]
            = internal::wiener7_integrate_cdf<T_partials_return, false, false,
                                              true>(
                  [&](auto&&... args) {
                    return internal::wiener4_cdf_grad_v<T_partials_return>(
                        args...);
                  },
                  hcubature_err, params, dim, xmin, xmax,
                  maximal_evaluations_hcubature, absolute_error_hcubature,
                  relative_error_hcubature / 2)
              / cdf;
      }
    }
    if (!is_constant_all<T_sw>::value) {
      if (sw_value == 0) {
        partials<6>(ops_partials)[i] = 0;
      } else {
        if (st0_value == 0 && sv_value == 0) {
          deriv = internal::estimate_with_err_check<T_partials_return, 5>(
              [&](auto&&... args) {
                return internal::conditionally_grad_sw_cdf<T_partials_return,
                                                           true>(
                    internal::wiener7_cdf_grad_sw<T_partials_return>, args...);
              },
              hcubature_err, y_value - t0_value, a_value, v_value, w_value,
              sv_value, sw_value, 0,
              log_error_absolute
                  - LOG_TWO);  // added wildcard 0. delete later somehow
          deriv = deriv / cdf - 1 / sw_value;
        } else {
          deriv = internal::wiener7_integrate_cdf<T_partials_return, false,
                                                  false, false, false, true>(
                      [&](auto&&... args) {
                        return internal::wiener4_cdf_grad_w<T_partials_return>(
                            args...);
                      },
                      hcubature_err, params, dim, xmin, xmax,
                      maximal_evaluations_hcubature, absolute_error_hcubature,
                      relative_error_hcubature / 2)
                  / cdf;
        }
        partials<6>(ops_partials)[i] = deriv;
      }
    }
    if (!is_constant_all<T_st0>::value) {
      if (st0_value == 0) {
        partials<7>(ops_partials)[i] = 0;
      } else if (y_value - (t0_value + st0_value) <= 0) {
        partials<7>(ops_partials)[i] = -1 / st0_value;
      } else {
        const T_partials_return t0_st0 = t0_value + st0_value;
        if (sw_value == 0 && sv_value == 0) {
          deriv = internal::estimate_with_err_check<T_partials_return, 5>(
              [&](auto&&... args) {
                return internal::wiener4_distribution<true, T_partials_return>(
                    args...);
              },
              lerror_bound + log(st0_value), y_value - t0_st0, a_value, v_value,
              w_value, 0,
              log_error_absolute
                  - LOG_TWO);  // added wildcard 0. delete later somehow
          deriv = deriv / st0_value / cdf - 1 / st0_value;
        } else {
          const int dim_st = (sv_value != 0) + (sw_value != 0);
          const T_partials_return new_error = log_error_absolute - LOG_TWO;
          const auto& params_st = std::make_tuple(
              y_value, a_value, v_value, w_value, t0_st0, sv_value, sw_value, 0,
              new_error);  // added wildcard 0. delete later somehow
          deriv = internal::wiener7_integrate_cdf<T_partials_return>(
              [&](auto&&... args) {
                return internal::wiener4_distribution<true, T_partials_return>(
                    args...);
              },
              hcubature_err, params_st, dim_st, xmin, xmax,
              maximal_evaluations_hcubature, absolute_error_hcubature,
              relative_error_hcubature / 2);
          deriv = deriv / st0_value / cdf - 1 / st0_value;
        }
        partials<7>(ops_partials)[i] = deriv;
      }
    }
  }
  return result + ops_partials.build(lcdf);
}
}  // namespace math
}  // namespace stan
#endif
