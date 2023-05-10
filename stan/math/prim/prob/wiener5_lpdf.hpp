#ifndef STAN_MATH_PRIM_PROB_WIENER5_LPDF_HPP
#define STAN_MATH_PRIM_PROB_WIENER5_LPDF_HPP

#include <stan/math/prim/fun.hpp>

namespace stan {
namespace math {
namespace internal {

enum class FunType { Density, GradT, GradA, GradV, GradW, GradSV, GradSW };

/**
 * Implementation function for calculating the density or gradient w.r.t. a
 * specified parameter of the 5-parameter Wiener distribution
 *
 * @tparam FunTypeEnum An enum value indicating the quantity required
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param vn The drift rate
 * @param wn The relative starting point
 * @param sv The inter-trial variability of the drift rate
 * @param err The error tolerance for estimation
 * @param use_log Whether to return calculation w.r.t. the log value
 *
 */
template <FunType FunTypeEnum>
inline double wiener5_helper(const double& y, const double& a, const double& vn,
                             const double& wn, const double& sv,
                             const double& err = log(1e-12),
                             bool use_log = false) {
  double w = 1.0 - wn;
  double v = -vn;
  double sv_sqr = square(sv);
  double one_plus_svsqr_y = 1 + sv_sqr * y;
  double two_avw = 2 * a * v * w;
  double two_log_a = 2 * log(a);
  double ans;
  double ld = 0;

  if (FunTypeEnum == FunType::GradSV || FunTypeEnum == FunType::GradV) {
    ld = wiener5_helper<FunType::Density>(y, a, vn, wn, sv, err);
    if (FunTypeEnum == FunType::GradV) {
      ans = (a * (1 - wn) - vn * y);
      if (sv != 0) {
        ans /= 1 + sv_sqr * y;
      }
    } else {
      double t1 = -y / one_plus_svsqr_y;
      double t2 = (square(a * w) + two_avw * y + square(v * y))
                  / square(one_plus_svsqr_y);
      ans = sv * (t1 + t2);
    }
  } else {
    double lg1;
    if (sv != 0) {
      lg1 = (sv_sqr * square(a * w) - two_avw - square(v) * y) / 2.0
                / one_plus_svsqr_y
            - two_log_a - 0.5 * log(one_plus_svsqr_y);
    } else {
      lg1 = (-two_avw - square(v) * y) / 2.0 - two_log_a;
    }
    double ans0 = 0;
    if (FunTypeEnum != FunType::Density) {
      double var_a = (FunTypeEnum == FunType::GradA) ? w : a;
      double var_b = (FunTypeEnum == FunType::GradA) ? a : w;
      if (sv != 0) {
        if (FunTypeEnum == FunType::GradT) {
          ans0 = -0.5
                 * (square(sv_sqr) * (y + square(a * w))
                    + sv_sqr * (1 - two_avw) + square(v))
                 / square(one_plus_svsqr_y);
        } else {
          ans0 = (-v * var_a + sv_sqr * square(var_a) * var_b)
                 / one_plus_svsqr_y;
        }
      } else {
        ans0 = (FunTypeEnum == FunType::GradT) ? -0.5 * square(v) : -v * var_a;
      }
    }

    double es = (err - lg1) + (FunTypeEnum == FunType::GradT) * two_log_a;
    double two_es = 2.0 * es;
    double y_asq = y / square(a);
    double log_y_asq = log(y) - two_log_a;

    double K1, u_eps, arg;
    double K1_mult = (FunTypeEnum == FunType::Density) ? 2 : 3;
    K1 = (sqrt(K1_mult * y_asq) + w) / 2.0;
    if (FunTypeEnum == FunType::Density || FunTypeEnum == FunType::GradW) {
      u_eps = fmin(-1.0, LOG_TWO + LOG_PI + 2.0 * log_y_asq + two_es);
    } else {
      u_eps = fmin(-3.0,
                   (log(8.0) - log(27.0) + LOG_PI + 4.0 * log_y_asq + two_es));
    }
    double arg_mult = (FunTypeEnum == FunType::Density) ? 1 : 3;
    arg = -arg_mult * y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));

    double K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
    size_t kss = static_cast<size_t>(ceil(fmax(K1, K2)));

    static const double PISQ = square(pi());  // pi*pi
    if (FunTypeEnum == FunType::Density) {
      K1 = 1.0 / (pi() * sqrt(y_asq));
      double two_log_piy = -2.0 * (LOG_PI + log_y_asq + es);
      K2 = (two_log_piy >= 0) ? sqrt(two_log_piy / (PISQ * y_asq)) : 0.0;
    } else {
      K1_mult = (FunTypeEnum == FunType::GradT) ? 2 : 3;
      K1 = sqrt(K1_mult / y_asq) / pi();
      double u_eps_arg
          = (FunTypeEnum == FunType::GradT)
                ? log(3.0) - log(5.0) + LOG_PI + 2.0 * log_y_asq + es
                : log(4.0) - log(9.0) + 2.0 * LOG_PI + 3.0 * log_y_asq + two_es;
      u_eps = fmin(-1, u_eps_arg);
      arg_mult = (FunTypeEnum == FunType::GradT) ? (2.0 / PISQ / y_asq) : 1;
      arg = -arg_mult * (u_eps - sqrt(-2.0 * u_eps - 2.0));
      K2 = (FunTypeEnum == FunType::GradT)
               ? (arg > 0) ? sqrt(arg) : K1
               : (arg > 0) ? sqrt(arg / y_asq) / pi() : K1;
    }

    size_t kll = static_cast<size_t>(ceil(fmax(K1, K2)));

    double ld_err = 0;
    if (FunTypeEnum != FunType::Density) {
      if (FunTypeEnum == FunType::GradT) {
        ld_err = log(max(fabs(ans0 - 1.5 / y), fabs(ans0)));
      } else if (FunTypeEnum == FunType::GradA) {
        ld_err = log(max(fabs(ans0 + 1.0 / a), fabs(ans0 - 2.0 / a)));
      } else {
        ld_err = log(fabs(ans0));
      }
      ld = wiener5_helper<FunType::Density>(y, a, vn, wn, sv, err - ld_err);
    }

    double erg = NEGATIVE_INFTY;
    int newsign = 1;
    double scaling = (2 * kss <= kll) ? inv(2.0 * y_asq) : y_asq / 2.0;
    if (2 * kss <= kll) {
      double mult = (FunTypeEnum == FunType::Density) ? 1 : 3;
      double offset = (FunTypeEnum == FunType::GradW) ? y_asq : 0;
      double sqrt_offset = sqrt(offset);
      for (size_t k = kss; k >= std::fmax(0, kss - 1e3); k--) {
        double wp2k = w + 2.0 * k;
        double wm2k = w - 2.0 * k;
        int wp2k_sign = (wp2k > sqrt_offset) ? 1 : -1;
        int wm2k_sign = (wm2k > sqrt_offset) ? 1 : -1;
        double wp2k_quant = log(wp2k_sign * (wp2k - offset))
                            - (square(wp2k) - offset) * scaling;
        double wm2k_quant = log(wm2k_sign * (wm2k - offset))
                            - (square(wm2k) - offset) * scaling;
        std::forward_as_tuple(erg, newsign)
            = log_sum_exp_signed(erg, newsign, mult * wp2k_quant, wp2k_sign);
        std::forward_as_tuple(erg, newsign) = log_sum_exp_signed(
            erg, newsign, mult * wm2k_quant, -1 * wm2k_sign);
      }
      if (FunTypeEnum == FunType::Density) {
        ans = lg1 - 0.5 * LOG_TWO - LOG_SQRT_PI - 1.5 * log_y_asq + erg;
      } else if (FunTypeEnum == FunType::GradT) {
        ans = ans0 - 1.5 / y
              + newsign
                    * exp(lg1 - two_log_a - 1.5 * LOG_TWO - LOG_SQRT_PI
                          - 3.5 * log_y_asq + erg - ld);
      } else if (FunTypeEnum == FunType::GradA) {
        ans = ans0 + 1.0 / a
              - newsign
                    * exp(-0.5 * LOG_TWO - LOG_SQRT_PI - 2.5 * log(y)
                          + 2.0 * two_log_a + lg1 + erg - ld);
      } else if (FunTypeEnum == FunType::GradW) {
        ans = -(ans0
                - newsign
                      * exp(erg - (ld - lg1) - 2.5 * log_y_asq - 0.5 * LOG_TWO
                            - 0.5 * LOG_PI));
      }
    } else {
      double mult = 3;
      if (FunTypeEnum == FunType::Density) {
        mult = 1;
      } else if (FunTypeEnum == FunType::GradW) {
        mult = 2;
      }
      for (size_t k = kll; k >= std::fmax(1, kll - 1e3); k--) {
        double pi_k = k * pi();
        double check
            = (FunTypeEnum == FunType::GradW) ? cos(pi_k * w) : sin(pi_k * w);
        int check_sign = sign(check);
        double kll_quant
            = mult * log(k) - square(pi_k) * scaling + log(fabs(check));
        std::forward_as_tuple(erg, newsign)
            = log_sum_exp_signed(erg, newsign, kll_quant, check_sign);
      }
      if (FunTypeEnum == FunType::Density) {
        ans = lg1 + erg + LOG_PI;
      } else if (FunTypeEnum == FunType::GradT) {
        ans = ans0
              - newsign
                    * exp(lg1 - two_log_a + 3.0 * LOG_PI - LOG_TWO + erg - ld);
      } else if (FunTypeEnum == FunType::GradA) {
        ans = ans0 - 2.0 / a
              + newsign * exp(log(y) + lg1 - 3 * (log(a) - LOG_PI) + erg - ld);
      } else if (FunTypeEnum == FunType::GradW) {
        ans = -(ans0 + newsign * exp(erg - (ld - lg1) + 2 * LOG_PI));
      }
    }
  }
  if (FunTypeEnum == FunType::Density && use_log) {
    return exp(ans);
  }
  return use_log ? ans * exp(ld) : ans;
}

/**
 * Utility function for replacing a value with a specified error value
 */
template <size_t NestedIndex>
inline void assign_err(double arg, double err) {
  arg = err;
}

/**
 * Utility function for replacing a value with a specified error value,
 * overload for use when the value is stored within a tuple.
 */
template <size_t NestedIndex, typename... TArgs>
inline void assign_err(std::tuple<TArgs...>& args_tuple, double err) {
  std::get<NestedIndex>(args_tuple) = err;
}

/**
 * Utility function for estimating a function with a given set of arguments,
 * checking the result against a provided error tolerance, and re-estimating
 * the function with the increased error if it fails.
 *
 * @tparam ErrIndex Position of the error argument in the provided arguments
 * @tparam NestedIndex Nested position if the error argument is within a tuple
 * @tparam F Type of functor
 * @tparam ArgsTupleT Type of tuple of arguments for functor
 *
 * @param functor Function to apply
 * @param err Error value to check against
 * @param args_tuple Tuple of arguments to pass to functor
 * @param log_result Whether the function result is already on the log-scale
 */
template <size_t ErrIndex, size_t NestedIndex = 0, typename F,
          typename ArgsTupleT>
double estimate_with_err_check(const F& functor, double err,
                               ArgsTupleT&& args_tuple,
                               bool log_result = true) {
  double result = math::apply([&](auto&&... args) { return functor(args...); },
                              args_tuple);
  double lfabs_result = log_result ? log(fabs(result)) : fabs(result);
  if (lfabs_result < err) {
    ArgsTupleT err_args_tuple = args_tuple;
    assign_err<NestedIndex>(std::get<ErrIndex>(err_args_tuple),
                            err + lfabs_result);
    result = math::apply([&](auto&&... args) { return functor(args...); },
                         err_args_tuple);
  }
  return result;
}
}  // namespace internal

/**
 * Log-density function for the 5-parameter Wiener density.
 * See 'wiener_full_lpdf' for more comprehensive documentation
 *
 * @tparam T_y type of scalar
 * @tparam T_a type of boundary
 * @tparam T_t0 type of non-decision time
 * @tparam T_w type of relative starting point
 * @tparam T_v type of drift rate
 * @tparam T_sv type of inter-trial variability of drift rate
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param t0 The non-decision time
 * @param w The relative starting point
 * @param v The drift rate
 * @param sv The inter-trial variability of the drift rate
 * @return The log of the Wiener first passage time density with
 *  the specified arguments for upper boundary responses
 */
template <bool propto, typename T_y, typename T_a, typename T_t0, typename T_w,
          typename T_v, typename T_sv>
inline return_type_t<T_y, T_a, T_t0, T_w, T_v, T_sv> wiener5_lpdf(
    const T_y& y, const T_a& a, const T_t0& t0, const T_w& w, const T_v& v,
    const T_sv& sv, const double& prec) {
  using T_y_ref = ref_type_t<T_y>;
  using T_a_ref = ref_type_t<T_a>;
  using T_t0_ref = ref_type_t<T_t0>;
  using T_w_ref = ref_type_t<T_w>;
  using T_v_ref = ref_type_t<T_v>;
  using T_sv_ref = ref_type_t<T_sv>;

  const char* function_name = "wiener5_lpdf";
  check_consistent_sizes(function_name, "Random variable", y,
                         "Boundary separation", a, "Drift rate", v,
                         "A-priori bias", w, "Nondecision time", t0,
                         "Inter-trial variability in drift rate", sv);
  check_consistent_size(function_name, "Random variable", y, 1);
  check_consistent_size(function_name, "Boundary separation", a, 1);
  check_consistent_size(function_name, "Nondecision time", t0, 1);
  check_consistent_size(function_name, "A-priori bias", w, 1);
  check_consistent_size(function_name, "Drift rate", v, 1);
  check_consistent_size(function_name, "Inter-trial variability in drift rate",
                        sv, 1);

  T_y_ref y_ref = y;
  T_a_ref a_ref = a;
  T_t0_ref t0_ref = t0;
  T_w_ref w_ref = w;
  T_v_ref v_ref = v;
  T_sv_ref sv_ref = sv;

  check_positive_finite(function_name, "Random variable", value_of(y_ref));
  check_positive_finite(function_name, "Boundary separation", value_of(a_ref));
  check_nonnegative(function_name, "Nondecision time", value_of(t0_ref));
  check_finite(function_name, "Nondecision time", value_of(t0_ref));
  check_less(function_name, "A-priori bias", value_of(w_ref), 1);
  check_greater(function_name, "A-priori bias", value_of(w_ref), 0);
  check_finite(function_name, "Drift rate", value_of(v_ref));
  check_nonnegative(function_name, "Inter-trial variability in drift rate",
                    value_of(sv_ref));
  check_finite(function_name, "Inter-trial variability in drift rate",
               value_of(sv_ref));

  if (size_zero(y, a, t0, w, v) || size_zero(sv)) {
    return 0;
  }
  size_t N = max_size(y, a, t0, w, v, sv);
  if (!N) {
    return 0;
  }
  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_a_ref> a_vec(a_ref);
  scalar_seq_view<T_t0_ref> t0_vec(t0_ref);
  scalar_seq_view<T_w_ref> w_vec(w_ref);
  scalar_seq_view<T_v_ref> v_vec(v_ref);
  scalar_seq_view<T_sv_ref> sv_vec(sv_ref);
  size_t N_y_t0 = max_size(y, t0);

  for (size_t i = 0; i < N_y_t0; ++i) {
    if (y_vec[i] <= t0_vec[i]) {
      std::stringstream msg;
      msg << ", but must be greater than nondecision time = " << t0_vec[i];
      std::string msg_str(msg.str());
      throw_domain_error(function_name, "Random variable", y_vec[i], " = ",
                         msg_str.c_str());
    }
  }

  if (!include_summand<propto, T_y, T_a, T_t0, T_w, T_v, T_sv>::value) {
    return 0;
  }

  double lerror_bound_dens = log(1e-6);  // precision for density
  double lerror_bound = log(prec);       // precision for derivatives
  double labstol_wiener5 = log(1e-12);   // eps_abs(wiener5)
  double dens = 0.0;
  double ld = 0.0;
  operands_and_partials<T_y_ref, T_a_ref, T_t0_ref, T_w_ref, T_v_ref, T_sv_ref>
      ops_partials(y_ref, a_ref, t0_ref, w_ref, v_ref, sv_ref);

  static constexpr double LOG_FOUR = LOG_TWO + LOG_TWO;

  // calculate density and partials
  for (size_t i = 0; i < N; i++) {
    // Calculate 4-parameter model without inter-trial variabilities (if
    // sv_vec[i] == 0) or 5-parameter model with inter-trial variability in
    // drift rate (if sv_vec[i] != 0)
    const double y_val = y_vec.val(i);
    const double a_val = a_vec.val(i);
    const double t0_val = t0_vec.val(i);
    const double w_val = w_vec.val(i);
    const double v_val = v_vec.val(i);
    const double sv_val = sv_vec.val(i);

    const auto params = std::make_tuple(y_val - t0_val, a_val, v_val, w_val,
                                        sv_val, labstol_wiener5, false);

    dens = internal::estimate_with_err_check<5>(
        internal::wiener5_helper<internal::FunType::Density>,
        lerror_bound_dens - LOG_TWO, params, true);
    ld += dens;

    double new_est_err = dens + lerror_bound - LOG_FOUR;

    // computation of derivative for t and precision check in order to give
    // the value as deriv_y to edge1 and as -deriv_y to edge5
    double deriv_y = internal::estimate_with_err_check<5>(
        internal::wiener5_helper<internal::FunType::GradT>, new_est_err,
        params);

    // computation of derivatives and precision checks
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[i] = deriv_y;
    }
    if (!is_constant_all<T_a>::value) {
      ops_partials.edge2_.partials_[i] = internal::estimate_with_err_check<5>(
          internal::wiener5_helper<internal::FunType::GradA>, new_est_err,
          params);
    }
    if (!is_constant_all<T_t0>::value) {
      ops_partials.edge3_.partials_[i] = -deriv_y;
    }
    if (!is_constant_all<T_w>::value) {
      ops_partials.edge4_.partials_[i] = internal::estimate_with_err_check<5>(
          internal::wiener5_helper<internal::FunType::GradW>, new_est_err,
          params);
    }
    if (!is_constant_all<T_v>::value) {
      ops_partials.edge5_.partials_[i]
          = internal::wiener5_helper<internal::FunType::GradV>(
              y_val - t0_val, a_val, v_val, w_val, sv_val);
    }
    if (!is_constant_all<T_sv>::value) {
      ops_partials.edge6_.partials_[i]
          = internal::wiener5_helper<internal::FunType::GradSV>(
              y_val - t0_val, a_val, v_val, w_val, sv_val);
    }
  }  // end for loop
  return ops_partials.build(ld);
}  // end wiener5_lpdf
}  // namespace math
}  // namespace stan
#endif
