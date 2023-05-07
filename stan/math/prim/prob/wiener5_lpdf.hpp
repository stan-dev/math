#ifndef STAN_MATH_PRIM_PROB_WIENER5_LPDF_HPP
#define STAN_MATH_PRIM_PROB_WIENER5_LPDF_HPP

#include <stan/math/prim/fun.hpp>

namespace stan {
namespace math {
namespace internal {

template <typename... TArgs>
inline void assign_err(std::tuple<TArgs...>& args_tuple, double err) {
  std::get<8>(args_tuple) = err;
}
inline void assign_err(double arg, double err) { arg = err; }

template <size_t ErrIndex, typename F, typename ArgsTupleT>
double estimate_with_err_check(const F& functor, double err,
                               ArgsTupleT&& args_tuple,
                               bool log_result = true) {
  double result = math::apply([&](auto&&... args) { return functor(args...); },
                              args_tuple);
  double lfabs_result = log_result ? log(fabs(result)) : fabs(result);
  if (lfabs_result < err) {
    ArgsTupleT err_args_tuple = args_tuple;
    assign_err(std::get<ErrIndex>(err_args_tuple), err + lfabs_result);
    result = math::apply([&](auto&&... args) { return functor(args...); },
                         err_args_tuple);
  }
  return result;
}

enum class FunType { Density, GradT, GradA, GradV, GradW, GradSV };

template <FunType FunTypeEnum, typename KSSFuncType, typename KLLFuncType>
inline auto wiener5_helper(const double& y, const double& a, const double& vn,
                           const double& wn, const double& sv,
                           const double& err, const KSSFuncType& kss_functor,
                           const KLLFuncType& kll_functor) {
  double w = 1.0 - wn;
  double v = -vn;
  double sv_sqr = square(sv);
  double one_plus_svsqr_y = 1 + sv_sqr * y;
  double lg1;
  if (sv != 0) {
    lg1 = (sv_sqr * square(a * w) - 2 * a * v * w - square(v) * y) / 2.0
              / one_plus_svsqr_y
          - 2 * log(a) - 0.5 * log(one_plus_svsqr_y);
  } else {
    lg1 = (-2 * a * v * w - square(v) * y) / 2.0 - 2 * log(a);
  }
  double ans0 = 0;
  if (FunTypeEnum != FunType::Density) {
    double var_a = (FunTypeEnum == FunType::GradA) ? w : a;
    double var_b = (FunTypeEnum == FunType::GradA) ? a : w;
    if (sv != 0) {
      if (FunTypeEnum == FunType::GradT) {
        ans0 = -0.5
               * (square(sv_sqr) * (y + square(a * w))
                  + sv_sqr * (1 - 2 * a * v * w) + square(v))
               / square(one_plus_svsqr_y);
      } else {
        ans0 = (-v * var_a + sv_sqr * square(var_a) * var_b) / one_plus_svsqr_y;
      }
    } else {
      ans0 = (FunTypeEnum == FunType::GradT) ? -0.5 * square(v) : -v * var_a;
    }
  }

  double es = (err - lg1);
  if (FunTypeEnum == FunType::GradT) {
    es += 2.0 * log(a);
  }
  double y_asq = y / square(a);

  double K1, u_eps, arg;
  if (FunTypeEnum == FunType::Density) {
    K1 = (sqrt(2.0 * y_asq) + w) / 2.0;
    u_eps = fmin(-1.0, LOG_TWO + LOG_PI + 2.0 * log(y_asq) + 2.0 * (es));
    arg = -y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  } else {
    K1 = (sqrt(3.0 * y_asq) + w) / 2.0;
    if (FunTypeEnum == FunType::GradW) {
      u_eps = fmin(-1.0, 2.0 * es + LOG_TWO + LOG_PI + 2.0 * log(y_asq));
    } else {
      u_eps = fmin(
          -1.0, (log(8.0 / 27.0) + LOG_PI + 4.0 * log(y_asq) + 2.0 * es) / 3.0);
    }
    arg = -3.0 * y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  }

  double K2 = (arg > 0) ? 0.5 * (sqrt(arg) - w) : K1;
  size_t kss = static_cast<size_t>(ceil(fmax(K1, K2)));

  static const double PISQ = square(pi());  // pi*pi
  if (FunTypeEnum == FunType::Density) {
    K1 = 1.0 / (pi() * sqrt(y_asq));
    double two_log_piy = -2.0 * (log(pi() * y_asq) + es);
    K2 = (two_log_piy >= 0) ? sqrt(two_log_piy / (PISQ * y_asq)) : 0.0;
  } else if (FunTypeEnum == FunType::GradT) {
    K1 = sqrt(3.0 / y_asq) / pi();
    u_eps = fmin(-1.0, es + log(0.6) + LOG_PI + 2.0 * log(y_asq));
    arg = -2.0 / PISQ / y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
    K2 = (arg > 0) ? sqrt(arg) : K1;
  } else if (FunTypeEnum == FunType::GradW) {
    K1 = sqrt(2.0 / y_asq) / pi();
    u_eps = fmin(-1.0,
                 log(4.0 / 9.0) + 2.0 * LOG_PI + 3.0 * log(y_asq) + 2.0 * es);
    arg = -(u_eps - sqrt(-2.0 * u_eps - 2.0));
    K2 = (arg > 0) ? 1.0 / pi() * sqrt(arg / y_asq) : K1;
  }
  size_t kll = static_cast<size_t>(ceil(fmax(K1, K2)));

  double erg;
  int newsign;
  double fplus = NEGATIVE_INFTY;
  double fminus = NEGATIVE_INFTY;
  double twoy = 2.0 * y_asq;
  if (2 * kss <= kll) {
    double mult = (FunTypeEnum == FunType::Density) ? 1 : 3;
    double offset = (FunTypeEnum == FunType::GradW) ? y_asq : 0;
    double sqrt_offset = sqrt(offset);
    if (kss > 0) {
      for (size_t k = kss; k >= 1; k--) {
        double wp2k = w + 2.0 * k;
        double wm2k = w - 2.0 * k;

        if (wp2k > sqrt_offset) {
          fplus = log_sum_exp(
            mult * log(wp2k - offset) - (square(wp2k) - offset)
              / twoy,
            fplus);
        } else {
          fminus = log_sum_exp(
            mult * log(-(wp2k - offset)) - (square(wp2k) - offset)
              / twoy,
            fminus);
        }

        if (wm2k > sqrt_offset) {
          fplus = log_sum_exp(
            mult * log(wm2k - offset) - (square(wm2k) - offset)
              / twoy,
            fplus);
        } else {
          fminus = log_sum_exp(
            mult * log(-(wm2k - offset)) - (square(wm2k) - offset)
              / twoy,
            fminus);
        }
      }
    }
    if (w > sqrt_offset) {
      fplus = log_sum_exp(mult * log(w - offset) - square(w) / twoy, fplus);
    } else {
      fminus = log_sum_exp(mult * log(-(w - offset)) - square(w) / twoy,
                          fminus);
    }
    erg = (fplus < fminus) ? log_diff_exp(fminus, fplus)
                                  : log_diff_exp(fplus, fminus);
    newsign = (fplus < fminus) ? -1 : 1;
    return kss_functor(erg, newsign, lg1, ans0);
  } else {
    double mult;
    if (FunTypeEnum == FunType::Density) {
      mult = 1;
    } else if (FunTypeEnum == FunType::GradW) {
      mult = 2;
    } else {
      mult = 3;
    }
    double halfy = y_asq / 2.0;
    for (size_t k = kll; k >= 1; k--) {
      double pi_k = k * pi();
      double check
          = (FunTypeEnum == FunType::GradW) ? cos(pi_k * w) : sin(pi_k * w);
      if (check > 0) {
        fplus = log_sum_exp(mult * log(k) - square(pi_k) * halfy + log(check),
                            fplus);
      } else {
        fminus = log_sum_exp(mult * log(k) - square(pi_k) * halfy + log(-check),
                            fminus);
      }
    }
    erg = (fplus < fminus) ? log_diff_exp(fminus, fplus)
                                  : log_diff_exp(fplus, fminus);
    newsign = (fplus < fminus) ? -1 : 1;
    return kll_functor(erg, newsign, lg1, ans0);
  }
}

// calculate density in log
inline double dwiener5(const double& y, const double& a, const double& vn,
                       const double& wn, const double& sv, const double& err) {
  const auto& kss_functor = [&](double erg, int newsign, double lg1,
                                double ans0) {
    return lg1 - 0.5 * LOG_TWO - LOG_SQRT_PI - 1.5 * log(y / square(a)) + erg;
  };
  const auto& kll_functor = [&](double erg, int newsign, double lg1,
                                double ans0) { return lg1 + erg + LOG_PI; };

  return wiener5_helper<FunType::Density>(y, a, vn, wn, sv, err, kss_functor,
                                          kll_functor);
}
//-----------------------------------------------

// d/dt DENSITY
// calculate derivative of density with respect to t (in log, ans =
// d/dt(log(f))=d/dt f'/f; ans*exp(ld)=f' on normal scale)
inline double dtdwiener5(const double& y, const double& a, const double& vn,
                         const double& wn, const double& sv,
                         const double& err) {
  const auto& kss_functor = [&](double erg, int newsign, double lg1,
                                double ans0) {
    double ld = dwiener5(y, a, vn, wn, sv,
                         err - log(max(fabs(ans0 - 1.5 / y), fabs(ans0))));
    return ans0 - 1.5 / y
           + newsign
                 * exp(lg1 - 2.0 * log(a) - 1.5 * LOG_TWO - LOG_SQRT_PI
                       - 3.5 * log(y / square(a)) + erg - ld);
  };
  const auto& kll_functor = [&](double erg, int newsign, double lg1,
                                double ans0) {
    double ld = dwiener5(y, a, vn, wn, sv,
                         err - log(max(fabs(ans0 - 1.5 / y), fabs(ans0))));
    return ans0
           - newsign
                 * exp(lg1 - 2.0 * log(a) + 3.0 * LOG_PI - LOG_TWO + erg - ld);
  };

  return wiener5_helper<FunType::GradT>(y, a, vn, wn, sv, err, kss_functor,
                                        kll_functor);
}
//-----------------------------------------------

// d/da DENSITY
// calculate derivative of density with respect to a (in log, ans =
// d/da(log(f))=d/da f'/f; ans*exp(ld)=f' on normal scale)
inline double dadwiener5(const double& y, const double& a, const double& vn,
                         const double& wn, const double& sv, const double& err,
                         const int& normal_or_log) {
  double la = log(a);
  double ly = log(y);

  const auto& kss_functor
      = [&](double erg, int newsign, double lg1, double ans0) {
          double factor = lg1 - 3 * la;
          double ld = dwiener5(
              y, a, vn, wn, sv,
              err - log(max(fabs(ans0 + 1.0 / a), fabs(ans0 - 2.0 / a))));
          double ans = ans0 + 1.0 / a
                       - newsign
                             * exp(-0.5 * LOG_TWO - LOG_SQRT_PI - 2.5 * ly
                                   + 4.0 * la + lg1 + erg - ld);
          return (normal_or_log == 1) ? ans * exp(ld) : ans;
        };
  const auto& kll_functor = [&](double erg, int newsign, double lg1,
                                double ans0) {
    double factor = lg1 - 3 * la;
    double ld
        = dwiener5(y, a, vn, wn, sv,
                   err - log(max(fabs(ans0 + 1.0 / a), fabs(ans0 - 2.0 / a))));
    double ans
        = ans0 - 2.0 / a + newsign * exp(ly + factor + 3.0 * LOG_PI + erg - ld);
    return (normal_or_log == 1) ? ans * exp(ld) : ans;
  };

  return wiener5_helper<FunType::GradA>(y, a, vn, wn, sv, err, kss_functor,
                                        kll_functor);
}
//-----------------------------------------------

// d/dv DENSITY
// calculate derivative of density with respect to v (in log, ans =
// d/dv(log(f))=d/dv f'/f; ans*exp(ld)=f' on normal scale)
inline double dvdwiener5(const double& y, const double& a, const double& vn,
                         const double& wn, const double& sv) {
  double ans = (a * (1 - wn) - vn * y);
  if (sv != 0) {
    ans /= 1 + square(sv) * y;
  }
  return ans;
}
//-----------------------------------------------

// d/dw DENSITY
// calculate derivative of density with respect to w (in log, ans =
// d/dw(log(f))=d/dw f'/f; ans*exp(ld)=f' on normal scale)
inline double dwdwiener5(const double& y, const double& a, const double& vn,
                         const double& wn, const double& sv, const double& err,
                         const int& normal_or_log) {
  const auto& kss_functor
      = [&](double erg, int newsign, double lg1, double ans0) {
          double ld = dwiener5(y, a, vn, wn, sv, err - log(fabs(ans0)));
          double ls = -lg1 + ld;
          double ans = ans0
                       - newsign
                             * exp(erg - ls - 2.5 * log(y / square(a))
                                   - 0.5 * LOG_TWO - 0.5 * LOG_PI);
          return (normal_or_log == 1) ? -ans * exp(ld) : -ans;
        };
  const auto& kll_functor
      = [&](double erg, int newsign, double lg1, double ans0) {
          double ld = dwiener5(y, a, vn, wn, sv, err - log(fabs(ans0)));
          double ll = -lg1 + ld;
          double ans = ans0 + newsign * exp(erg - ll + 2 * LOG_PI);
          return (normal_or_log == 1) ? -ans * exp(ld) : -ans;
        };

  return wiener5_helper<FunType::GradW>(y, a, vn, wn, sv, err, kss_functor,
                                        kll_functor);
}
//-----------------------------------------------

// d/dsv DENSITY
// calculate derivative of density with respect to sv (in log, ans =
// d/dsv(log(f))=d/dsv f'/f; ans*exp(ld)=f' on normal scale)
inline double dsvdwiener5(const double& y, const double& a, const double& vn,
                          const double& wn, const double& sv) {
  double v, w;
  v = -vn;
  w = 1 - wn;
  double one_sqrsv_y = 1 + square(sv) * y;
  double t1 = -y / one_sqrsv_y;
  double t2 = (square(a * w) + 2 * a * v * w * y + square(v * y))
              / square(one_sqrsv_y);
  return sv * (t1 + t2);
}
//-----------------------------------------------

}  // namespace internal

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

  double error_bound_dens = 1e-6;  // precision for density
  double lerror_bound_dens = log(error_bound_dens);
  double error_bound = prec;  // precision for
  // derivatives (controllable by user)
  double lerror_bound = log(error_bound);  // log(alpha)
  double abstol = 0.0;
  double reltol = .9 * error_bound;  // eps_rel(Integration)
  double abstol_wiener5 = 1e-12;     // eps_abs(wiener5)
  double labstol_wiener5 = log(abstol_wiener5);
  // log(eps_abs(wiener5)
  int Meval = 6000;
  double dens = 0.0;
  double ld = 0.0;
  operands_and_partials<T_y_ref, T_a_ref, T_t0_ref, T_w_ref, T_v_ref, T_sv_ref>
      ops_partials(y_ref, a_ref, t0_ref, w_ref, v_ref, sv_ref);

  static constexpr double LOG_FOUR = LOG_TWO + LOG_TWO;
  static constexpr double LOG_POINT1 = -1;

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
                                        sv_val, labstol_wiener5);

    dens = internal::estimate_with_err_check<5>(
        internal::dwiener5, lerror_bound_dens - LOG_TWO, params, true);
    ld += dens;

    // computation of derivative for t and precision check in order to give
    // the value as deriv_y to edge1 and as -deriv_y to edge5
    double deriv_y = internal::estimate_with_err_check<5>(
        internal::dtdwiener5, dens + lerror_bound - LOG_TWO, params);

    // computation of derivatives and precision checks
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[i] = deriv_y;
    }
    if (!is_constant_all<T_a>::value) {
      ops_partials.edge2_.partials_[i] = internal::estimate_with_err_check<5>(
          internal::dadwiener5, dens + lerror_bound - LOG_FOUR,
          std::tuple_cat(params, std::make_tuple(0)));
    }
    if (!is_constant_all<T_t0>::value) {
      ops_partials.edge3_.partials_[i] = -deriv_y;
    }
    if (!is_constant_all<T_w>::value) {
      ops_partials.edge4_.partials_[i] = internal::estimate_with_err_check<5>(
          internal::dwdwiener5, dens + lerror_bound - LOG_FOUR,
          std::tuple_cat(params, std::make_tuple(0)));
    }
    if (!is_constant_all<T_v>::value) {
      ops_partials.edge5_.partials_[i]
          = internal::dvdwiener5(y_val - t0_val, a_val, v_val, w_val, sv_val);
    }
    if (!is_constant_all<T_sv>::value) {
      ops_partials.edge6_.partials_[i]
          = internal::dsvdwiener5(y_val - t0_val, a_val, v_val, w_val, sv_val);
    }
  }  // end for loop
  return ops_partials.build(ld);
}  // end wiener5_lpdf
}  // namespace math
}  // namespace stan
#endif
