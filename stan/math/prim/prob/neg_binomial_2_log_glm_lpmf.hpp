#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LOG_GLM_LPMF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LOG_GLM_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/multiply.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma_diff.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <vector>
#include <cmath>

namespace stan {
namespace math {
namespace internal {

/**
 * Return f(y[n]) for every instance n, with f evaluated once per distinct
 * count. The counts must be in [0, y_max].
 *
 * @tparam R return type of f
 * @tparam T_y_arr type of the array of counts
 * @tparam F type of the function
 * @param y counts
 * @param y_max largest count
 * @param f function of a count
 * @return array of f(y[n])
 */
template <typename R, typename T_y_arr, typename F>
inline Eigen::Array<R, Eigen::Dynamic, 1> per_count_values(const T_y_arr& y,
                                                           int y_max,
                                                           const F& f) {
  std::vector<R> table(y_max + 1);
  std::vector<char> done(y_max + 1, 0);
  Eigen::Array<R, Eigen::Dynamic, 1> out(y.size());
  for (Eigen::Index n = 0; n < y.size(); ++n) {
    const int k = y(n);
    if (!done[k]) {
      table[k] = f(k);
      done[k] = 1;
    }
    out(n) = table[k];
  }
  return out;
}

}  // namespace internal

/** \ingroup multivar_dists
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with Negative-Binomial-2 distribution and log link function.
 * The idea is that neg_binomial_2_log_glm_lpmf(y, x, alpha, beta, phi) should
 * compute a more efficient version of
 * neg_binomial_2_log_lpmf(y, alpha + x * beta, phi) by using analytically
 * simplified gradients.
 * If containers are supplied, returns the log sum of the probabilities.
 *
 * @tparam T_y type of positive int vector of variates (labels);
 * this can also be a single positive integer value;
 * @tparam T_x_scalar type of a scalar in the matrix of independent variables
 * (features)
 * @tparam T_x_rows compile-time number of rows of `x`. It can be either
 * `Eigen::Dynamic` or 1.
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector;
 * this can also be a scalar;
 * @tparam T_precision type of the (positive) precision(s);
 * this can be a vector (of the same length as y, for heteroskedasticity)
 * or a scalar.
 *
 * @param y failures count scalar or vector parameter. If it is a scalar it
 * will be broadcast - used for all instances.
 * @param x design matrix or row vector. If it is a row vector it will be
 * broadcast - used for all instances.
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @param phi (vector of) precision parameter(s)
 * @return log probability or log sum of probabilities
 * @throw std::invalid_argument if container sizes mismatch.
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if phi is infinite or non-positive.
 * @throw std::domain_error if y is negative.
 */
template <bool propto, typename T_y, typename T_x, typename T_alpha,
          typename T_beta, typename T_precision,
          require_matrix_t<T_x>* = nullptr>
inline return_type_t<T_x, T_alpha, T_beta, T_precision>
neg_binomial_2_log_glm_lpmf(const T_y& y, const T_x& x, const T_alpha& alpha,
                            const T_beta& beta, const T_precision& phi) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::exp;
  using Eigen::log1p;
  using Eigen::Matrix;
  constexpr int T_x_rows = T_x::RowsAtCompileTime;
  using T_partials_return
      = partials_return_t<T_y, T_x, T_alpha, T_beta, T_precision>;
  using T_precision_val = typename std::conditional_t<
      is_vector<T_precision>::value,
      Eigen::Array<partials_return_t<T_precision>, -1, 1>,
      partials_return_t<T_precision>>;
  using T_sum_val = typename std::conditional_t<
      is_vector<T_y>::value || is_vector<T_precision>::value,
      Eigen::Array<partials_return_t<T_y, T_precision>, -1, 1>,
      partials_return_t<T_y, T_precision>>;
  using T_theta_tmp =
      typename std::conditional_t<T_x_rows == 1, T_partials_return,
                                  Array<T_partials_return, Dynamic, 1>>;
  using T_x_ref = ref_type_if_not_constant_t<T_x>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_alpha>;
  using T_beta_ref = ref_type_if_not_constant_t<T_beta>;
  using T_phi_ref = ref_type_if_not_constant_t<T_precision>;

  const size_t N_instances = T_x_rows == 1 ? stan::math::size(y) : x.rows();
  const size_t N_attributes = x.cols();

  static constexpr const char* function = "neg_binomial_2_log_glm_lpmf";
  check_consistent_size(function, "Vector of dependent variables", y,
                        N_instances);
  check_consistent_size(function, "Weight vector", beta, N_attributes);
  check_consistent_size(function, "Vector of precision parameters", phi,
                        N_instances);
  check_consistent_size(function, "Vector of intercepts", alpha, N_instances);
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  const auto& beta_val = value_of(beta_ref);
  const auto& alpha_val = value_of(alpha_ref);
  const auto& beta_val_vec = to_ref(as_column_vector_or_scalar(beta_val));
  const auto& alpha_val_vec = to_ref(as_column_vector_or_scalar(alpha_val));
  check_finite(function, "Weight vector", beta_val_vec);
  check_finite(function, "Intercept", alpha_val_vec);

  if (size_zero(y, phi)) {
    return 0;
  }

  const auto& y_ref = to_ref(y);
  T_phi_ref phi_ref = phi;

  const auto& y_val = value_of(y_ref);
  const auto& phi_val = value_of(phi_ref);

  const auto& y_val_vec = to_ref(as_column_vector_or_scalar(y_val));
  const auto& phi_val_vec = to_ref(as_column_vector_or_scalar(phi_val));
  check_nonnegative(function, "Failures variables", y_val_vec);
  check_positive_finite(function, "Precision parameter", phi_val_vec);

  if constexpr (!include_summand<propto, T_x, T_alpha, T_beta,
                                 T_precision>::value) {
    return 0;
  }

  T_x_ref x_ref = x;

  const auto& x_val = to_ref_if<is_autodiff_v<T_beta>>(value_of(x_ref));

  const auto& y_arr = as_array_or_scalar(y_val_vec);
  const auto& phi_arr = as_array_or_scalar(phi_val_vec);

  Array<T_partials_return, Dynamic, 1> theta(N_instances);
  if constexpr (T_x_rows == 1) {
    T_theta_tmp theta_tmp = (x_val * beta_val_vec)(0, 0);
    theta = theta_tmp + as_array_or_scalar(alpha_val_vec);
  } else {
    theta = multiply(x_val, beta_val_vec).array();
    theta += as_array_or_scalar(alpha_val_vec);
  }
  check_finite(function, "Matrix of independent variables", theta);
  T_precision_val log_phi = log(phi_arr);
  // With u = theta - log(phi) and e = exp(-|u|), log1p_exp(u) and
  // log1p_exp(-u) share log1p(e), and e also gives the partials below.
  Array<T_partials_return, Dynamic, 1> u = theta - log_phi;
  Array<T_partials_return, Dynamic, 1> e = (u > 0).select(-u, u).exp();
  Array<T_partials_return, Dynamic, 1> log1p_e = log1p(e);
  // log1p(mu / phi) with mu = exp(theta)
  Array<T_partials_return, Dynamic, 1> log1p_exp_theta_m_log_phi
      = (u > 0).select(u + log1p_e, log1p_e);

  T_sum_val y_plus_phi = y_arr + phi_arr;

  // Compute the log-density. The log pmf of one instance, with
  // mu = exp(theta), is
  //
  //   lchoose(y + phi - 1, y) - y log1p(phi / mu) - phi log1p(mu / phi),
  //
  // the form of neg_binomial_2_log_lpmf. Formed as phi log(phi) - lgamma(phi)
  // + lgamma(y + phi) - (y + phi) log(mu + phi) + y theta, the terms are of
  // the size of phi log(phi) and cancel for large phi.
  T_partials_return logp = -sum(y_arr * (u > 0).select(log1p_e, log1p_e - u)
                                + phi_arr * log1p_exp_theta_m_log_phi);
  // With a scalar phi and a vector y, the lchoose, lgamma(y + 1) and digamma
  // terms depend on y only. Counts repeat, so each is evaluated once per
  // distinct count (internal::per_count_values), if the largest count stays
  // proportional to the number of instances.
  constexpr bool per_count
      = is_vector<T_y>::value && !is_vector<T_precision>::value;
  int y_max = -1;
  if constexpr (per_count) {
    if (N_instances > 0) {
      const int y_max_all = y_arr.maxCoeff();
      if (static_cast<size_t>(y_max_all) <= 4 * N_instances + 64) {
        y_max = y_max_all;
      }
    }
  }
  const bool use_table = y_max >= 0;
  // lchoose(y + phi - 1, y) = -lbeta(phi, y + 1) - log(y + phi). The form
  // binomial_coefficient_log(y + phi - 1, y) rounds phi to the precision of
  // y in y + phi - 1, which loses phi where it is small compared to y.
  if constexpr (include_summand<propto, T_precision>::value) {
    bool lchoose_done = false;
    if constexpr (per_count) {
      if (use_table) {
        const auto lchoose_neg = internal::per_count_values<T_partials_return>(
            y_arr, y_max,
            [&](int k) { return lbeta(phi_val, k + 1.0) + log(k + phi_val); });
        for (size_t n = 0; n < N_instances; ++n) {
          logp -= lchoose_neg(n);
        }
        lchoose_done = true;
      }
    }
    if (!lchoose_done) {
      if constexpr (is_vector<T_y>::value || is_vector<T_precision>::value) {
        scalar_seq_view<decltype(y_val_vec)> y_vec(y_val_vec);
        scalar_seq_view<decltype(phi_val_vec)> phi_vec(phi_val_vec);
        for (size_t n = 0; n < N_instances; ++n) {
          logp
              -= lbeta(phi_vec[n], y_vec[n] + 1.0) + log(y_vec[n] + phi_vec[n]);
        }
      } else {
        logp -= N_instances
                * (lbeta(phi_val, y_val + 1.0) + log(y_val + phi_val));
      }
    }
  }
  // Under propto, remove -lgamma(y + 1) (in lchoose), phi log(phi) for data
  // phi, and y theta for data x, alpha and beta (in the log1p_exp terms).
  if constexpr (!include_summand<propto>::value) {
    if constexpr (include_summand<propto, T_precision>::value) {
      if constexpr (is_vector<T_y>::value) {
        bool lgamma_done = false;
        if constexpr (per_count) {
          if (use_table) {
            logp += sum(internal::per_count_values<double>(
                y_arr, y_max, [](int k) { return lgamma(k + 1.0); }));
            lgamma_done = true;
          }
        }
        if (!lgamma_done) {
          logp += sum(lgamma(y_arr + 1.0));
        }
      } else {
        logp += sum(lgamma(y_arr + 1.0)) * N_instances;
      }
    } else {
      if constexpr (is_vector<T_precision>::value) {
        logp -= sum(phi_arr * log_phi);
      } else {
        logp -= N_instances * multiply_log(phi_val, phi_val);
      }
    }
    if constexpr (!include_summand<propto, T_x, T_alpha, T_beta>::value) {
      logp -= sum(y_arr * theta);
    }
  }

  // Compute the necessary derivatives.
  auto ops_partials
      = make_partials_propagator(x_ref, alpha_ref, beta_ref, phi_ref);
  if constexpr (is_any_autodiff_v<T_x, T_beta, T_alpha, T_precision>) {
    // mu / (mu + phi) and phi / (mu + phi), mu = exp(theta), are
    // 1 / (1 + e) and e / (1 + e) with e = exp(-|u|), u = theta - log(phi),
    // in the order given by the sign of u. mu itself overflows from
    // theta = 709.8, and the derivative of 1 - mu / (mu + phi) cancels.
    Array<T_partials_return, Dynamic, 1> inv_1p_e = 1.0 / (1.0 + e);
    Array<T_partials_return, Dynamic, 1> e_over_1p_e = e * inv_1p_e;
    Array<T_partials_return, Dynamic, 1> mu_over_mu_p_phi
        = (u > 0).select(inv_1p_e, e_over_1p_e);
    if constexpr (is_any_autodiff_v<T_x, T_beta, T_alpha>) {
      Matrix<T_partials_return, Dynamic, 1> theta_derivative
          = y_arr - y_plus_phi * mu_over_mu_p_phi;
      if constexpr (is_autodiff_v<T_beta>) {
        if constexpr (T_x_rows == 1) {
          edge<2>(ops_partials).partials_ = theta_derivative.sum() * x_val;
        } else if constexpr (T_x::ColsAtCompileTime == 1) {
          edge<2>(ops_partials).partials_[0]
              = dot_product(x_val, theta_derivative);
        } else {
          edge<2>(ops_partials).partials_
              = multiply(x_val.transpose(), theta_derivative);
        }
      }
      if constexpr (is_autodiff_v<T_x>) {
        if constexpr (T_x_rows == 1) {
          edge<0>(ops_partials).partials_
              = beta_val_vec * theta_derivative.sum();
        } else {
          edge<0>(ops_partials).partials_
              = (beta_val_vec * theta_derivative.transpose()).transpose();
        }
      }
      if constexpr (is_autodiff_v<T_alpha>) {
        if constexpr (is_vector<T_alpha>::value) {
          partials<1>(ops_partials) = std::move(theta_derivative);
        } else {
          partials<1>(ops_partials)[0] = sum(theta_derivative);
        }
      }
    }
    if constexpr (is_autodiff_v<T_precision>) {
      // 1 - (y + phi) / (mu + phi) = mu / (mu + phi) - (y / phi) phi
      // / (mu + phi) and log(phi) - log(mu + phi) are formed without
      // cancellation, and digamma_diff replaces
      // digamma(y + phi) - digamma(phi)
      Array<T_partials_return, Dynamic, 1> phi_over_mu_p_phi
          = (u > 0).select(e_over_1p_e, inv_1p_e);
      if constexpr (is_vector<T_precision>::value) {
        edge<3>(ops_partials).partials_
            = mu_over_mu_p_phi - phi_over_mu_p_phi * y_arr / phi_arr
              - log1p_exp_theta_m_log_phi + digamma_diff(phi_arr, y_arr);
      } else {
        bool dd_done = false;
        if constexpr (per_count) {
          if (use_table) {
            const auto dd = internal::per_count_values<T_partials_return>(
                y_arr, y_max, [&](int k) { return digamma_diff(phi_val, k); });
            partials<3>(ops_partials)[0]
                = sum(mu_over_mu_p_phi - phi_over_mu_p_phi * y_arr / phi_arr
                      - log1p_exp_theta_m_log_phi + dd);
            dd_done = true;
          }
        }
        if (!dd_done) {
          partials<3>(ops_partials)[0]
              = sum(mu_over_mu_p_phi - phi_over_mu_p_phi * y_arr / phi_arr
                    - log1p_exp_theta_m_log_phi + digamma_diff(phi_arr, y_arr));
        }
      }
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_x, typename T_alpha, typename T_beta,
          typename T_precision>
inline return_type_t<T_x, T_alpha, T_beta, T_precision>
neg_binomial_2_log_glm_lpmf(const T_y& y, const T_x& x, const T_alpha& alpha,
                            const T_beta& beta, const T_precision& phi) {
  return neg_binomial_2_log_glm_lpmf<false>(y, x, alpha, beta, phi);
}
}  // namespace math
}  // namespace stan
#endif
