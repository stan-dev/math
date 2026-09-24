#ifndef STAN_MATH_PRIM_PROB_CATEGORICAL_LOGIT_LPMF_HPP
#define STAN_MATH_PRIM_PROB_CATEGORICAL_LOGIT_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log_softmax.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

// CategoricalLog(n|theta)  [0 < n <= N, theta unconstrained], no checking
template <bool propto, typename T_prob, require_col_vector_t<T_prob>* = nullptr>
inline return_type_t<T_prob> categorical_logit_lpmf(int n, const T_prob& beta) {
  static constexpr const char* function = "categorical_logit_lpmf";
  check_bounded(function, "categorical outcome out of support", n, 1,
                beta.size());
  ref_type_t<T_prob> beta_ref = beta;

  // Autodiff args must be finite, data can be +/-inf
  if constexpr (is_constant<T_prob>::value) {
    check_not_nan(function, "log odds parameter", beta_ref);
    check_greater(function, "log odds parameter", beta_ref.maxCoeff(),
                  NEGATIVE_INFTY);
  } else {
    check_finite(function, "log odds parameter", beta_ref);
  }

  if constexpr (!include_summand<propto, T_prob>::value) {
    return 0.0;
  }

  // INFTY case: uniform over the +inf entries, zero probability elsewhere
  if constexpr (is_constant<T_prob>::value) {
    int num_infty = (beta_ref.array() == INFTY).count();
    if (num_infty > 0) {
      return beta_ref.coeff(n - 1) == INFTY ? -std::log(num_infty)
                                            : NEGATIVE_INFTY;
    }
  }

  // FIXME:  wasteful vs. creating term (n-1) if not vectorized
  return beta_ref.coeff(n - 1)
         - log_sum_exp(beta_ref);  // == log_softmax(beta)(n-1);
}

template <bool propto, typename T_prob, require_col_vector_t<T_prob>* = nullptr>
inline return_type_t<T_prob> categorical_logit_lpmf(const std::vector<int>& ns,
                                                    const T_prob& beta) {
  static constexpr const char* function = "categorical_logit_lpmf";

  check_bounded(function, "categorical outcome out of support", ns, 1,
                beta.size());
  ref_type_t<T_prob> beta_ref = beta;

  // Autodiff args must be finite, data can be +/-inf
  if constexpr (is_constant<T_prob>::value) {
    check_not_nan(function, "log odds parameter", beta_ref);
    if (beta_ref.size() > 0) {
      check_greater(function, "log odds parameter", beta_ref.maxCoeff(),
                    NEGATIVE_INFTY);
    }
  } else {
    check_finite(function, "log odds parameter", beta_ref);
  }

  if constexpr (!include_summand<propto, T_prob>::value) {
    return 0.0;
  }

  if (ns.empty()) {
    return 0.0;
  }

  // INFTY case: uniform over the +inf entries, zero probability elsewhere
  if constexpr (is_constant<T_prob>::value) {
    int num_infty = (beta_ref.array() == INFTY).count();
    if (num_infty > 0) {
      double lp = 0.0;
      for (int n : ns) {
        lp += beta_ref.coeff(n - 1) == INFTY ? -std::log(num_infty)
                                             : NEGATIVE_INFTY;
      }
      return lp;
    }
  }

  auto log_softmax_beta = to_ref(log_softmax(beta_ref));

  // FIXME:  replace with more efficient sum()
  Eigen::Matrix<return_type_t<T_prob>, Eigen::Dynamic, 1> results(ns.size());
  for (size_t i = 0; i < ns.size(); ++i) {
    results[i] = log_softmax_beta(ns[i] - 1);
  }
  return sum(results);
}

template <typename T_n, typename T_prob, require_st_integral<T_n>* = nullptr,
          require_col_vector_t<T_prob>* = nullptr>
inline return_type_t<T_prob> categorical_logit_lpmf(const T_n& ns,
                                                    const T_prob& beta) {
  return categorical_logit_lpmf<false>(ns, beta);
}

}  // namespace math
}  // namespace stan
#endif
