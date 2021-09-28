#ifndef STAN_MATH_PRIM_PROB_CATEGORICAL_LPMF_HPP
#define STAN_MATH_PRIM_PROB_CATEGORICAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

// Categorical(n|theta)  [0 < n <= N;   0 <= theta[n] <= 1;  SUM theta = 1]
template <bool propto, typename T_prob,
          require_eigen_col_vector_t<T_prob>* = nullptr>
return_type_t<T_prob> categorical_lpmf(int n, const T_prob& theta) {
  static const char* function = "categorical_lpmf";
  using std::log;

  check_bounded(function, "Number of categories", n, 1, theta.size());
  ref_type_t<T_prob> theta_ref = theta;
  check_simplex(function, "Probabilities parameter", value_of(theta_ref));

  if (include_summand<propto, T_prob>::value) {
    return log(theta_ref.coeff(n - 1));
  }
  return 0.0;
}

template <bool propto, typename T_prob,
          require_eigen_col_vector_t<T_prob>* = nullptr>
return_type_t<T_prob> categorical_lpmf(const std::vector<int>& ns,
                                       const T_prob& theta) {
  static const char* function = "categorical_lpmf";

  check_bounded(function, "element of outcome array", ns, 1, theta.size());
  ref_type_t<T_prob> theta_ref = theta;
  check_simplex(function, "Probabilities parameter", value_of(theta_ref));

  if (!include_summand<propto, T_prob>::value) {
    return 0.0;
  }

  if (ns.size() == 0) {
    return 0.0;
  }

  Eigen::Matrix<value_type_t<T_prob>, Eigen::Dynamic, 1> log_theta
      = log(theta_ref);
  Eigen::Matrix<return_type_t<T_prob>, Eigen::Dynamic, 1> log_theta_ns(
      ns.size());
  for (size_t i = 0; i < ns.size(); ++i) {
    log_theta_ns(i) = log_theta(ns[i] - 1);
  }

  return sum(log_theta_ns);
}

template <typename T_n, typename T_prob, require_st_integral<T_n>* = nullptr,
          require_eigen_col_vector_t<T_prob>* = nullptr>
inline return_type_t<T_prob> categorical_lpmf(const T_n& ns,
                                              const T_prob& theta) {
  return categorical_lpmf<false>(ns, theta);
}

}  // namespace math
}  // namespace stan
#endif
