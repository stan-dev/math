#ifndef STAN_MATH_PRIM_MAT_PROB_MULTINOMIAL_LPMF_HPP
#define STAN_MATH_PRIM_MAT_PROB_MULTINOMIAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/err/check_simplex.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <vector>

namespace stan {
namespace math {
// Multinomial(ns|N, theta)   [0 <= n <= N;  SUM ns = N;
//                            0 <= theta[n] <= 1;  SUM theta = 1]
template <bool propto, typename T_prob>
return_type_t<T_prob> multinomial_lpmf(
    const std::vector<int>& ns,
    const Eigen::Matrix<T_prob, Eigen::Dynamic, 1>& theta) {
  static const char* function = "multinomial_lpmf";

  using boost::math::tools::promote_args;

  return_type_t<T_prob> lp(0.0);
  check_nonnegative(function, "Number of trials variable", ns);
  check_simplex(function, "Probabilities parameter", theta);
  check_size_match(function, "Size of number of trials variable", ns.size(),
                   "rows of probabilities parameter", theta.rows());

  if (include_summand<propto>::value) {
    double sum = 1.0;
    for (int n : ns) {
      sum += n;
    }
    lp += lgamma(sum);
    for (int n : ns) {
      lp -= lgamma(n + 1.0);
    }
  }
  if (include_summand<propto, T_prob>::value) {
    for (unsigned int i = 0; i < ns.size(); ++i) {
      lp += multiply_log(ns[i], theta[i]);
    }
  }
  return lp;
}

template <typename T_prob>
return_type_t<T_prob> multinomial_lpmf(
    const std::vector<int>& ns,
    const Eigen::Matrix<T_prob, Eigen::Dynamic, 1>& theta) {
  return multinomial_lpmf<false>(ns, theta);
}

}  // namespace math
}  // namespace stan
#endif
