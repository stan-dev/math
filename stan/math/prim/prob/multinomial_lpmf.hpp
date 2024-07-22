#ifndef STAN_MATH_PRIM_PROB_MULTINOMIAL_LPMF_HPP
#define STAN_MATH_PRIM_PROB_MULTINOMIAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <vector>

namespace stan {
namespace math {
// Multinomial(ns|N, theta)   [0 <= n <= N;  SUM ns = N;
//                            0 <= theta[n] <= 1;  SUM theta = 1]
template <bool propto, typename T_prob,
          require_eigen_col_vector_t<T_prob>* = nullptr>
return_type_t<T_prob> multinomial_lpmf(const std::vector<int>& ns,
                                       const T_prob& theta) {
  static constexpr const char* function = "multinomial_lpmf";
  check_size_match(function, "Size of number of trials variable", ns.size(),
                   "rows of probabilities parameter", theta.rows());
  check_nonnegative(function, "Number of trials variable", ns);
  const auto& theta_ref = to_ref(theta);
  check_simplex(function, "Probabilities parameter", theta_ref);

  return_type_t<T_prob> lp(0.0);

  if (include_summand<propto>::value) {
    double sum = 1.0;
    for (int n : ns) {
      sum += n;
      lp -= lgamma(n + 1.0);
    }
    lp += lgamma(sum);
  }
  if (include_summand<propto, T_prob>::value) {
    for (unsigned int i = 0; i < ns.size(); ++i) {
      lp += multiply_log(ns[i], theta_ref.coeff(i));
    }
  }
  return lp;
}

template <typename T_prob>
return_type_t<T_prob> multinomial_lpmf(const std::vector<int>& ns,
                                       const T_prob& theta) {
  return multinomial_lpmf<false>(ns, theta);
}

}  // namespace math
}  // namespace stan
#endif
