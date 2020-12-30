#ifndef STAN_MATH_PRIM_PROB_CATEGORICAL_LPMF_HPP
#define STAN_MATH_PRIM_PROB_CATEGORICAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

// Categorical(n|theta)  [0 < n <= N;   0 <= theta[n] <= 1;  SUM theta = 1]
template <bool propto, typename T_n, typename T_prob>
return_type_t<T_prob> categorical_lpmf(const T_n& n, const T_prob& theta) {
  static const char* function = "categorical_lpmf";

  decltype(auto) theta_ref = to_ref(theta);
  scalar_seq_view<T_n> n_vec(n);
  vector_seq_view<T_prob> theta_vec(theta_ref);

  size_t vec_size = std::max(stan::math::size(n), stan::math::size_mvt(theta_ref));

  for (size_t i = 0; i < vec_size; ++i) {
    check_bounded(function, "Number of categories", n, 1, theta_vec[i].size());
    check_simplex(function, "Probabilities parameter", theta_vec[i]);
  }

  if (!include_summand<propto, T_prob>::value || stan::math::size(n) == 0) {
    return 0.0;
  }

  using T_plain = plain_type_t<decltype(theta_ref)>;

  T_plain log_theta = log(theta_ref);
  vector_seq_view<T_plain> log_theta_vec(log_theta);

  return_type_t<T_prob> lp(0);
  for (size_t i = 0; i < vec_size; ++i) {
    lp += log_theta_vec[i].coeffRef(n_vec[i] - 1);
  }

  return lp;
}

template <typename T_n, typename T_prob>
inline return_type_t<T_prob> categorical_lpmf(const T_n& ns,
                                              const T_prob& theta) {
  return categorical_lpmf<false>(ns, theta);
}

}  // namespace math
}  // namespace stan
#endif
