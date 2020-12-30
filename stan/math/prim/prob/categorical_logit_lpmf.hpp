#ifndef STAN_MATH_PRIM_PROB_CATEGORICAL_LOGIT_LPMF_HPP
#define STAN_MATH_PRIM_PROB_CATEGORICAL_LOGIT_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/log_softmax.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <vector>

namespace stan {
namespace math {

// CategoricalLog(n|theta)  [0 < n <= N, theta unconstrained], no checking
template <bool propto, typename T_n, typename T_prob>
return_type_t<T_prob> categorical_logit_lpmf(const T_n& n,
                                             const T_prob& beta) {
  static const char* function = "categorical_logit_lpmf";

  scalar_seq_view<T_n> n_vec(n);
  vector_seq_view<T_prob> beta_vec(beta);

  size_t vec_size = std::max(stan::math::size(n), stan::math::size_mvt(beta));

  for (size_t i = 0; i < vec_size; ++i) {
    check_bounded(function, "categorical outcome out of support", n, 1,
                  beta_vec[i].size());
    check_finite(function, "log odds parameter", beta_vec[i]);
  }

  if (!include_summand<propto, T_prob>::value || stan::math::size(n) == 0) {
    return 0.0;
  }

  using T_plain = plain_type_t<decltype(beta)>;

  T_plain log_softmax_beta = log_softmax(beta);
  vector_seq_view<T_plain> log_softmax_beta_vec(log_softmax_beta);

  return_type_t<T_prob> lp(0);
  for (size_t i = 0; i < vec_size; ++i) {
    lp += log_softmax_beta_vec[i].coeffRef(n_vec[i] - 1);
  }
  return lp;
}

template <typename T_n, typename T_prob>
inline return_type_t<T_prob> categorical_logit_lpmf(const T_n& ns,
                                                    const T_prob& beta) {
  return categorical_logit_lpmf<false>(ns, beta);
}

}  // namespace math
}  // namespace stan
#endif
