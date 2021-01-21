#ifndef STAN_MATH_PRIM_PROB_CATEGORICAL_LPMF_HPP
#define STAN_MATH_PRIM_PROB_CATEGORICAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log PMF of the categorical distribution. If containers of
 * integers and/or probabilities are supplied, returns the log sum of the
 * PMF.
 *
 * @tparam T_n type of integer parameters
 * @tparam T_prob type probability vector(s)
 * @param n integer parameter(s)
 * @param theta probability vector(s)
 * @return log probability
 */
template <bool propto, typename T_n, typename T_prob>
return_type_t<T_prob> categorical_lpmf(const T_n& n, const T_prob& theta) {
  static const char* function = "categorical_lpmf";

  ref_type_t<T_prob> theta_ref = theta;
  scalar_seq_view<T_n> n_vec(n);
  vector_seq_view<ref_type_t<T_prob>> theta_vec(theta_ref);

  size_t vec_size = std::max(size(n), size_mvt(theta_ref));

  if (!include_summand<propto, T_prob>::value || size(n) == 0) {
    return 0.0;
  }

  if (size_mvt(theta_ref) > 1) {
    check_consistent_sizes(function, "Integer parameter", n,
                           "Probabilities parameter", theta_ref);
  }

  for (size_t i = 0; i < vec_size; ++i) {
    check_bounded(function, "Number of categories", n_vec[i], 1,
                  theta_vec[i].size());
    check_simplex(function, "Probabilities parameter", theta_vec[i]);
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
