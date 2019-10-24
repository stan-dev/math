#ifndef STAN_MATH_PRIM_MAT_PROB_DIRICHLET_LPMF_HPP
#define STAN_MATH_PRIM_MAT_PROB_DIRICHLET_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/mat/err/check_simplex.hpp>
#include <stan/math/prim/mat/fun/lgamma.hpp>
#include <stan/math/prim/mat/fun/digamma.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * The log of the Dirichlet density for the given theta and
 * a vector of prior sample sizes, alpha.
 * Each element of alpha must be greater than 0.
 * Each element of theta must be greater than or 0.
 * Theta sums to 1.
 *
 * \f[
 * \theta\sim\mbox{Dirichlet}(\alpha_1,\ldots,\alpha_k)\\
 * \log(p(\theta\,|\,\alpha_1,\ldots,\alpha_k))=\log\left(
 * \frac{\Gamma(\alpha_1+\cdots+\alpha_k)}{\Gamma(\alpha_1)+
 * \cdots+\Gamma(\alpha_k)}* \left(\theta_1^{\alpha_1-1}+
 * \cdots+\theta_k^{\alpha_k-1}\right)\right)\\
 * =\log(\Gamma(\alpha_1+\cdots+\alpha_k))-\left(
 * \log(\Gamma(\alpha_1))+\cdots+\log(\Gamma(\alpha_k))\right)+
 * (\alpha_1-1)\log(\theta_1)+\cdots+(\alpha_k-1)\log(\theta_k)
 * \f]
 *
 * \f[
 * \frac{\partial }{\partial
 * \theta_x}\log(p(\theta\,|\,\alpha_1,\ldots,\alpha_k))=
 * \frac{\alpha_x-1}{\theta_x}
 * \f]
 *
 * \f[
 * \frac{\partial}{\partial\alpha_x}\log(p(\theta\,|\,\alpha_1,\ldots,\alpha_k))
 * =\psi_{(0)}(\sum\alpha)-\psi_{(0)}(\alpha_x)+\log\theta_x
 * \f]
 *
 * @param theta A scalar vector.
 * @param alpha Prior sample sizes.
 * @return The log of the Dirichlet density.
 * @throw std::domain_error if any element of alpha is less than
 * or equal to 0.
 * @throw std::domain_error if any element of theta is less than 0.
 * @throw std::domain_error if the sum of theta is not 1.
 * @tparam T_prob Type of scalar.
 * @tparam T_prior_size Type of prior sample sizes.
 */
template <bool propto, typename T_prob, typename T_prior_size>
return_type_t<T_prob, T_prior_size> dirichlet_lpmf(const T_prob& theta,
                                                   const T_prior_size& alpha) {
  static const char* function = "dirichlet_lpmf";

  using T_partials_return = partials_return_t<T_prob, T_prior_size>;
  using T_partials_vec = typename Eigen::Matrix<T_partials_return, -1, 1>;
  using T_partials_mat = typename Eigen::Matrix<T_partials_return, -1, -1>;

  vector_seq_view<T_prob> theta_vec(theta);
  vector_seq_view<T_prior_size> alpha_vec(alpha);
  const size_t t_length = max_size_mvt(theta, alpha);

  check_consistent_sizes(function, "probabilities", theta_vec[0],
                         "prior sample sizes", alpha_vec[0]);

  for (size_t t = 0; t < t_length; t++) {
    check_positive(function, "prior sample sizes", alpha_vec[t]);
    check_simplex(function, "probabilities", theta_vec[t]);
  }

  const size_t t_size = theta_vec[0].size();

  T_partials_mat theta_dbl(t_size, t_length);
  for (size_t t = 0; t < t_length; t++)
    theta_dbl.col(t) = value_of(theta_vec[t]);

  T_partials_mat alpha_dbl(t_size, t_length);
  for (size_t t = 0; t < t_length; t++)
    alpha_dbl.col(t) = value_of(alpha_vec[t]);

  T_partials_return lp(0.0);

  if (include_summand<propto, T_prior_size>::value) {
    lp += (lgamma(alpha_dbl.colwise().sum())
           - lgamma(alpha_dbl).colwise().sum())
              .sum();
  }

  const auto& alpha_m_1 = (alpha_dbl.array() - 1.0);
  const auto& theta_log = theta_dbl.array().log();

  if (include_summand<propto, T_prob, T_prior_size>::value) {
    lp += (theta_log * alpha_m_1).sum();
  }

  operands_and_partials<T_prob, T_prior_size> ops_partials(theta, alpha);
  if (!is_constant_all<T_prob>::value) {
    const auto& theta_deriv = (alpha_m_1 / theta_dbl.array()).matrix();
    for (size_t t = 0; t < t_length; t++)
      ops_partials.edge1_.partials_vec_[t] += theta_deriv.col(t);
  }

  if (!is_constant_all<T_prior_size>::value) {
    for (size_t t = 0; t < t_length; t++)
      ops_partials.edge2_.partials_vec_[t]
          += (digamma(alpha_dbl.col(t).sum())
              - digamma(alpha_dbl).col(t).array() + theta_log.col(t))
                 .matrix();
  }
  return ops_partials.build(lp);
}

template <typename T_prob, typename T_prior_size>
return_type_t<T_prob, T_prior_size> dirichlet_lpmf(const T_prob& theta,
                                                   const T_prior_size& alpha) {
  return dirichlet_lpmf<false>(theta, alpha);
}

}  // namespace math
}  // namespace stan
#endif
