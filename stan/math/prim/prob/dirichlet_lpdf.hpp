#ifndef STAN_MATH_PRIM_PROB_DIRICHLET_LPDF_HPP
#define STAN_MATH_PRIM_PROB_DIRICHLET_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size_mvt.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
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
 * @tparam T_prob type of scalar
 * @tparam T_prior_size type of prior sample sizes
 * @param theta A scalar vector.
 * @param alpha Prior sample sizes.
 * @return The log of the Dirichlet density.
 * @throw std::domain_error if any element of alpha is less than
 * or equal to 0.
 * @throw std::domain_error if any element of theta is less than 0.
 * @throw std::domain_error if the sum of theta is not 1.
 */
template <bool propto, typename T_prob, typename T_prior_size,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_prob, T_prior_size>* = nullptr>
return_type_t<T_prob, T_prior_size> dirichlet_lpdf(const T_prob& theta,
                                                   const T_prior_size& alpha) {
  using T_partials_return = partials_return_t<T_prob, T_prior_size>;
  using T_partials_array = typename Eigen::Array<T_partials_return, -1, -1>;
  using T_theta_ref = ref_type_t<T_prob>;
  using T_alpha_ref = ref_type_t<T_prior_size>;
  static constexpr const char* function = "dirichlet_lpdf";

  T_theta_ref theta_ref = theta;
  T_alpha_ref alpha_ref = alpha;
  vector_seq_view<T_theta_ref> theta_vec(theta_ref);
  vector_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  const size_t t_length = max_size_mvt(theta, alpha);

  check_consistent_sizes(function, "probabilities", theta_vec[0],
                         "prior sample sizes", alpha_vec[0]);

  for (size_t t = 0; t < t_length; t++) {
    check_positive(function, "prior sample sizes", alpha_vec[t]);
    check_simplex(function, "probabilities", theta_vec[t]);
  }

  const size_t t_size = theta_vec[0].size();

  T_partials_array theta_dbl(t_size, t_length);
  for (size_t t = 0; t < t_length; t++) {
    theta_dbl.col(t) = theta_vec.val(t);
  }
  T_partials_array alpha_dbl(t_size, t_length);
  for (size_t t = 0; t < t_length; t++) {
    alpha_dbl.col(t) = alpha_vec.val(t);
  }

  T_partials_return lp(0.0);

  if (include_summand<propto, T_prior_size>::value) {
    lp += (lgamma(alpha_dbl.colwise().sum())
           - lgamma(alpha_dbl).colwise().sum())
              .sum();
  }

  const auto& alpha_m_1
      = to_ref_if<!is_constant_all<T_prob>::value>(alpha_dbl - 1.0);
  const auto& theta_log
      = to_ref_if<!is_constant_all<T_prior_size>::value>(theta_dbl.log());

  if (include_summand<propto, T_prob, T_prior_size>::value) {
    lp += (theta_log * alpha_m_1).sum();
  }

  auto ops_partials = make_partials_propagator(theta_ref, alpha_ref);
  if (!is_constant_all<T_prob>::value) {
    for (size_t t = 0; t < t_length; t++) {
      partials_vec<0>(ops_partials)[t]
          += (alpha_m_1.col(t) / theta_dbl.col(t)).matrix();
    }
  }

  if (!is_constant_all<T_prior_size>::value) {
    for (size_t t = 0; t < t_length; t++) {
      partials_vec<1>(ops_partials)[t]
          += (digamma(alpha_dbl.col(t).sum()) - digamma(alpha_dbl.col(t))
              + theta_log.col(t))
                 .matrix();
    }
  }
  return ops_partials.build(lp);
}

template <typename T_prob, typename T_prior_size>
return_type_t<T_prob, T_prior_size> dirichlet_lpdf(const T_prob& theta,
                                                   const T_prior_size& alpha) {
  return dirichlet_lpdf<false>(theta, alpha);
}

}  // namespace math
}  // namespace stan
#endif
