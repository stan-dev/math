#ifndef STAN_MATH_PRIM_PROB_DIRICHLET_MULTINOMIAL_LPMF_HPP
#define STAN_MATH_PRIM_PROB_DIRICHLET_MULTINOMIAL_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <vector>


namespace stan {
namespace math {

/** \ingroup multivar_dists
 * The log of the Dirichlet-Multinomial probability for
 * the given integer vector \f$n\f$ and a vector of prior sample
 * sizes, \f$\alpha\f$.
 * Each element of \f$\alpha\f$ must be greater than 0.
 * Each element of \f$n\f$ must be greater than or equal to 0.
 * 
 * Suppose that \f$n = (n_1, \ldots, n_K)\f$ has a Dirichlet-multinomial
 * distribution with prior sample size
 * \f$\alpha = (\alpha_1,\ldots,\alpha_K)\f$ (also called the intensity):
 *
 * \f[
 * (n_1,\ldots,n_K)\sim\mbox{DirMult}(\alpha_1,\ldots,\alpha_K)
 * \f]
 * 
 * Write \f$N = n_1 + \cdots + n_K\f$ and 
 * \f$\alpha_0 = \alpha_1 + \cdots + \alpha_K\f$,
 * then the log probability mass function is given by
 * 
 * \f[
 * \log(p(n_1, \ldots, n_K\,|\,\alpha_1,\ldots,\alpha_K))=\log\left(
 * \frac{\Gamma(\alpha_0)\Gamma(N+1)}{\Gamma(\alpha_0 + N)}
 * \prod_{k=1}^K \frac{\Gamma(n_k + \alpha_k)}{\Gamma(\alpha_k) 
 * \Gamma(n_k+1)} \right)\\
 * = \log(N) + \log(B(\alpha_0, N)) - 
 * \sum_{k : n_k > 0} \bigl(\log(n_k) + \log(B(\alpha_k, n_k))\bigr)
 * \f]
 * 
 * The second identity is only valid for \f$N > 0\f$. For \f$N=0\f$, 
 * we have \f$\log(p(n\,|\,\alpha)) = 0\f$.
 *
 * @tparam T_prior_size type of prior sample sizes
 * @param ns A vector of integers.
 * @param alpha Prior sample sizes (or intensity vector).
 * @return The log of the Dirichlet-Multinomial probability.
 * @throw std::domain_error if any element of alpha is less than
 * or equal to 0, or infinite.
 * @throw std::domain_error any element of ns is less than 0.
 */

template <bool propto, typename T_prior_size,
          require_eigen_col_vector_t<T_prior_size>* = nullptr>
return_type_t<T_prior_size> dirichlet_multinomial_lpmf(
    const std::vector<int>& ns, const T_prior_size& alpha) {
  static const char* function = "dirichlet_multinomial_lpmf";
  check_size_match(function, "Size of number of trials variable", ns.size(),
                   "rows of prior size parameter", alpha.rows());
  check_nonnegative(function, "Number of trials variable", ns);
  const auto& alpha_ref = to_ref(alpha);
  check_positive_finite(function, "Prior size parameter", alpha_ref);

  return_type_t<T_prior_size> lp(0.0);

  if (!include_summand<propto, T_prior_size>::value) return lp;

  int sum_ns = 0;

  for (size_t i = 0; i < ns.size(); ++i) {
    if (ns[i] > 0) {
      lp -= lbeta(alpha_ref.coeff(i), ns[i]);
      sum_ns += ns[i];
    }
  }

  if (sum_ns == 0) return lp;

  value_type_t<T_prior_size> sum_alpha = alpha_ref.sum();

  lp += lbeta(sum_alpha, sum_ns);

  if (include_summand<propto>::value) {
    for (size_t i = 0; i < ns.size(); ++i) {
      if (ns[i] > 0) {
        lp -= log(ns[i]);
      }
    }
    lp += log(sum_ns);
  }

  return lp;
}

template <typename T_prior_size>
return_type_t<T_prior_size> dirichlet_multinomial_lpmf(
    const std::vector<int>& ns, const T_prior_size& alpha) {
  return dirichlet_multinomial_lpmf<false>(ns, alpha);
}

}  // namespace math
}  // namespace stan
#endif
