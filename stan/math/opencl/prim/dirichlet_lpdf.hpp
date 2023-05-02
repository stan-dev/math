#ifndef STAN_MATH_OPENCL_PRIM_DIRICHLET_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_DIRICHLET_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
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
template <bool propto, typename T_prob_cl, typename T_prior_size_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_prob_cl, T_prior_size_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_prob_cl, T_prior_size_cl>* = nullptr>
inline return_type_t<T_prob_cl, T_prior_size_cl> dirichlet_lpdf(
    const T_prob_cl& theta, const T_prior_size_cl& alpha) {
  static const char* function = "dirichlet_lpdf(OpenCL)";

  check_consistent_sizes(function, "probabilities", theta, "prior sample sizes",
                         alpha);

  if (max_size(theta, alpha) == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_prob_cl, T_prior_size_cl>::value) {
    return 0.0;
  }

  const auto& theta_val = value_of(theta);
  const auto& alpha_val = value_of(alpha);

  auto check_alpha_positive
      = check_cl(function, "prior sample sizes", alpha_val, "positive");
  auto alpha_positive = alpha_val > 0.0;
  auto check_theta_nonnegative
      = check_cl(function, "probabilities", theta_val, "nonnegative");
  auto theta_nonnegative = theta_val >= 0.0;

  auto theta_csum = colwise_sum(theta_val);
  auto alpha_m_1 = rowwise_optional_broadcast(alpha_val - 1.0);
  auto theta_log = rowwise_optional_broadcast(log(theta_val));

  matrix_cl<double> theta_csum_cl;
  matrix_cl<double> alpha_csum_cl;
  matrix_cl<double> lgamma_alpha_csum_cl;
  matrix_cl<double> theta_log_alpha_m_1_sum_cl;
  matrix_cl<double> theta_deriv_cl;
  matrix_cl<double> alpha_deriv_cl;
  if (theta.cols() > alpha.cols()) {
    auto alpha_csum = colwise_sum(alpha_val);
    auto lgamma_alpha_csum = colwise_sum(lgamma(alpha_val));
    matrix_cl<double> digamma_alpha_cl(alpha.rows(), alpha.cols());
    results(check_alpha_positive, alpha_csum_cl, lgamma_alpha_csum_cl,
            digamma_alpha_cl)
        = expressions(
            alpha_positive,
            calc_if<include_summand<propto, T_prior_size_cl>::value>(
                alpha_csum),
            calc_if<include_summand<propto, T_prior_size_cl>::value>(
                lgamma_alpha_csum),
            calc_if<!is_constant<T_prior_size_cl>::value>(digamma(alpha_val)));

    auto theta_log_alpha_m_1_sum = sum_2d(elt_multiply(theta_log, alpha_m_1));

    auto theta_deriv
        = elt_divide(alpha_m_1, rowwise_optional_broadcast(theta_val));
    auto alpha_deriv = theta_log - rowwise_optional_broadcast(digamma_alpha_cl);

    results(check_theta_nonnegative, theta_csum_cl, theta_log_alpha_m_1_sum_cl,
            theta_deriv_cl, alpha_deriv_cl)
        = expressions(
            theta_nonnegative, theta_csum,
            calc_if<include_summand<propto, T_prob_cl, T_prior_size_cl>::value>(
                theta_log_alpha_m_1_sum),
            calc_if<!is_constant<T_prob_cl>::value>(theta_deriv),
            calc_if<!is_constant<T_prior_size_cl>::value>(alpha_deriv));

    if (include_summand<propto, T_prior_size_cl>::value) {
      matrix_cl<double> alpha_csum_cl2;
      matrix_cl<double> lgamma_alpha_csum_cl2;
      while (alpha_csum_cl.rows() > 1) {
        results(alpha_csum_cl2, lgamma_alpha_csum_cl2) = expressions(
            calc_if<include_summand<propto, T_prior_size_cl>::value>(
                colwise_sum(alpha_csum_cl)),
            calc_if<include_summand<propto, T_prior_size_cl>::value>(
                colwise_sum(lgamma_alpha_csum_cl)));
        alpha_csum_cl = std::move(alpha_csum_cl2);
        lgamma_alpha_csum_cl = std::move(lgamma_alpha_csum_cl2);
      }
    }
    while (theta_csum_cl.rows() > 1) {
      theta_csum_cl = colwise_sum(theta_csum_cl).eval();
    }
  } else {
    auto alpha_csum = colwise_sum(alpha_val);
    auto lgamma_alpha_csum = colwise_sum(lgamma(alpha_val));
    if (alpha.cols() > theta.cols()) {
      matrix_cl<double> log_theta_cl(theta.rows(), theta.cols());
      results(check_theta_nonnegative, theta_csum_cl, log_theta_cl)
          = expressions(
              theta_nonnegative, theta_csum,
              calc_if<
                  include_summand<propto, T_prob_cl, T_prior_size_cl>::value>(
                  theta_log));

      auto log_theta_bc = rowwise_optional_broadcast(log_theta_cl);
      auto theta_log_alpha_m_1_sum
          = sum_2d(elt_multiply(log_theta_bc, alpha_m_1));

      auto theta_deriv
          = elt_divide(alpha_m_1, rowwise_optional_broadcast(theta_val));
      auto alpha_deriv
          = log_theta_bc - rowwise_optional_broadcast(digamma(alpha_val));

      results(check_alpha_positive, alpha_csum_cl, lgamma_alpha_csum_cl,
              theta_log_alpha_m_1_sum_cl, theta_deriv_cl, alpha_deriv_cl)
          = expressions(
              alpha_positive,
              calc_if<include_summand<propto, T_prior_size_cl>::value>(
                  alpha_csum),
              calc_if<include_summand<propto, T_prior_size_cl>::value>(
                  lgamma_alpha_csum),
              calc_if<
                  include_summand<propto, T_prob_cl, T_prior_size_cl>::value>(
                  theta_log_alpha_m_1_sum),
              calc_if<!is_constant<T_prob_cl>::value>(theta_deriv),
              calc_if<!is_constant<T_prior_size_cl>::value>(alpha_deriv));

      while (alpha_csum_cl.rows() > 1) {
        matrix_cl<double> alpha_csum_cl2;
        matrix_cl<double> lgamma_alpha_csum_cl2;
        results(alpha_csum_cl2, lgamma_alpha_csum_cl2) = expressions(
            calc_if<include_summand<propto, T_prior_size_cl>::value>(
                colwise_sum(alpha_csum_cl)),
            calc_if<include_summand<propto, T_prior_size_cl>::value>(
                colwise_sum(lgamma_alpha_csum_cl)));
        if (include_summand<propto, T_prior_size_cl>::value) {
          alpha_csum_cl = std::move(alpha_csum_cl2);
          lgamma_alpha_csum_cl = std::move(lgamma_alpha_csum_cl2);
        }
      }
      double theta_sum = sum(from_matrix_cl(theta_csum));
      check_cl(function, "sum of probabilities", theta_sum, "equal to 1")
          = (fabs(theta_sum - 1.0) <= CONSTRAINT_TOLERANCE);
    } else {
      auto theta_log_alpha_m_1_sum = sum_2d(elt_multiply(theta_log, alpha_m_1));

      auto theta_deriv
          = elt_divide(alpha_m_1, rowwise_optional_broadcast(theta_val));
      auto alpha_deriv
          = theta_log - rowwise_optional_broadcast(digamma(alpha_val));

      results(check_alpha_positive, check_theta_nonnegative, theta_csum_cl,
              alpha_csum_cl, lgamma_alpha_csum_cl, theta_log_alpha_m_1_sum_cl,
              theta_deriv_cl, alpha_deriv_cl)
          = expressions(
              alpha_positive, theta_nonnegative, theta_csum,
              calc_if<include_summand<propto, T_prior_size_cl>::value>(
                  alpha_csum),
              calc_if<include_summand<propto, T_prior_size_cl>::value>(
                  lgamma_alpha_csum),
              calc_if<
                  include_summand<propto, T_prob_cl, T_prior_size_cl>::value>(
                  theta_log_alpha_m_1_sum),
              calc_if<!is_constant<T_prob_cl>::value>(theta_deriv),
              calc_if<!is_constant<T_prior_size_cl>::value>(alpha_deriv));

      while (theta_csum_cl.rows() > 1) {
        matrix_cl<double> theta_csum_cl2;
        matrix_cl<double> alpha_csum_cl2;
        matrix_cl<double> lgamma_alpha_csum_cl2;
        results(theta_csum_cl2, alpha_csum_cl2, lgamma_alpha_csum_cl2)
            = expressions(
                colwise_sum(theta_csum_cl),
                calc_if<include_summand<propto, T_prior_size_cl>::value>(
                    colwise_sum(alpha_csum_cl)),
                calc_if<include_summand<propto, T_prior_size_cl>::value>(
                    colwise_sum(lgamma_alpha_csum_cl)));
        theta_csum_cl = std::move(theta_csum_cl2);
        if (include_summand<propto, T_prior_size_cl>::value) {
          alpha_csum_cl = std::move(alpha_csum_cl2);
          lgamma_alpha_csum_cl = std::move(lgamma_alpha_csum_cl2);
        }
      }
    }
  }

  if (theta.cols() >= alpha.cols()) {
    // transpose is there just because working on col vectors is more efficient
    // than on row vectors with kernel generator
    check_cl(function, "sum of probabilities", transpose(theta_csum_cl),
             "equal to 1")
        = (fabs(transpose(theta_csum_cl) - 1.0) <= CONSTRAINT_TOLERANCE);
  }

  double lp = 0.0;

  if (include_summand<propto, T_prior_size_cl>::value) {
    if (theta.cols() > alpha.cols()) {
      lp += (lgamma(from_matrix_cl(alpha_csum_cl)) * theta.cols()
             - from_matrix_cl(lgamma_alpha_csum_cl) * theta.cols())
                .sum();
    } else {
      lp += (lgamma(from_matrix_cl(alpha_csum_cl))
             - from_matrix_cl(lgamma_alpha_csum_cl))
                .sum();
    }
  }
  if (include_summand<propto, T_prob_cl, T_prior_size_cl>::value) {
    lp += from_matrix_cl(theta_log_alpha_m_1_sum_cl).sum();
  }

  auto ops_partials = make_partials_propagator(theta, alpha);

  if (!is_constant<T_prob_cl>::value) {
    if (theta.cols() < alpha.cols()) {
      partials<0>(ops_partials) = rowwise_sum(theta_deriv_cl);
    } else {
      partials<0>(ops_partials) = std::move(theta_deriv_cl);
    }
  }
  if (!is_constant<T_prior_size_cl>::value) {
    if (theta.cols() > alpha.cols()) {
      matrix_cl<double> tmp_cl
          = digamma(alpha_csum_cl) * static_cast<double>(theta.cols());
      partials<1>(ops_partials)
          = colwise_broadcast(tmp_cl) + rowwise_sum(alpha_deriv_cl);
    } else {
      matrix_cl<double> tmp_cl = digamma(alpha_csum_cl);
      partials<1>(ops_partials) = colwise_broadcast(tmp_cl) + alpha_deriv_cl;
    }
  }
  return ops_partials.build(lp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
