#ifndef STAN_MATH_PRIM_MAT_PROB_BERNOULLI_LOGIT_GLM_RNG_HPP
#define STAN_MATH_PRIM_MAT_PROB_BERNOULLI_LOGIT_GLM_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_consistent_size.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/fun/inv_logit.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns a draw from the Generalized Linear Model (GLM)
 * with Bernoulli distribution and logit link function.
 *
 * This is a convenience wrapper around
 * <code>bernoulli_logit_rng(alpha + x * beta, rng)</code>.
 * @tparam T_x type of the matrix of independent variables (features);
 * this should be an Eigen::Matrix whose number of columns should
 * match the length of beta; the number of rows is the number of
 * samples.
 * @tparam T_alpha type of the intercept(s); this can be a vector of
 * the same length as y) of intercepts or a single value (for models
 * with constant intercept); if a vector its length should match x's
 * row-count;
 * @tparam T_beta type of the weight vector;
 * @tparam RNG Type of pseudo-random number generator.
 * @param x design matrix
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @param rng Pseudo-random number generator.
 * @return Bernoulli logit glm random variate
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <typename T_x, typename T_alpha, typename T_beta, class RNG>
inline typename VectorBuilder<true, int, T_alpha>::type bernoulli_logit_glm_rng(
    const T_x &x, const T_alpha &alpha, const T_beta &beta, RNG &rng) {
  using boost::bernoulli_distribution;
  using boost::variate_generator;

  static const char *function = "bernoulli_logit_glm_rng";

  const size_t N = x.row(0).size();
  const size_t M = x.col(0).size();

  check_finite(function, "Matrix of independent variables", x);
  check_finite(function, "Weight vector", beta);
  check_finite(function, "Intercept", alpha);
  check_consistent_size(function, "Weight vector", beta, N);
  check_consistent_size(function, "Vector of intercepts", alpha, M);

  scalar_seq_view<T_beta> beta_vec(beta);
  Eigen::VectorXd beta_vector(N);
  for (int i = 0; i < N; ++i) {
    beta_vector[i] = beta_vec[i];
  }

  Eigen::VectorXd x_beta = x * beta_vector;

  scalar_seq_view<T_alpha> alpha_vec(alpha);

  VectorBuilder<true, int, T_alpha> output(M);

  for (size_t m = 0; m < M; ++m) {
    double theta_m = alpha_vec[m] + x_beta(m);
    variate_generator<RNG &, bernoulli_distribution<>> bernoulli_rng(
        rng, bernoulli_distribution<>(inv_logit(theta_m)));
    output[m] = bernoulli_rng();
  }

  return output.data();
}
}  // namespace math
}  // namespace stan
#endif
