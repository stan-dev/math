#ifndef STAN_MATH_PRIM_MAT_PROB_BERNOULLI_LOGIT_GLM_RNG_HPP
#define STAN_MATH_PRIM_MAT_PROB_BERNOULLI_LOGIT_GLM_RNG_HPP

#include <stan/math/prim/scal/prob/bernoulli_logit_rng.hpp>
#include <stan/math/prim/scal/err/check_consistent_size.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <Eigen/Dense>
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
inline std::vector<int>
bernoulli_logit_glm_rng(const T_x &x,
                        const T_alpha &alpha,
                        const T_beta &beta,
                        RNG& rng) {
  static const char *function = "bernoulli_logit_glm_rng";

  const size_t N = x.col(0).size();
  const size_t M = x.row(0).size();

  check_finite(function, "Weight vector", beta);
  check_finite(function, "Intercept", alpha);
  check_consistent_size(function, "Weight vector", beta, M);
  if (is_vector<T_alpha>::value)
    check_consistent_size(function, "Vector of intercepts", alpha, N);

  Eigen::VectorXd beta_matrix(M);
  for (int i = 0; i < M; ++i)
    beta_matrix[i] = beta[i];
  Eigen::VectorXd alpha_matrix(N);
  scalar_seq_view<T_alpha> alpha_vec(alpha);
  for (int i = 0; i < N; ++i)
    alpha_matrix[i] = alpha_vec[i];

  Eigen::VectorXd theta = (x * beta_matrix) + alpha_matrix;

  check_finite(function, "Matrix of independent variables", theta);

  return bernoulli_logit_rng(theta, rng);
}
}  // namespace math
}  // namespace stan
#endif
