#ifndef STAN_MATH_PRIM_PROB_BERNOULLI_LOGIT_GLM_RNG_HPP
#define STAN_MATH_PRIM_PROB_BERNOULLI_LOGIT_GLM_RNG_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup multivar_dists
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
    const T_x& x, const T_alpha& alpha, const T_beta& beta, RNG& rng) {
  using boost::bernoulli_distribution;
  using boost::variate_generator;
  using T_x_ref = ref_type_t<T_x>;
  using T_alpha_ref = ref_type_t<T_alpha>;
  using T_beta_ref = ref_type_t<T_beta>;

  const size_t N = x.cols();
  const size_t M = x.rows();

  static constexpr const char* function = "bernoulli_logit_glm_rng";
  check_consistent_size(function, "Weight vector", beta, N);
  check_consistent_size(function, "Vector of intercepts", alpha, M);
  T_x_ref x_ref = x;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_finite(function, "Matrix of independent variables", x_ref);
  check_finite(function, "Weight vector", beta_ref);
  check_finite(function, "Intercept", alpha_ref);

  const auto& beta_vector = as_column_vector_or_scalar(beta_ref);

  Eigen::VectorXd x_beta;
  if (is_vector<T_beta>::value) {
    x_beta = x_ref * beta_vector;
  } else {
    x_beta = (x_ref.array() * forward_as<double>(beta_vector)).rowwise().sum();
  }

  scalar_seq_view<T_alpha> alpha_vec(alpha_ref);

  VectorBuilder<true, int, T_alpha> output(M);

  for (size_t m = 0; m < M; ++m) {
    double theta_m = alpha_vec[m] + x_beta(m);
    variate_generator<RNG&, bernoulli_distribution<>> bernoulli_rng(
        rng, bernoulli_distribution<>(inv_logit(theta_m)));
    output[m] = bernoulli_rng();
  }

  return output.data();
}
}  // namespace math
}  // namespace stan
#endif
