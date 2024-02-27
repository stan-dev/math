#ifndef STAN_MATH_PRIM_PROB_CATEGORICAL_LOGIT_RNG_HPP
#define STAN_MATH_PRIM_PROB_CATEGORICAL_LOGIT_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/cumulative_sum.hpp>
#include <stan/math/prim/fun/softmax.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Return a draw from a Categorical distribution given a
 * a vector of unnormalized log probabilities and a psuedo-random
 * number generator.
 *
 * This is a convenience wrapper around
 * <code>categorical_rng(softmax(beta), rng)</code>.
 *
 * @tparam RNG Type of pseudo-random number generator.
 * @param beta Vector of unnormalized log probabilities.
 * @param rng Pseudo-random number generator.
 * @return Categorical random variate
 */
template <class RNG>
inline int categorical_logit_rng(const Eigen::VectorXd& beta, RNG& rng) {
  using boost::uniform_01;
  using boost::variate_generator;
  static constexpr const char* function = "categorical_logit_rng";
  check_finite(function, "Log odds parameter", beta);

  variate_generator<RNG&, uniform_01<> > uniform01_rng(rng, uniform_01<>());
  Eigen::VectorXd theta = softmax(beta);
  Eigen::VectorXd index = cumulative_sum(theta);

  double c = uniform01_rng();
  int b = 0;
  while (c > index(b)) {
    b++;
  }
  return b + 1;
}
}  // namespace math
}  // namespace stan
#endif
