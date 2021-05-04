#ifndef STAN_MATH_PRIM_PROB_DIRICHLET_RNG_HPP
#define STAN_MATH_PRIM_PROB_DIRICHLET_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Return a draw from a Dirichlet distribution with specified
 * parameters and pseudo-random number generator.
 *
 * For prior counts greater than zero, the usual algorithm that
 * draws gamma variates and normalizes is used.
 *
 * For prior counts less than zero (i.e., parameters with value
 * less than one), a log-scale version of the following algorithm
 * is used to deal with underflow:
 *
 * <blockquote>
 * G. Marsaglia and W. Tsang. A simple method for generating gamma
 * variables. ACM Transactions on Mathematical Software.
 * 26(3):363--372, 2000.
 * </blockquote>
 *
 * @tparam RNG type of pseudo-random number generator
 * @param alpha Prior count (plus 1) parameter for Dirichlet.
 * @param rng Pseudo-random number generator.
 */
template <class RNG>
inline Eigen::VectorXd dirichlet_rng(
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha, RNG& rng) {
  using boost::gamma_distribution;
  using boost::variate_generator;
  using boost::random::uniform_real_distribution;
  using Eigen::VectorXd;
  using std::exp;
  using std::log;

  // separate algorithm if any parameter is less than 1
  if (alpha.minCoeff() < 1) {
    variate_generator<RNG&, uniform_real_distribution<> > uniform_rng(
        rng, uniform_real_distribution<>(0.0, 1.0));
    VectorXd log_y(alpha.size());
    for (int i = 0; i < alpha.size(); ++i) {
      variate_generator<RNG&, gamma_distribution<> > gamma_rng(
          rng, gamma_distribution<>(alpha(i) + 1, 1));
      double log_u = log(uniform_rng());
      log_y(i) = log(gamma_rng()) + log_u / alpha(i);
    }
    double log_sum_y = log_sum_exp(log_y);
    VectorXd theta(alpha.size());
    for (int i = 0; i < alpha.size(); ++i) {
      theta(i) = exp(log_y(i) - log_sum_y);
    }
    return theta;
  }

  // standard normalized gamma algorithm
  Eigen::VectorXd y(alpha.rows());
  for (int i = 0; i < alpha.rows(); i++) {
    variate_generator<RNG&, gamma_distribution<> > gamma_rng(
        rng, gamma_distribution<>(alpha(i, 0), 1e-7));
    y(i) = gamma_rng();
  }
  return y / y.sum();
}

}  // namespace math
}  // namespace stan
#endif
