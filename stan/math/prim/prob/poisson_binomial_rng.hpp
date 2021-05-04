#ifndef STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_POISSON_BINOMIAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a pseudorandom Poisson binomial random variable for the given vector
 * of success parameters using the specified random number
 * generator.
 *
 * @tparam RNG class of rng
 * @param theta (Sequence of) chance of success parameter(s)
 * @param rng random number generator
 * @return a Poisson binomial distribution random variable
 * @throw std::domain_error if theta is not a valid probability
 */
template <typename T_theta, typename RNG,
          require_eigen_vt<std::is_arithmetic, T_theta>* = nullptr>
inline int poisson_binomial_rng(const T_theta& theta, RNG& rng) {
  static const char* function = "poisson_binomial_rng";
  check_finite(function, "Probability parameters", theta);
  check_bounded(function, "Probability parameters", value_of(theta), 0.0, 1.0);

  int y = 0;
  for (size_t i = 0; i < theta.size(); ++i) {
    boost::variate_generator<RNG&, boost::bernoulli_distribution<> >
        bernoulli_rng(rng, boost::bernoulli_distribution<>(theta(i)));
    y += bernoulli_rng();
  }

  return y;
}

}  // namespace math
}  // namespace stan
#endif
