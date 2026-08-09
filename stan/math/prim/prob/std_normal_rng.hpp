#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <random>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a standard Normal random variate using the specified
 * random number generator.
 *
 * @tparam RNG type of random number generator
 * @param rng random number generator
 * @return A standard normal random variate
 */
template <class RNG>
inline double std_normal_rng(RNG& rng) {
  std::normal_distribution<> norm_rng(0, 1);

  return norm_rng(rng);
}

}  // namespace math
}  // namespace stan
#endif
