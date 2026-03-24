#ifndef STAN_MATH_PRIM_PROB_GEOMETRIC_RNG_HPP
#define STAN_MATH_PRIM_PROB_GEOMETRIC_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <cmath>
#include <utility>
#include <vector>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a geometric random variate with the given success probability
 * parameter, using the given random number generator.
 *
 * theta can be a scalar or a one-dimensional container.
 *
 * The geometric distribution models the number of failures before the first
 * success, with support on {0, 1, 2, ...}. This is the number-of-failures
 * parameterization, consistent with the geometric distribution being a
 * special case of the negative binomial with r = 1.
 *
 * Uses the inverse CDF method: n = floor(log(U) / log(1 - theta))
 * where U ~ Uniform(0, 1).
 *
 * @tparam T_prob type of success probability parameter
 * @tparam RNG type of random number generator
 *
 * @param theta (Sequence of) success probability parameter(s)
 * @param rng random number generator
 * @return (Sequence of) geometric random variate(s)
 * @throw std::domain_error if theta is not in (0, 1]
 */
template <typename T_prob, typename RNG>
inline auto geometric_rng(T_prob&& theta, RNG& rng) {
  static constexpr const char* function = "geometric_rng";
  decltype(auto) theta_ref = to_ref(std::forward<T_prob>(theta));
  check_positive_finite(function, "Success probability parameter",
                        value_of(theta_ref));
  check_less_or_equal(function, "Success probability parameter",
                      value_of(theta_ref), 1.0);

  boost::random::uniform_real_distribution<double> uniform(0.0, 1.0);

  if constexpr (is_stan_scalar_v<T_prob>) {
    if (theta_ref == 1.0) {
      return 0;
    }
    double u = uniform(rng);
    return static_cast<int>(std::floor(std::log(u) / log1m(theta_ref)));
  } else {
    auto theta_arr = as_array_or_scalar(theta_ref);
    std::vector<int> result(theta_arr.size());
    for (int i = 0; i < theta_arr.size(); i++) {
      if (theta_arr[i] == 1.0) {
        result[i] = 0;
      } else {
        double u_i = uniform(rng);
        result[i] = static_cast<int>(
            std::floor(std::log(u_i) / log1m(theta_arr[i])));
      }
    }
    return result;
  }
}

}  // namespace math
}  // namespace stan
#endif
