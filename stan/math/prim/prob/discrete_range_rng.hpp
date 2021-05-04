#ifndef STAN_MATH_PRIM_PROB_DISCRETE_RANGE_RNG_HPP
#define STAN_MATH_PRIM_PROB_DISCRETE_RANGE_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return an integer random variate between the given lower and upper bounds
 * (inclusive) using the specified random number generator.
 *
 * `lower` and `upper` can each be a scalar or a one-dimensional container.
 * Any non-scalar inputs must be the same size.
 *
 * @tparam T_lower type of lower bound, either int or std::vector<int>
 * @tparam T_upper type of upper bound, either int or std::vector<int>
 * @tparam RNG type of random number generator
 *
 * @param lower lower bound
 * @param upper upper bound
 * @param rng random number generator
 * @return A (sequence of) integer random variate(s) between `lower` and
 * `upper`, both bounds included.
 * @throw std::domain_error if upper is smaller than lower.
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes.
 */
template <typename T_lower, typename T_upper, class RNG>
inline typename VectorBuilder<true, int, T_lower, T_upper>::type
discrete_range_rng(const T_lower& lower, const T_upper& upper, RNG& rng) {
  static const char* function = "discrete_range_rng";
  using boost::variate_generator;
  using boost::random::uniform_int_distribution;
  check_consistent_sizes(function, "Lower bound parameter", lower,
                         "Upper bound parameter", upper);
  check_greater_or_equal(function, "Upper bound parameter", upper, lower);

  scalar_seq_view<T_lower> lower_vec(lower);
  scalar_seq_view<T_upper> upper_vec(upper);
  size_t N = max_size(lower, upper);
  VectorBuilder<true, int, T_lower, T_upper> output(N);

  for (size_t n = 0; n < N; ++n) {
    variate_generator<RNG&, uniform_int_distribution<>> discrete_range_rng(
        rng, uniform_int_distribution<>(lower_vec[n], upper_vec[n]));

    output[n] = discrete_range_rng();
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
