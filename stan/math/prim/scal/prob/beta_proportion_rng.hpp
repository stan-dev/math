#ifndef STAN_MATH_PRIM_SCAL_PROB_BETA_PROPORTION_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BETA_PROPORTION_RNG_HPP

#include <boost/random/gamma_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/fun/log_sum_exp.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/max_size.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>

namespace stan {
namespace math {

/**
 * Return a Beta random variate parameterized with location and
 * precision parameters using the given random number generator.
 *
 * p and c can each be a scalar or a one-dimensional container. Any
 * non-scalar inputs must be the same size.
 *
 * @tparam T_loc Type of location parameter
 * @tparam T_prec Type of precision parameter
 * @tparam RNG type of random number generator
 * @param p (Sequence of) location parameter(s) in (0, 1)
 * @param c (Sequence of) positive finite precision parameter(s)
 * @param rng random number generator
 * @return (Sequence of) beta random variate(s)
 * @throw std::domain_error if p is outside of (0, 1)
 * @throw std::domain_error if c is nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_loc, typename T_prec, class RNG>
inline typename VectorBuilder<true, double, T_loc, T_prec>::type
beta_proportion_rng(const T_loc &p, const T_prec &c, RNG &rng) {
  static const char *function = "beta_proportion_rng";

  check_positive(function, "Location parameter", p);
  check_less_or_equal(function, "Location parameter", p, 1.0);
  check_positive_finite(function, "Precision parameter", c);
  check_consistent_sizes(function, "Location parameter", p,
                         "Precision parameter", c);

  scalar_seq_view<T_loc> p_vec(p);
  scalar_seq_view<T_prec> c_vec(c);
  size_t N = max_size(p, c);
  VectorBuilder<true, double, T_loc, T_prec> output(N);

  for (size_t n = 0; n < N; ++n) {
    double alpha = p_vec[n] * c_vec[n];
    double beta = (1.0 - p_vec[n]) * c_vec[n];
    output[n] = beta_rng(alpha, beta, rng);
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
