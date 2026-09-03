#ifndef STAN_MATH_PRIM_PROB_INV_GAUSSIAN_RNG_HPP
#define STAN_MATH_PRIM_PROB_INV_GAUSSIAN_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return an inverse Gaussian random variate for the given mean and shape
 * using the specified random number generator.
 *
 * mu and lambda can each be a scalar or a one-dimensional container. Any
 * non-scalar inputs must be the same size.
 *
 * <p>The algorithm used in inv_gaussian_rng is the transformation in:
 *
 * Generating Random Variates Using Transformations with Multiple Roots
 * J. R. Michael, W. R. Schucany and R. W. Haas
 * The American Statistician, Vol. 30, No. 2 (1976), pp. 88-90
 *
 * <p>A chi-square variate \f$w\f$ on one degree of freedom gives a pair of
 * candidate roots and a Bernoulli draw selects between them. The method is
 * exact, with no rejection step.
 *
 * <p>The smaller root is computed in reciprocal form; since
 * \f$(1 + u/2)^2 - (u + u^2/4) = 1\f$ exactly, it subtracts nothing and
 * stays accurate for \f$u = \mu w / \lambda\f$ up to \f$10^{20}\f$.
 *
 * @tparam T_loc type of mean parameter
 * @tparam T_shape type of shape parameter
 * @tparam RNG type of random number generator
 * @param mu (Sequence of) mean parameter(s)
 * @param lambda (Sequence of) shape parameter(s)
 * @param rng random number generator
 * @return (Sequence of) inverse Gaussian random variate(s)
 * @throw std::domain_error if mu or lambda is not positive and finite
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_loc, typename T_shape, class RNG>
inline typename VectorBuilder<true, double, T_loc, T_shape>::type
inv_gaussian_rng(const T_loc& mu, const T_shape& lambda, RNG& rng) {
  using boost::variate_generator;
  using boost::random::normal_distribution;
  using boost::random::uniform_01;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_lambda_ref = ref_type_t<T_shape>;
  static constexpr const char* function = "inv_gaussian_rng";
  check_consistent_sizes(function, "Mean parameter", mu, "Shape parameter",
                         lambda);
  T_mu_ref mu_ref = mu;
  T_lambda_ref lambda_ref = lambda;
  check_positive_finite(function, "Mean parameter", mu_ref);
  check_positive_finite(function, "Shape parameter", lambda_ref);

  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_lambda_ref> lambda_vec(lambda_ref);
  size_t N = max_size(mu, lambda);
  VectorBuilder<true, double, T_loc, T_shape> output(N);

  variate_generator<RNG&, normal_distribution<> > norm_rng(
      rng, normal_distribution<>(0, 1));
  variate_generator<RNG&, uniform_01<> > uniform01_rng(rng, uniform_01<>());

  for (size_t n = 0; n < N; ++n) {
    const double mu_dbl = mu_vec[n];
    const double lambda_dbl = lambda_vec[n];
    const double nu = norm_rng();
    const double w = nu * nu;
    const double u = mu_dbl * w / lambda_dbl;
    const double x = mu_dbl / (1.0 + 0.5 * u + std::sqrt(u + 0.25 * u * u));
    output[n]
        = (uniform01_rng() <= mu_dbl / (mu_dbl + x)) ? x : mu_dbl * mu_dbl / x;
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
