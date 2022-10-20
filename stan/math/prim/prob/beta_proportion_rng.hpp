#ifndef STAN_MATH_PRIM_PROB_BETA_PROPORTION_RNG_HPP
#define STAN_MATH_PRIM_PROB_BETA_PROPORTION_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a Beta random variate specified probability, location, and
 * precision parameters: beta_proportion_rng(y | mu, kappa) =
 * beta_rng(y | mu * kappa, (1 - mu) * kappa).  Any arguments other
 * than scalars must be containers of the same size.  With non-scalar
 * arguments, the return is a container matching the size of the
 * arguments with scalars broadcast as necessary.
 *
 * @tparam T_loc Type of location parameter
 * @tparam T_prec Type of precision parameter
 * @tparam RNG type of random number generator
 * @param mu (Sequence of) location parameter(s) in (0, 1)
 * @param kappa (Sequence of) positive finite precision parameter(s)
 * @param rng random number generator
 * @return (Sequence of) beta random variate(s)
 * @throw std::domain_error if mu is outside of (0, 1)
 * @throw std::domain_error if kappa is nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_loc, typename T_prec, class RNG>
inline typename VectorBuilder<true, double, T_loc, T_prec>::type
beta_proportion_rng(const T_loc &mu, const T_prec &kappa, RNG &rng) {
  using T_mu_ref = ref_type_t<T_loc>;
  using T_kappa_ref = ref_type_t<T_prec>;
  static constexpr const char *function = "beta_proportion_rng";
  check_consistent_sizes(function, "Location parameter", mu,
                         "Precision parameter", kappa);
  T_mu_ref mu_ref = mu;
  T_kappa_ref kappa_ref = kappa;
  check_positive(function, "Location parameter", mu_ref);
  check_less(function, "Location parameter", mu_ref, 1.0);
  check_positive_finite(function, "Precision parameter", kappa_ref);

  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_kappa_ref> kappa_vec(kappa_ref);
  size_t N = max_size(mu, kappa);
  VectorBuilder<true, double, T_loc, T_prec> output(N);

  for (size_t n = 0; n < N; ++n) {
    double alpha = mu_vec[n] * kappa_vec[n];
    double beta = kappa_vec[n] - alpha;
    output[n] = beta_rng(alpha, beta, rng);
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
