#ifndef STAN_MATH_PRIM_PROB_CHI_SQUARE_RNG_HPP
#define STAN_MATH_PRIM_PROB_CHI_SQUARE_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <boost/random/chi_squared_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a chi squared random variate with nu degrees of
 * freedom using the specified random number generator.
 *
 * nu can be a scalar or a one-dimensional container.
 *
 * @tparam T_deg type of degrees of freedom parameter
 * @tparam RNG class of random number generator
 * @param nu (Sequence of) positive degrees of freedom parameter(s)
 * @param rng random number generator
 * @return (Sequence of) chi squared random variate(s)
 * @throw std::domain_error if nu is nonpositive
 */
template <typename T_deg, class RNG>
inline typename VectorBuilder<true, double, T_deg>::type chi_square_rng(
    const T_deg& nu, RNG& rng) {
  using boost::variate_generator;
  using boost::random::chi_squared_distribution;
  using T_nu_ref = ref_type_t<T_deg>;
  static const char* function = "chi_square_rng";
  T_nu_ref nu_ref = nu;
  check_positive_finite(function, "Degrees of freedom parameter", nu_ref);

  scalar_seq_view<T_nu_ref> nu_vec(nu_ref);
  size_t N = stan::math::size(nu);
  VectorBuilder<true, double, T_deg> output(N);

  for (size_t n = 0; n < N; ++n) {
    variate_generator<RNG&, chi_squared_distribution<> > chi_square_rng(
        rng, chi_squared_distribution<>(nu_vec[n]));
    output[n] = chi_square_rng();
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
