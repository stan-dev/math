#ifndef STAN_MATH_PRIM_PROB_STUDENT_T_RNG_HPP
#define STAN_MATH_PRIM_PROB_STUDENT_T_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <boost/random/student_t_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a student-t random variate for the given degrees of freedom,
 * location, and scale using the specified random number generator.
 *
 * nu, mu, and sigma can each be a scalar or a one-dimensional container. Any
 * non-scalar inputs must be the same size.
 *
 * @tparam T_deg type of degrees of freedom parameter
 * @tparam T_loc type of location parameter
 * @tparam T_scale type of scale parameter
 * @tparam RNG type of random number generator
 *
 * @param nu (Sequence of) degrees of freedom parameter(s)
 * @param mu (Sequence of) location parameter(s)
 * @param sigma (Sequence of) scale parameter(s)
 * @param rng random number generator
 * @return Student-t random variate
 * @throw std::domain_error if nu is nonpositive, mu is infinite, or sigma
 * is nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_deg, typename T_loc, typename T_scale, class RNG>
inline typename VectorBuilder<true, double, T_deg, T_loc, T_scale>::type
student_t_rng(const T_deg& nu, const T_loc& mu, const T_scale& sigma,
              RNG& rng) {
  using T_nu_ref = ref_type_t<T_deg>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_sigma_ref = ref_type_t<T_scale>;
  using boost::variate_generator;
  using boost::random::student_t_distribution;
  static constexpr const char* function = "student_t_rng";
  check_consistent_sizes(function, "Degrees of freedom parameter", nu,
                         "Location parameter", mu, "Scale Parameter", sigma);
  T_nu_ref nu_ref = nu;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  check_positive_finite(function, "Degrees of freedom parameter", nu_ref);
  check_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Scale parameter", sigma_ref);

  scalar_seq_view<T_nu_ref> nu_vec(nu_ref);
  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_sigma_ref> sigma_vec(sigma_ref);
  size_t N = max_size(nu, mu, sigma);
  VectorBuilder<true, double, T_deg, T_loc, T_scale> output(N);

  for (size_t n = 0; n < N; ++n) {
    variate_generator<RNG&, student_t_distribution<> > rng_unit_student_t(
        rng, student_t_distribution<>(nu_vec[n]));
    output[n] = mu_vec[n] + sigma_vec[n] * rng_unit_student_t();
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
