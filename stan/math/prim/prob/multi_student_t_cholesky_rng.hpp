#ifndef STAN_MATH_PRIM_PROB_MULTI_STUDENT_T_CHOLESKY_RNG_HPP
#define STAN_MATH_PRIM_PROB_MULTI_STUDENT_T_CHOLESKY_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <stan/math/prim/prob/inv_gamma_rng.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Return a multivariate student-t random variate with the given degrees of
 * freedom location and Cholesky factor the scale matrix
 * using the specified random number generator.
 *
 * mu can be either an Eigen::VectorXd, an Eigen::RowVectorXd, or a std::vector
 * of either of those types.
 *
 * @tparam t_loc Type of location parameter
 * @tparam rng type of pseudo-random number generator
 * @param nu a scalar indicating the degrees of freedom parameter
 * @param mu An Eigen::VectorXd, Eigen::RowVectorXd, or std::vector
 *    of location values for the multivariate student t
 * @param L cholesky factor of the scale matrix
 * @param rng random number generator
 * @return eigen vector of multivariate student t random variates
 *    with the given nu, mu, L
 * @throw std::domain_error if L is not a Cholesky factor, any value in mu is
 *    not finite, nu is not positive, or nu is NaN
 * @throw std::invalid_argument if the length of (each) mu is not equal to the
 *    number of rows and columns in L
 */
template <typename T_loc, class RNG>
inline typename StdVectorBuilder<true, Eigen::VectorXd, T_loc>::type
multi_student_t_cholesky_rng(double nu, const T_loc& mu,
                             const Eigen::MatrixXd& L, RNG& rng) {
  using boost::normal_distribution;
  using boost::variate_generator;
  using boost::random::gamma_distribution;

  static constexpr const char* function = "multi_student_t_cholesky_rng";
  check_not_nan(function, "Degrees of freedom parameter", nu);
  check_positive(function, "Degrees of freedom parameter", nu);
  check_positive(function, "Scale matrix rows", L.rows());
  vector_seq_view<T_loc> mu_vec(mu);

  size_t N = size_mvt(mu);
  for (size_t i = 1; i < N; i++) {
    check_size_match(function,
                     "Size of one of the vectors of "
                     "the location variable",
                     mu_vec[i].size(),
                     "Size of the first vector of the "
                     "location variable",
                     mu_vec[i - 1].size());
  }

  check_size_match(function, "Size of random variable", mu_vec[0].size(),
                   "rows of scale parameter", L.rows());

  for (size_t i = 0; i < N; i++) {
    check_finite(function, "Location parameter", mu_vec[i]);
  }
  const auto& L_ref = to_ref(L);
  check_cholesky_factor(function, "L matrix", L_ref);

  StdVectorBuilder<true, Eigen::VectorXd, T_loc> output(N);

  variate_generator<RNG&, normal_distribution<> > std_normal_rng(
      rng, normal_distribution<>(0, 1));

  double w = inv_gamma_rng(nu / 2, nu / 2, rng);
  for (size_t n = 0; n < N; ++n) {
    Eigen::VectorXd z(L.cols());
    for (int i = 0; i < L.cols(); i++) {
      z(i) = std_normal_rng();
    }
    z *= std::sqrt(w);
    output[n] = as_column_vector_or_scalar(mu_vec[n]) + L * z;
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
