#ifndef STAN_MATH_PRIM_PROB_MULTI_STUDENT_T_RNG_HPP
#define STAN_MATH_PRIM_PROB_MULTI_STUDENT_T_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Return a multivariate student-t random variate with the given degrees of
 * freedom location and covariance using the specified random number generator.
 *
 * mu can be either an Eigen::VectorXd, an Eigen::RowVectorXd, or a std::vector
 * of either of those types.
 *
 * @tparam T_loc Type of location parameter
 * @tparam RNG Type of pseudo-random number generator
 * @param nu degrees of freedom parameter
 * @param mu (Sequence of) location parameter(s)
 * @param S Covariance matrix
 * @param rng random number generator
 * @throw std::domain_error if S is not positive definite, any value in mu is
 *    not finite, nu is not positive, or nu is NaN
 * @throw std::invalid_argument if the length of (each) mu is not equal to the
 *    number of rows and columns in S
 */
template <typename T_loc, class RNG>
inline typename StdVectorBuilder<true, Eigen::VectorXd, T_loc>::type
multi_student_t_rng(
    double nu, const T_loc& mu,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& S, RNG& rng) {
  using boost::normal_distribution;
  using boost::variate_generator;
  using boost::random::gamma_distribution;

  static const char* function = "multi_student_t_rng";
  check_not_nan(function, "Degrees of freedom parameter", nu);
  check_positive(function, "Degrees of freedom parameter", nu);
  check_positive(function, "Covariance matrix rows", S.rows());
  vector_seq_view<T_loc> mu_vec(mu);
  size_t size_mu = mu_vec[0].size();

  size_t N = size_mvt(mu);
  for (size_t i = 1; i < N; i++) {
    check_size_match(function,
                     "Size of one of the vectors of "
                     "the location variable",
                     mu_vec[i].size(),
                     "Size of the first vector of the "
                     "location variable",
                     size_mu);
  }

  for (size_t i = 0; i < N; i++) {
    check_finite(function, "Location parameter", mu_vec[i]);
  }
  const auto& S_ref = to_ref(S);
  check_not_nan(function, "Covariance matrix", S_ref);
  check_symmetric(function, "Covariance matrix", S_ref);
  Eigen::LLT<Eigen::MatrixXd> llt_of_S = S_ref.llt();
  check_pos_definite(function, "covariance matrix argument", llt_of_S);

  StdVectorBuilder<true, Eigen::VectorXd, T_loc> output(N);

  variate_generator<RNG&, normal_distribution<> > std_normal_rng(
      rng, normal_distribution<>(0, 1));
  variate_generator<RNG&, gamma_distribution<> > gamma_rng(
      rng, gamma_distribution<>(nu / 2.0, 2.0 / nu));

  double w = 1.0 / gamma_rng();
  for (size_t n = 0; n < N; ++n) {
    Eigen::VectorXd z(S.cols());
    for (int i = 0; i < S.cols(); i++) {
      z(i) = std::sqrt(w) * std_normal_rng();
    }

    output[n] = as_column_vector_or_scalar(mu_vec[n]) + llt_of_S.matrixL() * z;
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
