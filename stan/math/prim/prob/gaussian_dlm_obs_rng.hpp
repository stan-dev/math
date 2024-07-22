#ifndef STAN_MATH_PRIM_PROB_GAUSSIAN_DLM_OBS_RNG_HPP
#define STAN_MATH_PRIM_PROB_GAUSSIAN_DLM_OBS_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {

/** \ingroup multivar_dists
 * Return a multivariate normal random variate with the given location
 * and covariance using the specified random number generator.
 *
 * No error checking or templating, takes the LDLT directly to avoid
 * recomputation. Can sample from semidefinite covariance matrices.
 *
 * @tparam RNG type of pseudo-random number generator
 * @param mu location parameter
 * @param S_ldlt Eigen::LDLT of covariance matrix, semidefinite is okay
 * @param rng random number generator
 *
 */
template <class RNG>
inline Eigen::VectorXd multi_normal_semidefinite_rng(
    const Eigen::VectorXd &mu, const Eigen::LDLT<Eigen::MatrixXd> &S_ldlt,
    RNG &rng) {
  using boost::normal_distribution;
  using boost::variate_generator;

  variate_generator<RNG &, normal_distribution<>> std_normal_rng(
      rng, normal_distribution<>(0, 1));

  Eigen::VectorXd stddev = S_ldlt.vectorD().array().sqrt().matrix();
  size_t M = S_ldlt.vectorD().size();
  Eigen::VectorXd z(M);
  for (int i = 0; i < M; i++) {
    z(i) = stddev(i) * std_normal_rng();
  }

  Eigen::VectorXd Y
      = mu + (S_ldlt.transpositionsP().transpose() * (S_ldlt.matrixL() * z));
  // The inner paranthesis matter, transpositionsP() gives a
  // permutation matrix from pivoting and matrixL() gives a lower
  // triangular matrix. The types cannot be directly multiplied.

  return Y;
}

}  // namespace internal

/** \ingroup multivar_dists
 * Simulate random draw from Gaussian dynamic linear model (GDLM).
 * This distribution is equivalent to, for \f$t = 1:T\f$,
 * \f{eqnarray*}{
 * y_t & \sim N(F' \theta_t, V) \\
 * \theta_t & \sim N(G \theta_{t-1}, W) \\
 * \theta_0 & \sim N(m_0, C_0)
 * \f}
 *
 * @tparam RNG Type of pseudo-random number generator.
 * @param F A n x r matrix. The design matrix.
 * @param G A n x n matrix. The transition matrix.
 * @param V A r x r matrix. The observation covariance matrix.
 * @param W A n x n matrix. The state covariance matrix.
 * @param m0 A n x 1 matrix. The mean vector of the distribution
 * of the initial state.
 * @param C0 A n x n matrix. The covariance matrix of the distribution
 * of the initial state.
 * @param T a positive integer, how many timesteps to simulate.
 * @param rng Pseudo-random number generator.
 * @return A r x T matrix of simulated observations. Rows are
 * variables, columns are observations. First column is the state
 * after the first transition. Last column is the state after the last
 * transition. Initial state not returned.
 * @throw std::domain_error if a matrix is not symmetric or not
 * positive semi-definite. Or throw std::invalid_argument if a size is
 * wrong or any input is NaN or non-finite, or if T is not
 * positive. Require C0 in particular to be strictly positive
 * definite. V and W can be semidefinite.
 *
 */
template <class RNG>
inline Eigen::MatrixXd gaussian_dlm_obs_rng(const Eigen::MatrixXd &F,
                                            const Eigen::MatrixXd &G,
                                            const Eigen::MatrixXd &V,
                                            const Eigen::MatrixXd &W,
                                            const Eigen::VectorXd &m0,
                                            const Eigen::MatrixXd &C0,
                                            const int T, RNG &rng) {
  static constexpr const char *function = "gaussian_dlm_obs_rng";

  int r = F.cols();  // number of variables
  int n = G.rows();  // number of states

  check_size_match(function, "rows of F", F.rows(), "rows of G", n);
  check_finite(function, "F", F);
  check_square(function, "G", G);
  check_finite(function, "G", G);
  check_size_match(function, "rows of V", V.rows(), "cols of F", r);
  check_finite(function, "V", V);
  check_positive(function, "V rows", V.rows());
  check_symmetric(function, "V", V);
  check_size_match(function, "rows of W", W.rows(), "rows of G", n);
  check_finite(function, "W", W);
  check_positive(function, "W rows", W.rows());
  check_symmetric(function, "W", W);
  check_size_match(function, "rows of W", W.rows(), "rows of G", n);
  check_size_match(function, "size of m0", m0.size(), "rows of G", n);
  check_finite(function, "m0", m0);
  check_size_match(function, "rows of C0", C0.rows(), "rows of G", n);
  check_finite(function, "C0", C0);
  check_positive(function, "C0 rows", C0.rows());
  check_symmetric(function, "C0", C0);
  check_positive(function, "T", T);

  Eigen::LDLT<Eigen::MatrixXd> V_ldlt = V.ldlt();
  check_pos_semidefinite(function, "V", V_ldlt);
  Eigen::LDLT<Eigen::MatrixXd> W_ldlt = W.ldlt();
  check_pos_semidefinite(function, "W", W_ldlt);
  Eigen::LDLT<Eigen::MatrixXd> C0_ldlt = C0.ldlt();
  check_pos_semidefinite(function, "C0", C0_ldlt);

  Eigen::MatrixXd y(r, T);
  Eigen::VectorXd theta_t
      = internal::multi_normal_semidefinite_rng(m0, C0_ldlt, rng);
  for (int t = 0; t < T; ++t) {
    theta_t = internal::multi_normal_semidefinite_rng(
        Eigen::VectorXd(G * theta_t), W_ldlt, rng);
    y.col(t) = internal::multi_normal_semidefinite_rng(
        Eigen::VectorXd(F.transpose() * theta_t), V_ldlt, rng);
  }
  return y;
}

template <class RNG>
inline Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
gaussian_dlm_obs_rng(const Eigen::MatrixXd &F, const Eigen::MatrixXd &G,
                     const Eigen::VectorXd &V, const Eigen::MatrixXd &W,
                     const Eigen::VectorXd &m0, const Eigen::MatrixXd &C0,
                     const int T, RNG &rng) {
  return gaussian_dlm_obs_rng(F, G, Eigen::MatrixXd(V.asDiagonal()), W, m0, C0,
                              T, rng);
}

}  // namespace math
}  // namespace stan
#endif
