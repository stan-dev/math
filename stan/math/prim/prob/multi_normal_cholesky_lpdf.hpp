#ifndef STAN_MATH_PRIM_PROB_MULTI_NORMAL_CHOLESKY_LPDF_HPP
#define STAN_MATH_PRIM_PROB_MULTI_NORMAL_CHOLESKY_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/mdivide_left_tri.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * The log of the multivariate normal density for the given y, mu, and
 * a Cholesky factor L of the variance matrix.
 * Sigma = LL', a square, semi-positive definite matrix.
 *
 * Analytic expressions taken from
 * http://qwone.com/~jason/writing/multivariateNormal.pdf
 * written by Jason D. M. Rennie.
 *
 * All expressions are adapted to avoid (most) inversions and maximal
 * reuse of intermediates.
 *
 * @param y A scalar vector
 * @param mu The mean vector of the multivariate normal distribution.
 * @param L The Cholesky decomposition of a variance matrix
 * of the multivariate normal distribution
 * @return The log of the multivariate normal density.
 * @throw std::domain_error if LL' is not square, not symmetric,
 * or not semi-positive definite.
 * @tparam T_y Type of scalar.
 * @tparam T_loc Type of location.
 * @tparam T_covar Type of scale.
 */
template <bool propto, typename T_y, typename T_loc, typename T_covar>
return_type_t<T_y, T_loc, T_covar> multi_normal_cholesky_lpdf(
    const T_y& y, const T_loc& mu, const T_covar& L) {
  static const char* function = "multi_normal_cholesky_lpdf";
  using T_covar_elem = typename scalar_type<T_covar>::type;
  using T_return = return_type_t<T_y, T_loc, T_covar>;
  using T_partials_return = partials_return_t<T_y, T_loc, T_covar>;
  using matrix_partials_t
      = Eigen::Matrix<T_partials_return, Eigen::Dynamic, Eigen::Dynamic>;
  using vector_partials_t = Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1>;
  using row_vector_partials_t
      = Eigen::Matrix<T_partials_return, 1, Eigen::Dynamic>;

  check_consistent_sizes_mvt(function, "y", y, "mu", mu);
  size_t number_of_y = size_mvt(y);
  size_t number_of_mu = size_mvt(mu);
  if (number_of_y == 0 || number_of_mu == 0) {
    return 0;
  }
  vector_seq_view<T_y> y_vec(y);
  vector_seq_view<T_loc> mu_vec(mu);
  const size_t size_vec = max_size_mvt(y, mu);

  const int size_y = y_vec[0].size();
  const int size_mu = mu_vec[0].size();
  if (likely(size_vec > 1)) {
    // check size consistency of all random variables y
    int size_y_old = size_y;
    for (size_t i = 1, size_mvt_y = size_mvt(y); i < size_mvt_y; i++) {
      int size_y_new = y_vec[i].size();
      check_size_match(function,
                       "Size of one of the vectors of "
                       "the random variable",
                       size_y_new,
                       "Size of another vector of the "
                       "random variable",
                       size_y_old);
      size_y_old = size_y_new;
    }
    // check size consistency of all means mu
    int size_mu_old = size_mu;
    for (size_t i = 1, size_mvt_mu = size_mvt(mu); i < size_mvt_mu; i++) {
      int size_mu_new = mu_vec[i].size();
      check_size_match(function,
                       "Size of one of the vectors of "
                       "the location variable",
                       size_mu_new,
                       "Size of another vector of the "
                       "location variable",
                       size_mu_old);
      size_mu_old = size_mu_new;
    }
  }

  check_size_match(function, "Size of random variable", size_y,
                   "size of location parameter", size_mu);
  check_size_match(function, "Size of random variable", size_y,
                   "rows of covariance parameter", L.rows());
  check_size_match(function, "Size of random variable", size_y,
                   "columns of covariance parameter", L.cols());

  for (size_t i = 0; i < size_vec; i++) {
    check_finite(function, "Location parameter", mu_vec[i]);
    check_not_nan(function, "Random variable", y_vec[i]);
  }

  if (unlikely(size_y == 0)) {
    return T_return(0);
  }

  T_partials_return logp(0);
  operands_and_partials<T_y, T_loc, T_covar> ops_partials(y, mu, L);

  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * size_y * size_vec;
  }

  const matrix_partials_t inv_L_dbl
      = mdivide_left_tri<Eigen::Lower>(value_of(L));

  if (include_summand<propto, T_y, T_loc, T_covar_elem>::value) {
    for (size_t i = 0; i < size_vec; i++) {
      vector_partials_t y_minus_mu_dbl(size_y);
      for (int j = 0; j < size_y; j++) {
        y_minus_mu_dbl(j) = value_of(y_vec[i](j)) - value_of(mu_vec[i](j));
      }

      const row_vector_partials_t half
          = (inv_L_dbl.template triangularView<Eigen::Lower>() * y_minus_mu_dbl)
                .transpose();
      const vector_partials_t scaled_diff
          = (half * inv_L_dbl.template triangularView<Eigen::Lower>())
                .transpose();

      logp -= 0.5 * dot_self(half);

      if (!is_constant_all<T_y>::value) {
        for (int j = 0; j < size_y; j++) {
          ops_partials.edge1_.partials_vec_[i](j) -= scaled_diff(j);
        }
      }
      if (!is_constant_all<T_loc>::value) {
        for (int j = 0; j < size_y; j++) {
          ops_partials.edge2_.partials_vec_[i](j) += scaled_diff(j);
        }
      }
      if (!is_constant_all<T_covar>::value) {
        ops_partials.edge3_.partials_ += scaled_diff * half;
      }
    }
  }

  if (include_summand<propto, T_covar_elem>::value) {
    logp += sum(log(inv_L_dbl.diagonal())) * size_vec;
    if (!is_constant_all<T_covar>::value) {
      ops_partials.edge3_.partials_ -= size_vec * inv_L_dbl.transpose();
    }
  }

  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_covar>
inline return_type_t<T_y, T_loc, T_covar> multi_normal_cholesky_lpdf(
    const T_y& y, const T_loc& mu, const T_covar& L) {
  return multi_normal_cholesky_lpdf<false>(y, mu, L);
}

}  // namespace math
}  // namespace stan
#endif
