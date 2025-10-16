#ifndef STAN_MATH_PRIM_PROB_ORDERED_PROBIT_LPMF_HPP
#define STAN_MATH_PRIM_PROB_ORDERED_PROBIT_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log_diff_exp.hpp>
#include <stan/math/prim/fun/log1m_exp.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/prob/std_normal_lcdf.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <vector>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Returns the (natural) log probability of the specified array
 * of integers given the vector of continuous locations and
 * array of specified cutpoints in an ordered probit model.
 *
 * <p>Typically the continuous location
 * will be the dot product of a vector of regression coefficients
 * and a vector of predictors for the outcome.
 *
 * @tparam propto True if calculating up to a proportion.
 * @tparam T_y Type of y variable - `int` or `std::vector<int>`.
 * @tparam T_loc Location type - Eigen vector or scalar.
 * @tparam T_cut Cut-point type - Eigen vector or a std vector of Eigen vectors.
 * @param y integer or Array of integers
 * @param lambda Location.
 * @param c Positive increasing vectors of cutpoints.
 * @return Log probability of outcome given location and
 * cutpoints.
 * @throw std::domain_error If the outcome is not between 1 and
 * the number of cutpoints plus 2; if the cutpoint vector contains a
 * non-positive or non-finite value; or if the cutpoint vector is not sorted in
 * ascending order.
 * @throw std::invalid_argument If y and lambda are different
 * lengths; if the cutpoint vector is empty; if y and the array of cutpoints are
 * of different lengths.
 */
template <bool propto, typename T_y, typename T_loc, typename T_cut>
return_type_t<T_loc, T_cut> ordered_probit_lpmf(const T_y& y,
                                                const T_loc& lambda,
                                                const T_cut& c) {
  using std::exp;
  using std::log;
  using T_lambda_ref = ref_type_t<T_loc>;
  static constexpr const char* function = "ordered_probit";

  check_nonzero_size(function, "Cut-points", c);
  int N = math::size(lambda);
  int C_l = size_mvt(c);
  vector_seq_view<T_cut> c_vec(c);
  int K = c_vec[0].size() + 1;

  check_consistent_sizes(function, "Integers", y, "Locations", lambda);
  if (C_l > 1) {
    check_size_match(function, "Length of location variables ", N,
                     "Number of cutpoint vectors ", C_l);
  }

  check_bounded(function, "Random variable", y, 1, K);
  scalar_seq_view<T_y> y_vec(y);
  check_nonzero_size(function, "First cutpoint set", c_vec[0]);
  for (int i = 0; i < size_mvt(c); ++i) {
    check_size_match(function, "One cutpoint set", K - 1, "First cutpoint set",
                     c_vec[i].size());
    check_ordered(function, "Cut-points", c_vec[i]);
    if (K > 2) {
      check_finite(function, "Final cut point", c_vec[i].coeff(K - 2));
    }
    check_finite(function, "First cut point", c_vec[i].coeff(0));
  }

  T_lambda_ref lambda_ref = lambda;
  check_finite(function, "Location parameter", lambda_ref);
  scalar_seq_view<T_lambda_ref> lambda_vec(lambda_ref);

  return_type_t<T_loc, T_cut> logp_n(0.0);

  for (int i = 0; i < N; ++i) {
    int K = c_vec[i].size() + 1;

    if (y_vec[i] == 1) {
      logp_n += std_normal_lcdf(c_vec[i].coeff(0) - lambda_vec[i]);
    } else if (y_vec[i] == K) {
      logp_n += std_normal_lcdf(lambda_vec[i] - c_vec[i].coeff(K - 2));
    } else {
      logp_n += log_diff_exp(
          std_normal_lcdf(c_vec[i].coeff(y_vec[i] - 1) - lambda_vec[i]),
          std_normal_lcdf(c_vec[i].coeff(y_vec[i] - 2) - lambda_vec[i]));
    }
  }
  return logp_n;
}

template <typename T_y, typename T_loc, typename T_cut>
return_type_t<T_loc, T_cut> ordered_probit_lpmf(const T_y& y,
                                                const T_loc& lambda,
                                                const T_cut& c) {
  return ordered_probit_lpmf<false>(y, lambda, c);
}
}  // namespace math
}  // namespace stan
#endif
