#ifndef STAN_MATH_PRIM_MAT_PROB_ORDERED_LOGISTIC_LPMF_HPP
#define STAN_MATH_PRIM_MAT_PROB_ORDERED_LOGISTIC_LPMF_HPP

#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/size.hpp>
#include <stan/math/prim/mat/meta/vector_seq_view.hpp>
#include <stan/math/prim/mat/err/check_ordered.hpp>
#include <stan/math/prim/scal/fun/inv_logit.hpp>
#include <stan/math/prim/scal/fun/log1p_exp.hpp>
#include <stan/math/prim/scal/fun/log_inv_logit_diff.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_greater.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T>
struct ordLog_helper {
  /**
   * Returns the (natural) log probability of the specified integer
   * outcome given the continuous location and specified cutpoints
   * in an ordered logistic model.
   *
   * Error-checking handled by main distribution functions, with this
   * function only called once inputs have been validated.
   *
   * @tparam T Type of location & cutpoint variables.
   * @param y Outcome.
   * @param K Number of categories.
   * @param lambda Location.
   * @param c Positive increasing vector of cutpoints.
   * @return Log probability of outcome given location and
   * cutpoints.
   */
  T logp(const int& y, const int& K, const T& lambda,
         const Eigen::Matrix<T, Eigen::Dynamic, 1>& c) {
    if (y == 1)
      return -log1p_exp(lambda - c[0]);
    else if (y == K)
      return -log1p_exp(c[K - 2] - lambda);
    else
      return log_inv_logit_diff(lambda - c[y - 2], lambda - c[y - 1]);
  }

  /**
   * Returns a vector with the gradients of the the continuous location
   * and ordered cutpoints in an ordered logistic model. The first element
   * of the vector contains the gradient for the location variable (lambda),
   * followed by the gradients for the ordered cutpoints (c).
   *
   * Error-checking handled by main distribution functions, with this
   * function only called once inputs have been validated.
   *
   * @tparam T Type of location & cutpoint variables.
   * @param y Outcome.
   * @param K Number of categories.
   * @param lambda Location.
   * @param c Positive increasing vector of cutpoints.
   * @return Vector of gradients.
   */
  Eigen::Matrix<T, Eigen::Dynamic, 1> deriv(
      const int& y, const int& K, const T& lambda,
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& c) {
    using std::exp;

    Eigen::Matrix<T, Eigen::Dynamic, 1> d(
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(K));

    if (y == 1) {
      d[0] -= inv_logit(lambda - c[0]);
      d[1] -= d[0];
      return d;
    } else if (y == K) {
      d[0] += inv_logit(c[K - 2] - lambda);
      d[K - 1] -= d[0];
      return d;
    } else {
      d[y - 1]
          += inv(1 - exp(c[y - 1] - c[y - 2])) - inv_logit(c[y - 2] - lambda);
      d[y] += inv(1 - exp(c[y - 2] - c[y - 1])) - inv_logit(c[y - 1] - lambda);
      d[0] -= d[y] + d[y - 1];
      return d;
    }
  }
};

/**
 * Returns the (natural) log probability of the specified integer
 * outcome given the continuous location and specified cutpoints
 * in an ordered logistic model.
 *
 * <p>Typically the continous location
 * will be the dot product of a vector of regression coefficients
 * and a vector of predictors for the outcome.
 *
 * @tparam propto True if calculating up to a proportion.
 * @tparam T_loc Location type.
 * @tparam T_cut Cut-point type.
 * @param y Outcome.
 * @param lambda Location.
 * @param c Positive increasing vector of cutpoints.
 * @return Log probability of outcome given location and
 * cutpoints.
 * @throw std::domain_error If the outcome is not between 1 and
 * the number of cutpoints plus 2; if the cutpoint vector is
 * empty; if the cutpoint vector contains a non-positive,
 * non-finite value; or if the cutpoint vector is not sorted in
 * ascending order.
 */
template <bool propto, typename T_loc, typename T_cut>
typename return_type<T_loc, T_cut>::type ordered_logistic_lpmf(
    int y, const T_loc& lambda,
    const Eigen::Matrix<T_cut, Eigen::Dynamic, 1>& c) {
  static const char* function = "ordered_logistic";

  typedef typename stan::partials_return_type<
      T_loc, Eigen::Matrix<T_cut, Eigen::Dynamic, 1>>::type T_partials_return;

  typedef typename Eigen::Matrix<T_partials_return, -1, 1> T_partials_vec;

  int K = c.size() + 1;

  check_bounded(function, "Random variable", y, 1, K);
  check_finite(function, "Location parameter", lambda);
  check_greater(function, "Size of cut points parameter", c.size(), 0);
  check_ordered(function, "Cut-points", c);
  check_finite(function, "Final cut-point", c(c.size() - 1));
  check_finite(function, "First cut-point", c(0));

  scalar_seq_view<T_loc> lam_vec(lambda);
  T_partials_return lam_dbl = value_of(lam_vec[0]);

  vector_seq_view<Eigen::Matrix<T_cut, Eigen::Dynamic, 1>> c_vec(c);
  T_partials_vec c_dbl((K - 1));
  c_dbl = value_of(c_vec[0]).template cast<T_partials_return>();

  ordLog_helper<T_partials_return> ordhelp;

  T_partials_vec d = ordhelp.deriv(y, K, lam_dbl, c_dbl);

  operands_and_partials<T_loc, Eigen::Matrix<T_cut, Eigen::Dynamic, 1>>
      ops_partials(lambda, c);

  if (!is_constant_struct<T_loc>::value)
    ops_partials.edge1_.partials_[0] = d[0];

  if (!is_constant_struct<Eigen::Matrix<T_cut, Eigen::Dynamic, 1>>::value)
    ops_partials.edge2_.partials_ = d.tail(K - 1);

  return ops_partials.build(ordhelp.logp(y, K, lam_dbl, c_dbl));
}

template <typename T_loc, typename T_cut>
typename return_type<T_loc, T_cut>::type ordered_logistic_lpmf(
    int y, const T_loc& lambda,
    const Eigen::Matrix<T_cut, Eigen::Dynamic, 1>& c) {
  return ordered_logistic_lpmf<false>(y, lambda, c);
}

/**
 * Returns the (natural) log probability of the specified array
 * of integers given the vector of continuous locations and
 * specified cutpoints in an ordered logistic model.
 *
 * <p>Typically the continous location
 * will be the dot product of a vector of regression coefficients
 * and a vector of predictors for the outcome.
 *
 * @tparam propto True if calculating up to a proportion.
 * @tparam T_loc Location type.
 * @tparam T_cut Cut-point type.
 * @param y Array of integers
 * @param lambda Vector of continuous location variables.
 * @param c Positive increasing vector of cutpoints.
 * @return Log probability of outcome given location and
 * cutpoints.
 * @throw std::domain_error If the outcome is not between 1 and
 * the number of cutpoints plus 2; if the cutpoint vector is
 * empty; if the cutpoint vector contains a non-positive,
 * non-finite value; or if the cutpoint vector is not sorted in
 * ascending order.
 * @throw std::invalid_argument If y and lambda are different
 * lengths.
 */
template <bool propto, typename T_loc, typename T_cut>
typename return_type<T_loc, T_cut>::type ordered_logistic_lpmf(
    const std::vector<int>& y,
    const Eigen::Matrix<T_loc, Eigen::Dynamic, 1>& lambda,
    const Eigen::Matrix<T_cut, Eigen::Dynamic, 1>& c) {
  static const char* function = "ordered_logistic";

  typedef typename stan::partials_return_type<
      Eigen::Matrix<T_loc, Eigen::Dynamic, 1>,
      Eigen::Matrix<T_cut, Eigen::Dynamic, 1>>::type T_partials_return;

  typedef typename Eigen::Matrix<T_partials_return, -1, 1> T_partials_vec;

  int N = lambda.size();
  int K = c.size() + 1;

  check_consistent_sizes(function, "Integers", y, "Locations", lambda);
  check_bounded(function, "Random variable", y, 1, K);
  check_finite(function, "Location parameter", lambda);
  check_ordered(function, "Cut-points", c);
  check_greater(function, "Size of cut points parameter", c.size(), 0);
  check_finite(function, "Final cut-point", c(c.size() - 1));
  check_finite(function, "First cut-point", c(0));

  vector_seq_view<Eigen::Matrix<T_loc, Eigen::Dynamic, 1>> lam_vec(lambda);
  T_partials_vec lam_dbl
      = value_of(lam_vec[0]).template cast<T_partials_return>();

  vector_seq_view<Eigen::Matrix<T_cut, Eigen::Dynamic, 1>> c_vec(c);
  T_partials_vec c_dbl((K - 1));
  c_dbl = value_of(c_vec[0]).template cast<T_partials_return>();

  ordLog_helper<T_partials_return> ordhelp;

  T_partials_return logp(0.0);
  T_partials_vec lam_deriv(N);
  T_partials_vec c_deriv(T_partials_vec::Zero(K - 1));

  for (int n = 0; n < N; ++n) {
    T_partials_vec d = ordhelp.deriv(y[n], K, lam_dbl[n], c_dbl);

    logp += ordhelp.logp(y[n], K, lam_dbl[n], c_dbl);
    lam_deriv[n] = d[0];
    c_deriv += d.tail(K - 1);
  }

  operands_and_partials<Eigen::Matrix<T_loc, Eigen::Dynamic, 1>,
                        Eigen::Matrix<T_cut, Eigen::Dynamic, 1>>
      ops_partials(lambda, c);

  if (!is_constant_struct<Eigen::Matrix<T_loc, Eigen::Dynamic, 1>>::value)
    ops_partials.edge1_.partials_ = lam_deriv;

  if (!is_constant_struct<Eigen::Matrix<T_cut, Eigen::Dynamic, 1>>::value)
    ops_partials.edge2_.partials_ = c_deriv;

  return ops_partials.build(logp);
}

template <typename T_loc, typename T_cut>
typename return_type<T_loc, T_cut>::type ordered_logistic_lpmf(
    const std::vector<int>& y,
    const Eigen::Matrix<T_loc, Eigen::Dynamic, 1>& lambda,
    const Eigen::Matrix<T_cut, Eigen::Dynamic, 1>& c) {
  return ordered_logistic_lpmf<false>(y, lambda, c);
}

/**
 * Returns the (natural) log probability of the specified array
 * of integers given the vector of continuous locations and
 * array of specified cutpoints in an ordered logistic model.
 *
 * <p>Typically the continous location
 * will be the dot product of a vector of regression coefficients
 * and a vector of predictors for the outcome.
 *
 * @tparam propto True if calculating up to a proportion.
 * @tparam T_y Type of y variable (should be std::vector<int>).
 * @tparam T_loc Location type.
 * @tparam T_cut Cut-point type.
 * @param y Array of integers
 * @param lambda Vector of continuous location variables.
 * @param c array of Positive increasing vectors of cutpoints.
 * @return Log probability of outcome given location and
 * cutpoints.
 * @throw std::domain_error If the outcome is not between 1 and
 * the number of cutpoints plus 2; if the cutpoint vector is
 * empty; if the cutpoint vector contains a non-positive,
 * non-finite value; or if the cutpoint vector is not sorted in
 * ascending order.
 * @throw std::invalid_argument If y and lambda are different
 * lengths, or if y and the array of cutpoints are of different
 * lengths.
 */
template <bool propto, typename T_loc, typename T_cut>
typename return_type<T_loc, T_cut>::type ordered_logistic_lpmf(
    const std::vector<int>& y,
    const Eigen::Matrix<T_loc, Eigen::Dynamic, 1>& lambda,
    const std::vector<Eigen::Matrix<T_cut, Eigen::Dynamic, 1>>& c) {
  static const char* function = "ordered_logistic";

  typedef typename stan::partials_return_type<
      Eigen::Matrix<T_loc, Eigen::Dynamic, 1>,
      std::vector<Eigen::Matrix<T_cut, Eigen::Dynamic, 1>>>::type
      T_partials_return;

  typedef typename Eigen::Matrix<T_partials_return, -1, 1> T_partials_vec;

  int N = lambda.size();
  int K = c[0].size() + 1;

  for (int n = 0; n < N; ++n) {
    check_bounded(function, "Random variable", y[n], 1, K);
    check_greater(function, "Size of cut points parameter", c[n].size(), 0);
    check_ordered(function, "Cut-points", c[n]);
  }
  check_consistent_sizes(function, "Integers", y, "Locations", lambda);
  check_consistent_sizes(function, "Integers", y, "Cut-points", c);
  check_finite(function, "Location parameter", lambda);
  check_finite(function, "Cut-points", c);

  vector_seq_view<Eigen::Matrix<T_loc, Eigen::Dynamic, 1>> lam_vec(lambda);
  T_partials_vec lam_dbl
      = value_of(lam_vec[0]).template cast<T_partials_return>();

  vector_seq_view<std::vector<Eigen::Matrix<T_cut, Eigen::Dynamic, 1>>> c_vec(
      c);
  std::vector<T_partials_vec> c_dbl(N);
  for (int n = 0; n < N; ++n)
    c_dbl[n] = value_of(c_vec[n]).template cast<T_partials_return>();

  ordLog_helper<T_partials_return> ordhelp;

  T_partials_return logp(0.0);
  T_partials_vec lam_deriv(N);
  std::vector<T_partials_vec> c_deriv(N);

  for (int n = 0; n < N; ++n) {
    T_partials_vec d = ordhelp.deriv(y[n], K, lam_dbl[n], c_dbl[n]);

    logp += ordhelp.logp(y[n], K, lam_dbl[n], c_dbl[n]);
    lam_deriv[n] = d[0];
    c_deriv[n] = d.tail(K - 1);
  }

  operands_and_partials<Eigen::Matrix<T_loc, Eigen::Dynamic, 1>,
                        std::vector<Eigen::Matrix<T_cut, Eigen::Dynamic, 1>>>
      ops_partials(lambda, c);

  if (!is_constant_struct<Eigen::Matrix<T_loc, Eigen::Dynamic, 1>>::value)
    ops_partials.edge1_.partials_ = lam_deriv;

  if (!is_constant_struct<
          std::vector<Eigen::Matrix<T_cut, Eigen::Dynamic, 1>>>::value) {
    for (int n = 0; n < N; ++n)
      ops_partials.edge2_.partials_vec_[n] = c_deriv[n];
  }

  return ops_partials.build(logp);
}

template <typename T_loc, typename T_cut>
typename return_type<T_loc, T_cut>::type ordered_logistic_lpmf(
    const std::vector<int>& y,
    const Eigen::Matrix<T_loc, Eigen::Dynamic, 1>& lambda,
    const std::vector<Eigen::Matrix<T_cut, Eigen::Dynamic, 1>>& c) {
  return ordered_logistic_lpmf<false>(y, lambda, c);
}
}  // namespace math
}  // namespace stan
#endif
