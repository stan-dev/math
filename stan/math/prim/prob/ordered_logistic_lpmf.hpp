#ifndef STAN_MATH_PRIM_PROB_ORDERED_LOGISTIC_LPMF_HPP
#define STAN_MATH_PRIM_PROB_ORDERED_LOGISTIC_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/is_integer.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <stan/math/prim/fun/log_inv_logit_diff.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Returns the (natural) log probability of the specified array
 * of integers given the vector of continuous locations and
 * specified cutpoints in an ordered logistic model.
 *
 * <p>Typically the continuous location
 * will be the dot product of a vector of regression coefficients
 * and a vector of predictors for the outcome
 *
  \f[
    \frac{\partial }{\partial \lambda} =
    \begin{cases}\\
    -\mathrm{logit}^{-1}(\lambda - c_1) & \mbox{if } k = 1,\\
    -(((1-e^{c_{k-1}-c_{k-2}})^{-1} - \mathrm{logit}^{-1}(c_{k-2}-\lambda)) +
    ((1-e^{c_{k-2}-c_{k-1}})^{-1} - \mathrm{logit}^{-1}(c_{k-1}-\lambda)))
    & \mathrm{if } 1 < k < K, \mathrm{and}\\
    \mathrm{logit}^{-1}(c_{K-2}-\lambda) & \mathrm{if } k = K.
    \end{cases}
  \f]

  \f[
    \frac{\partial }{\partial \lambda} =
    \begin{cases}
    -\mathrm{logit}^{-1}(\lambda - c_1) & \text{if } k = 1,\\
    -(((1-e^{c_{k-1}-c_{k-2}})^{-1} - \mathrm{logit}^{-1}(c_{k-2}-\lambda)) +
    ((1-e^{c_{k-2}-c_{k-1}})^{-1} - \mathrm{logit}^{-1}(c_{k-1}-\lambda)))
    & \text{if } 1 < k < K, \text{ and}\\
    \mathrm{logit}^{-1}(c_{K-2}-\lambda) & \text{if } k = K.
    \end{cases}
  \f]
 *
 * @tparam propto True if calculating up to a proportion.
 * @tparam T_y Y variable type (integer or array of integers).
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
template <bool propto, typename T_y, typename T_loc, typename T_cut,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_cut>* = nullptr>
return_type_t<T_loc, T_cut> ordered_logistic_lpmf(const T_y& y,
                                                  const T_loc& lambda,
                                                  const T_cut& c) {
  using T_partials_return = partials_return_t<T_loc, T_cut>;
  using T_cuts_val = partials_return_t<T_cut>;
  using T_y_ref = ref_type_t<T_y>;
  using T_lambda_ref = ref_type_if_t<!is_constant<T_loc>::value, T_loc>;
  using T_cut_ref = ref_type_if_t<!is_constant<T_cut>::value, T_cut>;
  using Eigen::Array;
  using Eigen::Dynamic;
  static const char* function = "ordered_logistic";

  T_cut_ref c_ref = c;
  vector_seq_view<T_cut_ref> c_vec(c_ref);
  int N = size(lambda);
  int C_l = size_mvt(c);

  check_consistent_sizes(function, "Integers", y, "Locations", lambda);
  if (C_l > 1) {
    check_size_match(function, "Length of location variables ", N,
                     "Number of cutpoint vectors ", C_l);
  }

  T_y_ref y_ref = y;
  T_lambda_ref lambda_ref = lambda;

  scalar_seq_view<T_y_ref> y_seq(y_ref);

  decltype(auto) lambda_val = to_ref(as_value_array_or_scalar(lambda_ref));

  check_finite(function, "Location parameter", lambda_val);
  if (C_l == 0 || N == 0) {
    return 0;
  }
  int K = c_vec[0].size() + 1;
  check_bounded(function, "Random variable", y_ref, 1, K);

  for (int i = 0; i < C_l; i++) {
    check_size_match(function, "One cutpoint set", c_vec[i].size(),
                     "First cutpoint set", K - 1);
    check_ordered(function, "Cut-points", c_vec[i]);
    if (K > 1) {
      if (K > 2) {
        check_finite(function, "Final cut-point", c_vec[i].coeff(K - 2));
      }
      check_finite(function, "First cut-point", c_vec[i].coeff(0));
    }
  }
  if (!include_summand<propto, T_loc, T_cut>::value) {
    return 0.0;
  }

  scalar_seq_view<decltype(lambda_val)> lam_vec(lambda_val);
  T_partials_return logp(0.0);
  Array<T_cuts_val, Dynamic, 1> cuts_y1(N), cuts_y2(N);
  for (int i = 0; i < N; i++) {
    int c = y_seq[i];
    if (c != K) {
      cuts_y1.coeffRef(i) = value_of(c_vec[i].coeff(c - 1));
    } else {
      cuts_y1.coeffRef(i) = INFINITY;
    }
    if (c != 1) {
      cuts_y2.coeffRef(i) = value_of(c_vec[i].coeff(c - 2));
    } else {
      cuts_y2.coeffRef(i) = -INFINITY;
    }
  }

  Array<T_partials_return, Dynamic, 1> cut2 = lambda_val - cuts_y2;
  Array<T_partials_return, Dynamic, 1> cut1 = lambda_val - cuts_y1;

  // Not immediately evaluating next two expressions benefits performance
  auto m_log_1p_exp_cut1
      = (cut1 > 0.0).select(-cut1, 0) - (-cut1.abs()).exp().log1p();
  auto m_log_1p_exp_m_cut2
      = (cut2 <= 0.0).select(cut2, 0) - (-cut2.abs()).exp().log1p();

  if (is_vector<T_y>::value) {
    Eigen::Map<const Eigen::Matrix<value_type_t<T_y>, Eigen::Dynamic, 1>> y_vec(
        y_seq.data(), y_seq.size());
    auto log1m_exp_cuts_diff = log1m_exp(cut1 - cut2);
    logp = y_vec.cwiseEqual(1)
               .select(m_log_1p_exp_cut1,
                       y_vec.cwiseEqual(K).select(m_log_1p_exp_m_cut2,
                                                  m_log_1p_exp_m_cut2
                                                      + log1m_exp_cuts_diff
                                                      + m_log_1p_exp_cut1))
               .sum();
  } else {
    if (y_seq[0] == 1) {
      logp = m_log_1p_exp_cut1.sum();
    } else if (y_seq[0] == K) {
      logp = m_log_1p_exp_m_cut2.sum();
    } else {
      logp = (m_log_1p_exp_m_cut2 + log1m_exp(cut1 - cut2).array()
              + m_log_1p_exp_cut1)
                 .sum();
    }
  }

  operands_and_partials<T_lambda_ref, T_cut_ref> ops_partials(lambda_ref,
                                                              c_ref);
  if (!is_constant_all<T_loc, T_cut>::value) {
    Array<T_partials_return, Dynamic, 1> exp_m_cut1 = exp(-cut1);
    Array<T_partials_return, Dynamic, 1> exp_m_cut2 = exp(-cut2);
    Array<T_partials_return, Dynamic, 1> exp_cuts_diff = exp(cuts_y2 - cuts_y1);
    Array<T_partials_return, Dynamic, 1> d1
        = (cut2 > 0).select(exp_m_cut2 / (1 + exp_m_cut2), 1 / (1 + exp(cut2)))
          - exp_cuts_diff / (exp_cuts_diff - 1);
    Array<T_partials_return, Dynamic, 1> d2
        = 1 / (1 - exp_cuts_diff)
          - (cut1 > 0).select(exp_m_cut1 / (1 + exp_m_cut1),
                              1 / (1 + exp(cut1)));
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge1_.partials_ = d1 - d2;
    }
    if (!is_constant_all<T_cut>::value) {
      for (int i = 0; i < N; i++) {
        int c = y_seq[i];
        if (c != K) {
          ops_partials.edge2_.partials_vec_[i][c - 1] += d2.coeff(i);
        }
        if (c != 1) {
          ops_partials.edge2_.partials_vec_[i][c - 2] -= d1.coeff(i);
        }
      }
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_cut>
return_type_t<T_loc, T_cut> ordered_logistic_lpmf(const T_y& y,
                                                  const T_loc& lambda,
                                                  const T_cut& c) {
  return ordered_logistic_lpmf<false>(y, lambda, c);
}

}  // namespace math
}  // namespace stan
#endif
