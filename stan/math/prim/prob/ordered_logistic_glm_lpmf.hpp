#ifndef STAN_MATH_PRIM_PROB_ORDERED_LOGISTIC_GLM_LPMF_HPP
#define STAN_MATH_PRIM_PROB_ORDERED_LOGISTIC_GLM_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/isfinite.hpp>
#include <stan/math/prim/fun/log1m_exp.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Returns the log PMF of the ordinal regression Generalized Linear Model
 * (GLM). This is equivalent to and faster than ordered_logistic_lpmf(y, x *
 * beta, cuts).
 *
 * @tparam T_y type of integer vector of classes. It can be either
 * `std::vector<int>` or `int`.
 * @tparam T_x type of the matrix of independent variables (features)
 * @tparam T_beta type of the vector of weights
 * @tparam T_cuts type of the vector of cutpoints
 *
 * @param y a scalar or vector of classes. If it is a scalar it will be
 * broadcast - used for all instances. Values should be between 1 and number
 * of classes, including endpoints.
 * @param x design matrix or row vector. If it is a row vector it will be
 * broadcast - used for all instances.
 * @param beta weight vector
 * @param cuts cutpoints vector
 * @return log probability
 * @throw std::domain_error If any class is not between 1 and
 * the number of cutpoints plus 2 or if the cutpoint vector is not sorted in
 * ascending order or any input is not finite
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y, typename T_x, typename T_beta,
          typename T_cuts, require_matrix_t<T_x>* = nullptr,
          require_all_col_vector_t<T_beta, T_cuts>* = nullptr>
return_type_t<T_x, T_beta, T_cuts> ordered_logistic_glm_lpmf(
    const T_y& y, const T_x& x, const T_beta& beta, const T_cuts& cuts) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;
  using std::exp;
  using std::isfinite;
  constexpr int T_x_rows = T_x::RowsAtCompileTime;
  using T_cuts_partials = partials_return_t<T_cuts>;
  using T_xbeta_partials = partials_return_t<T_x, T_beta>;
  using T_partials_return = partials_return_t<T_y, T_x, T_beta, T_cuts>;
  typedef typename std::conditional_t<
      T_x_rows == 1, T_xbeta_partials,
      Eigen::Matrix<T_xbeta_partials, Eigen::Dynamic, 1>>
      T_location;
  using T_y_ref = ref_type_t<T_y>;
  using T_x_ref = ref_type_if_t<!is_constant<T_x>::value, T_x>;
  using T_beta_ref = ref_type_if_t<!is_constant<T_beta>::value, T_beta>;
  using T_cuts_ref = ref_type_if_t<!is_constant<T_cuts>::value, T_cuts>;

  const size_t N_instances = T_x_rows == 1 ? stan::math::size(y) : x.rows();
  const size_t N_attributes = x.cols();
  const size_t N_classes = stan::math::size(cuts) + 1;

  static const char* function = "ordered_logistic_glm_lpmf";
  check_consistent_size(function, "Vector of dependent variables", y,
                        N_instances);
  check_consistent_size(function, "Weight vector", beta, N_attributes);
  T_y_ref y_ref = y;
  T_cuts_ref cuts_ref = cuts;
  const auto& cuts_val = value_of(cuts_ref);
  const auto& cuts_val_vec = to_ref(as_column_vector_or_scalar(cuts_val));
  check_bounded(function, "Vector of dependent variables", y_ref, 1, N_classes);
  check_ordered(function, "Cut-points", cuts_val_vec);
  if (N_classes > 1) {
    if (N_classes > 2) {
      check_finite(function, "Final cut-point", cuts_val_vec[N_classes - 2]);
    }
    check_finite(function, "First cut-point", cuts_val_vec[0]);
  }

  if (size_zero(y, cuts)) {
    return 0;
  }
  if (!include_summand<propto, T_x, T_beta, T_cuts>::value) {
    return 0;
  }

  T_x_ref x_ref = x;
  T_beta_ref beta_ref = beta;

  const auto& x_val = to_ref_if<!is_constant<T_beta>::value>(value_of(x_ref));
  const auto& beta_val = value_of(beta_ref);

  const auto& beta_val_vec = to_ref_if<!is_constant<T_x>::value>(
      as_column_vector_or_scalar(beta_val));

  scalar_seq_view<T_y_ref> y_seq(y_ref);
  Array<T_cuts_partials, Dynamic, 1> cuts_y1(N_instances), cuts_y2(N_instances);
  for (int i = 0; i < N_instances; i++) {
    int c = y_seq[i];
    if (c != N_classes) {
      cuts_y1.coeffRef(i) = cuts_val_vec.coeff(c - 1);
    } else {
      cuts_y1.coeffRef(i) = INFINITY;
    }
    if (c != 1) {
      cuts_y2.coeffRef(i) = cuts_val_vec.coeff(c - 2);
    } else {
      cuts_y2.coeffRef(i) = -INFINITY;
    }
  }

  T_location location = x_val * beta_val_vec;
  if (!isfinite(sum(location))) {
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Matrix of independent variables", x);
  }

  Array<T_partials_return, Dynamic, 1> cut2
      = as_array_or_scalar(location) - cuts_y2;
  Array<T_partials_return, Dynamic, 1> cut1
      = as_array_or_scalar(location) - cuts_y1;

  // Not immediately evaluating next two expressions benefits performance
  auto m_log_1p_exp_cut1
      = (cut1 > 0.0).select(-cut1, 0) - (-cut1.abs()).exp().log1p();
  auto m_log_1p_exp_m_cut2
      = (cut2 <= 0.0).select(cut2, 0) - (-cut2.abs()).exp().log1p();

  T_partials_return logp(0);
  if (is_vector<T_y>::value) {
    Eigen::Map<const Eigen::Matrix<int, Eigen::Dynamic, 1>> y_vec(y_seq.data(),
                                                                  y_seq.size());
    logp = y_vec.cwiseEqual(1)
               .select(m_log_1p_exp_cut1,
                       y_vec.cwiseEqual(N_classes).select(
                           m_log_1p_exp_m_cut2,
                           m_log_1p_exp_m_cut2 + log1m_exp(cut1 - cut2).array()
                               + m_log_1p_exp_cut1))
               .sum();
  } else {
    if (y_seq[0] == 1) {
      logp = m_log_1p_exp_cut1.sum();
    } else if (y_seq[0] == N_classes) {
      logp = m_log_1p_exp_m_cut2.sum();
    } else {
      logp = (m_log_1p_exp_m_cut2 + log1m_exp(cut1 - cut2).array()
              + m_log_1p_exp_cut1)
                 .sum();
    }
  }

  auto ops_partials = make_partials_propagator(x_ref, beta_ref, cuts_ref);
  if (!is_constant_all<T_x, T_beta, T_cuts>::value) {
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
    if (!is_constant_all<T_x, T_beta>::value) {
      Matrix<T_partials_return, 1, Dynamic> location_derivative = d1 - d2;
      if (!is_constant_all<T_x>::value) {
        if (T_x_rows == 1) {
          edge<0>(ops_partials).partials_
              = beta_val_vec * location_derivative.sum();
        } else {
          edge<0>(ops_partials).partials_
              = (beta_val_vec * location_derivative).transpose();
        }
      }
      if (!is_constant_all<T_beta>::value) {
        if (T_x_rows == 1) {
          edge<1>(ops_partials).partials_
              = (location_derivative * x_val.replicate(N_instances, 1))
                    .transpose();
        } else {
          edge<1>(ops_partials).partials_
              = (location_derivative * x_val).transpose();
        }
      }
    }
    if (!is_constant_all<T_cuts>::value) {
      for (int i = 0; i < N_instances; i++) {
        int c = y_seq[i];
        if (c != N_classes) {
          partials<2>(ops_partials)[c - 1] += d2.coeff(i);
        }
        if (c != 1) {
          partials<2>(ops_partials)[c - 2] -= d1.coeff(i);
        }
      }
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_x, typename T_beta, typename T_cuts>
return_type_t<T_x, T_beta, T_cuts> ordered_logistic_glm_lpmf(
    const T_y& y, const T_x& x, const T_beta& beta, const T_cuts& cuts) {
  return ordered_logistic_glm_lpmf<false>(y, x, beta, cuts);
}

}  // namespace math
}  // namespace stan

#endif
