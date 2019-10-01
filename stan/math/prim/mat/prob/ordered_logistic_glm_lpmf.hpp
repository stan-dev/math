#ifndef STAN_MATH_PRIM_MAT_PROB_ORDERED_LOGISTIC_GLM_LPMF_HPP
#define STAN_MATH_PRIM_MAT_PROB_ORDERED_LOGISTIC_GLM_LPMF_HPP

#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/arr/fun/value_of_rec.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/mat/err/check_ordered.hpp>
#include <stan/math/prim/arr/err/check_ordered.hpp>
#include <stan/math/prim/mat/fun/log1m_exp.hpp>
#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the log PMF of the ordinal regression Generalized Linear Model (GLM).
 * This is equivalent to and faster than ordered_logistic_lpmf(y, x * beta,
 * cuts).
 *
 * @tparam T_y type of integer vector of classes. It can be either
 * `std::vector<int>` or `int`.
 * @tparam T_x_scalar type of elements in the matrix of independent variables
 * (features)
 * @tparam T_x_rows compile-time number of rows of `x`. It can be either
 * `Eigen::Dynamic` or 1.
 * @tparam T_beta_scalar type of a scalar in the vector of weights
 * @tparam T_cuts_scalar type of a scalar in the vector of cutpoints
 * @param y a scalar or vector of classes. If it is a scalar it will be
 * broadcast - used for all instances. Values should be between 1 and number of
 * classes, including endpoints.
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
template <bool propto, typename T_y, typename T_x_scalar, int T_x_rows,
          typename T_beta_scalar, typename T_cuts_scalar>
typename stan::return_type_t<T_x_scalar, T_beta_scalar, T_cuts_scalar>
ordered_logistic_glm_lpmf(
    const T_y& y, const Eigen::Matrix<T_x_scalar, T_x_rows, Eigen::Dynamic>& x,
    const Eigen::Matrix<T_beta_scalar, Eigen::Dynamic, 1>& beta,
    const Eigen::Matrix<T_cuts_scalar, Eigen::Dynamic, 1>& cuts) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;
  using std::exp;
  using std::isfinite;

  typedef typename partials_return_type<T_y, T_x_scalar, T_beta_scalar,
                                        T_cuts_scalar>::type T_partials_return;
  typedef typename std::conditional_t<T_x_rows == 1, double,
                                      Array<double, Dynamic, 1>>
      T_location;

  static const char* function = "ordered_logistic_glm_lpmf";

  const size_t N_instances = T_x_rows == 1 ? length(y) : x.rows();
  const size_t N_attributes = x.cols();
  const size_t N_classes = length(cuts) + 1;

  if (is_vector<T_y>::value && T_x_rows != 1) {
    check_consistent_size(function, "Vector of dependent variables", y,
                          N_instances);
  }
  check_consistent_size(function, "Weight vector", beta, N_attributes);
  check_bounded(function, "Vector of dependent variables", y, 1, N_classes);
  check_ordered(function, "Cut-points", cuts);
  if (N_classes > 1) {
    if (N_classes > 2) {
      check_finite(function, "Final cut-point", cuts[N_classes - 2]);
    }
    check_finite(function, "First cut-point", cuts[0]);
  }

  if (size_zero(y, cuts))
    return 0;

  if (!include_summand<propto, T_x_scalar, T_beta_scalar, T_cuts_scalar>::value)
    return 0;

  const auto& x_val = value_of_rec(x);
  const auto& beta_val = value_of_rec(beta);
  const auto& cuts_val = value_of_rec(cuts);

  const auto& beta_val_vec = as_column_vector_or_scalar(beta_val);
  const auto& cuts_val_vec = as_column_vector_or_scalar(cuts_val);

  scalar_seq_view<T_y> y_seq(y);
  Array<double, Dynamic, 1> cuts_y1(N_instances), cuts_y2(N_instances);
  for (int i = 0; i < N_instances; i++) {
    int c = y_seq[i];
    if (c != N_classes) {
      cuts_y1[i] = cuts_val_vec[c - 1];
    } else {
      cuts_y1[i] = INFINITY;
    }
    if (c != 1) {
      cuts_y2[i] = cuts_val_vec[c - 2];
    } else {
      cuts_y2[i] = -INFINITY;
    }
  }

  T_location location = x_val * beta_val_vec;
  if (!isfinite(sum(location))) {
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Matrix of independent variables", x);
  }

  Array<double, Dynamic, 1> cut2 = location - cuts_y2;
  Array<double, Dynamic, 1> cut1 = location - cuts_y1;

  // Not immediately evaluating next two expressions benefits performance
  auto m_log_1p_exp_cut1
      = (cut1 > 0.0).select(-cut1, 0) - (-cut1.abs()).exp().log1p();
  auto m_log_1p_exp_m_cut2
      = (cut2 <= 0.0).select(cut2, 0) - (-cut2.abs()).exp().log1p();

  T_partials_return logp(0);
  if (is_vector<T_y>::value) {
    Eigen::Map<const Eigen::Matrix<int, Eigen::Dynamic, 1>> y_vec(&y_seq[0],
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

  operands_and_partials<Matrix<T_x_scalar, T_x_rows, Dynamic>,
                        Eigen::Matrix<T_beta_scalar, Eigen::Dynamic, 1>,
                        Eigen::Matrix<T_cuts_scalar, Eigen::Dynamic, 1>>
      ops_partials(x, beta, cuts);
  if (!is_constant_all<T_x_scalar, T_beta_scalar, T_cuts_scalar>::value) {
    Array<double, Dynamic, 1> exp_m_cut1 = exp(-cut1);
    Array<double, Dynamic, 1> exp_m_cut2 = exp(-cut2);
    Array<double, Dynamic, 1> exp_cuts_diff = exp(cuts_y2 - cuts_y1);
    Array<double, Dynamic, 1> d1
        = (cut2 > 0).select(exp_m_cut2 / (1 + exp_m_cut2), 1 / (1 + exp(cut2)))
          - exp_cuts_diff / (exp_cuts_diff - 1);
    Array<double, Dynamic, 1> d2
        = 1 / (1 - exp_cuts_diff)
          - (cut1 > 0).select(exp_m_cut1 / (1 + exp_m_cut1),
                              1 / (1 + exp(cut1)));
    if (!is_constant_all<T_x_scalar, T_beta_scalar>::value) {
      Matrix<double, 1, Dynamic> location_derivative = d1 - d2;
      if (!is_constant_all<T_x_scalar>::value) {
        if (T_x_rows == 1) {
          ops_partials.edge1_.partials_
              = beta_val_vec * location_derivative.sum();
        } else {
          ops_partials.edge1_.partials_
              = (beta_val_vec * location_derivative).transpose();
        }
      }
      if (!is_constant_all<T_beta_scalar>::value) {
        if (T_x_rows == 1) {
          ops_partials.edge2_.partials_
              = (location_derivative * x_val.replicate(N_instances, 1))
                    .transpose();
        } else {
          ops_partials.edge2_.partials_
              = (location_derivative * x_val).transpose();
        }
      }
    }
    if (!is_constant_all<T_cuts_scalar>::value) {
      for (int i = 0; i < N_instances; i++) {
        int c = y_seq[i];
        if (c != N_classes) {
          ops_partials.edge3_.partials_[c - 1] += d2[i];
        }
        if (c != 1) {
          ops_partials.edge3_.partials_[c - 2] -= d1[i];
        }
      }
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_x_scalar, int T_x_rows,
          typename T_beta_scalar, typename T_cuts_scalar>
typename stan::return_type_t<T_x_scalar, T_beta_scalar, T_cuts_scalar>
ordered_logistic_glm_lpmf(
    const T_y& y, const Eigen::Matrix<T_x_scalar, T_x_rows, Eigen::Dynamic>& x,
    const Eigen::Matrix<T_beta_scalar, Eigen::Dynamic, 1>& beta,
    const Eigen::Matrix<T_cuts_scalar, Eigen::Dynamic, 1>& cuts) {
  return ordered_logistic_glm_lpmf<false>(y, x, beta, cuts);
}

}  // namespace math
}  // namespace stan

#endif
