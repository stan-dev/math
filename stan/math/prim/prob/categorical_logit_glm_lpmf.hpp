#ifndef STAN_MATH_PRIM_PROB_CATEGORICAL_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_PRIM_PROB_CATEGORICAL_LOGIT_GLM_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/isfinite.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with categorical distribution and logit (softmax) link function.
 *
 * @tparam T_y type of classes. It can be either `std::vector<int>` or `int`.
 * @tparam T_x_scalar type of the matrix of independent variables (features)
 * @tparam T_alpha type of the intercept vector
 * @tparam T_beta type of the matrix of weights
 * @param y a scalar or vector of classes. If it is a scalar it will be
 * broadcast - used for all instances. Values should be between 1 and number
 * of classes, including endpoints.
 * @param x design matrix or row vector. If it is a row vector it will be
 * broadcast - used for all instances.
 * @param alpha intercept vector (in log odds)
 * @param beta weight matrix
 * @return log probability or log sum of probabilities
 * @throw std::domain_error x, beta or alpha is infinite or y is not within
 * bounds
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y, typename T_x, typename T_alpha,
          typename T_beta, require_matrix_t<T_x>* = nullptr,
          require_col_vector_t<T_alpha>* = nullptr,
          require_matrix_t<T_beta>* = nullptr>
return_type_t<T_x, T_alpha, T_beta> categorical_logit_glm_lpmf(
    const T_y& y, const T_x& x, const T_alpha& alpha, const T_beta& beta) {
  using T_partials_return = partials_return_t<T_x, T_alpha, T_beta>;
  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;
  using std::isfinite;
  using std::log;
  using T_y_ref = ref_type_t<T_y>;
  using T_x_ref = ref_type_if_not_constant_t<T_x>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_alpha>;
  using T_beta_ref = ref_type_if_not_constant_t<T_beta>;
  using T_beta_partials = partials_type_t<scalar_type_t<T_beta>>;
  constexpr int T_x_rows = T_x::RowsAtCompileTime;

  const size_t N_instances = T_x_rows == 1 ? stan::math::size(y) : x.rows();
  const size_t N_attributes = x.cols();
  const size_t N_classes = beta.cols();

  static const char* function = "categorical_logit_glm_lpmf";
  check_consistent_size(function, "Vector of dependent variables", y,
                        N_instances);
  check_consistent_size(function, "Intercept vector", alpha, N_classes);
  check_size_match(function, "x.cols()", N_attributes, "beta.rows()",
                   beta.rows());
  if (size_zero(y) || N_classes == 1) {
    return 0;
  }
  T_y_ref y_ref = y;
  check_bounded(function, "categorical outcome out of support", y_ref, 1,
                N_classes);

  if (!include_summand<propto, T_x, T_alpha, T_beta>::value) {
    return 0;
  }

  T_x_ref x_ref = x;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  const auto& x_val = to_ref_if<!is_constant<T_beta>::value>(value_of(x_ref));
  const auto& alpha_val = value_of(alpha_ref);
  const auto& beta_val
      = to_ref_if<!is_constant<T_x>::value>(value_of(beta_ref));

  const auto& alpha_val_vec = as_column_vector_or_scalar(alpha_val).transpose();

  Array<T_partials_return, T_x_rows, Dynamic> lin
      = (x_val * beta_val).rowwise() + alpha_val_vec;
  Array<T_partials_return, T_x_rows, 1> lin_max
      = lin.rowwise().maxCoeff();  // This is used to prevent overflow when
                                   // calculating softmax/log_sum_exp and
                                   // similar expressions
  Array<T_partials_return, T_x_rows, Dynamic> exp_lin
      = exp(lin.colwise() - lin_max);
  Array<T_partials_return, T_x_rows, 1> inv_sum_exp_lin
      = 1 / exp_lin.rowwise().sum();

  T_partials_return logp = log(inv_sum_exp_lin).sum() - lin_max.sum();
  if (T_x_rows == 1) {
    logp *= N_instances;
  }
  scalar_seq_view<T_y_ref> y_seq(y_ref);
  for (int i = 0; i < N_instances; i++) {
    if (T_x_rows == 1) {
      logp += lin(0, y_seq[i] - 1);
    } else {
      logp += lin(i, y_seq[i] - 1);
    }
  }
  // TODO(Tadej) maybe we can replace previous block with the following line
  // when we have newer Eigen  T_partials_return logp =
  // lin(Eigen::all,y-1).sum() + log(inv_sum_exp_lin).sum() - lin_max.sum();

  if (!isfinite(logp)) {
    check_finite(function, "Weight vector", beta_ref);
    check_finite(function, "Intercept", alpha_ref);
    check_finite(function, "Matrix of independent variables", x_ref);
  }

  // Compute the derivatives.
  auto ops_partials = make_partials_propagator(x_ref, alpha_ref, beta_ref);

  if (!is_constant_all<T_x>::value) {
    if (T_x_rows == 1) {
      Array<T_beta_partials, 1, Dynamic> beta_y = beta_val.col(y_seq[0] - 1);
      for (int i = 1; i < N_instances; i++) {
        beta_y += beta_val.col(y_seq[i] - 1).array();
      }
      edge<0>(ops_partials).partials_
          = beta_y
            - (exp_lin.matrix() * beta_val.transpose()).array().colwise()
                  * inv_sum_exp_lin * N_instances;
    } else {
      Array<T_beta_partials, Dynamic, Dynamic> beta_y(N_instances,
                                                      N_attributes);
      for (int i = 0; i < N_instances; i++) {
        beta_y.row(i) = beta_val.col(y_seq[i] - 1);
      }
      edge<0>(ops_partials).partials_
          = beta_y
            - (exp_lin.matrix() * beta_val.transpose()).array().colwise()
                  * inv_sum_exp_lin;
      // TODO(Tadej) maybe we can replace previous block with the following
      // line when we have newer Eigen  partials<0>(ops_partials) = beta_val(y
      // - 1, all) - (exp_lin.matrix() * beta.transpose()).colwise() *
      // inv_sum_exp_lin;
    }
  }
  if (!is_constant_all<T_alpha, T_beta>::value) {
    Array<T_partials_return, T_x_rows, Dynamic> neg_softmax_lin
        = exp_lin.colwise() * -inv_sum_exp_lin;
    if (!is_constant_all<T_alpha>::value) {
      if (T_x_rows == 1) {
        edge<1>(ops_partials).partials_
            = neg_softmax_lin.colwise().sum() * N_instances;
      } else {
        partials<1>(ops_partials) = neg_softmax_lin.colwise().sum();
      }
      for (int i = 0; i < N_instances; i++) {
        partials<1>(ops_partials)[y_seq[i] - 1] += 1;
      }
    }
    if (!is_constant_all<T_beta>::value) {
      Matrix<T_partials_return, Dynamic, Dynamic> beta_derivative
          = x_val.transpose().template cast<T_partials_return>()
            * neg_softmax_lin.matrix();
      if (T_x_rows == 1) {
        beta_derivative *= N_instances;
      }

      for (int i = 0; i < N_instances; i++) {
        if (T_x_rows == 1) {
          beta_derivative.col(y_seq[i] - 1) += x_val;
        } else {
          beta_derivative.col(y_seq[i] - 1) += x_val.row(i);
        }
      }
      // TODO(Tadej) maybe we can replace previous loop with the following
      // line when we have newer Eigen  partials<2>(ops_partials)(Eigen::all,
      // y
      // - 1) += x_val.colwise.sum().transpose();

      partials<2>(ops_partials) = std::move(beta_derivative);
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_x, typename T_alpha, typename T_beta>
return_type_t<T_x, T_alpha, T_beta> categorical_logit_glm_lpmf(
    const T_y& y, const T_x& x, const T_alpha& alpha, const T_beta& beta) {
  return categorical_logit_glm_lpmf<false>(y, x, alpha, beta);
}

}  // namespace math
}  // namespace stan

#endif
