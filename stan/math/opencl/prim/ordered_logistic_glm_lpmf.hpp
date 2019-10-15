#ifndef STAN_MATH_OPENCL_PRIM_ORDERED_LOGISTIC_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_ORDERED_LOGISTIC_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/arr/fun/value_of_rec.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/mat/err/check_ordered.hpp>
#include <stan/math/prim/arr/err/check_ordered.hpp>
#include <stan/math/prim/mat/fun/log1m_exp.hpp>
#include <stan/math/prim/mat/fun/sum.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/ordered_logistic_glm_lpmf.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the log PMF of the ordinal regression Generalized Linear Model (GLM).
 * This is equivalent to and faster than ordered_logistic_lpmf(y, x * beta,
 * cuts).
 * This is an overload of the GLM in
 * prim/mar/prob/ordered_logistic_glm_lpmf.hpp that is implemented in OpenCL.
 *
 * @tparam T_beta_scalar type of a scalar in the vector of weights
 * @tparam T_cuts_scalar type of a scalar in the vector of cutpoints
 * @param y_cl a scalar or vector of classes on OpenCL device. If it is a scalar
 * it will be broadcast - used for all instances. Values should be between 1 and
 * number of classes, including endpoints.
 * @param x_cl design matrix or row vector on OpenCL device. This overload does
 * not support broadcasting of a row vector x!
 * @param beta weight vector
 * @param cuts cutpoints vector
 * @return log probability
 * @throw std::domain_error If any class is not between 1 and
 * the number of cutpoints plus 2 or if the cutpoint vector is not sorted in
 * ascending order or any input is not finite
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_beta_scalar, typename T_cuts_scalar>
typename stan::return_type_t<T_beta_scalar, T_cuts_scalar>
ordered_logistic_glm_lpmf(
    const matrix_cl<int>& y_cl, const matrix_cl<double>& x_cl,
    const Eigen::Matrix<T_beta_scalar, Eigen::Dynamic, 1>& beta,
    const Eigen::Matrix<T_cuts_scalar, Eigen::Dynamic, 1>& cuts) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;
  using std::isfinite;
  using T_partials_return =
      typename partials_return_type<T_beta_scalar, T_cuts_scalar>::type;

  static const char* function = "ordered_logistic_glm_lpmf";

  const size_t N_instances = x_cl.rows();
  const size_t N_attributes = x_cl.cols();
  const size_t N_classes = length(cuts) + 1;

  if (y_cl.size() != 1) {
    check_size_match(function, "Rows of ", "x_cl", N_instances, "rows of ",
                     "y_cl", y_cl.size());
  }
  check_consistent_size(function, "Weight vector", beta, N_attributes);
  check_ordered(function, "Cut-points", cuts);
  if (N_classes > 1) {
    if (N_classes > 2) {
      check_finite(function, "Final cut-point", cuts[N_classes - 2]);
    }
    check_finite(function, "First cut-point", cuts[0]);
  }

  if (N_instances == 0 || N_classes == 1) {
    return 0;
  }

  if (!include_summand<propto, T_beta_scalar, T_cuts_scalar>::value) {
    return 0;
  }

  const auto& beta_val = value_of_rec(beta);
  const auto& cuts_val = value_of_rec(cuts);

  const auto& beta_val_vec = as_column_vector_or_scalar(beta_val);
  const auto& cuts_val_vec = as_column_vector_or_scalar(cuts_val);

  operands_and_partials<Eigen::Matrix<T_beta_scalar, Eigen::Dynamic, 1>,
                        Eigen::Matrix<T_cuts_scalar, Eigen::Dynamic, 1>>
      ops_partials(beta, cuts);
  const int local_size
      = opencl_kernels::ordered_logistic_glm.get_option("LOCAL_SIZE_");
  const int wgs = (N_instances + local_size - 1) / local_size;

  const matrix_cl<double> beta_cl(beta_val);
  const matrix_cl<double> cuts_cl(cuts_val);

  bool need_location_derivative = !is_constant_all<T_beta_scalar>::value;
  bool need_cuts_derivative = !is_constant_all<T_cuts_scalar>::value;
  matrix_cl<double> logp_cl(wgs, 1);
  matrix_cl<double> location_sum_cl(wgs, 1);
  matrix_cl<double> location_derivative_cl(need_location_derivative ? 1 : 0,
                                           N_instances);
  matrix_cl<double> cuts_derivative_cl(need_cuts_derivative ? wgs : 0,
                                       N_classes - 1);

  try {
    opencl_kernels::ordered_logistic_glm(
        cl::NDRange(local_size * wgs), cl::NDRange(local_size), location_sum_cl,
        logp_cl, location_derivative_cl, cuts_derivative_cl, y_cl, x_cl,
        beta_cl, cuts_cl, N_instances, N_attributes, N_classes,
        y_cl.size() != 1, need_location_derivative, need_cuts_derivative);
  } catch (const cl::Error& e) {
    check_opencl_error(function, e);
  }

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (!std::isfinite(sum(from_matrix_cl(location_sum_cl)))) {
    check_finite(function, "Weight vector", beta);
    check_bounded(function, "Vector of dependent variables",
                  from_matrix_cl(y_cl), 1, N_classes);
    check_finite(function, "Matrix of independent variables",
                 from_matrix_cl(x_cl));
  }

  if (!is_constant_all<T_beta_scalar>::value) {
    ops_partials.edge1_.partials_
        = from_matrix_cl<1, Dynamic>(location_derivative_cl * x_cl);
  }
  if (!is_constant_all<T_cuts_scalar>::value) {
    ops_partials.edge2_.partials_
        = from_matrix_cl(cuts_derivative_cl).colwise().sum();
  }
  return ops_partials.build(logp);
}

template <typename T_beta_scalar, typename T_cuts_scalar>
typename return_type<T_beta_scalar, T_cuts_scalar>::type
ordered_logistic_glm_lpmf(
    const matrix_cl<int>& y, const matrix_cl<double>& x,
    const Eigen::Matrix<T_beta_scalar, Eigen::Dynamic, 1>& beta,
    const Eigen::Matrix<T_cuts_scalar, Eigen::Dynamic, 1>& cuts) {
  return ordered_logistic_glm_lpmf<false>(y, x, beta, cuts);
}

}  // namespace math
}  // namespace stan

#endif
#endif
