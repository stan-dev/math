#ifndef STAN_MATH_OPENCL_PRIM_ORDERED_LOGISTIC_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_ORDERED_LOGISTIC_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/opencl/rev/operands_and_partials.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/opencl/kernels/ordered_logistic_glm_lpmf.hpp>

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log1m_exp.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the log PMF of the ordinal regression Generalized Linear Model (GLM).
 * This is equivalent to and faster than ordered_logistic_lpmf(y, x * beta,
 * cuts).
 * This is an overload of the GLM in
 * prim/prob/ordered_logistic_glm_lpmf.hpp that is implemented in OpenCL.
 *
 * @tparam T_beta type the vector of weights
 * @tparam T_cuts type the vector of cutpoints
 * @param y a scalar or vector of classes on OpenCL device. If it is a scalar
 * it will be broadcast - used for all instances. Values should be between 1 and
 * number of classes, including endpoints.
 * @param x design matrix or row vector on OpenCL device. This overload does
 * not support broadcasting of a row vector x!
 * @param beta weight vector
 * @param cuts cutpoints vector
 * @return log probability
 * @throw std::domain_error If any class is not between 1 and
 * the number of cutpoints plus 2 or if the cutpoint vector is not sorted in
 * ascending order or any input is not finite
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y, typename T_x, typename T_beta,
          typename T_cuts,
          require_all_prim_or_rev_kernel_expression_t<T_y, T_x, T_beta,
                                                      T_cuts>* = nullptr>
return_type_t<T_x, T_beta, T_cuts> ordered_logistic_glm_lpmf(
    const T_y& y, const T_x& x, const T_beta& beta, const T_cuts& cuts) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;
  using std::isfinite;
  using T_partials_return = partials_return_t<T_beta, T_cuts>;
  constexpr bool is_y_vector = !is_stan_scalar<T_y>::value;

  static const char* function = "ordered_logistic_glm_lpmf";

  const size_t N_instances = x.rows();
  const size_t N_attributes = x.cols();
  const size_t N_classes = math::size(cuts) + 1;

  if (is_y_vector) {
    check_size_match(function, "Rows of ", "x", N_instances, "rows of ", "y",
                     math::size(y));
  }
  check_size_match(function, "Columns of ", "x", N_attributes, "Size of",
                   "beta", math::size(beta));

  const auto& cuts_val = eval(value_of(cuts));
  if (N_classes >= 2) {
    auto cuts_head = block_zero_based(cuts_val, 0, 0, math::size(cuts) - 1, 1);
    auto cuts_tail = block_zero_based(cuts_val, 1, 0, math::size(cuts) - 1, 1);
    check_cl(function, "Cuts", cuts_head, "ordered and finite")
        = cuts_head < cuts_tail && isfinite(cuts_head) && isfinite(cuts_tail);
  } else {
    check_cl(function, "Cuts", cuts_val, "finite") = isfinite(cuts_val);
  }

  if (N_instances == 0 || N_classes == 1) {
    return 0;
  }
  if (!include_summand<propto, T_x, T_beta, T_cuts>::value) {
    return 0;
  }

  const auto& y_val = eval(value_of(y));
  const auto& x_val = eval(value_of(x));
  const auto& beta_val = eval(value_of(beta));

  const auto& y_val_cl = to_matrix_cl(y_val);

  const int local_size
      = opencl_kernels::ordered_logistic_glm.get_option("LOCAL_SIZE_");
  const int wgs = (N_instances + local_size - 1) / local_size;

  bool need_location_derivative = !is_constant_all<T_x, T_beta>::value;
  bool need_cuts_derivative = !is_constant_all<T_cuts>::value;
  matrix_cl<double> logp_cl(wgs, 1);
  matrix_cl<double> location_sum_cl(wgs, 1);
  matrix_cl<double> location_derivative_cl(need_location_derivative ? 1 : 0,
                                           N_instances);
  matrix_cl<double> cuts_derivative_cl(N_classes - 1,
                                       need_cuts_derivative ? wgs : 0);

  try {
    opencl_kernels::ordered_logistic_glm(
        cl::NDRange(local_size * wgs), cl::NDRange(local_size), location_sum_cl,
        logp_cl, location_derivative_cl, cuts_derivative_cl, y_val_cl, x_val,
        beta_val, cuts_val, N_instances, N_attributes, N_classes, is_y_vector,
        need_location_derivative, need_cuts_derivative);
  } catch (const cl::Error& e) {
    check_opencl_error(function, e);
  }

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (!std::isfinite(sum(from_matrix_cl(location_sum_cl)))) {
    check_cl(function, "Vector of dependent variables", y_val,
             "between 0 and number of classes")
        = y_val >= 1 && y_val <= static_cast<int>(N_classes);
    check_cl(function, "Design matrix", x_val, "finite") = isfinite(x_val);
    check_cl(function, "Weight vector", beta_val, "finite")
        = isfinite(beta_val);
  }

  auto ops_partials = make_partials_propagator(x, beta, cuts);
  if (!is_constant_all<T_x>::value) {
    partials<0>(ops_partials)
        = transpose(location_derivative_cl) * transpose(beta_val);
  }
  if (!is_constant_all<T_beta>::value) {
    matrix_cl<double> edge2_partials_transpose = location_derivative_cl * x_val;
    partials<1>(ops_partials) = matrix_cl<double>(
        edge2_partials_transpose.buffer(), edge2_partials_transpose.cols(),
        edge2_partials_transpose.rows());
    if (beta.rows() != 0) {
      edge<1>(ops_partials)
          .partials_.add_write_event(
              edge2_partials_transpose.write_events().back());
    }
  }
  if (!is_constant_all<T_cuts>::value) {
    if (wgs == 1) {
      partials<2>(ops_partials) = std::move(cuts_derivative_cl);
    } else {
      partials<2>(ops_partials) = rowwise_sum(cuts_derivative_cl);
    }
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan

#endif
#endif
