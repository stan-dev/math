#ifndef STAN_MATH_OPENCL_PRIM_ORDERED_LOGISTIC_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_ORDERED_LOGISTIC_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/kernels/add.hpp>
#include <stan/math/opencl/kernels/ordered_logistic_lpmf.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the (natural) log probability of the specified array
 * of integers given the vector of continuous locations and
 * specified cutpoints in an ordered logistic model.
 *
 * <p>Typically the continuous lambda
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
 * @tparam T_loc lambda type.
 * @tparam T_cut Cut-point type.
 * @param y Array of integers
 * @param lambda Vector of continuous lambda variables.
 * @param cuts Positive increasing vector of cutpoints.
 * @return Log probability of outcome given lambda and
 * cutpoints.
 * @throw std::domain_error If the outcome is not between 1 and
 * the number of cutpoints plus 2; if the cutpoint vector is
 * empty; if the cutpoint vector contains a non-positive,
 * non-finite value; or if the cutpoint vector is not sorted in
 * ascending order.
 * @throw std::invalid_argument If y and lambda are different
 * lengths.
 */
template <bool propto, typename T_y_cl, typename T_loc_cl, typename T_cuts_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                      T_cuts_cl>* = nullptr>
inline return_type_t<T_y_cl, T_loc_cl, T_cuts_cl> ordered_logistic_lpmf(
    const T_y_cl& y, const T_loc_cl& lambda, const T_cuts_cl& cuts) {
  constexpr bool is_y_vector = !is_stan_scalar<T_y_cl>::value;
  static const char* function = "ordered_logistic_lpmf(OpenCL)";

  if (size(y) != 1) {
    check_size_match(function, "Size of ", "y", math::size(y), "Size of",
                     "lambda", math::size(lambda));
  }

  int N_instances = max_size(y, lambda);
  int N_classes = cuts.rows() + 1;
  int N_cut_sets = cuts.cols();

  if (N_cut_sets > 1) {
    check_size_match(function, "Length of lambda variables ", N_instances,
                     "Number of cutpoint vectors ", N_cut_sets);
  }
  if (N_instances == 0 || N_classes == 1) {
    return 0.0;
  }
  const auto& cuts_val = eval(value_of(cuts));
  if (N_classes >= 2) {
    auto cuts_head
        = block_zero_based(cuts_val, 0, 0, cuts.rows() - 1, N_cut_sets);
    auto cuts_tail
        = block_zero_based(cuts_val, 1, 0, cuts.rows() - 1, N_cut_sets);
    check_cl(function, "Cuts", cuts_head, "ordered and finite")
        = cuts_head < cuts_tail && isfinite(cuts_head) && isfinite(cuts_tail);
  } else if (N_classes == 1) {
    check_cl(function, "Cuts", cuts_val, "finite") = isfinite(cuts_val);
  }

  if (!include_summand<propto, T_loc_cl, T_cuts_cl>::value) {
    return 0.0;
  }

  const auto& y_val = eval(value_of(y));
  const auto& lambda_val = eval(value_of(lambda));

  const auto& y_val_cl = to_matrix_cl(y_val);

  const int local_size
      = opencl_kernels::ordered_logistic.get_option("LOCAL_SIZE_");
  const int wgs = (N_instances + local_size - 1) / local_size;

  bool need_lambda_derivative = !is_constant_all<T_loc_cl>::value;
  bool need_cuts_derivative = !is_constant_all<T_cuts_cl>::value;
  bool need_broadcasting = N_cut_sets == 1 && N_instances != 1;
  matrix_cl<double> logp_cl(wgs, 1);
  matrix_cl<double> lambda_derivative_cl(N_instances,
                                         need_lambda_derivative ? 1 : 0);
  matrix_cl<double> cuts_derivative_cl(
      N_classes - 1,
      need_cuts_derivative ? (need_broadcasting ? wgs : N_cut_sets) : 0);

  try {
    opencl_kernels::ordered_logistic(
        cl::NDRange(local_size * wgs), cl::NDRange(local_size), logp_cl,
        lambda_derivative_cl, cuts_derivative_cl, y_val_cl, lambda_val,
        cuts_val, N_instances, N_classes, is_y_vector, !need_broadcasting,
        need_lambda_derivative, need_cuts_derivative);
  } catch (const cl::Error& e) {
    check_opencl_error(function, e);
  }

  double logp = sum(from_matrix_cl(logp_cl));

  if (!std::isfinite(logp)) {
    results(check_cl(function, "Vector of dependent variables", y_val,
                     "between 0 and number of classes"),
            check_cl(function, "lambda vector", lambda_val, "finite"))
        = expressions(y_val >= 1 && y_val <= static_cast<int>(N_classes),
                      isfinite(lambda_val));
  }
  operands_and_partials<T_loc_cl, T_cuts_cl> ops_partials(lambda, cuts);

  if (!is_constant_all<T_loc_cl>::value) {
    ops_partials.edge1_.partials_ = lambda_derivative_cl;
  }
  if (!is_constant_all<T_cuts_cl>::value) {
    if (need_broadcasting) {
      ops_partials.edge2_.partials_ = rowwise_sum(cuts_derivative_cl);
    } else {
      ops_partials.edge2_.partials_ = std::move(cuts_derivative_cl);
    }
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
