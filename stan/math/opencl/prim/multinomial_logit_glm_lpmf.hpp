#ifndef STAN_MATH_OPENCL_PRIM_MULTINOMIAL_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_MULTINOMIAL_LOGIT_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/opencl/rev/operands_and_partials.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/kernels/multinomial_logit_glm_lpmf.hpp>

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

#include <vector>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with multinomial distribution and softmax (logit) link function.
 * This is an OpenCL overload of
 * `prim/prob/multinomial_logit_glm_lpmf.hpp`.
 * Alpha can be either a shared 1×K row vector or an N×K per-instance matrix.
 *
 * @tparam T_x type of the design matrix (N×M kernel expression)
 * @tparam T_alpha type of the intercept (1×K or N×K kernel expression)
 * @tparam T_beta type of the weight matrix (M×K kernel expression)
 * @param y outcome count vectors: `y[n]` is a length-K vector of non-negative
 *   integer counts for instance n
 * @param x design matrix (N×M) on OpenCL device
 * @param alpha intercept: 1×K broadcast row or N×K per-instance matrix
 * @param beta weight matrix (M×K) on OpenCL device
 * @return log sum of multinomial log PMFs over all N instances
 * @throw std::domain_error if any element of x or beta is infinite or NaN,
 * or if alpha contains `+inf` or NaN (`-inf` forces the corresponding softmax
 * probability to zero and is allowed)
 * @throw std::domain_error if any count in y is negative
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_x, typename T_alpha, typename T_beta,
          require_all_prim_or_rev_kernel_expression_t<T_x, T_alpha,
                                                      T_beta>* = nullptr>
inline return_type_t<T_x, T_alpha, T_beta> multinomial_logit_glm_lpmf(
    const std::vector<std::vector<int>>& y, T_x&& x, T_alpha&& alpha,
    T_beta&& beta) {
  using T_partials_return = partials_return_t<T_x, T_alpha, T_beta>;
  static constexpr const char* function = "multinomial_logit_glm_lpmf";

  const int N_instances = x.rows();
  const int N_classes = beta.cols();

  check_size_match(function, "Rows of", "x", N_instances, "size of", "y",
                   y.size());
  check_size_match(function, "Columns of", "beta", N_classes, "columns of",
                   "alpha", alpha.cols());
  check_size_match(function, "Columns of", "x", x.cols(), "rows of", "beta",
                   beta.rows());

  const int alpha_rows = alpha.rows();
  const bool is_alpha_vector = alpha_rows == 1;
  if (!is_alpha_vector) {
    check_size_match(function, "Rows of", "alpha", alpha_rows, "rows of", "x",
                     N_instances);
  }

  if (N_instances == 0) {
    return 0;
  }
  for (int n = 0; n < N_instances; ++n) {
    check_size_match(function, "Size of outcome vector", y[n].size(),
                     "number of classes", N_classes);
    check_nonnegative(function, "outcome counts", y[n]);
  }

  if constexpr (!include_summand<propto, T_x, T_alpha, T_beta>::value) {
    return 0;
  }

  Eigen::MatrixXi y_mat(N_instances, N_classes);
  for (int n = 0; n < N_instances; ++n)
    for (int k = 0; k < N_classes; ++k)
      y_mat(n, k) = y[n][k];
  matrix_cl<int> y_cl(y_mat);

  auto x_val = eval(value_of(x));
  auto alpha_val = eval(value_of(alpha));
  auto beta_val = eval(value_of(beta));

  matrix_cl<double> x_beta_cl = x_val * beta_val;

  const int local_size
      = opencl_kernels::multinomial_logit_glm.get_option("LOCAL_SIZE_");
  const int wgs = (N_instances + local_size - 1) / local_size;

  constexpr bool need_delta = is_any_autodiff_v<T_x, T_alpha, T_beta>;

  matrix_cl<double> logp_cl(wgs, 1);
  matrix_cl<double> delta_cl(0, 0);
  if constexpr (need_delta)
    delta_cl = matrix_cl<double>(N_instances, N_classes);

  try {
    opencl_kernels::multinomial_logit_glm(
        cl::NDRange(local_size * wgs), cl::NDRange(local_size), logp_cl,
        delta_cl, y_cl, x_beta_cl, alpha_val, N_instances, N_classes,
        is_alpha_vector, need_delta, !propto);
  } catch (const cl::Error& e) {
    check_opencl_error(function, e);
  }

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (!std::isfinite(logp)) {
    check_cl(function, "Design matrix", x_val, "finite") = isfinite(x_val);
    check_cl(function, "Intercept", alpha_val, "finite") = isfinite(alpha_val);
    check_cl(function, "Weight matrix", beta_val, "finite")
        = isfinite(beta_val);
  }

  auto ops_partials = make_partials_propagator(std::forward<T_x>(x),
                                               std::forward<T_alpha>(alpha),
                                               std::forward<T_beta>(beta));
  if constexpr (need_delta) {
    if constexpr (is_autodiff_v<T_x>)
      partials<0>(ops_partials) = delta_cl * transpose(beta_val);
    if constexpr (is_autodiff_v<T_alpha>)
      partials<1>(ops_partials)
          = is_alpha_vector ? colwise_sum(delta_cl) : delta_cl;
    if constexpr (is_autodiff_v<T_beta>)
      partials<2>(ops_partials) = transpose(x_val) * delta_cl;
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan

#endif
#endif
