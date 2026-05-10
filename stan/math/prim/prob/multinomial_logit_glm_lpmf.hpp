#ifndef STAN_MATH_PRIM_PROB_MULTINOMIAL_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_PRIM_PROB_MULTINOMIAL_LOGIT_GLM_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/lmultiply.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_matrix.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with multinomial distribution and softmax (logit) link function.
 * Efficiently computes
 * \f$\sum_n \mathrm{multinomial\_logit\_lpmf}(y_n \mid \eta_n)\f$
 * with \f$\eta_n = x_n \beta + \alpha_n\f$ and analytically derived gradients.
 *
 * Writing \f$S_n = \sum_k y_{nk}\f$ for the per-instance count total and
 * \f$p_{nk} = \mathrm{softmax}(\eta_n)_k
 *   = \exp(\eta_{nk}) / \sum_{k'} \exp(\eta_{nk'})\f$
 * for the category probabilities, the log PMF is
 *
 * \f[
 * \log p(y\,|\,x,\alpha,\beta) = \sum_{n=1}^N \left[
 *   \log\Gamma(S_n + 1) - \sum_{k=1}^K \log\Gamma(y_{nk} + 1)
 *   + \sum_{k=1}^K y_{nk}\,\log p_{nk}
 * \right].
 * \f]
 *
 * Defining the residual
 * \f$\delta_{nk} = y_{nk} - S_n\,p_{nk}\f$
 * (so that \f$\partial\log p/\partial\eta_{nk} = \delta_{nk}\f$), the
 * gradients are
 *
 * \f[
 * \frac{\partial\log p}{\partial\alpha_{nk}} = \delta_{nk},
 * \f]
 *
 * \f[
 * \frac{\partial\log p}{\partial\beta_{mk}} = \sum_n x_{nm}\,\delta_{nk},
 * \f]
 *
 * \f[
 * \frac{\partial\log p}{\partial x_{nm}} = \sum_k \delta_{nk}\,\beta_{mk}.
 * \f]
 *
 * When \f$\alpha\f$ is broadcast (a single \f$1\times K\f$ row used for all
 * instances), \f$\partial\log p/\partial\alpha_k = \sum_n \delta_{nk}\f$.
 *
 * @tparam T_x type of the design matrix (N x M)
 * @tparam T_alpha type of the intercept; either a row vector (1 x K) broadcast
 * over all instances, or a matrix (N x K) with per-instance intercepts.
 * Analogous to a scalar vs vector alpha in univariate GLMs: the K dimension
 * corresponds to the class (output) space, not the instance dimension.
 * @tparam T_beta type of the weight matrix (M x K)
 * @param y outcome count vectors: `y[n]` is a length-K vector of non-negative
 * integer counts for instance n; the counts need not sum to a fixed total
 * @param x design matrix (N x M)
 * @param alpha intercept; row vector (1 x K) broadcast across all instances,
 * or matrix (N x K) with one bias row per instance
 * @param beta weight matrix (M x K)
 * @return log sum of multinomial log PMFs over all N instances
 * @throw std::domain_error if any element of x or beta is infinite or NaN, or
 * if alpha contains `+inf` or NaN (`-inf` forces the corresponding softmax
 * probability to zero and is allowed; if all classes have `-inf` alpha for
 * some instance but that instance has nonzero counts, the result is undefined),
 * or if any count in y is negative
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_x, typename T_alpha, typename T_beta,
          require_eigen_matrix_dynamic_t<T_x>* = nullptr,
          require_matrix_t<T_alpha>* = nullptr,
          require_matrix_t<T_beta>* = nullptr>
inline return_type_t<T_x, T_alpha, T_beta> multinomial_logit_glm_lpmf(
    const std::vector<std::vector<int>>& y, const T_x& x, const T_alpha& alpha,
    const T_beta& beta) {
  using T_partials_return = partials_return_t<T_x, T_alpha, T_beta>;
  using Eigen::Array;
  using Eigen::Dynamic;
  using T_x_ref = ref_type_if_not_constant_t<T_x>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_alpha>;
  using T_beta_ref = ref_type_if_not_constant_t<T_beta>;
  constexpr int T_alpha_rows = T_alpha::RowsAtCompileTime;
  constexpr bool need_delta = is_any_autodiff_v<T_x, T_alpha, T_beta>;

  const size_t N_instances = x.rows();
  const size_t N_classes = beta.cols();

  static constexpr const char* function = "multinomial_logit_glm_lpmf";
  check_size_match(function, "Rows of outcome vectors", y.size(),
                   "number of instances", N_instances);
  check_size_match(function, "Columns of intercept", alpha.cols(),
                   "number of classes", N_classes);
  if constexpr (T_alpha_rows != 1) {
    check_size_match(function, "Rows of intercept", alpha.rows(),
                     "rows of design matrix", x.rows());
  }
  check_size_match(function, "Columns of design matrix", x.cols(),
                   "rows of weight matrix", beta.rows());

  if (size_zero(y)) {
    return 0;
  }
  for (size_t n = 0; n < N_instances; ++n) {
    check_size_match(function, "Size of outcome vector", y[n].size(),
                     "number of classes", N_classes);
    check_nonnegative(function, "outcome counts", y[n]);
  }

  if constexpr (!include_summand<propto, T_x, T_alpha, T_beta>::value) {
    return 0;
  }

  T_x_ref x_ref = x;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  const auto& x_val = to_ref_if<is_autodiff_v<T_beta>>(value_of(x_ref));
  const auto& beta_val = to_ref_if<is_autodiff_v<T_x>>(value_of(beta_ref));

  const Array<T_partials_return, Dynamic, Dynamic> eta = [&]() {
    if constexpr (T_alpha_rows == 1) {
      return ((x_val * beta_val).rowwise() + value_of(alpha_ref)).array();
    } else {
      return (x_val * beta_val + value_of(alpha_ref)).array();
    }
  }();

  // Row-max shift for numerical stability; cancels in log-softmax.
  const Array<T_partials_return, Dynamic, Dynamic> shifted_eta
      = (eta.colwise() - eta.rowwise().maxCoeff()).eval();
  const Array<T_partials_return, Dynamic, Dynamic> exp_eta = exp(shifted_eta);
  const Array<T_partials_return, Dynamic, Dynamic> softmax_mat
      = exp_eta.colwise() / exp_eta.rowwise().sum();

  const Array<double, Dynamic, Dynamic> y_mat = to_matrix(y).array();
  const Array<double, Dynamic, 1> instance_totals = y_mat.rowwise().sum();

  // lmultiply implements 0*log(0)=0: classes with softmax=0 and y=0 contribute 0.
  T_partials_return logp = lmultiply(y_mat, softmax_mat).sum();
  if constexpr (include_summand<propto>::value) {
    logp += lgamma(instance_totals + 1.0).sum() - lgamma(y_mat + 1.0).sum();
  }

  if (!std::isfinite(logp)) {
    check_finite(function, "Weight matrix", beta_ref);
    check_finite(function, "Intercept", alpha_ref);
    check_finite(function, "Matrix of independent variables", x_ref);
  }

  auto ops_partials = make_partials_propagator(x_ref, alpha_ref, beta_ref);
  if constexpr (need_delta) {
    // δ[n,k] = y[n,k] - S_n·p[n,k]
    const Array<T_partials_return, Dynamic, Dynamic> delta
        = y_mat.template cast<T_partials_return>()
          - softmax_mat.colwise()
                * instance_totals.template cast<T_partials_return>();

    if constexpr (is_autodiff_v<T_alpha>) {
      if constexpr (T_alpha_rows == 1)
        partials<1>(ops_partials) = delta.colwise().sum();
      else
        partials<1>(ops_partials) = delta;
    }
    if constexpr (is_autodiff_v<T_beta>) {
      partials<2>(ops_partials)
          = x_val.transpose().template cast<T_partials_return>() * delta.matrix();
    }
    if constexpr (is_autodiff_v<T_x>) {
      edge<0>(ops_partials).partials_ = delta.matrix() * beta_val.transpose();
    }
  }
  return ops_partials.build(logp);
}

template <typename T_x, typename T_alpha, typename T_beta,
          require_eigen_matrix_dynamic_t<T_x>* = nullptr,
          require_matrix_t<T_alpha>* = nullptr,
          require_matrix_t<T_beta>* = nullptr>
inline return_type_t<T_x, T_alpha, T_beta> multinomial_logit_glm_lpmf(
    const std::vector<std::vector<int>>& y, const T_x& x, const T_alpha& alpha,
    const T_beta& beta) {
  return multinomial_logit_glm_lpmf<false>(y, x, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
