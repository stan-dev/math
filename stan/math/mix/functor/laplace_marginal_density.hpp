
#ifndef STAN_MATH_MIX_FUNCTOR_LAPLACE_MARGINAL_DENSITY_HPP
#define STAN_MATH_MIX_FUNCTOR_LAPLACE_MARGINAL_DENSITY_HPP
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/functor/wolfe_line_search.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/quad_form_diag.hpp>
#include <stan/math/prim/functor/iter_tuple_nested.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>

/**
 * @file
 * Reference for calculations of marginal and its gradients:
 * Margossian et al (2020), https://arxiv.org/abs/2004.12550
 * and Margossian (2023), https://arxiv.org/pdf/2306.14976
 */

namespace stan {
namespace math {

/**
 * Options for the laplace sampler
 */
struct laplace_options_base {
  /* Size of the blocks in block diagonal hessian*/
  int hessian_block_size{1};
  /**
   * Which Newton solver to use: (B matrix in equation 1 of
   * https://arxiv.org/pdf/2306.14976) (1) method using the cholesky
   * decomposition of `W` (the negative Hessian of log likelihood) (2) method
   * using the cholesky decomposition of `K` (the covariance matrix) (3) method
   * using an LU decomposition (more general, but slower)
   */
  int solver{1};
  /**
   * iterations end when difference in objective function is less than tolerance
   * Default is sqrt(machine_epsilon)
   */
  double tolerance{1.49012e-08};
  /* Maximum number of steps*/
  int max_num_steps{500};
  laplace_line_search_options line_search;
};

template <bool HasInitTheta>
struct laplace_options;

template <>
struct laplace_options<false> : public laplace_options_base {};

template <>
struct laplace_options<true> : public laplace_options_base {
  /* Value for user supplied initial theta  */
  Eigen::VectorXd theta_0{0};
};

using laplace_options_default = laplace_options<false>;
using laplace_options_user_supplied = laplace_options<true>;
namespace internal {

template <typename ThetaVec, typename WR, typename L_t, typename A_vec,
          typename ThetaGrad, typename LU_t, typename KRoot>
struct laplace_density_estimates {
  /* log marginal density */
  double lmd{std::numeric_limits<double>::infinity()};
  /* ThetaVec at the mode */
  ThetaVec theta;
  /* negative hessian or sqrt of negative hessian */
  WR W_r;
  /* Lower left of cholesky decomposition of stabilized inverse covariance */
  L_t L;
  /* inverse covariance times theta at the mode */
  A_vec a;
  /* the gradient of the log density with respect to theta */
  ThetaGrad theta_grad;
  /* LU matrix from solver 3 */
  LU_t LU;
  /* Cholesky of the covariance matrix */
  KRoot K_root;
  int solver_used{1};
  laplace_density_estimates(double lmd_, ThetaVec&& theta_, WR&& W_r_, L_t&& L_,
                            A_vec&& a_, ThetaGrad&& theta_grad_, LU_t&& LU_,
                            KRoot&& K_root_, int solver_used_)
      : lmd(lmd_),
        theta(std::move(theta_)),
        W_r(std::move(W_r_)),
        L(std::move(L_)),
        a(std::move(a_)),
        theta_grad(std::move(theta_grad_)),
        LU(std::move(LU_)),
        K_root(std::move(K_root_)),
        solver_used(solver_used_) {}
};

/**
 * Returns the principal square root of a block diagonal matrix.
 * @tparam WRootMat A type inheriting from `Eigen::EigenBase`.
 * @param W_root The output matrix to store the square root.
 * @param W The input block diagonal matrix.
 * @param block_size The size of each block in the block diagonal matrix.
 */
template <typename WRootMat>
inline void block_matrix_sqrt(WRootMat& W_root,
                              const Eigen::SparseMatrix<double>& W,
                              const Eigen::Index block_size) {
  int n_block = W.cols() / block_size;
  Eigen::MatrixXd local_block(block_size, block_size);
  Eigen::MatrixXd local_block_sqrt(block_size, block_size);
  Eigen::MatrixXd sqrt_t_mat = Eigen::MatrixXd::Zero(block_size, block_size);
  // No block operation available for sparse matrices, so we have to loop
  // See https://eigen.tuxfamily.org/dox/group__TutorialSparse.html#title7
  for (int i = 0; i < n_block; i++) {
    sqrt_t_mat.setZero();
    local_block
        = W.block(i * block_size, i * block_size, block_size, block_size);
    if (Eigen::isnan(local_block.array()).any()) {
      throw std::domain_error(
          std::string("Error in block_matrix_sqrt: "
                      "NaNs detected in block diagonal starting at (")
          + std::to_string(i) + ", " + std::to_string(i) + ")");
    }
    // Issue here, sqrt is done over T of the complex schur
    Eigen::RealSchur<Eigen::MatrixXd> schurOfA(local_block);
    // Compute Schur decomposition of arg
    const auto& t_mat = schurOfA.matrixT();
    const auto& u_mat = schurOfA.matrixU();
    // Check if diagonal of schur is not positive
    if ((t_mat.diagonal().array() < 0).any()) {
      throw std::domain_error(
          std::string("Error in block_matrix_sqrt: "
                      "values less than 0 detected in block diagonal's schur "
                      "decomposition starting at (")
          + std::to_string(i) + ", " + std::to_string(i) + ")");
    }
    try {
      // Compute square root of T
      Eigen::matrix_sqrt_quasi_triangular(t_mat, sqrt_t_mat);
      // Compute square root of arg
      local_block_sqrt = u_mat * sqrt_t_mat * u_mat.adjoint();
    } catch (const std::exception& e) {
      throw std::domain_error(
          "Error in block_matrix_sqrt: "
          "The matrix is not positive definite");
    }
    for (int k = 0; k < block_size; k++) {
      for (int j = 0; j < block_size; j++) {
        W_root.coeffRef(i * block_size + j, i * block_size + k)
            = local_block_sqrt(j, k);
      }
    }
  }
}

/**
 * @brief Performs a Cholesky decomposition on a block diagonal matrix.
 * @tparam WRootMat A type inheriting from `Eigen::EigenBase`.
 * @param W_root The output matrix to store the square root.
 * @param W The input block diagonal matrix.
 * @param block_size The size of each block in the block diagonal matrix.
 */
template <typename WRootMat>
inline void block_matrix_chol_L(WRootMat& W_root,
                                const Eigen::SparseMatrix<double>& W,
                                const Eigen::Index block_size) {
  int n_block = W.cols() / block_size;
  Eigen::MatrixXd local_block(block_size, block_size);
  Eigen::MatrixXd local_block_sqrt(block_size, block_size);
  Eigen::MatrixXd sqrt_t_mat = Eigen::MatrixXd::Zero(block_size, block_size);
  // No block operation available for sparse matrices, so we have to loop
  // See https://eigen.tuxfamily.org/dox/group__TutorialSparse.html#title7
  for (int i = 0; i < n_block; i++) {
    sqrt_t_mat.setZero();
    local_block
        = W.block(i * block_size, i * block_size, block_size, block_size);
    if (Eigen::isnan(local_block.array()).any()) {
      throw std::domain_error(
          std::string("Error in block_matrix_sqrt: "
                      "NaNs detected in block diagonal starting at (")
          + std::to_string(i) + ", " + std::to_string(i) + ")");
    }
    try {
      // Compute square root of T
      Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt(local_block);
      if (llt.info() != Eigen::Success) {
        throw std::runtime_error("Cholesky failed on block "
                                 + std::to_string(i));
      }
      const auto Lb = llt.matrixL();
      for (int k = 0; k < block_size; k++) {
        for (int j = k; j < block_size; j++) {
          W_root.coeffRef(i * block_size + j, i * block_size + k) = Lb(j, k);
        }
      }
    } catch (const std::exception& e) {
      // As a backup do the schur decomposition for this block diagonal
      local_block
          = W.block(i * block_size, i * block_size, block_size, block_size);
      // Issue here, sqrt is done over T of the complex schur
      Eigen::RealSchur<Eigen::MatrixXd> schurOfA(local_block);
      // Compute Schur decomposition of arg
      const auto& t_mat = schurOfA.matrixT();
      const auto& u_mat = schurOfA.matrixU();
      // Check if diagonal of schur is not positive
      if ((t_mat.diagonal().array() < 0).any()) {
        throw std::domain_error(
            std::string("Error in block_matrix_sqrt: "
                        "values less than 0 detected in block diagonal's schur "
                        "decomposition starting at (")
            + std::to_string(i) + ", " + std::to_string(i) + ")");
      }
      try {
        // Compute square root of T
        Eigen::matrix_sqrt_quasi_triangular(t_mat, sqrt_t_mat);
        // Compute square root of arg
        local_block_sqrt.noalias() = u_mat * sqrt_t_mat * u_mat.adjoint();
      } catch (const std::exception& e) {
        throw std::domain_error(
            "Error in block_matrix_sqrt: "
            "The matrix is not positive definite");
      }
      for (int k = 0; k < block_size; k++) {
        for (int j = 0; j < block_size; j++) {
          W_root.coeffRef(i * block_size + j, i * block_size + k)
              = local_block_sqrt(j, k);
        }
      }
      throw std::domain_error(
          "Error in block_matrix_sqrt: "
          "The matrix is not positive definite");
    }
  }
}

/**
 * Collect the adjoints from the input and add them to the output.
 * @tparam ZeroInput If true, the adjoints of the input will be set to zero
 * @tparam Output A tuple or type where all scalar types are `arithmetic` types
 * @tparam Input A tuple or type where all scalar types are `var` types
 * @param output The output to which the adjoints will be added
 * @param input The input from which the adjoints will be collected
 */
template <bool ZeroInput = false, typename Output, typename Input,
          require_t<is_all_arithmetic_scalar<Output>>* = nullptr,
          require_t<is_all_var_scalar<Input>>* = nullptr>
inline void collect_adjoints(Output& output, Input&& input) {
  return iter_tuple_nested(
      [](auto&& output_i, auto&& input_i) {
        using output_i_t = std::decay_t<decltype(output_i)>;
        if constexpr (is_std_vector_v<output_i_t>) {
          Eigen::Map<Eigen::Matrix<double, -1, 1>> output_map(output_i.data(),
                                                              output_i.size());
          Eigen::Map<Eigen::Matrix<var, -1, 1>> input_map(input_i.data(),
                                                          input_i.size());
          output_map.array() += input_map.adj().array();
          if constexpr (ZeroInput) {
            input_map.adj().setZero();
          }
        } else if constexpr (is_eigen_v<output_i_t>) {
          output_i.array() += input_i.adj().array();
          if constexpr (ZeroInput) {
            input_i.adj().setZero();
          }
        } else if constexpr (is_stan_scalar_v<output_i_t>) {
          output_i += input_i.adj();
          if constexpr (ZeroInput) {
            input_i.adj() = 0;
          }
        } else {
          static_assert(
              sizeof(std::decay_t<output_i_t>*) == 0,
              "INTERNAL ERROR:(laplace_marginal_lpdf) set_zero_adjoints was "
              "not able to deduce the actions needed for the given type. "
              "This is an internal error, please report it: "
              "https://github.com/stan-dev/math/issues");
        }
      },
      std::forward<Output>(output), std::forward<Input>(input));
}

/**
 * Throws an error if the parameter contains NaN or Inf values.
 * @tparam NameStr Type of the name string, e.g. `std::string` or `char*`.
 * @tparam ParamStr Type of the parameter string, e.g. `std::string` or `char*`.
 * @tparam Param Type of the parameter such as a vector, matrix, or scalar.
 * @param name_str Name of the function or context where the error occurred.
 * @param param_str Name of the parameter that contains NaN or Inf values.
 * @param param The parameter to check for NaN or Inf values.
 */
template <typename NameStr, typename ParamStr, typename Param>
inline STAN_COLD_PATH void throw_nan(NameStr&& name_str, ParamStr&& param_str,
                                     Param&& param) {
  std::string msg = std::string("Error in ") + name_str + ": "
                    + std::string(param_str) + " contains NaN values";
  if ((param.array().isNaN() || !param.array().isFinite()).all()) {
    msg += " for all values.";
    throw std::domain_error(msg);
  }
  msg += " at indices [";
  for (int i = 0; i < param.size(); ++i) {
    if (std::isnan(param(i)) || std::isinf(param(i))) {
      msg += std::to_string(i) + ", ";
    }
  }
  msg.pop_back();
  msg.pop_back();
  msg += "].";
  throw std::domain_error(msg);
}

/**
 * @brief  Curvature-aware Barzilai–Borwein (BB) step length with robust
 * safeguards.
 *
 * Given successive parameter displacements \f$s = x_k - x_{k-1}\f$ and
 * gradients \f$g_k\f$, \f$g_{k-1}\f$, this routine forms
 * \f$y = g_k - g_{k-1}\f$ and computes the two classical BB candidates
 *
 * \f{align*}{
 *   \alpha_{\text{BB1}} &= \frac{\langle s,s\rangle}{\langle s,y\rangle},\\
 *   \alpha_{\text{BB2}} &= \frac{\langle s,y\rangle}{\langle y,y\rangle},
 * \f}
 *
 * then chooses between them using the **spectral cosine**
 * \f$r = \cos^2\!\angle(s,y) = \dfrac{\langle s,y\rangle^2}
 *                                      {\langle s,s\rangle\,\langle
 * y,y\rangle}\in[0,1]\f$:
 *
 *  - if \f$r > 0.9\f$ (well-aligned curvature) and the previous line search
 *    did **≤ 1** backtrack, prefer the “long” step \f$\alpha_{\text{BB1}}\f$;
 *  - if \f$0.1 \le r \le 0.9\f$, take the neutral geometric mean
 *    \f$\sqrt{\alpha_{\text{BB1}}\alpha_{\text{BB2}}}\f$;
 *  - otherwise default to the “short” step \f$\alpha_{\text{BB2}}\f$.
 *
 * All candidates are clamped into \f$[\text{min\_alpha},\,\text{max\_alpha}]\f$
 * and must be finite and positive.
 * If the curvature scalars are ill-posed (non-finite or too small),
 * \f$\langle s,y\rangle \le \varepsilon\f$, or if `last_backtracks == 99`
 * (explicitly disabling BB for this iteration), the function falls back to a
 * **safe** step:
 * use `prev_step` when finite and positive, otherwise \f$1.0\f$, then clamp to
 * \f$[\text{min\_alpha},\,\text{max\_alpha}]\f$.
 *
 * @param s          Displacement between consecutive iterates
 *                   (\f$s = x_k - x_{k-1}\f$).
 * @param g_curr     Gradient at the current iterate \f$g_k\f$.
 * @param g_prev     Gradient at the previous iterate \f$g_{k-1}\f$.
 * @param prev_step  Previously accepted step length (used by the fallback).
 * @param last_backtracks
 *                   Number of backtracking contractions performed by the most
 *                   recent line search; set to 99 to force the safe fallback.
 * @param min_alpha  Lower bound for the returned step length.
 * @param max_alpha  Upper bound for the returned step length.
 *
 * @return A finite, positive BB-style step length \f$\alpha \in
 *         [\text{min\_alpha},\,\text{max\_alpha}]\f$ suitable for seeding a
 *         line search or as a spectral preconditioner scale.
 *
 * @note Uses \f$\varepsilon=10^{-16}\f$ to guard against division by very
 *       small curvature terms, and applies `std::abs` to BB ratios to avoid
 *       negative steps; descent is enforced by the line search.
 * @warning The vectors must have identical size. Non-finite inputs yield the
 *          safe fallback.
 */
inline double barzilai_borwein_step_size(const Eigen::VectorXd& s,
                                         const Eigen::VectorXd& g_curr,
                                         const Eigen::VectorXd& g_prev,
                                         double prev_step, int last_backtracks,
                                         double min_alpha, double max_alpha) {
  // Fallbacks
  auto safe_fallback = [&]() -> double {
    double a = std::clamp(
        prev_step > 0.0 && std::isfinite(prev_step) ? prev_step : 1.0,
        min_alpha, max_alpha);
    return a;
  };

  const Eigen::VectorXd y = g_curr - g_prev;
  const double sty = s.dot(y);
  const double sts = s.squaredNorm();
  const double yty = y.squaredNorm();

  // Basic validity checks
  constexpr double eps = 1e-16;
  if (!(std::isfinite(sty) && std::isfinite(sts) && std::isfinite(yty))
      || sts <= eps || yty <= eps || sty <= eps || last_backtracks == 99) {
    return safe_fallback();
  }

  // BB candidates
  double alpha_bb1 = std::clamp(std::abs(sts / sty), min_alpha, max_alpha);
  double alpha_bb2 = std::clamp(std::abs(sty / yty), min_alpha, max_alpha);

  // Safeguard candidates
  if (!std::isfinite(alpha_bb1) || !std::isfinite(alpha_bb2) || alpha_bb1 <= 0.0
      || alpha_bb2 <= 0.0) {
    return safe_fallback();
  }

  // Spectral cosine r = cos^2(angle(s, y)) in [0,1]
  const double r = (sty * sty) / (sts * yty);

  // Heuristic thresholds (robust defaults)
  constexpr double kLoose = 0.9;  // "nice" curvature
  constexpr double kTight = 0.1;  // "dodgy" curvature

  double alpha0 = alpha_bb2;  // default to short BB for robustness
  if (r > kLoose && last_backtracks <= 1) {
    // Spectrum looks friendly and line search was not harsh -> try long BB
    alpha0 = alpha_bb1;
  } else if (r >= kTight && r <= kLoose) {
    // Neither clearly friendly nor clearly dodgy -> neutral middle
    alpha0 = std::sqrt(alpha_bb1 * alpha_bb2);
  }  // else keep alpha_bb2

  // Clip to user bounds
  alpha0 = std::clamp(alpha0, min_alpha, max_alpha);

  if (!std::isfinite(alpha0) || alpha0 <= 0.0) {
    return safe_fallback();
  }
  return alpha0;
}

/**
 * For a latent Gaussian model with hyperparameters phi and
 * latent variables theta, and observations y, this function computes
 * an approximation of the log marginal density, p(y | phi).
 * This is done by marginalizing out theta, using a Laplace
 * approxmation. The latter is obtained by finding the mode,
 * via Newton's method, and computing the Hessian of the likelihood.
 *
 * The convergence criterion for the Newton is a small change in
 * log marginal density. The user controls the tolerance (i.e.
 * threshold under which change is deemed small enough) and
 * maximum number of steps.
 *
 * A description of this algorithm can be found in:
 *  - (2023) Margossian, "General Adjoint-Differentiated Laplace approximation",
 *    https://arxiv.org/pdf/2306.14976.
 * Additional references include:
 *  - (2020) Margossian et al, "HMC using an adjoint-differentiated Laplace...",
 *    NeurIPS, https://arxiv.org/abs/2004.12550.
 *  - (2006) Rasmussen and Williams, "Gaussian Processes for Machine Learning",
 *    second edition, MIT Press, algorithm 3.1.
 *
 * Variables needed for the gradient or generating quantities
 * are stored by reference.
 *
 * @tparam LLFun Type with a valid `operator(ThetaVec,  InnerLLTupleArgs)`
 * where `InnerLLTupleArgs` are the elements of `LLTupleArgs`
 * @tparam LLTupleArgs A tuple whose elements follow the types required for
 * `LLFun`
 * \laplace_common_template_args
 * @param[in] ll_fun A log likelihood functor
 * @param[in] ll_args Tuple containing parameters for `LLFun`
 * @param[in] covariance The covariance matrix for the latent Gaussian
 * \laplace_common_args
 * @param[in] options A set of options for tuning the solver
 * \msg_arg
 *
 * @return A struct containing
 * 1. lmd the log marginal density, p(y | phi)
 * 2. covariance the evaluated covariance function for the latent gaussian
 * variable
 * 3. theta a vector to store the mode
 * 4. W_r A sparse matrix containing the square root of the negative
 *    hessian, if solver 1 or 2 are used.
 * 5. L cholesky decomposition of stabilized inverse covariance
 * 6. a element in the Newton step
 * 7. l_grad the log density of the likelihood, evaluated at the mode
 *
 */
template <typename LLFun, typename LLTupleArgs, typename CovarMat,
          bool InitTheta,
          require_t<is_all_arithmetic_scalar<CovarMat, LLTupleArgs>>* = nullptr>
inline auto laplace_marginal_density_est(
    LLFun&& ll_fun, LLTupleArgs&& ll_args, CovarMat&& covariance,
    const laplace_options<InitTheta>& options, std::ostream* msgs) {
  if constexpr (InitTheta) {
    check_nonzero_size("laplace_marginal", "initial guess", options.theta_0);
    check_finite("laplace_marginal", "initial guess", options.theta_0);
  }
  check_nonnegative("laplace_marginal", "tolerance", options.tolerance);
  check_positive("laplace_marginal", "max_num_steps", options.max_num_steps);
  check_positive("laplace_marginal", "hessian_block_size",
                 options.hessian_block_size);
  check_square("laplace_marginal", "covariance", covariance);

  const Eigen::Index theta_size = covariance.rows();

  if (unlikely(theta_size % options.hessian_block_size != 0 ||
    theta_size < options.hessian_block_size)) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << "laplace_marginal_density: The hessian size (" << theta_size
          << ", " << theta_size << ")";
      if (theta_size % options.hessian_block_size != 0) {
        msg << " is not divisible by the hessian block size (";
      } else {
        msg << " is smaller than the hessian block size (";
      }
      msg << options.hessian_block_size
          << ")"
             ". Use a hessian block size such as [1, ";
      for (int i = 2; i < 12; ++i) {
        if (theta_size % i == 0) {
          msg << i << ", ";
        }
      }
      msg << "... " << theta_size;
      msg << "].";
      throw std::domain_error(msg.str());
    }();
  }

  auto throw_overstep = [](const auto max_num_steps) STAN_COLD_PATH {
    throw std::domain_error(
        std::string("laplace_marginal_density: max number of iterations: ")
        + std::to_string(max_num_steps) + " exceeded.");
  };
  // Wolfe optimizes over the latent 'a' space
  auto obj_fun = [&](const Eigen::VectorXd& a_val, auto&& theta_val) -> double {
    return -0.5 * a_val.dot(theta_val)
           + laplace_likelihood::log_likelihood(ll_fun, theta_val, ll_args,
                                                msgs);
  };
  auto theta_grad_f = [&ll_fun, &ll_args, &msgs](auto&& theta_val) {
    return laplace_likelihood::theta_grad(ll_fun, theta_val, ll_args, msgs);
  };
  internal::WolfeInfo wolfe_info(
      obj_fun, theta_size,
      [theta_size, &options]() -> decltype(auto) {
        if constexpr (InitTheta) {
          return options.theta_0;
        } else {
          return Eigen::VectorXd::Zero(theta_size);
        }
      }(),
      theta_grad_f);
  auto&& curr = wolfe_info.curr_;
  auto&& prev = wolfe_info.prev_;
  Eigen::MatrixXd B(theta_size, theta_size);
  Eigen::VectorXd b(theta_size);
  // 'a' gradient
  auto grad_fun = [&covariance](auto&& step) {
    return -covariance * step.a() + covariance * step.theta_grad();
  };
  auto update_step = [&covariance, &obj_fun, &theta_grad_f, &grad_fun](
                         auto& step_info, auto&& /* curr */, auto&& prev,
                         auto& eval_in, auto&& p) {
    step_info.a() = prev.a() + eval_in.alpha() * p;
    step_info.theta().noalias() = covariance * step_info.a();
    step_info.theta_grad() = theta_grad_f(step_info.theta());
    eval_in.obj() = obj_fun(step_info.a(), step_info.theta());
    eval_in.dir() = grad_fun(step_info).dot(p);
  };
  Eigen::VectorXd prev_g(theta_size);
  auto update_line_search = [&grad_fun, &covariance, &prev_g, &update_step,
                             &theta_grad_f, &options,
                             &msgs](auto&& wolfe_status, auto&& wolfe_info,
                                    auto&& curr, auto&& prev) {
    wolfe_info.p_ = curr.a() - prev.a();
    prev_g.noalias() = grad_fun(prev);
    wolfe_info.init_dir_ = prev_g.dot(wolfe_info.p_);
    // Flip direction if not ascending
    if (wolfe_info.init_dir_ < 0) {
      wolfe_info.p_ = -wolfe_info.p_;
      wolfe_info.init_dir_ = -wolfe_info.init_dir_;
    }
    auto scratch = wolfe_info.scratch_;
    scratch.alpha() = 1.0;
    while (scratch.alpha() > options.line_search.min_alpha) {
      try {
        update_step(scratch, curr, prev, scratch.eval_, wolfe_info.p_);
        if (std::isnan(scratch.eval_.obj()) || std::isinf(scratch.eval_.obj())
            || std::isnan(scratch.eval_.dir())
            || std::isinf(scratch.eval_.dir())) {
          scratch.alpha() *= options.line_search.tau;
          continue;
        }
      } catch (const std::exception& e) {
        scratch.alpha() *= options.line_search.tau;
        continue;
      }
      break;
    }
    if (scratch.alpha() <= options.line_search.min_alpha) {
      wolfe_status.success_ = false;
      return true;
    }
    if (options.line_search.max_iterations == 0) {
      if (scratch.alpha() > options.line_search.min_alpha) {
        curr.update(scratch);
        wolfe_status.success_ = true;
        return false;
      }
    } else {
      curr.alpha() = barzilai_borwein_step_size(
          wolfe_info.p_, grad_fun(scratch), prev_g, prev.alpha(),
          wolfe_status.num_backtracks_, options.line_search.min_alpha,
          options.line_search.max_alpha);
      wolfe_status = internal::wolfe_line_search(wolfe_info, update_step,
                                                 options.line_search, msgs);
    }
    return abs(curr.obj() - prev.obj()) < options.tolerance
           || (!wolfe_status.success_ && curr.obj() <= prev.obj());
  };
  auto set_next_iter = [&options](auto&& curr, auto&& prev) {
    prev.update(curr);
    curr.alpha() = std::clamp(curr.alpha(), 0.0, options.line_search.max_alpha);
  };
  /**
   * On the final loop if we found a better wolfe step, but we are going to
   * exit, we want to make sure all of our return values are with the most
   * recent wolfe step that was accepted. So we do one final loop to update
   * our return values.
   */
  bool allow_bounce = false;
  bool final_loop = false;
  bool finish_update = false;
  WolfeStatus wolfe_status;
  // Start with safe step size
  wolfe_status.num_backtracks_ = 99;
  Eigen::Index step_iter = 0;
  try {
    if (options.solver == 1) {
      if (options.hessian_block_size == 1) {
        //   std::cout << "Solver: 1Diag" << std::endl;
        Eigen::VectorXd W_r(theta_size);
        for (; step_iter <= options.max_num_steps; step_iter++) {
          auto W = laplace_likelihood::diagonal_hessian(ll_fun, prev.theta(),
                                                        ll_args, msgs);
          for (Eigen::Index j = 0; j < W.size(); j++) {
            if (W.coeff(j) < 0) {
              throw std::domain_error(
                  "laplace_marginal_density: Hessian matrix is not positive "
                  "definite");
            } else {
              W_r.coeffRef(j) = std::sqrt(W.coeff(j));
            }
          }
          B.noalias() = Eigen::MatrixXd::Identity(theta_size, theta_size)
                        + W_r.asDiagonal() * covariance * W_r.asDiagonal();
          Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt_B(B);
          if (llt_B.info() != Eigen::Success) {
            double jitter_try = 1e-10;
            for (; jitter_try < 1e-5; jitter_try *= 10) {
              B.diagonal().array() += jitter_try;
              llt_B.compute(B);
              if (llt_B.info() == Eigen::Success) {
                break;
              }
            }
            if (llt_B.info() != Eigen::Success) {
              throw std::domain_error(
                  "laplace_marginal_density: Cholesky failed in iteration "
                  + std::to_string(step_iter));
            }
          }
          auto L = llt_B.matrixL();
          auto LT = llt_B.matrixU();
          b.noalias()
              = (W.array() * prev.theta().array()).matrix() + prev.theta_grad();
          curr.a().noalias()
              = b
                - W_r.asDiagonal()
                      * LT.solve(L.solve(W_r.cwiseProduct(covariance * b)));
          if (!final_loop) {
            finish_update
                = update_line_search(wolfe_status, wolfe_info, curr, prev);
          }
          if (finish_update) {
            if (!final_loop && wolfe_status.success_) {
              // Do one final loop with exact wolfe conditions
              final_loop = true;
              // NOTE: Swapping here so we need to swap prev and curr later
              set_next_iter(curr, prev);
              continue;
            }
            const double B_log_determinant
                = 2.0 * llt_B.matrixLLT().diagonal().array().log().sum();
            /**
             * NOTE: At this point we are return prev because either
             * 1. Line search + tolerance passed and we swapped prev<->curr
             * 2. Line search failed and we want to return the previous step
             */
            return laplace_density_estimates{
                prev.obj() - 0.5 * B_log_determinant,
                std::move(prev.theta()),
                Eigen::SparseMatrix<double>(W_r.asDiagonal()),
                Eigen::MatrixXd(L),
                std::move(prev.a()),
                std::move(prev.theta_grad()),
                Eigen::PartialPivLU<Eigen::MatrixXd>{},
                Eigen::MatrixXd(0, 0),
                1};
          } else {
            set_next_iter(curr, prev);
          }
        }
      } else {
        Eigen::SparseMatrix<double> W_r(theta_size, theta_size);
        Eigen::Index block_size = options.hessian_block_size;
        W_r.reserve(Eigen::VectorXi::Constant(W_r.cols(), block_size));
        const Eigen::Index n_block = W_r.cols() / block_size;
        // Prefill W_r so we only make space once
        for (Eigen::Index ii = 0; ii < n_block; ii++) {
          for (Eigen::Index k = 0; k < block_size; k++) {
            for (Eigen::Index j = 0; j < block_size; j++) {
              W_r.insert(ii * block_size + j, ii * block_size + k) = 1.0;
            }
          }
        }
        W_r.makeCompressed();
        for (; step_iter <= options.max_num_steps; step_iter++) {
          auto W = laplace_likelihood::block_hessian(
              ll_fun, prev.theta(), options.hessian_block_size, ll_args, msgs);
          for (Eigen::Index j = 0; j < W.rows(); j++) {
            if (W.coeff(j, j) < 0) {
              throw std::domain_error(
                  "laplace_marginal_density: Hessian matrix is not positive "
                  "definite");
            }
          }
          block_matrix_sqrt(W_r, W, options.hessian_block_size);
          B.noalias() = Eigen::MatrixXd::Identity(theta_size, theta_size)
                        + W_r * (covariance * W_r);
          Eigen::LLT<Eigen::MatrixXd> llt_B(B);
          if (llt_B.info() != Eigen::Success) {
            double jitter_try = 1e-10;
            for (; jitter_try < 1e-5; jitter_try *= 10) {
              B.diagonal().array() += jitter_try;
              llt_B.compute(B);
              if (llt_B.info() == Eigen::Success) {
                break;
              }
            }
            if (llt_B.info() != Eigen::Success) {
              throw std::domain_error(
                  "laplace_marginal_density: Cholesky failed in iteration "
                  + std::to_string(step_iter));
            }
          }
          auto L = llt_B.matrixL();
          auto LT = llt_B.matrixU();
          b.noalias() = W * prev.theta() + prev.theta_grad();
          curr.a().noalias()
              = b - W_r * LT.solve(L.solve(W_r * (covariance * b)));
          if (!final_loop) {
            finish_update
                = update_line_search(wolfe_status, wolfe_info, curr, prev);
          }
          if (finish_update) {
            if (!final_loop && wolfe_status.success_) {
              // Do one final loop with exact wolfe conditions
              final_loop = true;
              set_next_iter(curr, prev);
              continue;
            }
            const double B_log_determinant
                = 2.0 * llt_B.matrixLLT().diagonal().array().log().sum();
            return laplace_density_estimates{
                prev.obj() - 0.5 * B_log_determinant,
                std::move(prev.theta()),
                std::move(W_r),
                Eigen::MatrixXd(L),
                std::move(prev.a()),
                std::move(prev.theta_grad()),
                Eigen::PartialPivLU<Eigen::MatrixXd>{},
                Eigen::MatrixXd(0, 0),
                1};
          } else {
            set_next_iter(curr, prev);
          }
        }
      }
      throw_overstep(options.max_num_steps);
    }
  } catch (const std::exception& e) {
    allow_bounce = true;
    if (msgs != nullptr) {
      (*msgs) << "Solver 1 failed at iteration " << step_iter
              << " with error: " << e.what() << std::endl;
      (*msgs) << "Attempting to switch to solver 2 (LLT decomposition)."
              << std::endl;
    }
  }
  try {
    if (options.solver == 2 || allow_bounce) {
      Eigen::MatrixXd K_root
          = covariance.template selfadjointView<Eigen::Lower>().llt().matrixL();
      for (; step_iter <= options.max_num_steps; step_iter++) {
        auto W = laplace_likelihood::block_hessian(
            ll_fun, prev.theta(), options.hessian_block_size, ll_args, msgs);
        B.noalias() = Eigen::MatrixXd::Identity(theta_size, theta_size)
                      + K_root.transpose() * W * K_root;
        Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt_B(B);
        if (llt_B.info() != Eigen::Success) {
          double jitter_try = 1e-10;
          for (; jitter_try < 1e-5; jitter_try *= 10) {
            B.diagonal().array() += jitter_try;
            llt_B.compute(B);
            if (llt_B.info() == Eigen::Success) {
              break;
            }
          }
          if (llt_B.info() != Eigen::Success) {
            throw std::domain_error(
                "laplace_marginal_density: Cholesky failed in iteration "
                + std::to_string(step_iter));
          }
        }
        auto L = llt_B.matrixL();
        auto LT = llt_B.matrixU();
        b.noalias() = W * prev.theta() + prev.theta_grad();
        curr.a().noalias()
            = K_root.transpose().template triangularView<Eigen::Upper>().solve(
                LT.solve(L.solve(K_root.transpose() * b)));
        if (!final_loop) {
          finish_update
              = update_line_search(wolfe_status, wolfe_info, curr, prev);
        }
        if (finish_update) {
          if (!final_loop && wolfe_status.success_) {
            // Do one final loop with exact wolfe conditions
            final_loop = true;
            // NOTE: Swapping here so we need to swap prev and curr later
            set_next_iter(curr, prev);
            continue;
          }
          const double B_log_determinant
              = 2.0 * llt_B.matrixLLT().diagonal().array().log().sum();
          return laplace_density_estimates{
              prev.obj() - 0.5 * B_log_determinant,
              std::move(prev.theta()),
              std::move(W),
              std::move(Eigen::MatrixXd(L)),
              std::move(prev.a()),
              std::move(prev.theta_grad()),
              Eigen::PartialPivLU<Eigen::MatrixXd>{},
              std::move(K_root),
              2};
        } else {
          set_next_iter(curr, prev);
        }
      }
      throw_overstep(options.max_num_steps);
    }
  } catch (const std::exception& e) {
    allow_bounce = true;
    if (msgs != nullptr) {
      (*msgs) << "Solver 2 failed at iteration " << step_iter
              << " with error: " << e.what() << std::endl;
      (*msgs) << "Attempting to switch to solver 3 (LU decomposition)."
              << std::endl;
    }
  }
  if (options.solver == 3 || allow_bounce) {
    //    std::cout << "Solver: 3" << std::endl;
    Eigen::PartialPivLU<Eigen::MatrixXd> LU(theta_size);
    for (; step_iter <= options.max_num_steps; step_iter++) {
      auto W = laplace_likelihood::block_hessian(
          ll_fun, prev.theta(), options.hessian_block_size, ll_args, msgs);
      LU.compute(Eigen::MatrixXd::Identity(theta_size, theta_size)
                 + covariance * W);
      // L on lower and U on upper triangular
      b.noalias() = W * prev.theta() + prev.theta_grad();
      curr.a().noalias() = b - W * LU.solve(covariance * b);
      if (!final_loop) {
        finish_update
            = update_line_search(wolfe_status, wolfe_info, curr, prev);
      }
      if (finish_update) {
        if (!final_loop && wolfe_status.success_) {
          // Do one final loop with exact wolfe conditions
          final_loop = true;
          // NOTE: Swapping here so we need to swap prev and curr later
          set_next_iter(curr, prev);
          continue;
        }
        const double B_log_determinant
            = LU.matrixLU().diagonal().array().log().sum();
        return laplace_density_estimates{prev.obj() - 0.5 * B_log_determinant,
                                         std::move(prev.theta()),
                                         std::move(W),
                                         Eigen::MatrixXd(0, 0),
                                         std::move(prev.a()),
                                         std::move(prev.theta_grad()),
                                         std::move(LU),
                                         Eigen::MatrixXd(0, 0),
                                         3};
      } else {
        set_next_iter(curr, prev);
      }
    }
    throw_overstep(options.max_num_steps);
  }
  throw std::domain_error(
      std::string("You chose a solver (") + std::to_string(options.solver)
      + ") that is not valid. Please choose either 1, 2, or 3.");
}
}  // namespace internal
/**
 * For a latent Gaussian model with global parameters phi, latent
 * variables theta, and observations y, this function computes
 * an approximation of the log marginal density, p(y | phi).
 * This is done by marginalizing out theta, using a Laplace
 * approxmation. The latter is obtained by finding the mode,
 * using a custom Newton method, and the Hessian of the likelihood.
 *
 * The convergence criterion for the Newton is a small change in
 * log marginal density. The user controls the tolerance (i.e.
 * threshold under which change is deemed small enough) and
 * maximum number of steps.
 *
 * Wrapper for when the hyperparameters are passed as a double.
 *
 * @tparam LLFun Type with a valid `operator(ThetaVec, InnerLLTupleArgs)`
 * where `InnerLLTupleArgs` are the elements of `LLTupleArgs`
 * @tparam LLTupleArgs A tuple whose elements follow the types required for
 * `LLFun`
 * \laplace_common_template_args
 * @param[in] ll_fun A log likelihood functor
 * @param[in] ll_args Tuple containing parameters for `LLFun`
 * \laplace_common_args
 * @param[in] options A set of options for tuning the solver
 * \msg_arg
 * @return the log maginal density, p(y | phi)
 */
template <
    typename LLFun, typename LLTupleArgs, typename CovarFun, typename CovarArgs,
    bool InitTheta,
    require_t<is_all_arithmetic_scalar<CovarArgs, LLTupleArgs>>* = nullptr>
inline double laplace_marginal_density(
    LLFun&& ll_fun, LLTupleArgs&& ll_args, CovarFun&& covariance_function,
    CovarArgs&& covar_args, const laplace_options<InitTheta>& options,
    std::ostream* msgs) {
  Eigen::MatrixXd covariance = stan::math::apply(
      [msgs, &covariance_function](auto&&... args) {
        return covariance_function(std::forward<decltype(args)>(args)..., msgs);
      },
      std::forward<CovarArgs>(covar_args));
  return internal::laplace_marginal_density_est(
             std::forward<LLFun>(ll_fun), std::forward<LLTupleArgs>(ll_args),
             std::move(covariance), options, msgs)
      .lmd;
}

namespace internal {

/**
 * Collects the adjoints from the input and adds them to the output.
 * @tparam Output A tuple or type where all scalar types are `arithmetic` types
 * @tparam Input A tuple or type where all scalar types are `arithmetic` types
 * @param output The output to which the adjoints will be added
 * @param input The input from which the adjoints will be collected
 */
template <typename Output, typename Input,
          require_t<is_all_arithmetic_scalar<Output>>* = nullptr,
          require_t<is_all_arithmetic_scalar<Input>>* = nullptr>
inline void collect_adjoints(Output&& output, Input&& input) {
  return iter_tuple_nested(
      [](auto&& output_i, auto&& input_i) {
        using output_i_t = std::decay_t<decltype(output_i)>;
        if constexpr (is_std_vector_v<output_i_t>) {
          Eigen::Map<Eigen::Matrix<double, -1, 1>> output_map(output_i.data(),
                                                              output_i.size());
          Eigen::Map<Eigen::Matrix<double, -1, 1>> input_map(input_i.data(),
                                                             input_i.size());
          output_map.array() += input_map.array();
        } else if constexpr (is_eigen_v<output_i_t>) {
          output_i.array() += input_i.array();
        } else if constexpr (is_stan_scalar_v<output_i_t>) {
          output_i += input_i;
        } else {
          static_assert(
              sizeof(std::decay_t<output_i_t>*) == 0,
              "INTERNAL ERROR:(laplace_marginal_lpdf) set_zero_adjoints was "
              "not able to deduce the actions needed for the given type. "
              "This is an internal error, please report it: "
              "https://github.com/stan-dev/math/issues");
        }
      },
      std::forward<Output>(output), std::forward<Input>(input));
}
/**
 * Base case for zero sized tuples
 */
template <bool ZeroInput = false>
inline void constexpr copy_compute_s2(const std::tuple<>& output,
                                      const std::tuple<>& input) noexcept {}

/**
 * Copies the adjoints from the input to the output, scaling them by 0.5.
 * @tparam ZeroInput If true, the adjoints of the input will be set to zero
 * @tparam Output A tuple or type where all scalar types are `arithmetic` types
 * @tparam Input A tuple or type where all scalar types are `var` types
 * @param output The output to which the adjoints will be added
 * @param input The input from which the adjoints will be collected
 */
template <bool ZeroInput = false, typename Output, typename Input,
          require_t<is_all_arithmetic_scalar<Output>>* = nullptr,
          require_t<is_any_var_scalar<Input>>* = nullptr>
inline void copy_compute_s2(Output&& output, Input&& input) {
  return iter_tuple_nested(
      [](auto&& output_i, auto&& input_i) {
        using output_i_t = std::decay_t<decltype(output_i)>;
        if constexpr (is_std_vector_v<output_i_t>) {
          Eigen::Map<Eigen::Matrix<double, -1, 1>> output_map(output_i.data(),
                                                              output_i.size());
          Eigen::Map<Eigen::Matrix<var, -1, 1>> input_map(input_i.data(),
                                                          input_i.size());
          output_map.array() += 0.5 * input_map.adj().array();
          if constexpr (ZeroInput) {
            input_map.adj().setZero();
          }
        } else if constexpr (is_eigen_v<output_i_t>) {
          output_i.array() += 0.5 * input_i.adj().array();
          if constexpr (ZeroInput) {
            input_i.adj().setZero();
          }
        } else if constexpr (is_stan_scalar_v<output_i_t>) {
          output_i += (0.5 * input_i.adj());
          if constexpr (ZeroInput) {
            input_i.adj() = 0;
          }
        } else {
          static_assert(
              sizeof(std::decay_t<output_i_t>*) == 0,
              "INTERNAL ERROR:(laplace_marginal_lpdf) set_zero_adjoints was "
              "not able to deduce the actions needed for the given type. "
              "This is an internal error, please report it: "
              "https://github.com/stan-dev/math/issues");
        }
      },
      std::forward<Output>(output), std::forward<Input>(input));
}

template <typename T>
inline constexpr decltype(auto) filter_var_scalar_types(T&& t) {
  return stan::math::filter_map<is_any_var_scalar>(
      [](auto&& arg) -> decltype(auto) {
        return std::forward<decltype(arg)>(arg);
      },
      std::forward<T>(t));
}
/**
 * Creates an arena type from the input with initialized with zeros
 * @tparam Input Possibly a tuple, std::vector, Eigen type, or scalar
 * @param input The input to be converted to an arena type
 */
template <typename Input>
inline constexpr auto make_zeroed_arena(Input&& input) {
  if constexpr (is_tuple_v<Input>) {
    return stan::math::filter_map<is_any_var_scalar>(
        [](auto&& output_i) { return make_zeroed_arena(output_i); }, input);
  } else if constexpr (is_std_vector_v<Input>) {
    if constexpr (!is_var_v<value_type_t<Input>>) {
      const auto output_size = input.size();
      arena_t<std::vector<decltype(make_zeroed_arena(input[0]))>> ret;
      ret.reserve(output_size);
      for (Eigen::Index i = 0; i < output_size; ++i) {
        ret.push_back(make_zeroed_arena(input[i]));
      }
      return ret;
    } else {
      return arena_t<std::vector<double>>(input.size(), 0.0);
    }
  } else if constexpr (is_eigen_v<Input>) {
    return arena_t<promote_scalar_t<double, Input>>(
        plain_type_t<promote_scalar_t<double, Input>>::Zero(input.rows(),
                                                            input.cols()));
  } else if constexpr (is_var<Input>::value) {
    return static_cast<double>(0.0);
  }
}

/**
 * Used in reverse pass to collect adjoints to the output
 * @tparam Output A tuple or type where all scalar types are `var` types
 * @tparam Input A tuple or type where all scalar types are `arithmetic` types
 * @param output The output to which the adjoints will be added
 * @param ret The vari object containing the adjoint to be added
 * @param input The input from which the adjoints will be collected
 */
template <typename Output, typename Input>
inline void collect_adjoints(Output&& output, const vari* ret, Input&& input) {
  if constexpr (is_tuple_v<Output>) {
    static_assert(sizeof(std::decay_t<Output>*) == 0,
                  "INTERNAL ERROR:(laplace_marginal_lpdf) "
                  "Accumulate Adjoints called on a tuple, but tuples cannot be "
                  "on the reverse mode stack! "
                  "This is an internal error, please report it: "
                  "https://github.com/stan-dev/math/issues");
  } else if constexpr (is_std_vector_v<Output>) {
    if constexpr (!is_var_v<value_type_t<Output>>) {
      const auto output_size = output.size();
      for (std::size_t i = 0; i < output_size; ++i) {
        collect_adjoints(output[i], ret, input[i]);
      }
    } else {
      Eigen::Map<Eigen::Matrix<var, -1, 1>> output_map(output.data(),
                                                       output.size());
      Eigen::Map<const Eigen::Matrix<double, -1, 1>> input_map(input.data(),
                                                               input.size());
      output_map.array().adj() += ret->adj_ * input_map.array();
    }
  } else if constexpr (is_eigen_v<Output>) {
    output.adj().array() += ret->adj_ * input.array();
  } else if constexpr (is_var_v<Output>) {
    output.adj() += ret->adj_ * input;
  }
}

/**
 * Collects adjoints from a tuple or std::vector of tuples
 * @tparam Output A tuple or std::vector of tuples where all scalar types are
 * `var` types
 * @tparam Input A tuple or std::vector of tuples where all scalar types are
 * `arithmetic` types
 * @param ret The vari object containing the adjoint to be added
 * @param output The output to which the adjoints will be added
 * @param input The input from which the adjoints will be collected
 */
template <typename Output, typename Input>
inline void reverse_pass_collect_adjoints(var ret, Output&& output,
                                          Input&& input) {
  if constexpr (is_tuple_v<Output>) {
    stan::math::for_each(
        [ret](auto&& inner_arg, auto&& inner_input) mutable {
          reverse_pass_collect_adjoints(
              ret, std::forward<decltype(inner_arg)>(inner_arg),
              std::forward<decltype(inner_input)>(inner_input));
        },
        std::forward<Output>(output), std::forward<Input>(input));
  } else if constexpr (is_std_vector_containing_tuple_v<Output>) {
    for (std::size_t i = 0; i < output.size(); ++i) {
      reverse_pass_collect_adjoints(ret, output[i], input[i]);
    }
  } else {
    reverse_pass_callback(
        [vi = ret.vi_, arg_arena = to_arena(std::forward<Output>(output)),
         input_arena = to_arena(std::forward<Input>(input))]() mutable {
          collect_adjoints(arg_arena, vi, input_arena);
        });
  }
}
}  // namespace internal
/**
 * For a latent Gaussian model with global parameters phi, latent
 * variables theta, and observations y, this function computes
 * an approximation of the log marginal density, p(y | phi).
 * This is done by marginalizing out theta, using a Laplace
 * approxmation. The latter is obtained by finding the mode,
 * using a custom Newton method, and the Hessian of the likelihood.
 *
 * The convergence criterion for the Newton is a small change in
 * the log marginal density. The user controls the tolerance (i.e.
 * threshold under which change is deemed small enough) and
 * maximum number of steps.
 *
 * Wrapper for when the global parameter is passed as a double.
 *
 * @tparam LLFun Type with a valid `operator(ThetaVec,  InnerLLTupleArgs)`
 * where `InnerLLTupleArgs` are the elements of `LLTupleArgs`
 * @tparam LLTupleArgs A tuple whose elements follow the types required for
 * `LLFun`
 * \laplace_common_template_args
 * @param[in] ll_fun A log likelihood functor
 * @param[in] ll_args Tuple containing parameters for `LLFun`
 * \laplace_common_args
 * @param[in] options A set of options for tuning the solver
 * \msg_arg
 * @return the log maginal density, p(y | phi)
 */
template <typename LLFun, typename LLTupleArgs, typename CovarFun,
          typename CovarArgs, bool InitTheta,
          require_t<is_any_var_scalar<LLTupleArgs, CovarArgs>>* = nullptr>
inline auto laplace_marginal_density(const LLFun& ll_fun, LLTupleArgs&& ll_args,
                                     CovarFun&& covariance_function,
                                     CovarArgs&& covar_args,
                                     const laplace_options<InitTheta>& options,
                                     std::ostream* msgs) {
  auto covar_args_refs = to_ref(std::forward<CovarArgs>(covar_args));
  auto ll_args_refs = to_ref(std::forward<LLTupleArgs>(ll_args));
  // Solver 1, 2, 3
  constexpr bool ll_args_contain_var = is_any_var_scalar<LLTupleArgs>::value;
  auto partial_parm = internal::make_zeroed_arena(ll_args_refs);
  auto covar_args_adj = internal::make_zeroed_arena(covar_args_refs);
  double lmd = 0.0;
  {
    nested_rev_autodiff nested;

    // Make one hard copy here
    using laplace_likelihood::internal::COPY_TYPE;
    using laplace_likelihood::internal::deep_copy_vargs;
    auto ll_args_copy = deep_copy_vargs<var>(ll_args_refs);
    auto covar_args_copy = deep_copy_vargs<var>(covar_args_refs);
    auto covariance = stan::math::apply(
        [&covariance_function, &msgs](auto&&... args) {
          if constexpr (is_any_var_scalar_v<decltype(args)...>) {
            return to_var_value(covariance_function(args..., msgs));
          } else {
            return covariance_function(args..., msgs);
          }
        },
        covar_args_copy);
    auto md_est = internal::laplace_marginal_density_est(
        ll_fun, value_of(ll_args_copy), value_of(covariance), options, msgs);
    if constexpr (ll_args_contain_var) {
      laplace_likelihood::ll_arg_grad(ll_fun, md_est.theta, ll_args_copy, msgs);
    }
    // Solver 1, 2
    arena_t<Eigen::MatrixXd> R(md_est.theta.size(), md_est.theta.size());
    // Solver 3
    arena_t<Eigen::MatrixXd> LU_solve_covariance(
        covariance.rows() * (md_est.solver_used == 3),
        covariance.cols() * (md_est.solver_used == 3));
    // Solver 1, 2, 3
    arena_t<Eigen::VectorXd> s2(md_est.theta.size());

    // Return references to var types
    auto ll_args_filter = internal::filter_var_scalar_types(ll_args_copy);
    stan::math::for_each(
        [](auto&& output_i, auto&& ll_arg_i) {
          if constexpr (is_any_var_scalar_v<decltype(ll_arg_i)>) {
            internal::collect_adjoints<true>(output_i, ll_arg_i);
          }
        },
        partial_parm, ll_args_filter);
    if (md_est.solver_used == 1) {
      if (options.hessian_block_size == 1) {
        arena_t<Eigen::MatrixXd> tmp = md_est.W_r.toDense();
        md_est.L.template triangularView<Eigen::Lower>().solveInPlace(tmp);
        R.noalias() = tmp.transpose() * tmp;
        arena_t<Eigen::MatrixXd> C
            = md_est.L.template triangularView<Eigen::Lower>().solve(
                md_est.W_r * value_of(covariance));
        if constexpr (!ll_args_contain_var) {
          s2.deep_copy(
              (0.5
               * (value_of(covariance).diagonal()
                  - (C.transpose() * C).diagonal())
                     .cwiseProduct(laplace_likelihood::third_diff(
                         ll_fun, md_est.theta, value_of(ll_args_copy), msgs))));
        } else {
          arena_t<Eigen::MatrixXd> A = value_of(covariance) - C.transpose() * C;
          auto s2_tmp = laplace_likelihood::compute_s2(
              ll_fun, md_est.theta, A, options.hessian_block_size, ll_args_copy,
              msgs);
          s2.deep_copy(s2_tmp);
          internal::copy_compute_s2<true>(partial_parm, ll_args_filter);
        }

      } else {
        arena_t<Eigen::MatrixXd> tmp = md_est.W_r.toDense();
        md_est.L.template triangularView<Eigen::Lower>().solveInPlace(tmp);
        R.noalias() = tmp.transpose() * tmp;
        arena_t<Eigen::MatrixXd> C
            = md_est.L.template triangularView<Eigen::Lower>().solve(
                md_est.W_r * value_of(covariance));
        arena_t<Eigen::MatrixXd> A = value_of(covariance) - C.transpose() * C;
        auto s2_tmp = laplace_likelihood::compute_s2(ll_fun, md_est.theta, A,
                                                     options.hessian_block_size,
                                                     ll_args_copy, msgs);
        s2.deep_copy(s2_tmp);
        internal::copy_compute_s2<true>(partial_parm, ll_args_filter);
      }
    } else if (md_est.solver_used == 2) {
      R = md_est.W_r
          - md_est.W_r * md_est.K_root
                * md_est.L.transpose()
                      .template triangularView<Eigen::Upper>()
                      .solve(
                          md_est.L.template triangularView<Eigen::Lower>()
                              .solve(md_est.K_root.transpose() * md_est.W_r));

      arena_t<Eigen::MatrixXd> C
          = md_est.L.template triangularView<Eigen::Lower>().solve(
              md_est.K_root.transpose());
      auto s2_tmp = laplace_likelihood::compute_s2(
          ll_fun, md_est.theta, (C.transpose() * C).eval(),
          options.hessian_block_size, ll_args_copy, msgs);
      s2.deep_copy(s2_tmp);
      internal::copy_compute_s2<true>(partial_parm, ll_args_filter);
    } else {  // options.solver with LU decomposition
      LU_solve_covariance = md_est.LU.solve(value_of(covariance));
      auto I_minus_BinvKW
          = Eigen::MatrixXd::Identity(md_est.W_r.rows(), md_est.W_r.cols())
            - LU_solve_covariance * md_est.W_r;
      R = md_est.W_r * I_minus_BinvKW;  // == W - W B^{-1} K W
      arena_t<Eigen::MatrixXd> A
          = value_of(covariance)
            - value_of(covariance) * md_est.W_r * LU_solve_covariance;
      auto s2_tmp = laplace_likelihood::compute_s2(ll_fun, md_est.theta, A,
                                                   options.hessian_block_size,
                                                   ll_args_copy, msgs);
      s2.deep_copy(s2_tmp);
      internal::copy_compute_s2<true>(partial_parm, ll_args_filter);
    }
    lmd = md_est.lmd;
    if constexpr (is_any_var_scalar_v<scalar_type_t<CovarArgs>>) {
      arena_t<Eigen::MatrixXd> K_adj_arena
          = 0.5 * md_est.a * md_est.a.transpose() - 0.5 * R
            + s2 * md_est.theta_grad.transpose()
            - (R * (covariance.val() * s2)) * md_est.theta_grad.transpose();
      var Z = make_callback_var(
          0.0, [covariance, K_adj_arena](auto&& vi) mutable {
            covariance.adj().array() += vi.adj() * K_adj_arena.array();
          });
      grad(Z.vi_);
      auto covar_args_filter
          = internal::filter_var_scalar_types(covar_args_copy);
      internal::collect_adjoints(covar_args_adj, covar_args_filter);
    }
    if constexpr (ll_args_contain_var) {
      arena_t<Eigen::VectorXd> v;
      if (md_est.solver_used == 1 || md_est.solver_used == 2) {
        v = value_of(covariance) * s2
            - value_of(covariance) * R * value_of(covariance) * s2;
      } else {
        v = LU_solve_covariance * s2;
      }
      laplace_likelihood::diff_eta_implicit(ll_fun, v, md_est.theta,
                                            ll_args_copy, msgs);
      internal::collect_adjoints<true>(partial_parm, ll_args_filter);
    }
  }
  var ret(lmd);
  if constexpr (is_any_var_scalar_v<CovarArgs>) {
    auto covar_args_filter = internal::filter_var_scalar_types(covar_args_refs);
    internal::reverse_pass_collect_adjoints(ret, covar_args_filter,
                                            covar_args_adj);
  }
  if constexpr (ll_args_contain_var) {
    auto ll_args_filter = internal::filter_var_scalar_types(ll_args_refs);
    internal::reverse_pass_collect_adjoints(ret, ll_args_filter, partial_parm);
  }
  return ret;
}

}  // namespace math
}  // namespace stan

#endif
