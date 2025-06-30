#ifndef STAN_MATH_MIX_FUNCTOR_LAPLACE_MARGINAL_DENSITY_HPP
#define STAN_MATH_MIX_FUNCTOR_LAPLACE_MARGINAL_DENSITY_HPP
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
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
#include <optional>

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
  /* Maximum number of steps in line search */
  int max_steps_line_search{0};
  /* iterations end when difference in objective function is less than tolerance
   */
  double tolerance{1e-6};
  /* Maximum number of steps*/
  int max_num_steps{100};
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

template <typename Covar, typename ThetaVec, typename WR, typename L_t,
          typename A_vec, typename ThetaGrad, typename LU_t, typename KRoot>
struct laplace_density_estimates {
  /* log marginal density */
  double lmd{std::numeric_limits<double>::infinity()};
  /* Evaluated covariance function for the latent gaussian variable */
  Covar covariance;
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
  laplace_density_estimates(double lmd_, Covar&& covariance_, ThetaVec&& theta_,
                            WR&& W_r_, L_t&& L_, A_vec&& a_,
                            ThetaGrad&& theta_grad_, LU_t&& LU_,
                            KRoot&& K_root_)
      : lmd(lmd_),
        covariance(std::move(covariance_)),
        theta(std::move(theta_)),
        W_r(std::move(W_r_)),
        L(std::move(L_)),
        a(std::move(a_)),
        theta_grad(std::move(theta_grad_)),
        LU(std::move(LU_)),
        K_root(std::move(K_root_)) {}
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
 * @brief Performs a simple line search
 *
 * @tparam AVec   Type of the parameter update vector (`a`), e.g.
 * Eigen::VectorXd.
 * @tparam APrev  Type of the previous parameter vector (`a_prev`), same shape
 * as AVec.
 * @tparam ThetaVec Type of the transformed vector (`theta`), e.g. Σ·a.
 * @tparam LLFun  Functor type for computing the log‐likelihood.
 * @tparam LLArgs Tuple or pack type forwarded to `ll_fun`.
 * @tparam Covar  Matrix type for the covariance Σ, e.g. Eigen::MatrixXd.
 * @tparam Msgs   Diagnostics container type for capturing warnings/errors.
 *
 * @param[in,out] objective_new On entry: objective at the full‐step `a` (must
 * satisfy objective_new < objective_old). On exit:  best objective found.
 * @param[in,out] a On entry: candidate parameter vector. On exit:  updated to
 * the step achieving the lowest objective.
 * @param[in,out] theta On entry: Σ·a for the initial candidate. On exit:  Σ·a
 * for the accepted best step.
 * @param[in,out] a_prev On entry: previous parameter vector, with objective
 * `objective_old`. On exit: rolled forward to each newly accepted step.
 * @param[in] ll_fun Callable that computes the log‐likelihood given `(theta,
 * ll_args, msgs)`.
 * @param[in] ll_args Arguments forwarded to `ll_fun` at each evaluation.
 * @param[in] covariance Covariance matrix Σ used to compute `theta = Σ·a`.
 * @param[in] max_steps_line_search Maximum number of iterations.
 * @param[in] objective_old Objective value at the initial `a_prev` (used as f₀
 * for the first pass).
 * @param[in] tolerance Minimum tolerance to accept a step
 * @param[in,out] msgs Pointer to a diagnostics container; may be used by
 * `ll_fun` to record warnings.
 */
template <typename AVec, typename APrev, typename ThetaVec, typename LLFun,
          typename LLArgs, typename Covar, typename Msgs>
inline void line_search(double& objective_new, AVec& a, ThetaVec& theta,
                        APrev& a_prev, LLFun&& ll_fun, LLArgs&& ll_args,
                        Covar&& covariance, const int max_steps_line_search,
                        const double objective_old, double tolerance,
                        Msgs* msgs) {
  Eigen::VectorXd a_tmp(a.size());
  double objective_new_tmp = 0.0;
  double objective_old_tmp = objective_old;
  Eigen::VectorXd theta_tmp(covariance.rows());
  for (int j = 0;
       j < max_steps_line_search && (objective_new < objective_old_tmp); ++j) {
    a_tmp.noalias() = a_prev + 0.5 * (a - a_prev);
    theta_tmp.noalias() = covariance * a_tmp;
    if (!theta_tmp.allFinite()) {
      break;
    } else {
      objective_new_tmp = -0.5 * a_tmp.dot(theta_tmp)
                          + laplace_likelihood::log_likelihood(
                              ll_fun, theta_tmp, ll_args, msgs);
      if (objective_new_tmp < objective_new) {
        a_prev.swap(a);
        a.swap(a_tmp);
        theta.swap(theta_tmp);
        objective_old_tmp = objective_new;
        objective_new = objective_new_tmp;
      } else {
        break;
      }
    }
  }
}

/**
 * Set all adjoints of the output to zero.
 */
template <typename Output>
inline void set_zero_adjoint(Output&& output) {
  if constexpr (is_all_arithmetic_scalar_v<Output>) {
    return;
  } else {
    return iter_tuple_nested(
        [](auto&& output_i) {
          using output_i_t = std::decay_t<decltype(output_i)>;
          if constexpr (is_all_arithmetic_scalar_v<output_i_t>) {
            return;
          } else if constexpr (is_std_vector<output_i_t>::value) {
            for (Eigen::Index i = 0; i < output_i.size(); ++i) {
              output_i[i].adj() = 0;
            }
          } else if constexpr (is_eigen_v<output_i_t>) {
            output_i.adj().setZero();
          } else if constexpr (is_stan_scalar_v<output_i_t>) {
            output_i.adj() = 0;
          } else {
            static_assert(
                sizeof(std::decay_t<output_i_t>*) == 0,
                "INTERNAL ERROR:(laplace_marginal_lpdf) set_zero_adjoints was "
                "not able to deduce the actions needed for the given type. "
                "This is an internal error, please report it: "
                "https://github.com/stan-dev/math/issues");
          }
        },
        std::forward<Output>(output));
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
  if ((Eigen::isnan(param.array()) || Eigen::isinf(param.array())).all()) {
    msg += " for all values.";
    throw std::domain_error(msg);
  }
  msg += " at indices [";
  for (int i = 0; i < param.size(); ++i) {
    if (std::isnan(param(i) || std::isinf(param(i)))) {
      msg += std::to_string(i) + ", ";
    }
  }
  msg.pop_back();
  msg.pop_back();
  msg += "].";
  throw std::domain_error(msg);
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
template <typename LLFun, typename LLTupleArgs, typename CovarFun,
          typename CovarArgs, bool InitTheta,
          require_t<is_all_arithmetic_scalar<CovarArgs>>* = nullptr>
inline auto laplace_marginal_density_est(
    LLFun&& ll_fun, LLTupleArgs&& ll_args, CovarFun&& covariance_function,
    CovarArgs&& covar_args, const laplace_options<InitTheta>& options,
    std::ostream* msgs) {
  using Eigen::MatrixXd;
  using Eigen::SparseMatrix;
  using Eigen::VectorXd;
  if constexpr (InitTheta) {
    check_nonzero_size("laplace_marginal", "initial guess", options.theta_0);
    check_finite("laplace_marginal", "initial guess", options.theta_0);
  }
  check_nonnegative("laplace_marginal", "tolerance", options.tolerance);
  check_positive("laplace_marginal", "max_num_steps", options.max_num_steps);
  check_positive("laplace_marginal", "hessian_block_size",
                 options.hessian_block_size);
  check_nonnegative("laplace_marginal", "max_steps_line_search",
                    options.max_steps_line_search);

  Eigen::MatrixXd covariance = stan::math::apply(
      [msgs, &covariance_function](auto&&... args) {
        return covariance_function(args..., msgs);
      },
      covar_args);
  check_square("laplace_marginal", "covariance", covariance);

  const Eigen::Index theta_size = covariance.rows();

  if (unlikely(theta_size % options.hessian_block_size != 0)) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << "laplace_marginal_density: The hessian size (" << theta_size
          << ", " << theta_size
          << ") is not divisible by the hessian block size ("
          << options.hessian_block_size
          << ")"
             ". Try a hessian block size such as [1, ";
      for (int i = 2; i < 12; ++i) {
        if (theta_size % i == 0) {
          msg << i << ", ";
        }
      }
      msg.str().pop_back();
      msg.str().pop_back();
      msg << "].";
      throw std::domain_error(msg.str());
    }();
  }

  auto throw_overstep = [](const auto max_num_steps) STAN_COLD_PATH {
    throw std::domain_error(
        std::string("laplace_marginal_density: max number of iterations: ")
        + std::to_string(max_num_steps) + " exceeded.");
  };
  auto ll_args_vals = value_of(ll_args);
  Eigen::VectorXd theta = [theta_size, &options]() {
    if constexpr (InitTheta) {
      return options.theta_0;
    } else {
      return Eigen::VectorXd::Zero(theta_size);
    }
  }();
  double objective_old = std::numeric_limits<double>::lowest();
  double objective_new = std::numeric_limits<double>::lowest() + 1;
  Eigen::VectorXd a_prev = Eigen::VectorXd::Zero(theta_size);
  Eigen::MatrixXd B(theta_size, theta_size);
  Eigen::VectorXd a(theta_size);
  Eigen::VectorXd b(theta_size);
  if (options.solver == 1) {
    if (options.hessian_block_size == 1) {
      for (Eigen::Index i = 0; i <= options.max_num_steps; i++) {
        auto [theta_grad, W] = laplace_likelihood::diff(
            ll_fun, theta, options.hessian_block_size, ll_args, msgs);
        Eigen::VectorXd W_r(W.rows());
        // Compute matrix square-root of W. If all elements of W are positive,
        // do an element wise square-root. Else try a matrix square-root
        for (Eigen::Index i = 0; i < W.rows(); i++) {
          if (W.coeff(i, i) < 0) {
            throw std::domain_error(
                "laplace_marginal_density: Hessian matrix is not positive "
                "definite");
          } else {
            W_r.coeffRef(i) = std::sqrt(W.coeff(i, i));
          }
        }
        B.noalias() = MatrixXd::Identity(theta_size, theta_size)
                      + W_r.asDiagonal() * covariance * W_r.asDiagonal();
        Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt_B(B);
        auto L = llt_B.matrixL();
        auto LT = llt_B.matrixU();
        b.noalias() = W.diagonal().cwiseProduct(theta) + theta_grad;
        a.noalias()
            = b
              - W_r.asDiagonal()
                    * LT.solve(L.solve(W_r.cwiseProduct(covariance * b)));
        // Simple Newton step
        theta.noalias() = covariance * a;
        objective_old = objective_new;
        if (unlikely(
                (Eigen::isinf(theta.array()) || Eigen::isnan(theta.array()))
                    .any())) {
          throw_nan("laplace_marginal_density", "theta", theta);
        }
        objective_new = -0.5 * a.dot(theta)
                        + laplace_likelihood::log_likelihood(
                            ll_fun, theta, ll_args_vals, msgs);
        if (options.max_steps_line_search) {
          line_search(objective_new, a, theta, a_prev, ll_fun, ll_args_vals,
                      covariance, options.max_steps_line_search, objective_old,
                      options.tolerance, msgs);
        }
        // Check for convergence
        if (abs(objective_new - objective_old) < options.tolerance) {
          const double B_log_determinant
              = 2.0 * llt_B.matrixLLT().diagonal().array().log().sum();
          // Overwrite W instead of making a new sparse matrix
          W.diagonal() = W_r;
          return laplace_density_estimates{
              objective_new - 0.5 * B_log_determinant,
              std::move(covariance),
              std::move(theta),
              std::move(W),
              Eigen::MatrixXd(L),
              std::move(a),
              std::move(theta_grad),
              Eigen::PartialPivLU<Eigen::MatrixXd>{},
              Eigen::MatrixXd(0, 0)};
        } else {
          a_prev = std::move(a);
          set_zero_adjoint(ll_args);
        }
      }
    } else {
      Eigen::SparseMatrix<double> W_r(theta_size, theta_size);
      Eigen::Index block_size = options.hessian_block_size;
      W_r.reserve(Eigen::VectorXi::Constant(W_r.cols(), block_size));
      const Eigen::Index n_block = W_r.cols() / block_size;
      // Prefill W_r so we only make space once
      for (Eigen::Index i = 0; i < n_block; i++) {
        for (Eigen::Index k = 0; k < block_size; k++) {
          for (Eigen::Index j = 0; j < block_size; j++) {
            W_r.insert(i * block_size + j, i * block_size + k) = 1.0;
          }
        }
      }
      W_r.makeCompressed();
      for (Eigen::Index i = 0; i <= options.max_num_steps; i++) {
        auto [theta_grad, W] = laplace_likelihood::diff(
            ll_fun, theta, options.hessian_block_size, ll_args, msgs);
        for (Eigen::Index i = 0; i < W.rows(); i++) {
          if (W.coeff(i, i) < 0) {
            throw std::domain_error(
                "laplace_marginal_density: Hessian matrix is not positive "
                "definite");
          }
        }
        block_matrix_chol_L(W_r, W, options.hessian_block_size);
        B.noalias() = MatrixXd::Identity(theta_size, theta_size)
                      + W_r * (covariance * W_r);
        Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt_B(B);
        auto L = llt_B.matrixL();
        auto LT = llt_B.matrixU();
        b.noalias() = W * theta + theta_grad;
        a.noalias() = b - W_r * LT.solve(L.solve(W_r * (covariance * b)));
        // Simple Newton step
        theta.noalias() = covariance * a;
        objective_old = objective_new;
        if (unlikely(
                (Eigen::isinf(theta.array()) || Eigen::isnan(theta.array()))
                    .any())) {
          throw_nan("laplace_marginal_density", "theta", theta);
        }
        objective_new = -0.5 * a.dot(value_of(theta))
                        + laplace_likelihood::log_likelihood(
                            ll_fun, value_of(theta), ll_args_vals, msgs);
        if (options.max_steps_line_search > 0) {
          line_search(objective_new, a, theta, a_prev, ll_fun, ll_args_vals,
                      covariance, options.max_steps_line_search, objective_old,
                      options.tolerance, msgs);
        }
        // Check for convergence
        if (abs(objective_new - objective_old) < options.tolerance) {
          const double B_log_determinant
              = 2.0 * llt_B.matrixLLT().diagonal().array().log().sum();
          return laplace_density_estimates{
              objective_new - 0.5 * B_log_determinant,
              std::move(covariance),
              std::move(theta),
              std::move(W_r),
              Eigen::MatrixXd(L),
              std::move(a),
              std::move(theta_grad),
              Eigen::PartialPivLU<Eigen::MatrixXd>{},
              Eigen::MatrixXd(0, 0)};
        } else {
          a_prev = a;
          set_zero_adjoint(ll_args);
        }
      }
    }
    throw_overstep(options.max_num_steps);
  } else if (options.solver == 2) {
    Eigen::MatrixXd K_root
        = covariance.template selfadjointView<Eigen::Lower>().llt().matrixL();
    for (Eigen::Index i = 0; i <= options.max_num_steps; i++) {
      auto [theta_grad, W] = laplace_likelihood::diff(
          ll_fun, theta, options.hessian_block_size, ll_args, msgs);
      B.noalias() = MatrixXd::Identity(theta_size, theta_size)
                    + K_root.transpose() * W * K_root;
      Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt_B(B);
      auto L = llt_B.matrixL();
      auto LT = llt_B.matrixU();
      b.noalias() = W * theta + theta_grad;
      a.noalias()
          = K_root.transpose().template triangularView<Eigen::Upper>().solve(
              LT.solve(L.solve(K_root.transpose() * b)));
      // Simple Newton step
      theta.noalias() = covariance * a;
      objective_old = objective_new;
      if (unlikely((Eigen::isinf(theta.array()) || Eigen::isnan(theta.array()))
                       .any())) {
        throw_nan("laplace_marginal_density", "theta", theta);
      }
      objective_new = -0.5 * a.dot(theta)
                      + laplace_likelihood::log_likelihood(ll_fun, theta,
                                                           ll_args_vals, msgs);
      // linesearch
      if (options.max_steps_line_search > 0) {
        line_search(objective_new, a, theta, a_prev, ll_fun, ll_args_vals,
                    covariance, options.max_steps_line_search, objective_old,
                    options.tolerance, msgs);
      }
      // Check for convergence
      if (abs(objective_new - objective_old) < options.tolerance) {
        const double B_log_determinant
            = 2.0 * llt_B.matrixLLT().diagonal().array().log().sum();
        return laplace_density_estimates{
            objective_new - 0.5 * B_log_determinant,
            std::move(covariance),
            std::move(theta),
            std::move(W),
            std::move(Eigen::MatrixXd(L)),
            std::move(a),
            std::move(theta_grad),
            Eigen::PartialPivLU<Eigen::MatrixXd>{},
            std::move(K_root)};
      } else {
        a_prev = a;
        set_zero_adjoint(ll_args);
      }
    }
    throw_overstep(options.max_num_steps);
  } else if (options.solver == 3) {
    for (Eigen::Index i = 0; i <= options.max_num_steps; i++) {
      auto [theta_grad, W] = laplace_likelihood::diff(
          ll_fun, theta, options.hessian_block_size, ll_args, msgs);
      Eigen::PartialPivLU<Eigen::MatrixXd> LU(
          MatrixXd::Identity(theta_size, theta_size) + covariance * W);
      // L on upper and U on lower triangular
      b.noalias() = W * theta + theta_grad;
      a.noalias() = b - W * LU.solve(covariance * b);
      // Simple Newton step
      theta.noalias() = covariance * a;
      objective_old = objective_new;
      if (((Eigen::isinf(theta.array()) || Eigen::isnan(theta.array()))
               .any())) {
        throw_nan("laplace_marginal_density", "theta", theta);
      }
      objective_new = -0.5 * a.dot(value_of(theta))
                      + laplace_likelihood::log_likelihood(
                          ll_fun, value_of(theta), ll_args_vals, msgs);

      if (options.max_steps_line_search > 0) {
        line_search(objective_new, a, theta, a_prev, ll_fun, ll_args_vals,
                    covariance, options.max_steps_line_search, objective_old,
                    options.tolerance, msgs);
      }
      if (abs(objective_new - objective_old) < options.tolerance) {
        // TODO(Charles): There has to be a simple trick for this
        const double B_log_determinant = log(LU.determinant());
        return laplace_density_estimates{
            objective_new - 0.5 * B_log_determinant,
            std::move(covariance),
            std::move(theta),
            std::move(W),
            Eigen::MatrixXd(0, 0),
            std::move(a),
            std::move(theta_grad),
            std::move(LU),
            Eigen::MatrixXd(0, 0)};
      } else {
        a_prev = a;
        set_zero_adjoint(ll_args);
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
  return internal::laplace_marginal_density_est(
             std::forward<LLFun>(ll_fun), std::forward<LLTupleArgs>(ll_args),
             std::forward<CovarFun>(covariance_function),
             std::forward<CovarArgs>(covar_args), options, msgs)
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
    using laplace_likelihood::internal::conditional_copy_and_promote;
    using laplace_likelihood::internal::COPY_TYPE;
    auto ll_args_copy
        = conditional_copy_and_promote<is_any_var_scalar, var, COPY_TYPE::DEEP>(
            ll_args_refs);

    auto md_est = internal::laplace_marginal_density_est(
        ll_fun, ll_args_copy, covariance_function, value_of(covar_args_refs),
        options, msgs);

    // Solver 1, 2
    arena_t<Eigen::MatrixXd> R(md_est.theta.size(), md_est.theta.size());
    // Solver 3
    arena_t<Eigen::MatrixXd> LU_solve_covariance;
    // Solver 1, 2, 3
    arena_t<Eigen::VectorXd> s2(md_est.theta.size());

    // Return references to var types
    auto ll_args_filter = internal::filter_var_scalar_types(ll_args_copy);
    stan::math::for_each(
        [](auto&& output_i, auto&& ll_arg_i) {
          if (is_any_var_scalar_v<decltype(ll_arg_i)>) {
            internal::collect_adjoints<true>(output_i, ll_arg_i);
          }
        },
        partial_parm, ll_args_filter);
    if (options.solver == 1) {
      if (options.hessian_block_size == 1) {
        // TODO(Steve): Solve without casting from sparse to dense
        Eigen::MatrixXd tmp
            = md_est.L.template triangularView<Eigen::Lower>().solve(
                md_est.W_r.toDense());
        R = tmp.transpose() * tmp;
        arena_t<Eigen::MatrixXd> C
            = md_est.L.template triangularView<Eigen::Lower>().solve(
                md_est.W_r * md_est.covariance);
        if constexpr (!ll_args_contain_var) {
          s2.deep_copy(
              (0.5
               * (md_est.covariance.diagonal() - (C.transpose() * C).diagonal())
                     .cwiseProduct(laplace_likelihood::third_diff(
                         ll_fun, md_est.theta, value_of(ll_args_copy), msgs))));
        } else {
          arena_t<Eigen::MatrixXd> A = md_est.covariance - C.transpose() * C;
          auto s2_tmp = laplace_likelihood::compute_s2(
              ll_fun, md_est.theta, A, options.hessian_block_size, ll_args_copy,
              msgs);
          s2.deep_copy(s2_tmp);
          internal::copy_compute_s2<true>(partial_parm, ll_args_filter);
        }

      } else {
        Eigen::MatrixXd tmp
            = md_est.L.template triangularView<Eigen::Lower>().solve(
                md_est.W_r.toDense());
        R = tmp.transpose() * tmp;
        arena_t<Eigen::MatrixXd> C
            = md_est.L.template triangularView<Eigen::Lower>().solve(
                md_est.W_r * md_est.covariance);
        arena_t<Eigen::MatrixXd> A = md_est.covariance - C.transpose() * C;
        auto s2_tmp = laplace_likelihood::compute_s2(ll_fun, md_est.theta, A,
                                                     options.hessian_block_size,
                                                     ll_args_copy, msgs);
        s2.deep_copy(s2_tmp);
        internal::copy_compute_s2<true>(partial_parm, ll_args_filter);
      }
    } else if (options.solver == 2) {
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
      LU_solve_covariance = md_est.LU.solve(md_est.covariance);
      R = md_est.W_r - md_est.W_r * LU_solve_covariance * md_est.W_r;
      arena_t<Eigen::MatrixXd> A
          = md_est.covariance
            - md_est.covariance * md_est.W_r * LU_solve_covariance;
      auto s2_tmp = laplace_likelihood::compute_s2(ll_fun, md_est.theta, A,
                                                   options.hessian_block_size,
                                                   ll_args_copy, msgs);
      s2.deep_copy(s2_tmp);
      internal::copy_compute_s2<true>(partial_parm, ll_args_filter);
    }
    lmd = md_est.lmd;
    if constexpr (is_any_var_scalar_v<scalar_type_t<CovarArgs>>) {
      [&covar_args_refs, &covar_args_adj, &md_est, &R, &s2,
       &covariance_function, &msgs]() mutable {
        const nested_rev_autodiff nested;
        auto covar_args_copy
            = laplace_likelihood::internal::conditional_copy_and_promote<
                is_any_var_scalar, var,
                laplace_likelihood::internal::COPY_TYPE::DEEP>(covar_args_refs);

        var_value<Eigen::MatrixXd> K_var = to_var_value(stan::math::apply(
            [&covariance_function, &msgs](auto&&... args) {
              return covariance_function(args..., msgs);
            },
            covar_args_copy));
        arena_t<Eigen::MatrixXd> K_adj_arena
            = 0.5 * md_est.a * md_est.a.transpose() - 0.5 * R
              + s2 * md_est.theta_grad.transpose()
              - (R * (K_var.val() * s2)) * md_est.theta_grad.transpose();
        var Z = make_callback_var(0.0, [K_var, K_adj_arena](auto&& vi) mutable {
          K_var.adj().array() += vi.adj() * K_adj_arena.array();
        });
        grad(Z.vi_);
        auto covar_args_filter
            = internal::filter_var_scalar_types(covar_args_copy);
        internal::collect_adjoints(covar_args_adj, covar_args_filter);
      }();
    }
    if constexpr (ll_args_contain_var) {
      arena_t<Eigen::VectorXd> v;
      if (options.solver == 1 || options.solver == 2) {
        v = md_est.covariance * s2
            - md_est.covariance * R * md_est.covariance * s2;
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
