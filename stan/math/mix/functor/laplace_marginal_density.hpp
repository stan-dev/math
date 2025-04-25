#ifndef STAN_MATH_MIX_FUNCTOR_LAPLACE_MARGINAL_DENSITY_HPP
#define STAN_MATH_MIX_FUNCTOR_LAPLACE_MARGINAL_DENSITY_HPP
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/quad_form_diag.hpp>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>

#include <cmath>

// Reference for calculations of marginal and its gradients:
// Margossian et al (2020), https://arxiv.org/abs/2004.12550
// and Margossian (2022), https://doi.org/10.7916/0wsc-kz90

// TODO(Charles) -- either use Eigen's .solve() or mdivide_left_tri
// The code needs to be more consistent

namespace stan {
namespace math {

/**
 * Options for the laplace sampler
 */
struct laplace_options {
  /* Size of the blocks in block diagonal hessian*/
  int hessian_block_size{1};
  /**
   * Which Newton solver to use:
   * (1) method using the root of W
   * (2) method using the root of the covariance
   * (3) method using an LU decomposition
   */
  int solver{1};
  /* Maximum number of steps in line search*/
  int max_steps_line_search{0};
  /* iterations end when difference in objective function is less than tolerance
   */
  double tolerance{1e-12};
  /* Maximum number of steps*/
  int64_t max_num_steps{100};
};

template <typename Covar, typename Theta, typename WR, typename L_t,
          typename A_vec, typename ThetaGrad, typename LU_t, typename KRoot>
struct laplace_density_estimates {
  // log marginal density
  double lmd{std::numeric_limits<double>::infinity()};
  // Evaluated covariance function for the latent gaussian variable
  Covar covariance;
  // Mode
  Theta theta;
  // the square root of the negative Hessian or the negative Hessian, depending
  // on which solver we use
  WR W_r;
  // cholesky decomposition of stabilized inverse covariance
  L_t L;
  // element in the Newton step
  A_vec a;
  // the gradient of the log density with respect to theta
  ThetaGrad theta_grad;
  // LU matrix
  LU_t LU;
  // Cholesky of the covariance matrix
  KRoot K_root;
  laplace_density_estimates(double lmd_, Covar&& covariance_, Theta&& theta_,
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
 * Function to compute the pseudo target, $\tilde Z$,
 * with a custom derivative method
 * NOTE: we actually don't need to compute the pseudo-target, only its
 * derivative
 * @tparam Kmat Type inheriting from `Eigen::EigenBase` with dynamic rows and
 * columns
 * @tparam AVec Type of matrix of initial tangents
 * @tparam RMat Type of the stable R matrix
 * @tparam LGradVec Type of the gradient of the log likelihood
 * @tparam S2Vec Type of the s2 vector
 */
template <
    typename KMat, typename AVec, typename RMat, typename LGradVec,
    typename S2Vec,
    require_eigen_matrix_dynamic_vt<std::is_floating_point, KMat>* = nullptr>
inline constexpr double laplace_pseudo_target(KMat&& /* K */, AVec&& /* a */,
                                              RMat&& /* R */,
                                              LGradVec&& /* l_grad */,
                                              S2Vec&& /* s2 */) {
  return static_cast<double>(0.0);
}

/**
 * Overload function for case where K is passed as a matrix of var
 * @tparam Kmat Type inheriting from `Eigen::EigenBase` with dynamic rows and
 * columns
 * @tparam AVec Type inheriting from `Eigen::EigenBase` with dynamic columns and
 * a single row
 * @tparam RMat Type inheriting from `Eigen::EigenBase` with dynamic rows and
 * columns
 * @tparam LGradVec Type inheriting from `Eigen::EigenBase` with dynamic rows
 * and a single column
 * @tparam S2Vec Type of s2 vector
 * @param K Covariance matrix
 * @param a Saved a vector from Newton solver
 * @param R Stable R matrix
 * @param l_grad Saved gradient of log likelihood
 * @param s2 Gradient of log determinant w.r.t latent Gaussian variable
 */
template <typename KMat, typename AVec, typename RMat, typename LGradVec,
          typename S2Vec,
          require_eigen_matrix_dynamic_vt<is_var, KMat>* = nullptr>
inline auto laplace_pseudo_target(KMat&& K, AVec&& a, RMat&& R,
                                  LGradVec&& l_grad, S2Vec&& s2) {
  const Eigen::Index dim_theta = K.rows();
  auto K_arena = to_arena(std::forward<KMat>(K));
  auto&& a_ref = to_ref(std::forward<AVec>(a));
  auto&& R_ref = to_ref(std::forward<RMat>(R));
  auto&& s2_ref = to_ref(std::forward<S2Vec>(s2));
  auto&& l_grad_ref = to_ref(std::forward<LGradVec>(l_grad));
  arena_matrix<Eigen::MatrixXd> K_adj_arena
      = 0.5 * a_ref * a_ref.transpose() - 0.5 * R_ref
        + s2_ref * l_grad_ref.transpose()
        - (R_ref * (value_of(K_arena) * s2_ref)) * l_grad_ref.transpose();
  return make_callback_var(0.0, [K_arena, K_adj_arena](auto&& vi) mutable {
    K_arena.adj().array() += vi.adj() * K_adj_arena.array();
  });
}

/**
 * Return the matrix square-root for a block diagonal matrix
 * @param W Block-diagonal matrix
 * @param block_size Size of blocks in W
 */
inline Eigen::SparseMatrix<double> block_matrix_sqrt(
    const Eigen::SparseMatrix<double>& W, const Eigen::Index block_size) {
  int n_block = W.cols() / block_size;
  Eigen::MatrixXd local_block(block_size, block_size);
  Eigen::MatrixXd local_block_sqrt(block_size, block_size);
  Eigen::MatrixXd sqrt_t_mat = Eigen::MatrixXd::Zero(block_size, block_size);
  Eigen::SparseMatrix<double> W_root(W.rows(), W.cols());
  W_root.reserve(Eigen::VectorXi::Constant(W_root.cols(), block_size));

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
        W_root.insert(i * block_size + j, i * block_size + k)
            = local_block_sqrt(j, k);
      }
    }
  }
  W_root.makeCompressed();
  return W_root;
}
template <typename AVec, typename APrev, typename ThetaVec, typename LLFun,
          typename LLArgs, typename Covar, typename Msgs>
inline auto line_search(double& objective_new, AVec&& a, APrev& a_prev,
                        ThetaVec&& theta, LLFun&& ll_fun, LLArgs&& ll_args,
                        Covar&& covariance, const int max_steps_line_search,
                        const double objective_old, Msgs* msgs) {
  Eigen::VectorXd a_tmp(a.size());
  double objective_new_tmp = 0.0;
  double objective_old_tmp = objective_old;
  Eigen::VectorXd theta_tmp(covariance.rows());
  int j = 0;
  for (; j < max_steps_line_search && (objective_new < objective_old_tmp);
       ++j) {
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
  return std::make_tuple(objective_new, std::move(a), std::move(theta));
}

template <bool ZeroInput = false, typename Output, typename Input1,
          require_t<is_any_var_scalar<Input1>>* = nullptr>
inline void collect_adjoints(Output& output, Input1&& precalc) {
  if constexpr (is_tuple<Output>::value) {
    stan::math::for_each(
        [](auto& output_i, auto&& precalc_i) {
          collect_adjoints(output_i, precalc_i);
        },
        output, precalc);
  } else if constexpr (is_std_vector<Output>::value) {
    if constexpr (is_stan_scalar<value_type_t<Output>>::value) {
      Eigen::Map<Eigen::Matrix<double, -1, 1>> output_map(output.data(),
                                                          output.size());
      Eigen::Map<Eigen::Matrix<var, -1, 1>> precalc_map(precalc.data(),
                                                        precalc.size());
      output_map.array() += precalc_map.adj().array();
      if constexpr (ZeroInput) {
        precalc_map.adj().setZero();
      }
    } else {
      const auto output_size = output.size();
      for (std::size_t i = 0; i < output_size; ++i) {
        collect_adjoints(output[i], precalc[i]);
      }
    }
  } else if constexpr (is_eigen<Output>::value) {
    output.array() += precalc.adj().array();
    if constexpr (ZeroInput) {
      precalc.adj().setZero();
    }
  } else if constexpr (is_stan_scalar_v<Output>) {
    output += precalc.adj();
    if constexpr (ZeroInput) {
      precalc.adj() = 0;
    }
  } else {
    static_assert(1, "We missed!!!");
  }
}

template <typename NameStr, typename ParamStr, typename Param>
STAN_COLD_PATH void throw_nan(NameStr&& name_str, ParamStr&& param_str,
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
 * TODO(Charles): add more robust convergence criterion.
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
 * @tparam LLFun Type with a valid `operator(Theta,  InnerLLTupleArgs)`
 * where `InnerLLTupleArgs` are the elements of `LLTupleArgs`
 * @tparam LLTupleArgs A tuple whose elements follow the types required for
 * `LLFun`
 * @tparam CovarFun Type with a valid `operator(InnerCovarArgs)` where
 * `InnerCovarArgs` are a parameter pack of the element types of
 * `CovarTupleArgs`
 * @tparam Theta Type derived from `Eigen::EigenBase` with dynamic rows and a
 * single column
 * @tparam CovarTupleArgs A tuple whose elements follow the types required for
 * `CovarFun`
 * @param[in] ll_fun A log likelihood functor
 * @param[in] ll_args Tuple containing parameters for `LLFun`
 * @param[in] covariance_function Functor for the covariance function
 * @param[in] theta_0 the initial guess for the Laplace optimization
 * @param[in,out] msgs stream for messages from likelihood and covariance
 * @param[in] options A set of options for tuning the solver
 * @param[in] covar_args Tuple of arguments to pass to the covariance matrix
 * functor
 *
 * @return A struct containing
 * 1. lmd the log marginal density, p(y | phi)
 * 2. covariance the evaluated covariance function for the latent gaussian
 * variable
 * 3. theta a vector to store the mode
 * 4. W_r a vector to store the square root of the
 *                 negative Hessian or the negative Hessian, depending
 *                 on which solver we use
 * 5. L cholesky decomposition of stabilized inverse covariance
 * 6. a element in the Newton step
 * 7. l_grad the log density of the likelihood, evaluated at the mode
 *
 */
template <typename LLFun, typename LLTupleArgs, typename CovarFun,
          typename Theta, typename CovarTupleArgs,
          require_t<is_all_arithmetic_scalar<Theta, CovarTupleArgs>>* = nullptr,
          require_eigen_vector_t<Theta>* = nullptr>
inline auto laplace_marginal_density_est(LLFun&& ll_fun, LLTupleArgs&& ll_args,
                                         Theta&& theta_0,
                                         CovarFun&& covariance_function,
                                         CovarTupleArgs&& covar_args,
                                         const laplace_options& options,
                                         std::ostream* msgs) {
  using Eigen::MatrixXd;
  using Eigen::SparseMatrix;
  using Eigen::VectorXd;
  check_nonzero_size("laplace_marginal", "initial guess", theta_0);
  check_finite("laplace_marginal", "initial guess", theta_0);
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
  auto throw_overstep = [](const auto max_num_steps) STAN_COLD_PATH {
    throw std::domain_error(
        std::string("laplace_marginal_density: max number of iterations: ")
        + std::to_string(max_num_steps) + " exceeded.");
  };
  auto ll_args_vals = value_of(ll_args);
  const Eigen::Index theta_size = theta_0.size();
  std::decay_t<Theta> theta = theta_0;
  double objective_old = std::numeric_limits<double>::lowest();
  double objective_new = std::numeric_limits<double>::lowest() + 1;
  Eigen::VectorXd a_prev = Eigen::VectorXd::Zero(theta_size);
  Eigen::MatrixXd B(theta_size, theta_size);
  Eigen::VectorXd a(theta_size);
  Eigen::VectorXd b(theta_size);
  if (options.solver == 1 && options.hessian_block_size == 1) {
    for (Eigen::Index i = 0; i <= options.max_num_steps; i++) {
      auto [theta_grad, W] = laplace_likelihood::diff(
          ll_fun, theta, options.hessian_block_size, ll_args, msgs);

      // Compute matrix square-root of W. If all elements of W are positive,
      // do an element wise square-root. Else try a matrix square-root
      for (Eigen::Index i = 0; i < W.rows(); i++) {
        if (W.coeff(i, i) < 0) {
          throw std::domain_error(
              "laplace_marginal_density: Hessian matrix is not positive "
              "definite");
        }
      }
      Eigen::SparseMatrix<double> W_r = W.cwiseSqrt();
      // TODO(Charles): Need better way to handle negative diagonals
      /*
      if (W_is_spd) {
        W_r = W.cwiseSqrt();
      } else {
        W_r = block_matrix_sqrt(W, options.hessian_block_size);
      }
      */
      // TODO(Steve): Memory can be made once out of the loop
      // This is our main cost
      B = MatrixXd::Identity(theta_size, theta_size)
          + W_r.diagonal().asDiagonal() * covariance
                * W_r.diagonal().asDiagonal();
      Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt_B(B);
      auto L = llt_B.matrixL();
      auto LT = llt_B.matrixU();
      b.noalias() = W.diagonal().cwiseProduct(theta) + theta_grad;
      a.noalias() = b
                    - W_r
                          * LT.solve(L.solve(
                              W_r.diagonal().cwiseProduct(covariance * b)));
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
      if (options.max_steps_line_search) {
        std::tie(objective_new, a, theta)
            = line_search(objective_new, std::move(a), a_prev, std::move(theta),
                          ll_fun, ll_args_vals, covariance,
                          options.max_steps_line_search, objective_old, msgs);
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
            std::move(Eigen::MatrixXd(L)),
            std::move(a),
            std::move(theta_grad),
            Eigen::PartialPivLU<Eigen::MatrixXd>{},
            Eigen::MatrixXd(0, 0)};
      } else {
        a_prev = std::move(a);
        laplace_likelihood::internal::set_zero_adjoint(ll_args);
      }
    }
    throw_overstep(options.max_num_steps);
  } else if (options.solver == 1 && !(options.hessian_block_size == 1)) {
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
      Eigen::SparseMatrix<double> W_r
          = block_matrix_sqrt(W, options.hessian_block_size);
      B.noalias() = MatrixXd::Identity(theta_size, theta_size)
                    + W_r * (covariance * W_r);
      Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt_B(B);
      auto L = llt_B.matrixL();
      auto LT = llt_B.matrixU();
      b.noalias() = W * theta + theta_grad;
      a.noalias() = b - W_r * LT.solve(L.solve(W_r * (covariance * b)));
      // Simple Newton step
      theta = covariance * a;
      objective_old = objective_new;
      if (unlikely((Eigen::isinf(theta.array()) || Eigen::isnan(theta.array()))
                       .any())) {
        throw_nan("laplace_marginal_density", "theta", theta);
      }
      objective_new = -0.5 * a.dot(value_of(theta))
                      + laplace_likelihood::log_likelihood(
                          ll_fun, value_of(theta), ll_args_vals, msgs);
      if (options.max_steps_line_search > 0) {
        std::tie(objective_new, a, theta)
            = line_search(objective_new, std::move(a), a_prev, std::move(theta),
                          ll_fun, ll_args_vals, covariance,
                          options.max_steps_line_search, objective_old, msgs);
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
            std::move(Eigen::MatrixXd(L)),
            std::move(a),
            std::move(theta_grad),
            Eigen::PartialPivLU<Eigen::MatrixXd>{},
            Eigen::MatrixXd(0, 0)};
      } else {
        a_prev = a;
        laplace_likelihood::internal::set_zero_adjoint(ll_args);
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
      // TODO(Charles) Throw if theta is not finite?
      if (unlikely((Eigen::isinf(theta.array()) || Eigen::isnan(theta.array()))
                       .any())) {
        throw_nan("laplace_marginal_density", "theta", theta);
      }
      objective_new = -0.5 * a.dot(theta)
                      + laplace_likelihood::log_likelihood(ll_fun, theta,
                                                           ll_args_vals, msgs);
      // linesearch
      if (options.max_steps_line_search > 0) {
        std::tie(objective_new, a, theta)
            = line_search(objective_new, std::move(a), a_prev, std::move(theta),
                          ll_fun, ll_args_vals, covariance,
                          options.max_steps_line_search, objective_old, msgs);
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
        laplace_likelihood::internal::set_zero_adjoint(ll_args);
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
      theta = covariance * a;
      objective_old = objective_new;
      if (((Eigen::isinf(theta.array()) || Eigen::isnan(theta.array()))
               .any())) {
        throw_nan("laplace_marginal_density", "theta", theta);
      }
      objective_new = -0.5 * a.dot(value_of(theta))
                      + laplace_likelihood::log_likelihood(
                          ll_fun, value_of(theta), ll_args_vals, msgs);

      // TODO(Charles): How do we handle NA values in theta?
      // linesearch
      // CHECK -- does linesearch work for options.solver 2?
      if (options.max_steps_line_search > 0) {
        std::tie(objective_new, a, theta)
            = line_search(objective_new, std::move(a), a_prev, std::move(theta),
                          ll_fun, ll_args_vals, covariance,
                          options.max_steps_line_search, objective_old, msgs);
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
        laplace_likelihood::internal::set_zero_adjoint(ll_args);
      }
    }
    throw_overstep(options.max_num_steps);
  }
  throw std::domain_error(
      std::string("You chose a solver (") + std::to_string(options.solver)
      + ") that is not valid. Please choose either 1, 2, or 3.");
}

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
 * @tparam LLFun Type with a valid `operator(Theta, InnerLLTupleArgs)`
 * where `InnerLLTupleArgs` are the elements of `LLTupleArgs`
 * @tparam LLTupleArgs A tuple whose elements follow the types required for
 * `LLFun`
 * @tparam CovarFun Type with a valid `operator(InnerCovarTupleArgs)` where
 * `InnerCovarTupleArgs` are a parameter pack of the element types of
 * `CovarTupleArgs`
 * @tparam Theta Type derived from `Eigen::EigenBase` with dynamic rows and a
 * single column
 * @tparam CovarTupleArgs A tuple whose elements follow the types required for
 * `CovarFun`
 * @param[in] ll_fun A log likelihood functor
 * @param[in] ll_args Tuple containing parameters for `LLFun`
 * @param[in] covariance_function Functor for the covariance function
 * @param[in] theta_0 the initial guess for the Laplace optimization
 * @param[in,out] msgs stream for messages from likelihood and covariance
 * @param[in] options A set of options for tuning the solver
 * @param[in] covar_args Tuple of arguments to pass to the covariance matrix
 * functor
 * @return the log maginal density, p(y | phi)
 */
template <typename LLFun, typename LLTupleArgs, typename CovarFun,
          typename Theta, typename CovarTupleArgs,
          require_t<is_all_arithmetic_scalar<Theta, CovarTupleArgs,
                                             LLTupleArgs>>* = nullptr,
          require_eigen_vector_t<Theta>* = nullptr>
inline double laplace_marginal_density(LLFun&& ll_fun, LLTupleArgs&& ll_args,
                                       Theta&& theta_0,
                                       CovarFun&& covariance_function,
                                       CovarTupleArgs&& covar_args,
                                       const laplace_options& options,
                                       std::ostream* msgs) {
  return laplace_marginal_density_est(
             std::forward<LLFun>(ll_fun), std::forward<LLTupleArgs>(ll_args),
             std::forward<Theta>(theta_0),
             std::forward<CovarFun>(covariance_function),
             std::forward<CovarTupleArgs>(covar_args), options, msgs)
      .lmd;
}
template <typename Output, typename Input1>
inline void collect_adjoints(Output&& output, const vari* ret,
                             Input1&& precalc) {
  if constexpr (is_tuple<Output>::value) {
    static_assert(1,
                  "Accumulate Adjoints called on a tuple, but tuples cannot be "
                  "on the reverse mode stack!");
  } else if constexpr (is_std_vector<Output>::value) {
    if constexpr (!is_var<value_type_t<Output>>::value) {
      const auto output_size = output.size();
      for (std::size_t i = 0; i < output_size; ++i) {
        collect_adjoints(output[i], ret, precalc[i]);
      }
    } else {
      Eigen::Map<Eigen::Matrix<var, -1, 1>> output_map(output.data(),
                                                       output.size());
      Eigen::Map<Eigen::Matrix<double, -1, 1>> precalc_map(precalc.data(),
                                                           precalc.size());
      output_map.array().adj() += ret->adj_ * precalc_map.array();
    }
  } else if constexpr (is_eigen<Output>::value) {
    output.adj().array() += ret->adj_ * precalc.array();
  } else if constexpr (is_var<Output>::value) {
    output.adj() += ret->adj_ * precalc;
  }
}

template <typename Output, typename Input1,
          require_t<is_all_arithmetic_scalar<Input1>>* = nullptr>
inline void collect_adjoints(Output&& output, Input1&& precalc) {
  if constexpr (is_tuple<Output>::value) {
    stan::math::for_each(
        [](auto&& output_i, auto&& precalc_i) {
          collect_adjoints(output_i, precalc_i);
        },
        output, precalc);
  } else if constexpr (is_std_vector<Output>::value) {
    const auto output_size = output.size();
    for (std::size_t i = 0; i < output_size; ++i) {
      collect_adjoints(output[i], precalc[i]);
    }
    if constexpr (!is_stan_scalar<value_type_t<Output>>::value) {
    } else {
      Eigen::Map<Eigen::Matrix<double, -1, 1>> output_map(output.data(),
                                                          output.size());
      Eigen::Map<Eigen::Matrix<double, -1, 1>> precalc_map(precalc.data(),
                                                           precalc.size());
      output_map.array() += precalc_map.array();
    }
  } else if constexpr (is_eigen<Output>::value) {
    output.array() += precalc.array();
  } else if constexpr (is_stan_scalar_v<Output>) {
    output += precalc;
  } else {
    static_assert(1, "We missed!!!");
  }
}

inline void constexpr copy_compute_s2(const std::tuple<>& output,
                                      const std::tuple<>& precalc) noexcept {}

template <typename Output, typename Input1,
          require_t<is_any_var_scalar<Input1>>* = nullptr>
inline void copy_compute_s2(Output&& output, Input1&& precalc) {
  if constexpr (is_tuple_v<Output>) {
    stan::math::for_each(
        [](auto& output_i, auto&& precalc_i) {
          if constexpr (is_any_var_scalar<Input1>::value) {
            copy_compute_s2(output_i, precalc_i);
          }
        },
        output, precalc);
  } else if constexpr (is_std_vector<Output>::value) {
    if constexpr (!is_stan_scalar<value_type_t<Output>>::value) {
      const auto output_size = output.size();
      for (std::size_t i = 0; i < output_size; ++i) {
        copy_compute_s2(output[i], precalc[i]);
      }
    } else {
      Eigen::Map<Eigen::Matrix<double, -1, 1>> output_map(output.data(),
                                                          output.size());
      Eigen::Map<Eigen::Matrix<var, -1, 1>> precalc_map(precalc.data(),
                                                        precalc.size());
      output_map.array() += 0.5 * precalc_map.adj().array();
    }
  } else if constexpr (is_eigen<Output>::value) {
    output.array() += 0.5 * precalc.adj().array();
  } else if constexpr (is_stan_scalar_v<Output>) {
    output += (0.5 * precalc.adj());
  } else {
    static_assert(1, "We missed!!!");
  }
}

template <typename Output, typename Input1, typename Input2,
          require_t<is_all_arithmetic_scalar<Input1, Input2>>* = nullptr>
inline void collect_adjoints(Output&& output, Input1&& precalc1,
                             Input2&& precalc2) {
  if constexpr (is_tuple_v<Output>) {
    stan::math::for_each(
        [](auto&& output_i, auto&& precalc1_i, auto&& precalc2_i) {
          collect_adjoints(output_i, precalc1_i, precalc2_i);
        },
        output, precalc1, precalc2);
  } else if constexpr (is_std_vector<Output>::value) {
    if constexpr (!is_stan_scalar<value_type_t<Output>>::value) {
      const auto output_size = output.size();
      for (std::size_t i = 0; i < output_size; ++i) {
        collect_adjoints(output[i], precalc1[i], precalc2[i]);
      }
    } else {
      Eigen::Map<Eigen::Matrix<double, -1, 1>> output_map(output.data(),
                                                          output.size());
      Eigen::Map<Eigen::Matrix<double, -1, 1>> precalc1_map(precalc1.data(),
                                                            precalc1.size());
      Eigen::Map<Eigen::Matrix<double, -1, 1>> precalc2_map(precalc2.data(),
                                                            precalc2.size());
      output_map.array() += precalc1_map.array() + precalc2_map.array();
    }
  } else if constexpr (is_eigen<Output>::value) {
    output.array() += precalc1.array() + precalc2.array();
  } else if constexpr (is_stan_scalar_v<Output>) {
    output += precalc1 + precalc2;
  } else {
    static_assert(1, "Collect adjoints missed!!!");
  }
}

template <typename T>
static constexpr bool is_dbl_nothrow_constructible_v
    = std::is_nothrow_constructible<
        promote_scalar_t<double, std::decay_t<T>>>::value;

template <typename Output>
inline constexpr auto make_zero(Output&& output) {
  if constexpr (is_tuple<Output>::value) {
    return stan::math::filter_map<is_any_var_scalar>(
        [](auto&& output_i) {
          return make_zero(output_i);
        },
        output);
  } else if constexpr (is_std_vector<Output>::value) {
    if constexpr (!is_var<value_type_t<Output>>::value) {
      const auto output_size = output.size();
      arena_t<promote_scalar_t<double, Output>> ret;
      ret.reserve(output_size);
      for (Eigen::Index i = 0; i < output_size; ++i) {
        ret.push_back(make_zero(output[i]));
      }
      return ret;
    } else {
      return arena_t<std::vector<double>>(output.size(), 0.0);
    }
  } else if constexpr (is_eigen<Output>::value) {
    return arena_t<promote_scalar_t<double, Output>>(
        plain_type_t<promote_scalar_t<double, Output>>::Zero(output.rows(),
                                                             output.cols()));
  } else if constexpr (is_var<Output>::value) {
    return static_cast<double>(0.0);
  }
}

template <typename Output,
          require_t<is_all_arithmetic_scalar<Output>>* = nullptr>
inline void print_adjoint(Output&& output) {
  if constexpr (is_tuple<Output>::value) {
    std::cout << "tuple adj\n";
    return stan::math::for_each(
        [](auto&& output_i) {
          return print_adjoint(output_i);
        },
        output);
  } else if constexpr (is_std_vector<Output>::value) {
    if constexpr (is_var<value_type_t<Output>>::value) {
      Eigen::Map<const Eigen::Matrix<double, -1, -1>> map_x(output.data(),
                                                            output.size());
      std::cout << "eigen adj: \n" << map_x << std::endl;
    } else {
      std::cout << "stdvec adjoint\n";
      for (int i = 0; i < output.size(); ++i) {
        print_adjoint(output[i]);
      }
    }
  } else if constexpr (is_eigen<Output>::value) {
    std::cout << "adj: \n" << output << std::endl;
  } else if constexpr (is_stan_scalar_v<Output>) {
    std::cout << "adj: " << output << std::endl;
  } else {
    static_assert(1, "print missed!!!");
  }
}

template <typename Output, require_t<is_any_var_scalar<Output>>* = nullptr>
inline void print_adjoint(Output&& output) {
  if constexpr (is_tuple<Output>::value) {
    std::cout << "tuple adj\n";
    return stan::math::for_each(
        [](auto&& output_i) {
          return print_adjoint(output_i);
        },
        output);
  } else if constexpr (is_std_vector<Output>::value) {
    if constexpr (is_var<value_type_t<Output>>::value) {
      Eigen::Map<const Eigen::Matrix<var, -1, -1>> map_x(output.data(),
                                                         output.size());
      std::cout << "eigen adj: \n" << map_x.adj() << std::endl;
    } else {
      std::cout << "stdvec adjoint\n";
      for (int i = 0; i < output.size(); ++i) {
        print_adjoint(output[i]);
      }
    }
  } else if constexpr (is_eigen<Output>::value) {
    std::cout << "adj: \n" << output.adj() << std::endl;
  } else if constexpr (is_stan_scalar_v<Output>) {
    std::cout << "adj: " << output.adj() << std::endl;
  } else {
    static_assert(1, "print missed!!!");
  }
}

template <typename Arg, typename Precalc>
inline void laplace_tuple_collect_adjoints(var ret, Arg&& arg,
                                           Precalc&& precalc) {
  if constexpr (is_tuple_v<Arg>) {
    stan::math::for_each(
        [ret](auto&& inner_arg, auto&& inner_precalc) mutable {
          laplace_tuple_collect_adjoints(ret,
           std::forward<decltype(inner_arg)>(inner_arg),
           std::forward<decltype(inner_precalc)>(inner_precalc));
        },
        std::forward<Arg>(arg), std::forward<Precalc>(precalc));
  } else {
    reverse_pass_callback(
        [vi = ret.vi_, arg_arena = to_arena(std::forward<Arg>(arg)),
        precalc_arena = to_arena(std::forward<Precalc>(precalc))]() mutable {
          collect_adjoints(arg_arena, vi, precalc_arena);
        });
  }
}

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
 * @tparam LLFun Type with a valid `operator(Theta,  InnerLLTupleArgs)`
 * where `InnerLLTupleArgs` are the elements of `LLTupleArgs`
 * @tparam LLTupleArgs A tuple whose elements follow the types required for
 * `LLFun`
 * @tparam CovarFun Type with a valid `operator(InnerCovarTupleArgs)` where
 * `InnerCovarTupleArgs` are a parameter pack of the element types of
 * `CovarTupleArgs`
 * @tparam Theta Type derived from `Eigen::EigenBase` with dynamic rows and a
 * single column
 * @tparam CovarTupleArgs A tuple whose elements follow the types required for
 * `CovarFun`
 * @param[in] ll_fun A log likelihood functor
 * @param[in] ll_args Tuple containing parameters for `LLFun`
 * @param[in] covariance_function Functor for the covariance function
 * @param[in] theta_0 the initial guess for the Laplace optimization
 * @param[in,out] msgs stream for messages from likelihood and covariance
 * @param[in] options A set of options for tuning the solver
 * @param[in] covar_args Tuple of arguments to pass to the covariance matrix
 * functor
 * @return the log maginal density, p(y | phi)
 */
template <
    typename LLFun, typename LLTupleArgs, typename CovarFun, typename Theta,
    typename CovarTupleArgs,
    require_t<is_any_var_scalar<Theta, LLTupleArgs, CovarTupleArgs>>* = nullptr,
    require_eigen_vector_t<Theta>* = nullptr>
inline auto laplace_marginal_density(const LLFun& ll_fun, LLTupleArgs&& ll_args,
                                     Theta&& theta_0,
                                     CovarFun&& covariance_function,
                                     CovarTupleArgs&& covar_args,
                                     const laplace_options& options,
                                     std::ostream* msgs) {
  auto covar_args_refs = to_ref(std::forward<CovarTupleArgs>(covar_args));
  auto ll_args_refs = to_ref(std::forward<LLTupleArgs>(ll_args));
  // Solver 1, 2, 3
  constexpr bool ll_args_contain_var = is_any_var_scalar<LLTupleArgs>::value;
  auto partial_parm = make_zero(ll_args_refs);
  auto covar_args_adj = make_zero(covar_args_refs);
  double lmd = 0.0;
  {
    nested_rev_autodiff nested;
    // Solver 1, 2
    arena_t<Eigen::MatrixXd> R;
    // Solver 3
    arena_t<Eigen::MatrixXd> LU_solve_covariance;
    // Solver 1, 2, 3
    arena_t<std::decay_t<Theta>> s2(theta_0.size());
    // Make one hard copy here
    using laplace_likelihood::internal::conditional_copy_and_promote;
    using laplace_likelihood::internal::COPY_TYPE;
    auto ll_args_copy
        = conditional_copy_and_promote<is_any_var_scalar, var, COPY_TYPE::DEEP>(
            ll_args_refs);

    auto md_est = laplace_marginal_density_est(
        ll_fun, ll_args_copy, value_of(theta_0), covariance_function,
        value_of(covar_args_refs), options, msgs);
    // Return references to var types
    auto ll_args_filter = stan::math::filter_map<is_any_var_scalar>(
        [](auto&& arg) -> decltype(auto) {
          return std::forward<decltype(arg)>(arg);
        },
        ll_args_copy);
    stan::math::for_each(
        [](auto&& output_i, auto&& ll_arg_i) {
          if (is_any_var_scalar_v<decltype(ll_arg_i)>) {
            collect_adjoints(output_i, ll_arg_i);
            laplace_likelihood::internal::set_zero_adjoint(ll_arg_i);
          }
        },
        partial_parm, ll_args_filter);

    if (options.solver == 1) {
      // TODO(Steve): Solve without casting from sparse to dense
      Eigen::MatrixXd tmp
          = md_est.L.template triangularView<Eigen::Lower>().solve(
              md_est.W_r.toDense());
      R = tmp.transpose() * tmp;
      arena_t<Eigen::MatrixXd> C
          = md_est.L.template triangularView<Eigen::Lower>().solve(
              md_est.W_r * md_est.covariance);
      if (!ll_args_contain_var && options.hessian_block_size == 1) {
        auto s2_tmp
            = (0.5
               * (md_est.covariance.diagonal() - (C.transpose() * C).diagonal())
                     .cwiseProduct(laplace_likelihood::third_diff(
                         ll_fun, md_est.theta, value_of(ll_args_copy), msgs)))
                  .eval();
        s2.deep_copy(s2_tmp);
      } else {
        arena_t<Eigen::MatrixXd> A = md_est.covariance - C.transpose() * C;
        auto s2_tmp = laplace_likelihood::compute_s2(ll_fun, md_est.theta, A,
                                                     options.hessian_block_size,
                                                     ll_args_copy, msgs);
        s2.deep_copy(s2_tmp);
        copy_compute_s2(partial_parm, ll_args_filter);
        laplace_likelihood::internal::set_zero_adjoint(ll_args_filter);
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
      copy_compute_s2(partial_parm, ll_args_filter);
      laplace_likelihood::internal::set_zero_adjoint(ll_args_filter);
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
      copy_compute_s2(partial_parm, ll_args_filter);
      laplace_likelihood::internal::set_zero_adjoint(ll_args_filter);
    }
    lmd = md_est.lmd;
    if constexpr (is_any_var_scalar_v<scalar_type_t<CovarTupleArgs>>) {
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
        //      var Z = laplace_pseudo_target(K_var, md_est.a, R,
        //      md_est.theta_grad, s2);
        arena_t<Eigen::MatrixXd> K_adj_arena
            = 0.5 * md_est.a * md_est.a.transpose() - 0.5 * R
              + s2 * md_est.theta_grad.transpose()
              - (R * (K_var.val() * s2)) * md_est.theta_grad.transpose();
        var Z = make_callback_var(0.0, [K_var, K_adj_arena](auto&& vi) mutable {
          K_var.adj().array() += vi.adj() * K_adj_arena.array();
        });
        grad(Z.vi_);
        auto covar_args_filter = stan::math::filter_map<is_any_var_scalar>(
            [](auto&& arg) -> decltype(auto) {
              return std::forward<decltype(arg)>(arg);
            },
            covar_args_copy);
        //      std::cout << "\ncovar args: " << std::endl;
        //      print_adjoint(covar_args_filter);
        collect_adjoints(covar_args_adj, covar_args_filter);
        //      std::cout << "\n______________\n";
        //      std::cout << "covar args adj: " << std::endl;
        //      print_adjoint(covar_args_adj);
        //      std::cout << "\n==============\n";
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
      collect_adjoints(partial_parm, ll_args_filter);
      laplace_likelihood::internal::set_zero_adjoint(ll_args_filter);
    }
  }
  var ret(lmd);
  if constexpr (is_any_var_scalar_v<CovarTupleArgs>) {
    auto covar_args_arena = stan::math::filter_map<is_any_var_scalar>(
        [](auto&& arg) { return to_arena(arg); }, covar_args_refs);
    laplace_tuple_collect_adjoints(ret, covar_args_arena, covar_args_adj);
  }
  if constexpr (ll_args_contain_var) {
    auto ll_args_filter = stan::math::filter_map<is_any_var_scalar>(
        [](auto&& arg) { return to_arena(arg); }, ll_args_refs);
    laplace_tuple_collect_adjoints(ret, ll_args_filter, partial_parm);
  }
  return ret;
}

}  // namespace math
}  // namespace stan

#endif
