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
  int hessian_block_size;
  /**
   * Which Newton solver to use:
   * (1) method using the root of W
   * (2) method using the root of the covariance
   * (3) method using an LU decomposition
   */
  int solver;
  /* Maximum number of steps in line search*/
  int max_steps_line_search;
  /* iterations end when difference in objective function is less than tolerance
   */
  double tolerance;
  /* Maximum number of steps*/
  int64_t max_num_steps;
};
template <typename Covar, typename Theta, typename WR, typename L_t,
          typename A_vec, typename ThetaGrad, typename EtaGrad, typename LU_t,
          typename KRoot>
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
  // the gradient of the log density with respect to eta
  EtaGrad eta_grad;
  // LU matrix
  LU_t LU;
  // Cholesky of the covariance matrix
  KRoot K_root;
  laplace_density_estimates(double lmd_, Covar&& covariance_, Theta&& theta_,
                            WR&& W_r_, L_t&& L_, A_vec&& a_,
                            ThetaGrad&& theta_grad_, EtaGrad&& eta_grad_,
                            LU_t&& LU_, KRoot&& K_root_) noexcept
      : lmd(lmd_),
        covariance(std::move(covariance_)),
        theta(std::move(theta_)),
        W_r(std::move(W_r_)),
        L(std::move(L_)),
        a(std::move(a_)),
        theta_grad(std::move(theta_grad_)),
        eta_grad(std::move(eta_grad_)),
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
  return 0;
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
  Eigen::SparseMatrix<double> W_root(W.rows(), W.cols());
  W_root.reserve(Eigen::VectorXi::Constant(W_root.cols(), block_size));

  // No block operation available for sparse matrices, so we have to loop
  // See https://eigen.tuxfamily.org/dox/group__TutorialSparse.html#title7
  for (int i = 0; i < n_block; i++) {
    for (int k = 0; k < block_size; k++) {
      for (int j = 0; j < block_size; j++) {
        local_block.coeffRef(j, k)
            = W.coeff(i * block_size + j, i * block_size + k);
      }
    }
    local_block_sqrt = local_block.sqrt();
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

/**
 * For a latent Gaussian model with hyperparameters phi and eta,
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
 * @tparam LLFun Type with a valid `operator(Theta, Eta, InnerLLTupleArgs)`
 * where `InnerLLTupleArgs` are the elements of `LLTupleArgs`
 * @tparam LLTupleArgs A tuple whose elements follow the types required for
 * `LLFun`
 * @tparam CovarFun Type with a valid `operator(InnerCovarArgs)` where
 * `InnerCovarArgs` are a parameter pack of the element types of `CovarArgs`
 * @tparam Theta Type derived from `Eigen::EigenBase` with dynamic rows and a
 * single column
 * @tparam Eta Type derived from `Eigen::EigenBase` with dynamic rows and a
 * single column
 * @tparam CovarArgs A tuple whose elements follow the types required for
 * `CovarFun`
 * @param[in] ll_fun A log likelihood functor
 * @param[in] ll_args Tuple containing parameters for `LLFun`
 * @param[in] covariance_function Functor for the covariance function
 * @param[in] eta hyperparameter (input for likelihood)
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
          typename Theta, typename Eta, typename CovarArgs,
          require_t<is_all_arithmetic_scalar<Eta, Theta, CovarArgs>>* = nullptr,
          require_eigen_vector_t<Theta>* = nullptr>
inline auto laplace_marginal_density_est(LLFun&& ll_fun, LLTupleArgs&& ll_args,
                                         Eta&& eta, Theta&& theta_0,
                                         CovarFun&& covariance_function,
                                         CovarArgs&& covar_args,
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
  auto line_search = [](auto& objective_new, auto& a, auto& theta,
                        const auto& a_old, const auto& covariance,
                        const auto& ll_fun, auto&& ll_args, const auto& eta,
                        const auto max_steps_line_search,
                        const auto objective_old, auto* msgs) mutable {
    for (int j = 0;
         j < max_steps_line_search && (objective_new < objective_old); ++j) {
      a = (a + a_old) * 0.5;  // TODO(Charles) -- generalize for any factor
      theta = covariance * a;
      if (Eigen::isfinite(theta.array()).sum()) {
        objective_new = -0.5 * a.dot(theta)
                        + laplace_likelihood::log_likelihood(ll_fun, theta, eta,
                                                             ll_args, msgs);
      } else {
        break;
      }
    }
  };
  promote_scalar_t<double, Theta> theta_grad;
  promote_scalar_t<double, Eta> eta_grad;
  const Eigen::Index theta_size = theta_0.size();
  std::decay_t<Theta> theta = theta_0;
  double objective_old = -1e+10;  // CHECK -- what value to use?
  double objective_new = -1e+10;
  Eigen::VectorXd a_old;
  if (options.solver == 1 && options.hessian_block_size == 1) {
    SparseMatrix<double> W;
    for (Eigen::Index i = 0; i <= options.max_num_steps; i++) {
      std::tie(theta_grad, eta_grad, W) = laplace_likelihood::diff(
          ll_fun, theta, eta, options.hessian_block_size, ll_args, msgs);

      // Compute matrix square-root of W. If all elements of W are positive,
      // do an element wise square-root. Else try a matrix square-root
      bool W_is_spd = true;
      for (Eigen::Index i = 0; i < theta_0.size(); i++) {
        W_is_spd &= (W.coeff(i, i) < 0);
      }
      Eigen::SparseMatrix<double> W_r;
      if (W_is_spd) {
        W_r = W.cwiseSqrt();
      } else {
        W_r = block_matrix_sqrt(W, options.hessian_block_size);
      }
      auto B = MatrixXd::Identity(theta_size, theta_size)
               + W_r.diagonal().asDiagonal() * covariance
                     * W_r.diagonal().asDiagonal();
      Eigen::MatrixXd L = std::move(B)
                              .template selfadjointView<Eigen::Lower>()
                              .llt()
                              .matrixL();
      const double B_log_determinant = 2.0 * L.diagonal().array().log().sum();
      VectorXd b = W.diagonal().cwiseProduct(theta) + theta_grad;
      Eigen::VectorXd a
          = b
            - W_r
                  * mdivide_left_tri<Eigen::Upper>(
                      L.transpose(),
                      mdivide_left_tri<Eigen::Lower>(
                          L, W_r.diagonal().cwiseProduct(covariance * b)));

      // Simple Newton step
      theta = covariance * a;
      objective_old = objective_new;
      if (!(Eigen::isinf(theta.array()).any())) {
        objective_new = -0.5 * a.dot(theta)
                        + laplace_likelihood::log_likelihood(ll_fun, theta, eta,
                                                             ll_args, msgs);
      }
      if (options.max_steps_line_search && i != 0) {
        line_search(objective_new, a, theta, a_old, covariance, ll_fun, ll_args,
                    eta, options.max_steps_line_search, objective_old, msgs);
      }
      a_old = a;
      // Check for convergence
      if (abs(objective_new - objective_old) < options.tolerance) {
        return laplace_density_estimates{
            objective_new - 0.5 * B_log_determinant,
            std::move(covariance),
            std::move(theta),
            std::move(W_r),
            std::move(L),
            std::move(a),
            std::move(theta_grad),
            std::move(eta_grad),
            Eigen::PartialPivLU<Eigen::MatrixXd>{},
            Eigen::MatrixXd(0, 0)};
      }
    }
    throw_overstep(options.max_num_steps);
  } else if (options.solver == 1 && !(options.hessian_block_size == 1)) {
    SparseMatrix<double> W;
    for (Eigen::Index i = 0; i <= options.max_num_steps; i++) {
      std::tie(theta_grad, eta_grad, W) = laplace_likelihood::diff(
          ll_fun, theta, eta, options.hessian_block_size, ll_args, msgs);
      Eigen::SparseMatrix<double> W_r
          = block_matrix_sqrt(W, options.hessian_block_size);
      auto B = MatrixXd::Identity(theta_size, theta_size)
               + W_r * (covariance * W_r);
      Eigen::MatrixXd L = std::move(B)
                              .template selfadjointView<Eigen::Lower>()
                              .llt()
                              .matrixL();
      const double B_log_determinant = 2.0 * L.diagonal().array().log().sum();
      VectorXd b = W * theta + theta_grad;
      Eigen::VectorXd a
          = b
            - W_r
                  * mdivide_left_tri<Eigen::Upper>(
                      transpose(L), mdivide_left_tri<Eigen::Lower>(
                                        L, W_r * (covariance * b)));
      // Simple Newton step
      theta = covariance * a;
      objective_old = objective_new;
      if (std::isfinite(theta.sum())) {
        objective_new = -0.5 * a.dot(theta)
                        + laplace_likelihood::log_likelihood(ll_fun, theta, eta,
                                                             ll_args, msgs);
      }
      if (options.max_steps_line_search > 0 && i != 0) {
        line_search(objective_new, a, theta, a_old, covariance, ll_fun, ll_args,
                    eta, options.max_steps_line_search, objective_old, msgs);
      }
      a_old = a;
      // Check for convergence
      if (abs(objective_new - objective_old) < options.tolerance) {
        return laplace_density_estimates{
            objective_new - 0.5 * B_log_determinant,
            std::move(covariance),
            std::move(theta),
            std::move(W_r),
            std::move(L),
            std::move(a),
            std::move(theta_grad),
            std::move(eta_grad),
            Eigen::PartialPivLU<Eigen::MatrixXd>{},
            Eigen::MatrixXd(0, 0)};
      }
    }
    throw_overstep(options.max_num_steps);
  } else if (options.solver == 2) {
    SparseMatrix<double> W;
    for (Eigen::Index i = 0; i <= options.max_num_steps; i++) {
      std::tie(theta_grad, eta_grad, W) = laplace_likelihood::diff(
          ll_fun, theta, eta, options.hessian_block_size, ll_args, msgs);
      Eigen::MatrixXd K_root
          = covariance.template selfadjointView<Eigen::Lower>().llt().matrixL();
      auto B = MatrixXd::Identity(theta_size, theta_size)
               + K_root.transpose() * W * K_root;
      Eigen::MatrixXd L = std::move(B)
                              .template selfadjointView<Eigen::Lower>()
                              .llt()
                              .matrixL();
      const double B_log_determinant = 2.0 * L.diagonal().array().log().sum();
      VectorXd b = W * theta + theta_grad;
      Eigen::VectorXd a = mdivide_left_tri<Eigen::Upper>(
          K_root.transpose(),
          mdivide_left_tri<Eigen::Upper>(
              L.transpose(),
              mdivide_left_tri<Eigen::Lower>(L, K_root.transpose() * b)));
      // Simple Newton step
      theta = covariance * a;
      objective_old = objective_new;
      if (std::isfinite(theta.sum())) {
        objective_new = -0.5 * a.dot(theta)
                        + laplace_likelihood::log_likelihood(ll_fun, theta, eta,
                                                             ll_args, msgs);
      }
      // linesearch
      if (options.max_steps_line_search > 0 && i != 0) {
        line_search(objective_new, a, theta, a_old, covariance, ll_fun, ll_args,
                    eta, options.max_steps_line_search, objective_old, msgs);
      }
      a_old = a;
      // Check for convergence
      if (abs(objective_new - objective_old) < options.tolerance) {
        return laplace_density_estimates{
            objective_new - 0.5 * B_log_determinant,
            std::move(covariance),
            std::move(theta),
            std::move(W),
            std::move(L),
            std::move(a),
            std::move(theta_grad),
            std::move(eta_grad),
            Eigen::PartialPivLU<Eigen::MatrixXd>{},
            std::move(K_root)};
      }
    }
    throw_overstep(options.max_num_steps);
  } else if (options.solver == 3) {
    SparseMatrix<double> W;
    for (Eigen::Index i = 0; i <= options.max_num_steps; i++) {
      std::tie(theta_grad, eta_grad, W) = laplace_likelihood::diff(
          ll_fun, theta, eta, options.hessian_block_size, ll_args, msgs);
      auto B = MatrixXd::Identity(theta_size, theta_size) + covariance * W;
      Eigen::PartialPivLU<Eigen::MatrixXd> LU
          = Eigen::PartialPivLU<Eigen::MatrixXd>(std::move(B));
      // L on upper and U on lower triangular
      auto&& U = LU.matrixLU();
      // Compute log-determinant (Charles: Verify this is correct)
      double B_log_determinant = 0.0;
      int signDet = LU.permutationP().determinant();  // +1 or -1
      for (Eigen::Index i = 0; i < U.rows(); ++i) {
        B_log_determinant += std::log(std::abs(U.coeff(i, i)));
        signDet *= (U.coeff(i, i) >= 0) ? 1 : -1;
      }
      B_log_determinant *= signDet;
      //      const double B_log_determinant = log(LU.determinant());
      VectorXd b = W * theta + theta_grad;
      Eigen::VectorXd a = b - W * LU.solve(covariance * b);
      // Simple Newton step
      theta = covariance * a;
      objective_old = objective_new;

      if (std::isfinite(theta.sum())) {
        objective_new = -0.5 * a.dot(theta)
                        + laplace_likelihood::log_likelihood(ll_fun, theta, eta,
                                                             ll_args, msgs);
      }

      // linesearch
      // CHECK -- does linesearch work for options.solver 2?
      if (options.max_steps_line_search > 0 && i != 0) {
        line_search(objective_new, a, theta, a_old, covariance, ll_fun, ll_args,
                    eta, options.max_steps_line_search, objective_old, msgs);
      }
      a_old = a;
      // Check for convergence
      if (abs(objective_new - objective_old) < options.tolerance) {
        return laplace_density_estimates{
            objective_new - 0.5 * B_log_determinant,
            std::move(covariance),
            std::move(theta),
            std::move(W),
            Eigen::MatrixXd(0, 0),
            std::move(a),
            std::move(theta_grad),
            std::move(eta_grad),
            std::move(LU),
            Eigen::MatrixXd(0, 0)};
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
 * @tparam LLFun Type with a valid `operator(Theta, Eta, InnerLLTupleArgs)`
 * where `InnerLLTupleArgs` are the elements of `LLTupleArgs`
 * @tparam LLTupleArgs A tuple whose elements follow the types required for
 * `LLFun`
 * @tparam CovarFun Type with a valid `operator(InnerCovarArgs)` where
 * `InnerCovarArgs` are a parameter pack of the element types of `CovarArgs`
 * @tparam Theta Type derived from `Eigen::EigenBase` with dynamic rows and a
 * single column
 * @tparam Eta Type derived from `Eigen::EigenBase` with dynamic rows and a
 * single column
 * @tparam CovarArgs A tuple whose elements follow the types required for
 * `CovarFun`
 * @param[in] ll_fun A log likelihood functor
 * @param[in] ll_args Tuple containing parameters for `LLFun`
 * @param[in] covariance_function Functor for the covariance function
 * @param[in] eta hyperparameter (input for likelihood)
 * @param[in] theta_0 the initial guess for the Laplace optimization
 * @param[in,out] msgs stream for messages from likelihood and covariance
 * @param[in] options A set of options for tuning the solver
 * @param[in] covar_args Tuple of arguments to pass to the covariance matrix
 * functor
 * @return the log maginal density, p(y | phi)
 */
template <typename LLFun, typename LLArgs, typename CovarFun, typename Eta,
          typename Theta, typename CovarArgs, require_eigen_t<Eta>* = nullptr,
          require_arithmetic_t<return_type_t<Eta>>* = nullptr,
          require_t<is_all_arithmetic_scalar<CovarArgs>>* = nullptr,
          require_eigen_vector_t<Theta>* = nullptr>
inline double laplace_marginal_density(LLFun&& ll_fun, LLArgs&& ll_args,
                                       Eta&& eta, Theta&& theta_0,
                                       CovarFun&& covariance_function,
                                       CovarArgs&& covar_args,
                                       const laplace_options& options,
                                       std::ostream* msgs) {
  return laplace_marginal_density_est(
             std::forward<LLFun>(ll_fun), std::forward<LLArgs>(ll_args),
             std::forward<Eta>(eta), std::forward<Theta>(theta_0),
             std::forward<CovarFun>(covariance_function),
             std::forward<CovarArgs>(covar_args), options, msgs)
      .lmd;
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
 * @tparam LLFun Type with a valid `operator(Theta, Eta, InnerLLTupleArgs)`
 * where `InnerLLTupleArgs` are the elements of `LLTupleArgs`
 * @tparam LLTupleArgs A tuple whose elements follow the types required for
 * `LLFun`
 * @tparam CovarFun Type with a valid `operator(InnerCovarArgs)` where
 * `InnerCovarArgs` are a parameter pack of the element types of `CovarArgs`
 * @tparam Theta Type derived from `Eigen::EigenBase` with dynamic rows and a
 * single column
 * @tparam Eta Type derived from `Eigen::EigenBase` with dynamic rows and a
 * single column
 * @tparam CovarArgs A tuple whose elements follow the types required for
 * `CovarFun`
 * @param[in] ll_fun A log likelihood functor
 * @param[in] ll_args Tuple containing parameters for `LLFun`
 * @param[in] covariance_function Functor for the covariance function
 * @param[in] eta hyperparameter (input for likelihood)
 * @param[in] theta_0 the initial guess for the Laplace optimization
 * @param[in,out] msgs stream for messages from likelihood and covariance
 * @param[in] options A set of options for tuning the solver
 * @param[in] covar_args Tuple of arguments to pass to the covariance matrix
 * functor
 * @return the log maginal density, p(y | phi)
 */
template <typename LLFun, typename LLTupleArgs, typename CovarFun,
          typename Theta, typename Eta, typename CovarArgs,
          require_eigen_t<Eta>* = nullptr,
          require_t<is_any_var_scalar<Theta, Eta, CovarArgs>>* = nullptr,
          require_eigen_vector_t<Theta>* = nullptr>
inline auto laplace_marginal_density(const LLFun& ll_fun, LLTupleArgs&& ll_args,
                                     Eta&& eta, Theta&& theta_0,
                                     CovarFun&& covariance_function,
                                     CovarArgs&& covar_args,
                                     const laplace_options& options,
                                     std::ostream* msgs) {
  auto args_refs = to_ref(std::forward<CovarArgs>(covar_args));
  auto eta_arena = to_arena(eta);
  auto md_est = laplace_marginal_density_est(
      ll_fun, ll_args, value_of(eta_arena), value_of(theta_0),
      covariance_function, value_of(args_refs), options, msgs);
  const Eigen::Index theta_size = md_est.theta.size();
  const Eigen::Index eta_size = eta_arena.size();
  using Eta_t = std::decay_t<Eta>;
  // Solver 1, 2
  arena_t<Eigen::MatrixXd> R;
  // Solver 3
  arena_t<Eigen::MatrixXd> LU_solve_covariance;
  // Solver 1, 2, 3
  arena_t<promote_scalar_t<double, Eta_t>> partial_parm;
  // Solver 1, 2, 3
  arena_t<std::decay_t<Theta>> s2;
  if (options.solver == 1) {
    // TODO(Steve): Solve without casting from sparse to dense
    Eigen::MatrixXd tmp
        = md_est.L.template triangularView<Eigen::Lower>().solve(
            md_est.W_r.toDense());
    R = tmp.transpose() * tmp;
    arena_t<Eigen::MatrixXd> C = mdivide_left_tri<Eigen::Lower>(
        md_est.L, md_est.W_r * md_est.covariance);
    if (options.hessian_block_size == 1 && eta_size == 0) {
      s2 = 0.5
           * (md_est.covariance.diagonal() - (C.transpose() * C).diagonal())
                 .cwiseProduct(laplace_likelihood::third_diff(
                     ll_fun, md_est.theta, value_of(eta_arena), ll_args, msgs));
    } else {
      // int block_size = (hessian_block_size == 0) ? hessian_block_size + 1
      //                                            : hessian_block_size;
      arena_t<Eigen::MatrixXd> A = md_est.covariance - C.transpose() * C;
      std::tie(s2, partial_parm) = laplace_likelihood::compute_s2(
          ll_fun, md_est.theta, value_of(eta_arena), A,
          options.hessian_block_size, ll_args, msgs);
    }
  } else if (options.solver == 2) {
    // TODO(Charles) -- use triangularView for K_root
    R = md_est.W_r
        - md_est.W_r * md_est.K_root
              * md_est.L.transpose()
                    .template triangularView<Eigen::Upper>()
                    .solve(
                        md_est.L.template triangularView<Eigen::Lower>().solve(
                            md_est.K_root.transpose() * md_est.W_r));

    arena_t<Eigen::MatrixXd> C
        = md_est.L.template triangularView<Eigen::Lower>().solve(
            md_est.K_root.transpose());
    std::tie(s2, partial_parm) = laplace_likelihood::compute_s2(
        ll_fun, md_est.theta, value_of(eta_arena), C.transpose() * C,
        options.hessian_block_size, ll_args, msgs);
  } else {  // options.solver with LU decomposition
    LU_solve_covariance = md_est.LU.solve(md_est.covariance);
    R = md_est.W_r - md_est.W_r * LU_solve_covariance * md_est.W_r;

    arena_t<Eigen::MatrixXd> A
        = md_est.covariance
          - md_est.covariance * md_est.W_r * LU_solve_covariance;
    std::tie(s2, partial_parm) = laplace_likelihood::compute_s2(
        ll_fun, md_est.theta, value_of(eta_arena), A,
        options.hessian_block_size, ll_args, msgs);
  }
  auto args_arena = to_arena(args_refs);
  if constexpr (is_any_var_scalar_v<CovarArgs> && !is_constant<Eta>::value) {
    if (eta_size != 0) {
      {
        const nested_rev_autodiff nested;
        Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> K_var
            = stan::math::apply(
                [&covariance_function, &msgs](auto&&... args) {
                  return covariance_function(args..., msgs);
                },
                args_arena);
        var Z
            = laplace_pseudo_target(K_var, md_est.a, R, md_est.theta_grad, s2);
        set_zero_all_adjoints_nested();
        grad(Z.vi_);
      }
      auto arg_adj_arena = stan::math::filter_map<has_var_scalar_type>(
          [](auto&& arg) { return to_arena(get_adj(arg)); }, args_arena);
      stan::math::for_each([](auto&& arg) { zero_adjoints(arg); }, args_arena);

      arena_t<std::decay_t<decltype(md_est.eta_grad)>> diff_eta
          = md_est.eta_grad;

      arena_t<Eigen::VectorXd> v;
      if (options.solver == 1 || options.solver == 2) {
        v = md_est.covariance * s2
            - md_est.covariance * R * md_est.covariance * s2;
      } else {
        v = LU_solve_covariance * s2;
      }
      if constexpr (Eta_t::RowsAtCompileTime != 0
                    && Eta_t::ColsAtCompileTime != 0) {
        arena_matrix<promote_scalar_t<double, Eta_t>> eta_adj_arena
            = md_est.eta_grad + partial_parm
              + laplace_likelihood::diff_eta_implicit(
                  ll_fun, v, md_est.theta, value_of(eta_arena), ll_args, msgs);
        return make_callback_var(
            md_est.lmd, [arg_adj_arena, args_arena, eta_arena,
                         eta_adj_arena](const auto& vi) mutable {
              stan::math::for_each(
                  [&vi](auto&& arg, auto&& arg_adj) {
                    if constexpr (is_var<scalar_type_t<decltype(arg)>>::value) {
                      internal::update_adjoints(arg, arg_adj, vi);
                    }
                  },
                  args_arena, arg_adj_arena);
              internal::update_adjoints(eta_arena, eta_adj_arena, vi);
            });
      }
    }
  } else if constexpr (is_any_var_scalar_v<scalar_type_t<CovarArgs>>) {
    {
      const nested_rev_autodiff nested;
      arena_t<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> K_var
          = stan::math::apply(
              [&covariance_function, &msgs](auto&&... args) {
                return covariance_function(args..., msgs);
              },
              args_arena);
      //  = covariance_function(x, phi_v, delta, delta_int, msgs);
      var Z = laplace_pseudo_target(K_var, md_est.a, R, md_est.theta_grad, s2);
      set_zero_all_adjoints_nested();
      grad(Z.vi_);
    }
    auto arg_adj_arena = stan::math::filter_map<has_var_scalar_type>(
        [](auto&& arg) { return to_arena(get_adj(arg)); }, args_arena);
    stan::math::for_each([](auto&& arg) { zero_adjoints(arg); }, args_arena);
    return make_callback_var(
        md_est.lmd, [arg_adj_arena, args_arena](const auto& vi) mutable {
          stan::math::for_each(
              [&vi](auto&& arg, auto&& arg_adj) {
                if constexpr (is_var<scalar_type_t<decltype(arg)>>::value) {
                  internal::update_adjoints(arg, arg_adj, vi);
                }
              },
              args_arena, arg_adj_arena);
        });
  } else if (!is_constant<Eta>::value && eta_size != 0) {
    if constexpr (Eta_t::RowsAtCompileTime != 0
                  && Eta_t::ColsAtCompileTime != 0) {
      arena_t<Eigen::VectorXd> diff_eta = md_est.eta_grad;

      arena_t<Eigen::VectorXd> v;
      if (options.solver == 1 || options.solver == 2) {
        v = md_est.covariance * s2
            - md_est.covariance * R * md_est.covariance * s2;
      } else {
        v = LU_solve_covariance * s2;
      }

      arena_matrix<Eigen::VectorXd> eta_adj_arena
          = md_est.eta_grad + partial_parm
            + laplace_likelihood::diff_eta_implicit(
                ll_fun, v, md_est.theta, value_of(eta_arena), ll_args, msgs);

      return make_callback_var(
          md_est.lmd, [eta_arena, eta_adj_arena](const auto& vi) mutable {
            internal::update_adjoints(eta_arena, eta_adj_arena, vi);
          });
    }
  }
  return var(0);
}

}  // namespace math
}  // namespace stan

#endif
