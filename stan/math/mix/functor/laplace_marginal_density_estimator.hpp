#ifndef STAN_MATH_MIX_FUNCTOR_LAPLACE_MARGINAL_DENSITY_ESTIMATOR_HPP
#define STAN_MATH_MIX_FUNCTOR_LAPLACE_MARGINAL_DENSITY_ESTIMATOR_HPP
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/generate_laplace_options.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/functor/wolfe_line_search.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun.hpp>
#include <stan/math/rev/functor.hpp>
#include <stan/math/mix/functor/barzilai_borwein_step_size.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/quad_form_diag.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/iter_tuple_nested.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
#include <mutex>

/**
 * @file
 * Reference for calculations of marginal and its gradients:
 * Margossian et al (2020), https://arxiv.org/abs/2004.12550
 * and Margossian (2023), https://arxiv.org/pdf/2306.14976
 */

namespace stan {
namespace math {

/**
 * Options for the Laplace approximation.
 */
struct laplace_options_base {
  /* Size of the blocks in block diagonal hessian*/
  int hessian_block_size{internal::laplace_default_hessian_block_size};  // 0
  /**
   * Which linear solver to use inside the Newton step.
   *
   * This selects how the system matrix `B` is formed and factorized.
   * For details, see equation 1 in:
   *   https://arxiv.org/pdf/2306.14976
   *
   * 1. Hessian-root Cholesky: form `B = I + W_r * Sigma * W_r` where
   *    `W_r = sqrt(W)` and `W` is the negative Hessian of the log-likelihood
   *    w.r.t. `theta`. This uses either a diagonal (`hessian_block_size == 1`)
   *    or block-diagonal (`hessian_block_size > 1`) approximation of `W`.
   * 2. Covariance-root Cholesky: precompute `K_root` such that
   *    `Sigma = K_root * K_root^T` and form `B = I + K_root^T * W * K_root`.
   * 3. General LU: form `B = I + Sigma * W` and factorize with LU.
   */
  int solver{internal::laplace_default_solver};  // 1
  /**
   * Iterations end when the absolute change in the optimization objective
   * is less than this tolerance.
   *
   * Note: the objective used for convergence is the one optimized by the
   * Newton/Wolfe loop (not the final Laplace-corrected log marginal density).
   */
  double tolerance{internal::laplace_default_tolerance};  // 2
  /* Maximum number of steps*/
  int max_num_steps{internal::laplace_default_max_num_steps};          // 3
  int allow_fallthrough{internal::laplace_default_allow_fallthrough};  // 4
  laplace_line_search_options line_search{
      internal::laplace_default_max_steps_line_search};  // 5
  laplace_options_base() = default;
  laplace_options_base(int hessian_block_size_, int solver_, double tolerance_,
                       int max_num_steps_, bool allow_fallthrough_,
                       int max_steps_line_search_)
      : hessian_block_size(hessian_block_size_),
        solver(solver_),
        tolerance(tolerance_),
        max_num_steps(max_num_steps_),
        allow_fallthrough(allow_fallthrough_),
        line_search(max_steps_line_search_) {}
};

template <bool HasInitTheta>
struct laplace_options;

template <>
struct laplace_options<false> : public laplace_options_base {
  laplace_options() = default;

  explicit laplace_options(int hessian_block_size_) {
    hessian_block_size = hessian_block_size_;
  }
};

template <>
struct laplace_options<true> : public laplace_options_base {
  /* Value for user supplied initial theta  */
  Eigen::VectorXd theta_0{0};  // 6

  template <typename ThetaVec>
  laplace_options(ThetaVec&& theta_0_, double tolerance_, int max_num_steps_,
                  int hessian_block_size_, int solver_,
                  int max_steps_line_search_, bool allow_fallthrough_)
      : laplace_options_base(hessian_block_size_, solver_, tolerance_,
                             max_num_steps_, allow_fallthrough_,
                             max_steps_line_search_),
        theta_0(value_of(std::forward<ThetaVec>(theta_0_))) {}
};

using laplace_options_default = laplace_options<false>;
using laplace_options_user_supplied = laplace_options<true>;

namespace internal {

template <typename Options>
inline constexpr auto tuple_to_laplace_options(Options&& ops) {
  using Ops = std::decay_t<Options>;
  if constexpr (is_tuple_v<Ops>) {
    if constexpr (!is_eigen_v<std::tuple_element_t<0, std::decay_t<Ops>>>) {
      static_assert(
          sizeof(std::decay_t<Ops>*) == 0,
          "ERROR:(laplace_marginal_lpdf) The first laplace argument is "
          "expected to be an Eigen vector of dynamic size representing the "
          "initial theta_0.");
    }
    if constexpr (!stan::is_inner_tuple_type_v<1, Ops, double>) {
      static_assert(
          sizeof(std::decay_t<Ops>*) == 0,
          "ERROR:(laplace_marginal_lpdf) The second laplace argument is "
          "expected to be a double representing the tolerance.");
    }
    if constexpr (!stan::is_inner_tuple_type_v<2, Ops, int>) {
      static_assert(
          sizeof(std::decay_t<Ops>*) == 0,
          "ERROR:(laplace_marginal_lpdf) The third laplace argument is "
          "expected to be an int representing the maximum number of steps for "
          "the laplace approximation.");
    }
    if constexpr (!stan::is_inner_tuple_type_v<3, Ops, int>) {
      static_assert(
          sizeof(std::decay_t<Ops>*) == 0,
          "ERROR:(laplace_marginal_lpdf) The fourth laplace argument is "
          "expected to be an int representing the solver.");
    }
    if constexpr (!stan::is_inner_tuple_type_v<4, Ops, int>) {
      static_assert(
          sizeof(std::decay_t<Ops>*) == 0,
          "ERROR:(laplace_marginal_lpdf) The fifth laplace argument is "
          "expected to be an int representing the max steps for the laplace "
          "approximaton's wolfe line search.");
    }
    constexpr bool is_fallthrough
        = stan::is_inner_tuple_type_v<
              5, Ops, int> || stan::is_inner_tuple_type_v<5, Ops, bool>;
    if constexpr (!is_fallthrough) {
      static_assert(
          sizeof(std::decay_t<Ops>*) == 0,
          "ERROR:(laplace_marginal_lpdf) The sixth laplace argument is "
          "expected to be an int representing allow fallthrough (0/1).");
    }
    auto defaults = laplace_options_default{};
    return laplace_options_user_supplied{
        value_of(std::get<0>(std::forward<Ops>(ops))),
        std::get<1>(ops),
        std::get<2>(ops),
        defaults.hessian_block_size,
        std::get<3>(ops),
        std::get<4>(ops),
        (std::get<5>(ops) > 0) ? true : false,
    };
  } else {
    return std::forward<Ops>(ops);
  }
}

template <typename ThetaVec, typename WR, typename L_t, typename A_vec,
          typename ThetaGrad, typename LU_t, typename KRoot>
struct laplace_density_estimates {
  /* log marginal density */
  double lmd{std::numeric_limits<double>::infinity()};
  /* ThetaVec at the mode */
  ThetaVec theta;
  /**
   * Solver-dependent Hessian quantity.
   *
   * - solver 1: sparse square root `W_r = sqrt(W)` where
   *   `W = -d^2/dtheta^2 log_likelihood(theta)` (diagonal or block-diagonal).
   * - solvers 2/3: sparse `W` itself.
   */
  WR W_r;
  /**
   * Solver-dependent factorization of the system matrix `B`.
   *
   * - solvers 1/2: lower Cholesky factor `L` such that `B = L * L^T`.
   * - solver 3: empty (0x0).
   */
  L_t L;
  /**
   * Mode in the `a` parameterization, where `theta = covariance * a`.
   * Equivalently, `a = covariance^{-1} * theta` when the inverse exists.
   */
  A_vec a;
  /** Gradient of the log-likelihood with respect to `theta` at the mode. */
  ThetaGrad theta_grad;
  /* LU matrix from solver 3 */
  LU_t LU;
  /**
   * Lower Cholesky factor of the covariance matrix.
   *
   * - solver 2: `covariance = K_root * K_root^T`.
   * - other solvers: empty (0x0).
   */
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
    if (!local_block.array().isFinite().any()) {
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
 * Validates the options for the Laplace approximation.
 *
 * @tparam InitTheta Whether an initial theta is provided
 * @tparam CovarMat Type of the covariance matrix
 * @param frame_name Name of the calling function (for error messages)
 * @param options The options to validate
 * @param covariance The covariance matrix (for size checks)
 */
template <bool InitTheta, typename CovarMat>
inline void validate_laplace_options(const char* frame_name,
                                     const laplace_options<InitTheta>& options,
                                     const CovarMat& covariance) {
  if constexpr (InitTheta) {
    check_nonzero_size(frame_name, "initial guess", options.theta_0);
    check_finite(frame_name, "initial guess", options.theta_0);
    if (unlikely(options.theta_0.size() != covariance.rows())) {
      std::stringstream msg;
      msg << frame_name << ": The size of the initial theta ("
          << options.theta_0.size()
          << ") vector must match the rows and columns of the covariance "
             "matrix ("
          << covariance.rows() << ", " << covariance.cols() << ").";
      throw std::domain_error(msg.str());
    }
  }
  check_nonnegative(frame_name, "tolerance", options.tolerance);
  check_positive(frame_name, "max_num_steps", options.max_num_steps);
  check_positive(frame_name, "hessian_block_size", options.hessian_block_size);
  check_square(frame_name, "covariance", covariance);

  const Eigen::Index theta_size = covariance.rows();
  if (unlikely(theta_size % options.hessian_block_size != 0
               || theta_size < options.hessian_block_size)) {
    throw std::domain_error(
        "laplace_marginal_density: Hessian block size mismatch.");
  }

  if (unlikely(options.solver < 1 || options.solver > 3)) {
    throw std::domain_error(
        "laplace_marginal_density: solver must be 1, 2, or 3. Got: "
        + std::to_string(options.solver));
  }
}

/**
 * @brief Holds the state for the Newton-Raphson optimization loop.
 *
 * This struct centralizes all state needed during the Newton iteration,
 * including the Wolfe line search state, workspace matrices, and convergence
 * flags. It is shared across different solver policies to avoid re-allocation
 * and to maintain progress (e.g., line search history) across solver fallbacks.
 *
 * @tparam ObjFun Type of the objective function callable
 * @tparam ThetaGradFun Type of the theta gradient function callable
 */
struct NewtonState {
  /** @brief Wolfe line search state including current/previous steps */
  WolfeInfo wolfe_info;

  /** @brief Status of the most recent Wolfe line search */
  WolfeStatus wolfe_status;

  /** @brief Cached proposal evaluated before the Wolfe line search. */
  WolfeData proposal;

  /** @brief Workspace vector: b = W * theta + grad(log_lik) */
  Eigen::VectorXd b;

  /** @brief Workspace matrix: B = I + W_r * Sigma * W_r (or similar) */
  Eigen::MatrixXd B;

  /** @brief Previous gradient for Barzilai-Borwein step calculation */
  Eigen::VectorXd prev_g;
  /**
   * On the final loop if we found a better wolfe step, but we are going to
   * exit, we want to make sure all of our return values are with the most
   * recent wolfe step that was accepted. So we do one final loop to update
   * our return values. This flag tells us if we are on the final loop and
   * need to update the values one more time.
   */
  bool final_loop = false;

  /**
   * @brief Constructs Newton state with a consistent (a_init, theta_init) pair.
   *
   * When the caller supplies a non-zero theta_init, a_init = Sigma^{-1} *
   * theta_init must be provided to maintain the invariant theta = Sigma * a.
   * @tparam ObjFun A callable type for the objective function
   * @tparam ThetaGradFun A callable type for the theta gradient function
   * @tparam CovarianceT A matrix type for the covariance (must support LLT solve)
   * @tparam ThetaInitializer A type for the initial theta (e.g., Eigen vector)
   * @param theta_size Dimension of the latent space
   * @param obj_fun Objective function: (a, theta) -> double
   * @param theta_grad_f Gradient function: theta -> grad
   * @param covariance Covariance matrix for the latent variables
   * @param a_init Initial a value consistent with theta_init
   * @param theta_init Initial theta value
   */
  template <typename ObjFun, typename ThetaGradFun, typename CovarianceT,
            typename ThetaInitializer>
  NewtonState(int theta_size, ObjFun&& obj_fun, ThetaGradFun&& theta_grad_f,
              CovarianceT&& covariance, ThetaInitializer&& theta_init)
      : wolfe_info(std::forward<ObjFun>(obj_fun),
                   covariance.llt().solve(theta_init),
                   std::forward<ThetaInitializer>(theta_init),
                   std::forward<ThetaGradFun>(theta_grad_f)),
        proposal(theta_size),
        b(theta_size),
        B(theta_size, theta_size),
        prev_g(theta_size) {
    wolfe_status.num_backtracks_ = -1;  // Safe initial value for BB step
  }

  /**
   * @brief Access the current step state (mutable).
   * @return Reference to current WolfeStep
   */
  auto& curr() & { return wolfe_info.curr_; }

  /**
   * @brief Access the current step state (const).
   * @return Const reference to current WolfeStep
   */
  const auto& curr() const& { return wolfe_info.curr_; }
  auto&& curr() && { return std::move(wolfe_info).curr(); }
  /**
   * @brief Access the previous step state (mutable).
   * @return Reference to previous WolfeStep
   */
  auto& prev() & { return wolfe_info.prev_; }

  /**
   * @brief Access the previous step state (const).
   * @return Const reference to previous WolfeStep
   */
  const auto& prev() const& { return wolfe_info.prev_; }
  auto&& prev() && { return std::move(wolfe_info).prev(); }
  auto& proposal_step() & { return proposal; }
  const auto& proposal_step() const& { return proposal; }
  auto&& proposal_step() && { return std::move(proposal); }
  template <typename Options>
  inline void update_next_step(const Options& options) {
    this->prev().swap(this->curr());
    this->curr().alpha()
        = std::clamp(this->curr().alpha(), 0.0, options.line_search.max_alpha);
  }
};

/**
 * @brief Factorize B with jittering fallback
 * @tparam LLT Cholesky factorization type
 * @tparam B_t Matrix type
 * @param llt_B Cholesky factorization object
 * @param B Matrix to factorize
 * @param min_jitter Minimum jitter to try
 * @param max_jitter Maximum jitter to try
 */
template <typename LLT, typename B_t>
inline void llt_with_jitter(LLT& llt_B, B_t& B, double min_jitter = 1e-10,
                            double max_jitter = 1e-5) {
  llt_B.compute(B);
  if (llt_B.info() != Eigen::Success) {
    double prev_jitter = 0.0;
    double jitter_try = min_jitter;
    for (; jitter_try < max_jitter; jitter_try *= 10) {
      // Remove previously added jitter before adding the new (larger) amount,
      // so that the total diagonal perturbation is exactly jitter_try.
      B.diagonal().array() += (jitter_try - prev_jitter);
      prev_jitter = jitter_try;
      llt_B.compute(B);
      if (llt_B.info() == Eigen::Success) {
        break;
      }
    }
    if (llt_B.info() != Eigen::Success) {
      throw std::domain_error(
          "laplace_marginal_density: Cholesky failed after adding jitter up to "
          + std::to_string(jitter_try));
    }
  }
}

/**
 * @brief Solver Policy 1 (Diagonal): Cholesky decomposition using W.
 *
 * This solver is used when `hessian_block_size == 1`. It computes the
 * diagonal of the negative Hessian of the log-likelihood, takes its
 * square root, and forms the system matrix B = I + W_r * Sigma * W_r
 * for Cholesky factorization.
 *
 * The solver is the fastest option but only valid when the Hessian
 * is truly diagonal (no cross-terms between latent variables).
 *
 * @note This solver corresponds to `solver == 1` in the original code.
 */
struct CholeskyWSolverDiag {
  /** @brief Square root of diagonal Hessian: W_r[j] = sqrt(W[j]) */
  Eigen::VectorXd W_r_diag;

  /** @brief Diagonal Hessian values from the likelihood */
  Eigen::VectorXd W_diag;

  /** @brief Cholesky factorization of B = I + W_r * Sigma * W_r */
  Eigen::LLT<Eigen::MatrixXd> llt_B;

  template <typename NewtonStateT, typename CovarMat>
  CholeskyWSolverDiag(const NewtonStateT& state, const CovarMat& covariance)
      : W_r_diag(Eigen::VectorXd::Zero(state.b.size())), W_diag(0), llt_B() {}
  /**
   * @brief Perform one Newton step using diagonal Hessian solver.
   *
   * Computes the diagonal Hessian, forms B = I + W_r * Sigma * W_r,
   * performs Cholesky factorization, and solves for the new `a` vector.
   *
   * @tparam NewtonStateT Type of the Newton state
   * @tparam LLFun Type of the log-likelihood functor
   * @tparam LLTupleArgs Type of the likelihood arguments tuple
   * @tparam CovarMat Type of the covariance matrix
   * @param[in,out] state Shared Newton state (modified: B, b,
   * proposal_step().a())
   * @param[in] ll_fun Log-likelihood functor
   * @param[in,out] ll_args Additional arguments for the likelihood
   * @param[in] covariance Prior covariance matrix Sigma
   * @param[in] hessian_block_size Ignored for diagonal solver
   * @param[in,out] msgs Output stream for diagnostic messages (may be nullptr)
   * @throws std::domain_error If Hessian is not positive definite or Cholesky
   * fails
   */
  template <typename NewtonStateT, typename LLFun, typename LLTupleArgs,
            typename CovarMat>
  void solve_step(NewtonStateT& state, const LLFun& ll_fun,
                  const LLTupleArgs& ll_args, const CovarMat& covariance,
                  int /*hessian_block_size*/, std::ostream* msgs) {
    const Eigen::Index theta_size = state.b.size();

    // 1. Compute diagonal Hessian
    W_diag = laplace_likelihood::diagonal_hessian(ll_fun, state.prev().theta(),
                                                  ll_args, msgs);
    for (Eigen::Index j = 0; j < W_diag.size(); j++) {
      if (W_diag.coeff(j) < 0 || !std::isfinite(W_diag.coeff(j))) {
        throw std::domain_error(
            "laplace_marginal_density: Hessian matrix is not positive "
            "definite");
      } else {
        W_r_diag.coeffRef(j) = std::sqrt(W_diag.coeff(j));
      }
    }

    // 2. Formulate B = I + W_r * Sigma * W_r
    state.B.noalias()
        = Eigen::MatrixXd::Identity(theta_size, theta_size)
          + W_r_diag.asDiagonal() * covariance * W_r_diag.asDiagonal();

    // 3. Factorize B with jittering fallback
    llt_with_jitter(llt_B, state.B);
    // 4. Solve for the raw Newton proposal in a-space.
    state.b.noalias() = (W_diag.array() * state.prev().theta().array()).matrix()
                        + state.prev().theta_grad();
    auto L = llt_B.matrixL();
    auto LT = llt_B.matrixU();
    state.proposal_step().a().noalias()
        = state.b
          - W_r_diag.asDiagonal()
                * LT.solve(
                    L.solve(W_r_diag.cwiseProduct(covariance * state.b)));
  }

  /**
   * @brief Compute log determinant of B from Cholesky factor.
   * @return log(det(B)) = 2 * sum(log(diag(L)))
   */
  double compute_log_determinant() const {
    return 2.0 * llt_B.matrixLLT().diagonal().array().log().sum();
  }

  /**
   * @brief Build the final result structure.
   *
   * @tparam NewtonStateT Type of the Newton state
   * @param[in] state Newton state containing converged values
   * @param[in] log_det Log determinant from compute_log_determinant()
   * @return laplace_density_estimates with solver_used = 1
   */
  template <typename NewtonStateT>
  auto build_result(NewtonStateT& state, double log_det) {
    return laplace_density_estimates{
        state.prev().obj() - 0.5 * log_det,
        std::move(state).prev().theta(),
        Eigen::SparseMatrix<double>(W_r_diag.asDiagonal()),
        Eigen::MatrixXd(llt_B.matrixL()),
        std::move(state).prev().a(),
        std::move(state).prev().theta_grad(),
        Eigen::PartialPivLU<Eigen::MatrixXd>{},
        Eigen::MatrixXd(0, 0),
        1};
  }
};

/**
 * @brief Solver Policy 1 (Block): Cholesky decomposition using block W.
 *
 * This solver is used when `hessian_block_size > 1`. It computes
 * the block-diagonal negative Hessian of the log-likelihood, computes
 * its principal square root via `block_matrix_sqrt`, and forms the
 * system matrix B = I + W_r * Sigma * W_r for Cholesky factorization.
 *
 * The sparse structure of `W_r` is initialized in the constructor to match the
 * problem dimensions.
 *
 * @note This solver corresponds to `solver == 1` with `hessian_block_size > 1`.
 */
struct CholeskyWSolverBlock {
  /** @brief Sparse square root of block Hessian */
  Eigen::SparseMatrix<double> W_r;

  /** @brief Sparse block-diagonal Hessian from likelihood */
  Eigen::SparseMatrix<double> W_block;

  /** @brief Cholesky factorization of B = I + W_r * Sigma * W_r */
  Eigen::LLT<Eigen::MatrixXd> llt_B;

  template <typename NewtonStateT>
  CholeskyWSolverBlock(const NewtonStateT& state, int hessian_block_size)
      : W_r(state.b.size(), state.b.size()) {
    const Eigen::Index theta_size = state.b.size();
    W_r.reserve(Eigen::VectorXi::Constant(theta_size, hessian_block_size));
    const Eigen::Index n_block = theta_size / hessian_block_size;
    for (Eigen::Index ii = 0; ii < n_block; ii++) {
      for (Eigen::Index k = 0; k < hessian_block_size; k++) {
        for (Eigen::Index j = 0; j < hessian_block_size; j++) {
          W_r.insert(ii * hessian_block_size + j, ii * hessian_block_size + k)
              = 1.0;
        }
      }
    }
    W_r.makeCompressed();
  }

  /**
   * @brief Perform one Newton step using block-diagonal Hessian solver.
   *
   * Computes the block Hessian, its square root via Schur decomposition,
   * forms B = I + W_r * Sigma * W_r, performs Cholesky factorization,
   * and solves for the new `a` vector.
   *
   * @tparam NewtonStateT Type of the Newton state
   * @tparam LLFun Type of the log-likelihood functor
   * @tparam LLTupleArgs Type of the likelihood arguments tuple
   * @tparam CovarMat Type of the covariance matrix
   * @param[in,out] state Shared Newton state (modified: B, b,
   * proposal_step().a())
   * @param[in] ll_fun Log-likelihood functor
   * @param[in,out] ll_args Additional arguments for the likelihood
   * @param[in] covariance Prior covariance matrix Sigma
   * @param[in] hessian_block_size Size of each Hessian block (must divide
   * theta_size)
   * @param[in,out] msgs Output stream for diagnostic messages (may be nullptr)
   * @throws std::domain_error If Hessian is not positive definite or Cholesky
   * fails
   */
  template <typename NewtonStateT, typename LLFun, typename LLTupleArgs,
            typename CovarMat>
  void solve_step(NewtonStateT& state, const LLFun& ll_fun,
                  const LLTupleArgs& ll_args, const CovarMat& covariance,
                  int hessian_block_size, std::ostream* msgs) {
    const Eigen::Index theta_size = state.b.size();
    // 1. Compute block Hessian
    W_block = laplace_likelihood::block_hessian(
        ll_fun, state.prev().theta(), hessian_block_size, ll_args, msgs);

    for (Eigen::Index j = 0; j < W_block.rows(); j++) {
      if (W_block.coeff(j, j) < 0 || !std::isfinite(W_block.coeff(j, j))) {
        throw std::domain_error(
            "laplace_marginal_density: Hessian matrix is not positive "
            "definite");
      }
    }

    // 2. Compute W_r = sqrt(W)
    block_matrix_sqrt(W_r, W_block, hessian_block_size);

    // 3. Formulate B = I + W_r * Sigma * W_r
    state.B.noalias() = Eigen::MatrixXd::Identity(theta_size, theta_size)
                        + W_r * (covariance * W_r);

    // 4. Factorize B with jittering fallback
    llt_with_jitter(llt_B, state.B);

    // 5. Solve for the raw Newton proposal in a-space.
    state.b.noalias()
        = W_block * state.prev().theta() + state.prev().theta_grad();
    auto L = llt_B.matrixL();
    auto LT = llt_B.matrixU();
    state.proposal_step().a().noalias()
        = state.b - W_r * LT.solve(L.solve(W_r * (covariance * state.b)));
  }

  /**
   * @brief Compute log determinant of B from Cholesky factor.
   * @return log(det(B)) = 2 * sum(log(diag(L)))
   */
  double compute_log_determinant() const {
    return 2.0 * llt_B.matrixLLT().diagonal().array().log().sum();
  }

  /**
   * @brief Build the final result structure.
   *
   * @tparam NewtonStateT Type of the Newton state
   * @param[in] state Newton state containing converged values
   * @param[in] log_det Log determinant from compute_log_determinant()
   * @return laplace_density_estimates with solver_used = 1
   */
  template <typename NewtonStateT>
  auto build_result(NewtonStateT& state, double log_det) {
    return laplace_density_estimates{state.prev().obj() - 0.5 * log_det,
                                     std::move(state).prev().theta(),
                                     std::move(W_r),
                                     Eigen::MatrixXd(llt_B.matrixL()),
                                     std::move(state).prev().a(),
                                     std::move(state).prev().theta_grad(),
                                     Eigen::PartialPivLU<Eigen::MatrixXd>{},
                                     Eigen::MatrixXd(0, 0),
                                     1};
  }
};

/**
 * @brief Solver Policy 2: Cholesky decomposition of K (Covariance).
 *
 * This solver pre-computes the Cholesky factorization of the prior
 * covariance matrix Sigma = K_root * K_root^T. The system matrix becomes
 * B = I + K_root^T * W * K_root, which is factorized in each iteration.
 *
 * This approach is numerically more stable than Solver 1 when the
 * covariance matrix has a more complex structure.
 *
 * @note This solver corresponds to `solver == 2` in the original code.
 */
struct CholeskyKSolver {
  /** @brief Lower Cholesky factor of covariance: Sigma = K_root * K_root^T */
  Eigen::MatrixXd K_root;

  /** @brief Full (block) Hessian matrix from likelihood */
  Eigen::SparseMatrix<double> W_full;

  /** @brief Cholesky factorization of B = I + K_root^T * W * K_root */
  Eigen::LLT<Eigen::MatrixXd> llt_B;

  template <typename NewtonStateT, typename CovarMat>
  CholeskyKSolver(const NewtonStateT& state, const CovarMat& covariance)
      : K_root(0, 0), W_full(0, 0), llt_B() {
    auto K_root_llt = covariance.template selfadjointView<Eigen::Lower>().llt();
    if (K_root_llt.info() != Eigen::Success) {
      throw std::domain_error(
          "laplace_marginal_density: Cholesky of covariance failed at start");
    }
    K_root = std::move(K_root_llt.matrixL());
  }

  /**
   * @brief Perform one Newton step using covariance Cholesky solver.
   *
   * Computes the block diagonal Hessian, forms B = I + K_root^T * W * K_root,
   * performs Cholesky factorization, and solves for the new `a` vector
   * using triangular solves.
   *
   * @tparam NewtonStateT Type of the Newton state
   * @tparam LLFun Type of the log-likelihood functor
   * @tparam LLTupleArgs Type of the likelihood arguments tuple
   * @tparam CovarMat Type of the covariance matrix
   * @param[in] state Shared Newton state (modified: B, b, proposal_step().a())
   * @param[in] ll_fun Log-likelihood functor
   * @param[in] ll_args Additional arguments for the likelihood
   * @param[in] covariance Prior covariance matrix Sigma
   * @param[in] hessian_block_size Size of each Hessian block
   * @param[in,out] msgs Output stream for diagnostic messages (may be nullptr)
   * @throws std::domain_error If Cholesky factorization fails
   */

  template <typename NewtonStateT, typename LLFun, typename LLTupleArgs,
            typename CovarMat>
  void solve_step(NewtonStateT& state, const LLFun& ll_fun,
                  const LLTupleArgs& ll_args, const CovarMat& covariance,
                  int hessian_block_size, std::ostream* msgs) {
    const Eigen::Index theta_size = state.b.size();

    // 1. Compute Hessian
    W_full = laplace_likelihood::block_hessian(
        ll_fun, state.prev().theta(), hessian_block_size, ll_args, msgs);

    // 2. Formulate B = I + K^T * W * K
    state.B.noalias() = Eigen::MatrixXd::Identity(theta_size, theta_size)
                        + K_root.transpose() * (W_full * K_root);

    // 3. Factorize B with jittering fallback
    llt_with_jitter(llt_B, state.B);

    // 4. Solve for the raw Newton proposal in a-space.
    state.b.noalias()
        = W_full * state.prev().theta() + state.prev().theta_grad();
    auto L = llt_B.matrixL();
    auto LT = llt_B.matrixU();
    state.proposal_step().a().noalias()
        = K_root.transpose().template triangularView<Eigen::Upper>().solve(
            LT.solve(L.solve(K_root.transpose() * state.b)));
  }

  /**
   * @brief Compute log determinant of B from Cholesky factor.
   * @return log(det(B)) = 2 * sum(log(diag(L)))
   */
  double compute_log_determinant() const {
    return 2.0 * llt_B.matrixLLT().diagonal().array().log().sum();
  }

  /**
   * @brief Build the final result structure.
   *
   * @tparam NewtonStateT Type of the Newton state
   * @param state Newton state containing converged values
   * @param log_det Log determinant from compute_log_determinant()
   * @return laplace_density_estimates with solver_used = 2
   */
  template <typename NewtonStateT>
  auto build_result(NewtonStateT& state, double log_det) {
    return laplace_density_estimates{state.prev().obj() - 0.5 * log_det,
                                     std::move(state.prev().theta()),
                                     std::move(W_full),
                                     Eigen::MatrixXd(llt_B.matrixL()),
                                     std::move(state.prev().a()),
                                     std::move(state.prev().theta_grad()),
                                     Eigen::PartialPivLU<Eigen::MatrixXd>{},
                                     std::move(K_root),
                                     2};
  }
};

/**
 * @brief Solver Policy 3: LU Decomposition.
 *
 * This solver uses LU decomposition with partial pivoting to factorize
 * the system matrix B = I + Sigma * W. This is the most numerically
 * robust solver but also the slowest, as it does not exploit the
 * positive definiteness of B that the Cholesky solvers assume.
 *
 * This solver is used as a fallback when Solvers 1 and 2 fail due to
 * numerical issues (e.g., near-singular Hessians).
 *
 * @note This solver corresponds to `solver == 3` in the original code.
 */
struct LUSolver {
  /** @brief LU factorization of B = I + Sigma * W */
  Eigen::PartialPivLU<Eigen::MatrixXd> lu;

  /** @brief Full Hessian matrix from likelihood */
  Eigen::SparseMatrix<double> W_full;

  /**
   * @brief Perform one Newton step using LU decomposition solver.
   *
   * Computes the full Hessian, forms B = I + Sigma * W, performs
   * LU decomposition, and solves for the new `a` vector.
   *
   * @tparam NewtonStateT Type of the Newton state
   * @tparam LLFun Type of the log-likelihood functor
   * @tparam LLTupleArgs Type of the likelihood arguments tuple
   * @tparam CovarMat Type of the covariance matrix
   * @param[in,out] state Shared Newton state (modified: b, proposal_step().a())
   * @param[in] ll_fun Log-likelihood functor
   * @param[in,out] ll_args Additional arguments for the likelihood
   * @param[in] covariance Prior covariance matrix Sigma
   * @param[in] hessian_block_size Size of each Hessian block
   * @param[in,out] msgs Output stream for diagnostic messages (may be nullptr)
   */
  template <typename NewtonStateT, typename LLFun, typename LLTupleArgs,
            typename CovarMat>
  void solve_step(NewtonStateT& state, const LLFun& ll_fun,
                  const LLTupleArgs& ll_args, const CovarMat& covariance,
                  int hessian_block_size, std::ostream* msgs) {
    const Eigen::Index theta_size = state.b.size();

    // 1. Compute Hessian
    W_full = laplace_likelihood::block_hessian(
        ll_fun, state.prev().theta(), hessian_block_size, ll_args, msgs);

    // 2. Factorize B = I + Sigma * W
    lu.compute(Eigen::MatrixXd::Identity(theta_size, theta_size)
               + covariance * W_full);

    // 3. Solve for the raw Newton proposal in a-space.
    state.b.noalias()
        = W_full * state.prev().theta() + state.prev().theta_grad();
    state.proposal_step().a().noalias()
        = state.b - W_full * lu.solve(covariance * state.b);
  }

  /**
   * @brief Compute log determinant from LU factorization.
   *
   * @note This uses the diagonal of the combined LU matrix produced by Eigen
   * (equivalently the diagonal of U). It does not account for the sign of the
   * permutation; callers assume `det(B) > 0` in the Laplace correction.
   *
   * @return Sum of log of the LU diagonal entries.
   */
  double compute_log_determinant() const {
    return lu.matrixLU().diagonal().array().log().sum();
  }

  /**
   * @brief Build the final result structure.
   *
   * @tparam NewtonStateT Type of the Newton state
   * @param state Newton state containing converged values
   * @param log_det Log determinant from compute_log_determinant()
   * @return laplace_density_estimates with solver_used = 3
   */
  template <typename NewtonStateT>
  auto build_result(NewtonStateT& state, double log_det) {
    return laplace_density_estimates{state.prev().obj() - 0.5 * log_det,
                                     std::move(state).prev().theta(),
                                     std::move(W_full),
                                     Eigen::MatrixXd(0, 0),
                                     std::move(state).prev().a(),
                                     std::move(state).prev().theta_grad(),
                                     std::move(lu),
                                     Eigen::MatrixXd(0, 0),
                                     3};
  }
};

/**
 * @brief Run a Newton loop with a solver policy, updating the shared state.
 *
 * This helper centralizes the iteration/line-search logic shared across
 * solver policies while preserving the step counter and fallback behavior.
 * @tparam SolverPolicy Either Solver 1, 2, or 3 policy struct
 * @tparam NewtonStateT Type holding state with solver
 * @tparam OptionsT struct of options for laplace
 * @tparam LLFunT Callable accepting (ThetaVec, LLTupleArgs...)
 * @tparam LLTupleArgsT tuple of additional args for LLFunT
 * @tparam CovarMatT Type of the covariance matrix
 * @param[in,out] solver The solver policy instance
 * @param[in,out] state Shared Newton optimization state
 * @param[in] options Solver options
 * @param[in,out] step_iter Current iteration number (modified)
 * @param[in] ll_fun Log-likelihood functor
 * @param[in] ll_args Additional arguments for the likelihood
 * @param[in] covariance Prior covariance matrix Sigma
 * @param[in] update_fun Callable to update proposal state in line search
 * @param[in,out] msgs Output stream for diagnostic messages (may be nullptr)
 * @return laplace_density_estimates built from solver policy
 */
template <typename SolverPolicy, typename NewtonStateT, typename OptionsT,
          typename LLFunT, typename LLTupleArgsT, typename CovarMatT,
          typename UpdateFun>
inline auto run_newton_loop(SolverPolicy& solver, NewtonStateT& state,
                            const OptionsT& options, Eigen::Index& step_iter,
                            const LLFunT& ll_fun, const LLTupleArgsT& ll_args,
                            const CovarMatT& covariance, UpdateFun&& update_fun,
                            std::ostream* msgs) {
  bool finish_update = false;
  for (; step_iter <= options.max_num_steps; step_iter++) {
    solver.solve_step(state, ll_fun, ll_args, covariance,
                      options.hessian_block_size, msgs);
    if (!state.final_loop) {
      auto&& proposal = state.proposal_step();
      state.wolfe_info.p_ = proposal.a() - state.prev().a();
      state.prev_g.noalias() = -covariance * state.prev().a()
                               + covariance * state.prev().theta_grad();
      state.wolfe_info.init_dir_ = state.prev_g.dot(state.wolfe_info.p_);
      // Flip direction if not ascending
      state.wolfe_info.flip_direction();
      auto&& scratch = state.wolfe_info.scratch_;
      proposal.eval_.alpha() = 1.0;
      const bool proposal_valid
          = update_fun(proposal, state.curr(), state.prev(), proposal.eval_,
                       state.wolfe_info.p_);
      const bool cached_proposal_ok
          = proposal_valid && std::isfinite(proposal.obj())
            && std::isfinite(proposal.dir())
            && proposal.alpha() > options.line_search.min_alpha;
      if (!cached_proposal_ok) {
        state.wolfe_status
            = WolfeStatus{WolfeReturn::StepTooSmall, 1, 0, false};
      } else if (options.line_search.max_iterations == 0) {
        state.curr().update(proposal);
        state.wolfe_status = WolfeStatus{WolfeReturn::Continue, 1, 0, true};
      } else {
        Eigen::VectorXd s = proposal.a() - state.prev().a();
        auto full_step_grad
            = (-covariance * proposal.a() + covariance * proposal.theta_grad())
                  .eval();
        state.curr().alpha() = barzilai_borwein_step_size(
            s, full_step_grad, state.prev_g, state.prev().alpha(),
            state.wolfe_status.num_backtracks_, options.line_search.min_alpha,
            options.line_search.max_alpha);
        state.wolfe_status = internal::wolfe_line_search(
            state.wolfe_info, update_fun, options.line_search, msgs);
      }
      bool search_failed = !state.wolfe_status.accept_;
      const bool proposal_armijo_ok
          = cached_proposal_ok
            && internal::check_armijo(
                proposal.obj(), state.prev().obj(), proposal.alpha(),
                state.wolfe_info.init_dir_, options.line_search);
      if (search_failed && proposal_armijo_ok) {
        state.curr().update(proposal);
        state.wolfe_status
            = WolfeStatus{WolfeReturn::Armijo, state.wolfe_status.num_evals_,
                          state.wolfe_status.num_backtracks_, true};
        search_failed = false;
      }
      bool objective_converged
          = state.wolfe_status.accept_
            && std::abs(state.curr().obj() - state.prev().obj())
                   < options.tolerance;
      finish_update = objective_converged || search_failed;
    }
    if (finish_update) {
      if (!state.final_loop && state.wolfe_status.accept_) {
        // Do one final loop with exact wolfe conditions
        state.final_loop = true;
        state.update_next_step(options);
        continue;
      }
      return solver.build_result(state, solver.compute_log_determinant());
    } else {
      state.update_next_step(options);
    }
  }
  if (msgs) {
    (*msgs)
        << std::string(
               "WARNING(laplace_marginal_density): max number of iterations: ")
               + std::to_string(options.max_num_steps) + " exceeded.";
  }
  return solver.build_result(state, solver.compute_log_determinant());
}

/**
 * @brief Log a solver fallback event to the provided stream.
 * @param[in] allow_fallthrough If false, throw instead of logging
 * @param[in,out] msgs Output stream (may be nullptr)
 * @param[in] context Context string for the log
 * @param[in] iter Current iteration number
 * @param[in] failed_solver Name of the solver that failed
 * @param[in] next_solver Name of the solver being attempted next
 * @param[in] e Exception that caused the fallback
 */
inline void log_solver_fallback(const bool allow_fallthrough,
                                std::ostream* msgs, std::string_view context,
                                Eigen::Index iter,
                                std::string_view failed_solver,
                                std::string_view next_solver,
                                const std::exception& e) {
  // Build once so we don't interleave with other logs.
  std::ostringstream os;
  std::string msg_type = allow_fallthrough ? "WARNING" : "ERROR";
  os << "[" << context << "] " << msg_type << ": solver fallback\n"
     << "  " << std::left << std::setw(12) << "iteration:" << iter << "\n"
     << "  " << std::left << std::setw(12) << "failed:" << failed_solver << "\n"
     << "  " << std::left << std::setw(12) << "reason:" << e.what() << "\n"
     << "  " << std::left << std::setw(12) << "action:"
     << "trying " << next_solver << "\n"
     << "note: this warning message will only be displayed once."
     << "\n";
  if (allow_fallthrough && msgs) {
    (*msgs) << os.str();
  } else {
    throw std::domain_error(std::string("[") + std::string(context) + "]");
  }
}

template <bool InitTheta, typename Opts>
inline decltype(auto) theta_init_impl(Eigen::Index theta_size, Opts&& options) {
  if constexpr (InitTheta) {
    // If requested, use the prior mean as the initial value
    return std::decay_t<decltype(options)>(options).theta_0;
  } else {
    return Eigen::MatrixXd::Zero(theta_size, 1);
  }
}

/**
 * @brief Create the update function for the line search, capturing necessary
 * references.
 * @tparam ObjFun Callable type for the objective function (accepting (a,
 * theta))
 * @tparam ThetaGradFun Callable type for the theta gradient function (accepting
 * theta)
 * @tparam Covariance Type of the covariance matrix
 * @tparam Options Type of the options struct containing line search parameters
 * @param[in] obj_fun Objective function functor
 * @param[in] theta_grad_f Theta gradient functor
 * @param[in] covariance Prior covariance matrix Sigma
 * @param[in] options Options struct containing line search parameters
 * @return A callable update function for the line search, with signature:
 * ```
 *  bool update_fun(proposal, curr, prev, eval_in, p)
 * ```
 */
template <typename ObjFun, typename ThetaGradFun, typename Covariance,
          typename Options>
inline auto create_update_fun(ObjFun&& obj_fun, ThetaGradFun&& theta_grad_f,
                              Covariance&& covariance, Options&& options) {
  auto update_step = [&covariance, &obj_fun, &theta_grad_f](
                         auto& proposal, auto&& /* curr */, auto&& prev,
                         auto& eval_in, auto&& p) {
    try {
      proposal.a() = prev.a() + eval_in.alpha() * p;
      proposal.theta().noalias() = covariance * proposal.a();
      proposal.theta_grad() = theta_grad_f(proposal.theta());
      eval_in.obj() = obj_fun(proposal.a(), proposal.theta());
      eval_in.dir()
          = (-covariance * proposal.a() + covariance * proposal.theta_grad())
                .dot(p);
      return std::isfinite(eval_in.obj()) && std::isfinite(eval_in.dir());
    } catch (const std::exception&) {
      return false;
    }
  };
  auto backoff = [&options](auto& eval) {
    eval.alpha() *= options.line_search.tau;
    return eval.alpha() > options.line_search.min_alpha;
  };
  return
      [update_step_ = std::move(update_step), backoff_ = std::move(backoff)](
          auto& proposal, auto&& curr, auto&& prev, auto& eval_in, auto&& p) {
        return internal::retry_evaluate(update_step_, proposal, curr, prev,
                                        eval_in, p, backoff_);
      };
}

static STAN_THREADS_DEF std::once_flag fallback_warning;
/**
 * For a latent Gaussian model with hyperparameters phi and
 * latent variables theta, and observations y, this function computes
 * an approximation of the log marginal density, p(y | phi).
 * This is done by marginalizing out theta, using a Laplace
 * approximation. The latter is obtained by finding the mode,
 * via Newton's method, and computing the Hessian of the likelihood.
 *
 * The convergence criterion for the Newton/Wolfe loop is a small change in the
 * optimization objective (stored in the Wolfe step state). The user controls
 * the tolerance (i.e. threshold under which the change is deemed small enough)
 * and maximum number of steps.
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
 * This function returns intermediate quantities needed by
 * `laplace_marginal_density` (including its reverse-mode implementation) to
 * compute gradients and/or generate quantities. Which fields are populated
 * depends on the solver that converged (see `solver_used`).
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
 * @return A `laplace_density_estimates` containing:
 * 1. `lmd`: the Laplace approximation of the log marginal density, p(y | phi)
 * 2. `theta`: the mode of the latent variables
 * 3. `a`: the mode in the `a` parameterization (`theta = covariance * a`)
 * 4. `theta_grad`: gradient of the log-likelihood w.r.t. `theta` at the mode
 * 5. `W_r`: solver-dependent Hessian quantity (see struct docs)
 * 6. `L`: solver-dependent factorization of `B` (Cholesky for solvers 1/2)
 * 7. `LU`: LU factorization of `B` (solver 3 only)
 * 8. `K_root`: Cholesky factor of the covariance (solver 2 only)
 * 9. `solver_used`: the solver policy that produced the result
 *
 */
template <typename LLFun, typename LLTupleArgs, typename CovarMat,
          bool InitTheta,
          require_t<is_all_arithmetic_scalar<CovarMat, LLTupleArgs>>* = nullptr>
inline auto laplace_marginal_density_est(
    LLFun&& ll_fun, LLTupleArgs&& ll_args, CovarMat&& covariance,
    const laplace_options<InitTheta>& options, std::ostream* msgs) {
  internal::validate_laplace_options("laplace_marginal_density", options,
                                     covariance);
  const Eigen::Index theta_size = covariance.rows();
  // Wolfe optimizes over the latent 'a' space
  auto obj_fun = [&ll_fun, &ll_args, &msgs](const Eigen::VectorXd& a_val,
                                            auto&& theta_val) -> double {
    return -0.5 * a_val.dot(theta_val)
           + laplace_likelihood::log_likelihood(ll_fun, theta_val, ll_args,
                                                msgs);
  };
  auto theta_grad_f = [&ll_fun, &ll_args, &msgs](auto&& theta_val) {
    return laplace_likelihood::theta_grad(ll_fun, theta_val, ll_args, msgs);
  };
  decltype(auto) theta_init = theta_init_impl<InitTheta>(theta_size, options);
  // When the user supplies a non-zero theta_init, we must initialise a
  // consistently so that the invariant theta = Sigma * a holds.  Otherwise
  // the prior term -0.5 * a'*theta vanishes (a=0 while theta!=0), inflating
  // the initial objective and causing the Wolfe line search to reject the
  // first Newton step.
  auto state
      = NewtonState(theta_size, obj_fun, theta_grad_f, covariance, theta_init);
  // Start with safe step size
  auto update_fun = create_update_fun(
      std::move(obj_fun), std::move(theta_grad_f), covariance, options);
  Eigen::Index step_iter = 0;
  try {
    if (options.solver == 1) {
      if (options.hessian_block_size == 1) {
        CholeskyWSolverDiag solver(state, covariance);
        return run_newton_loop(solver, state, options, step_iter, ll_fun,
                               ll_args, covariance, update_fun, msgs);
      } else {
        CholeskyWSolverBlock solver(state, options.hessian_block_size);
        return run_newton_loop(solver, state, options, step_iter, ll_fun,
                               ll_args, covariance, update_fun, msgs);
      }
    }
  } catch (const std::exception& e) {
    const std::string solver_type
        = (options.hessian_block_size == 1) ? "Diagonal" : "Block";
    std::string failed = "solver 1 (" + solver_type + " Hessian-root Cholesky)";
    std::call_once(
        fallback_warning,
        [](auto&&... args) {
          log_solver_fallback(std::forward<decltype(args)>(args)...);
        },
        options.allow_fallthrough, msgs, "laplace_marginal_density", step_iter,
        std::move(failed), "solver 2 (Covariance-root Cholesky)", e);
  }
  try {
    if (options.solver == 2 || options.allow_fallthrough) {
      CholeskyKSolver solver(state, covariance);
      return run_newton_loop(solver, state, options, step_iter, ll_fun, ll_args,
                             covariance, update_fun, msgs);
    }
  } catch (const std::exception& e) {
    std::call_once(
        fallback_warning,
        [](auto&&... args) {
          log_solver_fallback(std::forward<decltype(args)>(args)...);
        },
        options.allow_fallthrough, msgs, "laplace_marginal_density", step_iter,
        "solver 2 (Covariance-root Cholesky)", "solver 3 (General LU solver)",
        e);
  }
  if (options.solver == 3 || options.allow_fallthrough) {
    LUSolver solver;
    return run_newton_loop(solver, state, options, step_iter, ll_fun, ll_args,
                           covariance, update_fun, msgs);
  }
  throw std::domain_error(
      std::string("You chose a solver (") + std::to_string(options.solver)
      + ") that is not valid. Please choose either 1, 2, or 3.");
}
}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
