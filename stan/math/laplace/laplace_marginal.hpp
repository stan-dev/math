#ifndef STAN_MATH_LAPLACE_LAPLACE_MARGINAL_HPP
#define STAN_MATH_LAPLACE_LAPLACE_MARGINAL_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/quad_form_diag.hpp>
#include <stan/math/prim/fun/diag_pre_multiply.hpp>
#include <stan/math/prim/fun/diag_post_multiply.hpp>
#include <stan/math/prim/fun/cholesky_decompose.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/rev/fun/cholesky_decompose.hpp>
#include <stan/math/laplace/laplace_pseudo_target.hpp>
#include <stan/math/laplace/block_matrix_sqrt.hpp>

#include <Eigen/Sparse>
#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>

#include <cmath>

// Reference for calculations of marginal and its gradients:
// Margossian et al, 2020, https://arxiv.org/abs/2004.12550

// TODO -- either use Eigen's .solve() or mdivide_left_tri.
// The code needs to be more consistent.

namespace stan {
namespace math {
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
 * TO DO: add more robust convergence criterion.
 *
 * This algorithm is adapted from Rasmussen and Williams,
 * "Gaussian Processes for Machine Learning", second edition,
 * MIT Press 2006, algorithm 3.1.
 *
 * Variables needed for the gradient or generating quantities
 * are stored by reference.
 *
 * @tparam D structure type for the likelihood object.
 * @tparam K structure type for the covariance object.
 * @tparam Tx type of x, which can in Stan be passed as a matrix or
 *            an array of vectors.
 * @param[in] D structure to compute and differentiate the log likelihood.
 * @param[in] K structure to compute the covariance function.
 * @param[in] phi hyperparameter (input for the covariance function).
 * @param[in] eta hyperparameter (input for likelihood).
 * @param[in] x fixed spatial data (input for the covariance function).
 * @param[in] delta additional fixed real data (input for covariance
 *            function).
 * @param[in] delta_int additional fixed integer data (input for covariance
 *            function).
 * @param[in, out] covariance the evaluated covariance function for the
 *                 latent gaussian variable.
 * @param[in, out] theta a vector to store the mode.
 * @param[in, out] W_r a vector to store the square root of the
 *                 negative Hessian or the negative Hessian, depending
 *                 on which solver we use.
 * @param[in, out] L cholesky decomposition of stabilized inverse covariance.
 * @param[in, out] a element in the Newton step
 * @param[in, out] l_grad the log density of the likelihood.
 * @param[in] theta_0 the initial guess for the mode.
 * @param[in] tolerance the convergence criterion for the Newton solver.
 * @param[in] max_num_steps maximum number of steps for the Newton solver.
 * @param[in] hessian_block_size the size of the block, where we assume
 *              the Hessian is block-diagonal.
 * @param[in] solver which Newton solver to use:
 *                     (1) method using the root of W.
 *                     (2) method using the root of the covariance.
 *                     (3) method using an LU decomposition.
 *
 * @return the log marginal density, p(y | phi).
 */
template <typename D, typename CovarFun, typename Tx>
inline double laplace_marginal_density(
    D&& diff_likelihood, CovarFun&& covariance_function,
    const Eigen::VectorXd& phi, const Eigen::VectorXd& eta, const Tx& x,
    const std::vector<double>& delta, const std::vector<int>& delta_int,
    /* In outs start here */
    Eigen::MatrixXd& covariance, Eigen::VectorXd& theta,
    Eigen::SparseMatrix<double>& W_r, Eigen::MatrixXd& L, Eigen::VectorXd& a,
    Eigen::VectorXd& l_grad, Eigen::PartialPivLU<Eigen::MatrixXd>& LU,
    Eigen::MatrixXd& K_root,
    /* In outs end here */
    const Eigen::VectorXd& theta_0,
    std::ostream* msgs = nullptr, const double tolerance = 1e-6,
    const long int max_num_steps = 100, const int hessian_block_size = 0,
    const int solver = 1, const int do_line_search = 0,
    const int max_steps_line_search = 10) {
  using Eigen::MatrixXd;
  using Eigen::SparseMatrix;
  using Eigen::VectorXd;

  // TODO: Figure out whether this should be on / off or what
  constexpr int diagonal_covariance = 0;
  // solver = 1;
  // hessian_block_size = 1;

  covariance = covariance_function(phi, x, delta, delta_int, msgs);

  if (diagonal_covariance) {
    covariance = covariance.diagonal().asDiagonal();
  }

  const Eigen::Index theta_size = theta_0.size();
  theta = theta_0;
  double objective_old = -1e+10;  // CHECK -- what value to use?
  double objective_inter = -1e+10;
  double objective_new;
  double B_log_determinant;
  Eigen::VectorXd a_old;
  Eigen::VectorXd a_new;
  Eigen::VectorXd theta_new;

  const bool is_hessian_block_size_zero = hessian_block_size == 0;
  if (is_hessian_block_size_zero && solver != 1) {
    constexpr const char* msg
        = "laplace_marginal_density: if treating the Hessian as diagonal"
          " we assume its matrix square-root can be computed."
          " If you don't want to compute the matrix square-root,"
          " set hessian_block_size to 1.";
    throw std::domain_error(std::string(msg));
  }
  Eigen::Index block_size = is_hessian_block_size_zero ? hessian_block_size + 1
                                                       : hessian_block_size;

  for (Eigen::Index i = 0; i <= max_num_steps; i++) {
    if (i == max_num_steps) {
      throw std::domain_error(
          std::string("laplace_marginal_density: max number of iterations: ")
          + std::to_string(max_num_steps) + " exceeded.");
    }

    SparseMatrix<double> W = -diff_likelihood.diff(theta, eta, l_grad, block_size);
    VectorXd b;
    {
      MatrixXd B;
      if (solver == 1) {
        if (is_hessian_block_size_zero) {
          W_r = W.cwiseSqrt();
          B = MatrixXd::Identity(theta_size, theta_size)
              + quad_form_diag(covariance, W_r.diagonal());
        } else {
          W_r = block_matrix_sqrt(W, block_size);
          B = MatrixXd::Identity(theta_size, theta_size)
              + W_r * (covariance * W_r);
        }

        L = cholesky_decompose(B);
        B_log_determinant = 2 * sum(L.diagonal().array().log());

        if (is_hessian_block_size_zero) {
          b = W.diagonal().cwiseProduct(theta) + l_grad.head(theta_size);
          a = b
              - W_r
                    * mdivide_left_tri<Eigen::Upper>(
                        transpose(L),
                        mdivide_left_tri<Eigen::Lower>(
                            L, W_r.diagonal().cwiseProduct(covariance * b)));
        } else {
          b = W * theta + l_grad.head(theta_size);
          a = b
              - W_r
                    * mdivide_left_tri<Eigen::Upper>(
                        transpose(L), mdivide_left_tri<Eigen::Lower>(
                                          L, W_r * (covariance * b)));
        }
      } else if (solver == 2) {
        // TODO -- use triangularView for K_root.
        W_r = W;

        if (diagonal_covariance) {
          K_root = covariance.cwiseSqrt();
        } else {
          K_root = cholesky_decompose(covariance);
        }
        B = MatrixXd::Identity(theta_size, theta_size)
            + K_root.transpose() * W * K_root;
        L = cholesky_decompose(B);
        B_log_determinant = 2 * sum(L.diagonal().array().log());
        b = W * theta + l_grad.head(theta_size);
        a = mdivide_left_tri<Eigen::Upper>(
            K_root.transpose(),
            mdivide_left_tri<Eigen::Upper>(
                L.transpose(),
                mdivide_left_tri<Eigen::Lower>(L, K_root.transpose() * b)));
      } else {
        W_r = W;
        B = MatrixXd::Identity(theta_size, theta_size) + covariance * W;
        LU = Eigen::PartialPivLU<Eigen::MatrixXd>(B);

        // TODO: compute log determinant directly.
        B_log_determinant = log(LU.determinant());

        b = W * theta + l_grad.head(theta_size);
        a = b - W * LU.solve(covariance * b);
      }
    }

    // Simple Newton step
    theta = covariance * a;

    if (i != 0)
      objective_old = objective_new;

    if (std::isfinite(theta.sum())) {
      objective_new
          = -0.5 * a.dot(theta) + diff_likelihood.log_likelihood(theta, eta);
    }

    // linesearch
    // CHECK -- does linesearch work for solver 2?
    int j = 0;
    if (do_line_search && i != 0) {
      while (
          j < max_steps_line_search
          && (objective_new < objective_old || !std::isfinite(theta.sum()))) {
        a = (a + a_old) * 0.5;  // TODO -- generalize for any factor.
        theta = covariance * a;
        if (std::isfinite(theta.sum())) {
          objective_new = -0.5 * a.dot(theta)
                          + diff_likelihood.log_likelihood(theta, eta);
        }

        j++;
      }
    }

    a_old = a;

    // Check for convergence.
    double objective_diff = abs(objective_new - objective_old);

    // if (i % 500 == 0) std::cout << "obj: " << objective_new << std::endl;
    // if (objective_diff < tolerance) std::cout << "iter: " << i << std::endl;

    if (objective_diff < tolerance)
      break;
  }

  return objective_new - 0.5 * B_log_determinant;
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
 * Wrapper for when the hyperparameters passed as a double.
 *
 * @tparam T type of the initial guess.
 * @tparam D structure type for the likelihood object.
 * @tparam K structure type for the covariance object.
 * @tparam Tx type of spatial data for covariance: in Stan, this can
 *            either be a matrix or an array of vectors.
 * @param[in] D structure to compute and differentiate the log likelihood.
 *            The object stores the sufficient stats for the observations.
 * @param[in] K structure to compute the covariance function.
 * @param[in] phi the global parameter (input for the covariance function).
 * @param[in] x data for the covariance function.
 * @param[in] delta additional fixed real data (input for covariance
 *            function).
 * @param[in] delta_int additional fixed integer data (input for covariance
 *            function).
 * @param[in] theta_0 the initial guess for the mode.
 * @param[in] tolerance the convergence criterion for the Newton solver.
 * @param[in] max_num_steps maximum number of steps for the Newton solver.
 * @return the log maginal density, p(y | phi).
 */
// TODO: Operands and partials version of this.
template <typename T, typename D, typename CovarFun, typename Tx>
inline double laplace_marginal_density(
    const D& diff_likelihood, CovarFun&& covariance_function,
    const Eigen::VectorXd& phi, const Eigen::VectorXd& eta, const Tx& x,
    const std::vector<double>& delta, const std::vector<int>& delta_int,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta_0,
    std::ostream* msgs = nullptr, const double tolerance = 1e-6,
    const long int max_num_steps = 100, const int hessian_block_size = 0,
    const int solver = 1, const int do_line_search = 0,
    const int max_steps_line_search = 10) {
  Eigen::VectorXd theta, a, l_grad;
  Eigen::MatrixXd L, covariance, K_root;
  Eigen::SparseMatrix<double> W_r;
  Eigen::PartialPivLU<Eigen::MatrixXd> LU;
  return laplace_marginal_density(
      diff_likelihood, covariance_function, phi, eta, x, delta, delta_int,
      covariance, theta, W_r, L, a, l_grad, LU, K_root, value_of(theta_0), msgs,
      tolerance, max_num_steps, hessian_block_size, solver, do_line_search,
      max_steps_line_search);
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
 * @tparam T0 type of the initial guess.
 * @tparam T1 type of the global parameters.
 * @tparam D structure type for the likelihood object.
 * @tparam K structure type for the covariance object.
 *@tparam Tx type for the spatial data passed to the covariance.
 * @param[in] D structure to compute and differentiate the log likelihood.
 *            The object stores the sufficient stats for the observations.
 * @param[in] K structure to compute the covariance function.
 * @param[in] phi the global parameter (input for the covariance function).
 * @param[in] x data for the covariance function.
 * @param[in] delta addition real data for covariance function.
 * @param[in] delta_int additional interger data for covariance function.
 * @param[in] theta_0 the initial guess for the mode.
 * @param[in] tolerance the convergence criterion for the Newton solver.
 * @param[in] max_num_steps maximum number of steps for the Newton solver.
 * @return the log maginal density, p(y | phi).
 */
template <typename T0, typename T1, typename T2, typename D, typename CovarFun,
          typename Tx>
inline auto laplace_marginal_density(
    const D& diff_likelihood, CovarFun&& covariance_function,
    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
    const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta, const Tx& x,
    const std::vector<double>& delta, const std::vector<int>& delta_int,
    const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
    std::ostream* msgs = nullptr, const double tolerance = 1e-6,
    const long int max_num_steps = 100, const int hessian_block_size = 0,
    const int solver = 1, const int do_line_search = 0,
    const int max_steps_line_search = 10) {
  Eigen::VectorXd theta;
  Eigen::VectorXd a;
  Eigen::VectorXd l_grad;
  Eigen::SparseMatrix<double> W_root;
  Eigen::MatrixXd L;
  Eigen::MatrixXd K_root;
  Eigen::MatrixXd covariance;
  double marginal_density_dbl;
  Eigen::PartialPivLU<Eigen::MatrixXd> LU;

  auto phi_arena = to_arena(phi);
  auto eta_arena = to_arena(eta);

  marginal_density_dbl = laplace_marginal_density(
      diff_likelihood, covariance_function, value_of(phi_arena),
      value_of(eta_arena), x, delta, delta_int, covariance, theta, W_root, L, a,
      l_grad, LU, K_root, value_of(theta_0), msgs, tolerance, max_num_steps,
      hessian_block_size, solver);
  const Eigen::Index theta_size = theta.size();
  const Eigen::Index phi_size_ = phi_arena.size();
  const Eigen::Index eta_size_ = eta_arena.size();
  /*
  for (Eigen::Index i = 0; i < phi_size_; i++) phi_[i] = phi(i).vi_;
  for (Eigen::Index i = 0; i < eta_size_; i++) eta_[i] = eta(i).vi_;
  */
  // CHECK -- is there a cleaner way of doing this?
  //  marginal_density_[0] = this;
  //  marginal_density_[0] = new vari(marginal_density, false);

  Eigen::MatrixXd R;
  Eigen::MatrixXd LU_solve_covariance;
  Eigen::VectorXd eta_dbl = value_of(eta_arena);
  Eigen::VectorXd partial_parm;
  Eigen::VectorXd s2;

  if (solver == 1) {
    Eigen::MatrixXd W_root_diag = W_root;
    R = W_root
        * L.transpose().triangularView<Eigen::Upper>().solve(
            L.triangularView<Eigen::Lower>().solve(W_root_diag));

    Eigen::MatrixXd C = mdivide_left_tri<Eigen::Lower>(L, W_root * covariance);
    if (hessian_block_size == 0 && eta_size_ == 0) {
      s2 = 0.5
           * (covariance.diagonal() - (C.transpose() * C).diagonal())
                 .cwiseProduct(diff_likelihood.third_diff(theta, eta_dbl));
    } else {
      int block_size = (hessian_block_size == 0) ? hessian_block_size + 1
                                                 : hessian_block_size;
      Eigen::MatrixXd A = covariance - C.transpose() * C;
      partial_parm = diff_likelihood.compute_s2(theta, eta_dbl, A, block_size);
      s2 = partial_parm.head(theta_size);
    }
  } else if (solver == 2) {
    // TODO -- use triangularView for K_root.
    R = W_root
        - W_root * K_root
              * L.transpose().triangularView<Eigen::Upper>().solve(
                  L.triangularView<Eigen::Lower>().solve(K_root.transpose()
                                                         * W_root));

    Eigen::MatrixXd C
        = L.triangularView<Eigen::Lower>().solve(K_root.transpose());
    Eigen::MatrixXd A = C.transpose() * C;
    partial_parm
        = diff_likelihood.compute_s2(theta, eta_dbl, A, hessian_block_size);
    s2 = partial_parm.head(theta_size);
  } else {  // solver with LU decomposition
    LU_solve_covariance = LU.solve(covariance);
    R = W_root - W_root * LU_solve_covariance * W_root;

    Eigen::MatrixXd A = covariance - covariance * W_root * LU_solve_covariance;
    // Eigen::MatrixXd A = covariance - covariance * R * covariance;
    partial_parm
        = diff_likelihood.compute_s2(theta, eta_dbl, A, hessian_block_size);
    s2 = partial_parm.head(theta_size);
  }

  arena_matrix<Eigen::VectorXd> phi_adj_arena;
  if (!is_constant<T1>::value && phi_size_ != 0) {
    const nested_rev_autodiff nested;
    Eigen::Matrix<var, Eigen::Dynamic, 1> phi_v = value_of(phi_arena);
    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> K_var
        = covariance_function(phi_v, x, delta, delta_int, msgs);
    Eigen::VectorXd l_grad_theta = l_grad.head(theta_size);
    var Z = laplace_pseudo_target(K_var, a, R, l_grad_theta, s2);
    set_zero_all_adjoints_nested();
    grad(Z.vi_);
    phi_adj_arena = phi_v.adj();
  }

  arena_matrix<Eigen::VectorXd> eta_adj_arena;
  if (!is_constant<T2>::value
      && eta_size_ != 0) {  // TODO: instead, check if eta contains var.
    Eigen::VectorXd diff_eta = l_grad.tail(eta_size_);

    Eigen::VectorXd v;
    if (solver == 1) {
      Eigen::MatrixXd W
          = W_root * W_root;  // CHECK -- store W from Newton step?
      v = covariance * s2 - covariance * R * covariance * s2;
    } else if (solver == 2) {
      v = covariance * s2 - covariance * R * covariance * s2;
    } else {
      v = LU_solve_covariance * s2;
    }

    eta_adj_arena = l_grad.tail(eta_size_) + partial_parm.tail(eta_size_)
                    + diff_likelihood.diff_eta_implicit(v, theta, eta_dbl);
  }

  return make_callback_var(
      marginal_density_dbl, [phi_arena, eta_arena, phi_adj_arena,
                             eta_adj_arena](const auto& vi) mutable {
        if (!is_constant<T1>::value && phi_adj_arena.size() != 0) {
          phi_arena.array().adj() += vi.adj() * phi_adj_arena.array();
        }
        if (!is_constant<T2>::value && eta_adj_arena.size() != 0) {
          eta_arena.array().adj() += vi.adj() * eta_adj_arena.array();
        }
      });
}

}  // namespace math
}  // namespace stan

#endif
