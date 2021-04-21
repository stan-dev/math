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

#include <iostream>
#include <istream>  // CHECK -- do we need this?
#include <fstream>  // CHECK -- do we need this?

// Reference for calculations of marginal and its gradients:
// Margossian et al, 2020, https://arxiv.org/abs/2004.12550


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
   * @param[in, out] W_root a vector to store the square root of the
   *                 diagonal negative Hessian.
   * @param[in, out] L cholesky decomposition of stabilized inverse covariance.
   * @param[in, out] a element in the Newton step
   * @param[in, out] l_grad the log density of the likelihood.
   * @param[in] theta_0 the initial guess for the mode.
   * @param[in] tolerance the convergence criterion for the Newton solver.
   * @param[in] max_num_steps maximum number of steps for the Newton solver.
   * @return the log marginal density, p(y | phi).
   */
  template <typename D, typename K, typename Tx>
  double
  laplace_marginal_density (const D& diff_likelihood,
                            const K& covariance_function,
                            const Eigen::VectorXd& phi,
                            const Eigen::VectorXd& eta,
                            const Tx& x,
                            const std::vector<double>& delta,
                            const std::vector<int>& delta_int,
                            Eigen::MatrixXd& covariance,
                            Eigen::VectorXd& theta,
                            Eigen::SparseMatrix<double>& W_r,
                            Eigen::MatrixXd& L,
                            Eigen::VectorXd& a,
                            Eigen::VectorXd& l_grad,
                            Eigen::PartialPivLU<Eigen::MatrixXd>& LU,
                            const Eigen::VectorXd& theta_0,
                            std::ostream* msgs = nullptr,
                            double tolerance = 1e-6,
                            long int max_num_steps = 100,
                            int hessian_block_size = 0,
                            int compute_W_root = 1) {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::SparseMatrix;

    int theta_size = theta_0.size();
    covariance = covariance_function(phi, x, delta, delta_int, msgs);
    theta = theta_0;
    double objective_old = - 1e+10;  // CHECK -- what value to use?
    double objective_inter = - 1e+10;
    double objective_new;
    double B_log_determinant;
    Eigen::VectorXd a_old;
    int j;

    if (hessian_block_size == 0 && compute_W_root == 0) {
      std::ostringstream message;
      message << "laplace_marginal_density: if treating the Hessian as diagonal"
        << " we assume its matrix square-root can be computed."
        << " If you don't want to compute the matrix square-root,"
        << " set hessian_block_size to 1.";
      throw boost::math::evaluation_error(message.str());
    }

    int block_size = (hessian_block_size == 0) ? hessian_block_size + 1
                      : hessian_block_size;

    for (int i = 0; i <= max_num_steps; i++) {
      if (i == max_num_steps) {
        std::ostringstream message;
        message << "laplace_marginal_density: max number of iterations:"
                << max_num_steps << " exceeded.";
        throw boost::math::evaluation_error(message.str());
      }

      SparseMatrix<double> W;
      diff_likelihood.diff(theta, eta, l_grad, W, block_size);
      W = - W;

      VectorXd b;
      {
        MatrixXd B;
        if (compute_W_root) {
          if (hessian_block_size == 0) {
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

          if (hessian_block_size == 0) {
            b = W.diagonal().cwiseProduct(theta) + l_grad.head(theta_size);
            a = b - W_r
              * mdivide_left_tri<Eigen::Upper>(transpose(L),
                 mdivide_left_tri<Eigen::Lower>(L,
                   W_r.diagonal().cwiseProduct(covariance * b)));
          } else {
            b = W * theta + l_grad.head(theta_size);
            a = b - W_r
              * mdivide_left_tri<Eigen::Upper>(transpose(L),
                  mdivide_left_tri<Eigen::Lower>(L,
                  W_r * (covariance * b)));
          }
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

      if (i != 0) objective_old = objective_new;
      objective_new = -0.5 * a.dot(theta)
        + diff_likelihood.log_likelihood(theta, eta);

      // linesearch method
      int do_line_search = 1;
      int max_steps_line_search = 10;
      if (do_line_search && i != 0) {  // CHECK -- no line search at first step?
        j = 0;

        // CHECK -- should we use a different convergence criterion?
        // while (j <= max_steps_line_search || objective_new < objective_old) {
        while (j <= max_steps_line_search || objective_new > objective_inter) {

          a = (a + a_old) * 0.5;  // CHECK -- generalize this for any reduction?
          theta = covariance * a;

          objective_inter = objective_new;
          objective_new = - 0.5 * a.dot(theta)
            + diff_likelihood.log_likelihood(theta, eta);

          j += 1;
        }
      }

      a_old = a;

      // Check for convergence.
      double objective_diff = abs(objective_new - objective_old);

       if (i % 500 == 0) std::cout << "obj: " << objective_new << std::endl;

      // if (objective_diff < tolerance) std::cout << "iter: " << i << std::endl;
      if (objective_diff < tolerance) break;
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
  template <typename T, typename D, typename K, typename Tx>
  double
  laplace_marginal_density (const D& diff_likelihood,
                            const K& covariance_function,
                            const Eigen::VectorXd& phi,
                            const Eigen::VectorXd& eta,
                            const Tx& x,
                            const std::vector<double>& delta,
                            const std::vector<int>& delta_int,
                            const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta_0,
                            std::ostream* msgs = nullptr,
                            double tolerance = 1e-6,
                            long int max_num_steps = 100,
                            int hessian_block_size = 0,
                            int compute_W_root = 1) {
    Eigen::VectorXd theta, a, l_grad;
    Eigen::MatrixXd L, covariance;
    Eigen::SparseMatrix<double> W_r;
    Eigen::PartialPivLU<Eigen::MatrixXd> LU;
    return laplace_marginal_density(diff_likelihood, covariance_function,
                                    phi, eta, x, delta, delta_int,
                                    covariance,
                                    theta, W_r, L, a, l_grad, LU,
                                    value_of(theta_0), msgs,
                                    tolerance, max_num_steps,
                                    hessian_block_size,
                                    compute_W_root);
  }

  /**
   * The vari class for the laplace marginal density.
   * The method is adapted from algorithm 5.1 in Rasmussen & Williams,
   * "Gaussian Processes for Machine Learning"
   * with modifications described in my (Charles Margossian)
   * thesis proposal.
   *
   * To make computation efficient, variables produced during the
   * Newton step are stored and reused. To avoid storing these variables
   * for too long, the sensitivies are computed in the constructor, and
   * stored for the chain method. Hence, we store a single small vector,
   * instead of multiple large matrices.
   */
  struct laplace_marginal_density_vari : public vari {
    /* dimension of hyperparameters. */
    int phi_size_;
    /* hyperparameters for covariance K. */
    vari** phi_;
    /* dimension of hyperparameters for likelihood. */
    int eta_size_;
    /* hyperparameters for likelihood. */
    vari** eta_;
    /* the marginal density of the observation, conditional on the
     * globl parameters. */
    vari** marginal_density_;
    /* An object to store the sensitivities of phi. */
    Eigen::VectorXd phi_adj_;
    /* An object to store the sensitivities of eta. */
    Eigen::VectorXd eta_adj_;

    template <typename T1, typename T2, typename K, typename D, typename Tx>
    laplace_marginal_density_vari
      (const D& diff_likelihood,
       const K& covariance_function,
       const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
       const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta,
       const Tx& x,
       const std::vector<double>& delta,
       const std::vector<int>& delta_int,
       double marginal_density,
       const Eigen::MatrixXd& covariance,
       const Eigen::VectorXd& theta,
       // const Eigen::MatrixXd& W_root,
       const Eigen::SparseMatrix<double>& W_r,
       const Eigen::MatrixXd& L,
       const Eigen::VectorXd& a,
       const Eigen::VectorXd& l_grad,
       const Eigen::PartialPivLU<Eigen::MatrixXd> LU,
       std::ostream* msgs = nullptr,
       int hessian_block_size = 0,
       int compute_W_root = 1)
      : vari(marginal_density),
        phi_size_(phi.size()),
        phi_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
	        phi.size())),
        eta_size_(eta.size()),
        eta_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
          eta.size())),
        marginal_density_(
          ChainableStack::instance_->memalloc_.alloc_array<vari*>(1)) {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using Eigen::MatrixXd;
      using Eigen::VectorXd;
      using Eigen::SparseMatrix;

      int theta_size = theta.size();
      for (int i = 0; i < phi_size_; i++) phi_[i] = phi(i).vi_;
      for (int i = 0; i < eta_size_; i++) eta_[i] = eta(i).vi_;

      // CHECK -- is there a cleaner way of doing this?
      marginal_density_[0] = this;
      marginal_density_[0] = new vari(marginal_density, false);

      MatrixXd R;
      Eigen::MatrixXd LU_solve_covariance;
      Eigen::VectorXd eta_dbl = value_of(eta);
      Eigen::VectorXd partial_parm;
      Eigen::VectorXd s2;

      if (compute_W_root == 1) {
        MatrixXd W_root_diag = W_r;
        R = W_r * L.transpose().triangularView<Eigen::Upper>()
                                      .solve(L.triangularView<Eigen::Lower>()
                                        .solve(W_root_diag));

        Eigen::MatrixXd C = mdivide_left_tri<Eigen::Lower>(L, W_r * covariance);
        if (hessian_block_size == 0 && eta_size_ == 0) {
          s2 = 0.5 * (covariance.diagonal()
                 - (C.transpose() * C).diagonal())
                  .cwiseProduct(diff_likelihood.third_diff(theta, eta_dbl));
        } else {
          int block_size = (hessian_block_size == 0) ? hessian_block_size + 1
                            : hessian_block_size;
          Eigen::MatrixXd A = covariance - C.transpose() * C;
          partial_parm
            = diff_likelihood.compute_s2(theta, eta_dbl, A, block_size);
          s2 = partial_parm.head(theta_size);
        }
      } else {  // we have not computed W_root.
        LU_solve_covariance = LU.solve(covariance);
        R = W_r - W_r * LU_solve_covariance * W_r;

        Eigen::MatrixXd A = covariance - covariance * W_r * LU_solve_covariance;
        // Eigen::MatrixXd A = covariance - covariance * R * covariance;
        partial_parm
          = diff_likelihood.compute_s2(theta, eta_dbl, A, hessian_block_size);
        s2 = partial_parm.head(theta_size);
      }

     phi_adj_ = Eigen::VectorXd(phi_size_);
     start_nested();
     try {
       Matrix<var, Dynamic, 1> phi_v = value_of(phi);
       Matrix<var, Dynamic, Dynamic>
         K_var = covariance_function(phi_v, x, delta, delta_int, msgs);
       Eigen::VectorXd l_grad_theta = l_grad.head(theta_size);
       var Z = laplace_pseudo_target(K_var, a, R, l_grad_theta, s2);

       set_zero_all_adjoints_nested();
       grad(Z.vi_);

       for (int j = 0; j < phi_size_; j++) phi_adj_[j] = phi_v(j).adj();

    } catch (const std::exception& e) {
      recover_memory_nested();
      throw;
    }
    recover_memory_nested();

    eta_adj_ = Eigen::VectorXd(eta_size_);
    if (eta_size_ != 0) {  // TODO: instead, check if eta contains var.
      VectorXd diff_eta = l_grad.tail(eta_size_);

      Eigen::VectorXd v;
      if (compute_W_root == 1) {
        Eigen::MatrixXd W = W_r * W_r;  // NOTE: store W from Newton step?
        v = covariance * s2
          - covariance * R * covariance * s2;
          // - covariance * W
          // * L.transpose().triangularView<Eigen::Upper>()
          //     . solve(L.triangularView<Eigen::Lower>()
          //       .solve(covariance * (covariance * s2)));
      } else {
        v = LU_solve_covariance * s2;
      }

      eta_adj_ = l_grad.tail(eta_size_) + partial_parm.tail(eta_size_)
        + diff_likelihood.diff_eta_implicit(v, theta, eta_dbl);
    }
  }

    void chain() {
      for (int j = 0; j < phi_size_; j++)
        phi_[j]->adj_ += marginal_density_[0]->adj_ * phi_adj_[j];

      for (int l = 0; l < eta_size_; l++)
        eta_[l]->adj_ += marginal_density_[0]->adj_ * eta_adj_[l];
    }
  };

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
  template <typename T0, typename T1, typename T2, typename D, typename K,
            typename Tx>
  T1 laplace_marginal_density
      (const D& diff_likelihood,
       const K& covariance_function,
       const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
       const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta,
       const Tx& x,
       const std::vector<double>& delta,
       const std::vector<int>& delta_int,
       const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
       std::ostream* msgs = nullptr,
       double tolerance = 1e-6,
       long int max_num_steps = 100,
       int hessian_block_size = 0,
       int compute_W_root = 1) {
    Eigen::VectorXd theta, a, l_grad;
    Eigen::SparseMatrix<double> W_root;
    Eigen::MatrixXd L;
    double marginal_density_dbl;
    Eigen::MatrixXd covariance;
    Eigen::PartialPivLU<Eigen::MatrixXd> LU;


    marginal_density_dbl
      = laplace_marginal_density(diff_likelihood,
                                 covariance_function,
                                 value_of(phi), value_of(eta),
                                 x, delta, delta_int, covariance,
                                 theta, W_root, L, a, l_grad, LU,
                                 value_of(theta_0),
                                 msgs,
                                 tolerance, max_num_steps,
                                 hessian_block_size,
                                 compute_W_root);

    // construct vari
    laplace_marginal_density_vari* vi0
      = new laplace_marginal_density_vari(diff_likelihood,
                                          covariance_function,
                                          phi, eta, x, delta, delta_int,
                                          marginal_density_dbl,
                                          covariance,
                                          theta, W_root, L, a, l_grad, LU,
                                          msgs, hessian_block_size,
                                          compute_W_root);

    var marginal_density = var(vi0->marginal_density_[0]);

    return marginal_density;
  }

}  // namespace math
}  // namespace stan

#endif
