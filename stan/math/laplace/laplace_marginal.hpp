#ifndef STAN_MATH_LAPLACE_LAPLACE_MARGINAL_HPP
#define STAN_MATH_LAPLACE_LAPLACE_MARGINAL_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/quad_form_diag.hpp>
#include <stan/math/prim/mat/fun/diag_pre_multiply.hpp>
#include <stan/math/prim/mat/fun/diag_post_multiply.hpp>
#include <stan/math/prim/mat/fun/cholesky_decompose.hpp>
#include <stan/math/prim/mat/fun/sqrt.hpp>
#include <stan/math/rev/mat/fun/cholesky_decompose.hpp>
#include <stan/math/laplace/laplace_likelihood.hpp>


// Reference for calculations of marginal and its gradients:
// Rasmussen and Williams,
// "Gaussian Processes for Machine Learning", first edition,
// Algorithms 3.1 and 5.1.
// The MIT Press, 2006.
// Note: where I didn't conflict with my own notation, I used their notation,
// which significantly helps when debuging the code.

namespace stan {
namespace math {
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
   * TO DO: add more robust convergence criterion.
   *
   * This algorithm is adapted from Rasmussen and Williams,
   * "Gaussian Processes for Machine Learning", second edition,
   * MIT Press 2006, algorithm 3.1.
   *
   * Variables needed for the gradient are stored by reference.
   *
   * @tparam D structure type for the likelihood object.
   * @tparam K structure type for the covariance object.
   * @param[in] theta_0 the initial guess for the mode.
   * @param[in] phi the global parameter (input for the covariance function).
   * @param[in] x data for the covariance function.
   * @param[in] D structure to compute and differentiate the log likelihood.
   *            The object stores the sufficient stats for the observations.
   * @param[in] K structure to compute the covariance function.
   * @param[in, out] theta a vector to store the mode.
   * @param[in, out] W_root a vector to store the square root of the 
   *                 diagonal negative Hessian.
   * @param[in, out] L cholesky decomposition of stabilized inverse covariance.
   * @param[in, out] a element in the Newton step
   * @param[in, out] l_grad the log density of the likelihood.
   * @param[in] tolerance the convergence criterion for the Newton solver.
   * @param[in] max_num_steps maximum number of steps for the Newton solver.
   * @return the log marginal density, p(y | phi).
   */
  template <typename D, typename K>
  double
  laplace_marginal_density (const Eigen::VectorXd& theta_0,
                            const Eigen::VectorXd& phi,
                            const std::vector<Eigen::VectorXd>& x,
                            const D& diff_likelihood,
                            const K& covariance_function,
                            Eigen::VectorXd& theta,
                            Eigen::VectorXd& W_root,
                            Eigen::MatrixXd& L,
                            Eigen::VectorXd& a,
                            Eigen::VectorXd& l_grad,
                            double tolerance = 1e-6,
                            long int max_num_steps = 100) {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    int group_size = theta_0.size();
    MatrixXd covariance = covariance_function(phi, x, group_size);
    theta = theta_0;
    double objective_old = - 1e+10;  // CHECK -- what value to use?
    double objective_new;

    for (int i = 0; i <= max_num_steps; i++) {
      if (i == max_num_steps) {
        std::ostringstream message;
        message << "laplace_marginal_density: max number of iterations:"
                << max_num_steps << " exceeded.";
        throw boost::math::evaluation_error(message.str());
      }

      // Compute variable a.
      VectorXd hessian;
      diff_likelihood.diff(theta, l_grad, hessian);
      VectorXd W = - hessian;
      W_root = sqrt(W);
      {
        MatrixXd B = MatrixXd::Identity(group_size, group_size)
          + quad_form_diag(covariance, W_root);
        L = cholesky_decompose(B);
      }
      VectorXd b = W.cwiseProduct(theta) + l_grad;
      a = b - W_root.asDiagonal() * mdivide_left_tri<Eigen::Upper>(transpose(L),
           mdivide_left_tri<Eigen::Lower>(L,
           diag_pre_multiply(W_root, multiply(covariance, b))));

      // Simple Newton step
      theta = covariance * a;

      // Check for convergence.
      if (i != 0) objective_old = objective_new;
      objective_new = -0.5 * a.dot(theta)
        + diff_likelihood.log_likelihood(theta);
      double objective_diff = abs(objective_new - objective_old);
      if (objective_diff < tolerance) break;
    }

    return objective_new - sum(L.diagonal().array().log());
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
   * Wrapper for when the global parameter is passed as a double..
   *
   * @tparam T type of the initial guess.
   * @tparam D structure type for the likelihood object.
   * @tparam K structure type for the covariance object.
   * @param[in] theta_0 the initial guess for the mode.
   * @param[in] phi the global parameter (input for the covariance function).
   * @param[in] x data for the covariance function.
   * @param[in] D structure to compute and differentiate the log likelihood.
   *            The object stores the sufficient stats for the observations.
   * @param[in] K structure to compute the covariance function.
   * @param[in] tolerance the convergence criterion for the Newton solver.
   * @param[in] max_num_steps maximum number of steps for the Newton solver.
   * @return the log maginal density, p(y | phi).
   */
  template <typename T, typename D, typename K>
  double
  laplace_marginal_density (const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta_0,
                            const Eigen::VectorXd& phi,
                            const std::vector<Eigen::VectorXd>& x,
                            const D& diff_likelihood,
                            const K& covariance_function,
                            double tolerance = 1e-6,
                            long int max_num_steps = 100) {
    Eigen::VectorXd theta, W_root, a, l_grad;
    Eigen::MatrixXd L;
    return laplace_marginal_density(value_of(theta_0), phi, x,
              diff_likelihood, covariance_function,
              theta, W_root, L, a, l_grad,
              tolerance, max_num_steps);
  }

  /**
   * A structure to the compute sensitivities of the covariance
   * function using forward mode autodiff. The functor is formatted
   * so that it can be passed to Jacobian(). This requires one input
   * vector and one output vector.
   *
   * TO DO: make this structure no templated. See comment by @SteveBronder.
   */
  template <typename K>
  struct covariance_sensitivities {
    /* input data for the covariance function. */
    std::vector<Eigen::VectorXd> x_;
    /* number of latent variables. */
    int theta_size_;
    /* structure to compute the covariance function. */
    K covariance_function_;
    
    covariance_sensitivities (const std::vector<Eigen::VectorXd>& x,
                              int theta_size,
                              const K& covariance_function) :
    x_(x), theta_size_(theta_size),
    covariance_function_(covariance_function) { }

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    operator() (const Eigen::Matrix<T, Eigen::Dynamic, 1>& phi) const {
      return to_vector(covariance_function_(phi, x_, theta_size_));
    }
  };

  /**
   * The vari class for the laplace marginal density.
   * The method is adapted from algorithm 5.1 in Rasmussen & Williams,
   * "Gaussian Processes for Machine Learning".
   * The covariance function is differentiated using forward autodiff.
   *
   * To make computation efficient, variables produced during the
   * Newton step are stored and reused. To avoid storing these variables
   * for too long, the sensitivies are computed in the constructor, and
   * store for the chain method. Hence, we store a single small vector,
   * instead of multiple large matrices.
   */
  struct laplace_marginal_density_vari : public vari {
    /* dimension of the global parameters. */
    int phi_size_;
    /* global parameters. */
    vari** phi_;
    /* the marginal density of the observation, conditional on the 
     * globl parameters. */
    vari** marginal_density_;
    /* An object to store the sensitivities of phi. */
    Eigen::VectorXd phi_adj_;

    template <typename T, typename K, typename D>
    laplace_marginal_density_vari
      (double marginal_density,
       const Eigen::Matrix<T, Eigen::Dynamic, 1>& phi,
       const std::vector<Eigen::VectorXd>& x,
       const D& diff_likelihood,
       const K& covariance_function,
       const Eigen::VectorXd& theta,
       const Eigen::VectorXd& W_root,
       const Eigen::MatrixXd& L,
       const Eigen::VectorXd& a,
       const Eigen::VectorXd& l_grad)
      : vari(marginal_density),
        phi_size_(phi.size()),
        phi_(ChainableStack::instance().memalloc_.alloc_array<vari*>(
	        phi.size())),
        // theta_size_(theta.size()),
        marginal_density_(
          ChainableStack::instance().memalloc_.alloc_array<vari*>(1)) {
      int theta_size = theta.size();
      for (int i = 0; i < phi_size_; i++) phi_[i] = phi(i).vi_;

      // CHECK -- is there a cleaner way of doing this?
      marginal_density_[0] = this;
      marginal_density_[0] = new vari(marginal_density, false);

      // compute derivatives of covariance matrix with respect to phi.
      covariance_sensitivities<K> f(x, theta_size, covariance_function);
      Eigen::MatrixXd diff_cov;
      Eigen::MatrixXd covariance;
      {
        Eigen::VectorXd covariance_vector;
        jacobian_fwd(f, value_of(phi), covariance_vector, diff_cov);
        covariance = to_matrix(covariance_vector, theta_size, theta_size);
      }

      // Now compute the full gradient (using algorithm 5.1 of R & W)
      // Check: is there an efficient way to solve / divide a diagonal matrix?
      Eigen::MatrixXd Z;
      {
        Eigen::MatrixXd W_root_diag = W_root.asDiagonal();
        Z = W_root_diag *
              L.transpose().triangularView<Eigen::Upper>()
               .solve(L.triangularView<Eigen::Lower>()
                 .solve(W_root_diag));
      }

      Eigen::MatrixXd
        C = mdivide_left_tri<Eigen::Lower>(L,
                  diag_pre_multiply(W_root, covariance));

      // CHECK -- should there be a minus sign here?
      Eigen::VectorXd s2 = 0.5 * (covariance.diagonal()
                 - (C.transpose() * C).diagonal())
                 .cwiseProduct(diff_likelihood.third_diff(theta));

      phi_adj_ = Eigen::VectorXd(phi_size_);
      for (int j = 0; j < phi_size_; j++) {
        Eigen::VectorXd j_col = diff_cov.col(j);
        C = to_matrix(j_col, theta_size, theta_size);
        double s1 = 0.5 * quad_form(C, a) - 0.5 * sum((Z * C).diagonal());
        Eigen::VectorXd b = C * l_grad;
        Eigen::VectorXd s3 = b - covariance * (Z * b);
        phi_adj_[j] = s1 + s2.dot(s3);
      }
    }

    void chain() { 
      for (int j = 0; j < phi_size_; j++)
        phi_[j]->adj_ += marginal_density_[0]->adj_ * phi_adj_[j];
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
   * @param[in] theta_0 the initial guess for the mode.
   * @param[in] phi the global parameter (input for the covariance function).
   * @param[in] x data for the covariance function.
   * @param[in] D structure to compute and differentiate the log likelihood.
   *            The object stores the sufficient stats for the observations.
   * @param[in] K structure to compute the covariance function.
   * @param[in] tolerance the convergence criterion for the Newton solver.
   * @param[in] max_num_steps maximum number of steps for the Newton solver.
   * @return the log maginal density, p(y | phi).
   */
  template <typename T0, typename T1, typename D, typename K>
  T1 laplace_marginal_density
      (const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta_0,
       const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
       const std::vector<Eigen::VectorXd>& x,
       const D& diff_likelihood,
       const K& covariance_function,
       double tolerance = 1e-6,
       long int max_num_steps = 100) {
    Eigen::VectorXd theta, W_root, a, l_grad;
    Eigen::MatrixXd L;
    double marginal_density_dbl
      = laplace_marginal_density(value_of(theta_0),
                                 value_of(phi),
                                 x, diff_likelihood,
                                 covariance_function,
                                 theta, W_root, L, a, l_grad,
                                 tolerance, max_num_steps);

    // construct vari
    laplace_marginal_density_vari* vi0
      = new laplace_marginal_density_vari(marginal_density_dbl,
                                          phi, x, diff_likelihood,
                                          covariance_function,
                                          theta, W_root, L, a, l_grad);

    var marginal_density = var(vi0->marginal_density_[0]);

    return marginal_density;
  }
          
  
}  // namespace math
}  // namespace stan

#endif
