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
// "Gaussian Processes for Machine Learning",
// Algorithms 3.1 and 5.1.
// The MIT Press, 2006.
// Note: where I didn't conflict with my own notation, I used their notation.

namespace stan {
namespace math {
  /**
   * Function to compute the log marginal density.
   * Variables needed for the gradient are stored by reference.
   * 
   * theta_0: initial guess for latent gaussian variable.
   * phi: global parameter (input in covariance function)
   * x: fixed data (input in covariance function)
   * D: structure to compute the conditional log likelihood and its
   *    derivatives.
   * K: a functor to compute the covariance function specified by the user.
   * theta: mode of the target.
   * W_root: squared root of the negative diagnonal Hessian.
   * L: Cholesky decomposition of B (stabilized covariance inverse).
   * a: component in the Newton step.
   * tolerance: function tolerance defined in terms of objective function.
   * max_num_steps: number of Newton iterations after which the algorithm stops.
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
        message << "gp_newton_solver: max number of iterations:"
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
      VectorXd b = elt_multiply(W, theta) + l_grad;
      a = b - diag_pre_multiply(W_root,
            mdivide_left_tri<Eigen::Upper>(transpose(L),
              mdivide_left_tri_low(L,
                diag_pre_multiply(W_root, multiply(covariance, b)))));

      // Simple Newton step
      theta = multiply(covariance, a);

      // Check for convergence.
      if (i != 0) objective_old = objective_new;
      objective_new = -0.5 * dot_product(a, theta)
        + diff_likelihood.log_likelihood(theta);
      double objective_diff = abs(objective_new - objective_old);
      if (objective_diff < tolerance)
        std::cout << "iterations: " << i << std::endl;
      if (objective_diff < tolerance) break;
    }

    return objective_new - sum(log(diagonal(L)));
  }

  // A wrapper that doesn't intermediate variables used to
  // compute derivatives.
  template <typename D, typename K>
  double
  laplace_marginal_density (const Eigen::VectorXd& theta_0,
                            const Eigen::VectorXd& phi,
                            const std::vector<Eigen::VectorXd>& x,
                            const D& diff_likelihood,
                            const K& covariance_function,
                            double tolerance = 1e-6,
                            long int max_num_steps = 100) {
    Eigen::VectorXd theta, W_root, a, l_grad;
    Eigen::MatrixXd L;
    return laplace_marginal_density(theta_0, phi, x,
              diff_likelihood, covariance_function,
              theta, W_root, L, a, l_grad,
              tolerance, max_num_steps);
  }

  // a structure to compute the sensitivities of the covariance
  // matrix using forward autodiff.
  template <typename K>
  struct covariance_sensitivities {
    std::vector<Eigen::VectorXd> x_;
    int theta_size_;
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

  // The vari class for the laplace marginal density.
  // The method for computing gradients couples built-in forward autodiff
  // with algorithm 5.1 in R & W.
  struct laplace_marginal_density_vari : public vari {
    int phi_size_;
    vari** phi_;
    int theta_size_;
    vari** marginal_density_;
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
            phi_size_)),
        theta_size_(theta.size()),
        marginal_density_(
          ChainableStack::instance().memalloc_.alloc_array<vari*>(1)) {
      for (int i = 0; i < phi_size_; i++) phi_[i] = phi(i).vi_;

      // CHECK -- is there a cleaner way of doing this?
      marginal_density_[0] = this;
      marginal_density_[0] = new vari(marginal_density, false);

      // compute derivatives of covariance matrix with respect to phi.
      // TO DO -- more efficient approach using covariance function.
      // EXP -- this seems to return the right derivative.
      covariance_sensitivities<K> f(x, theta_size_, covariance_function);
      Eigen::MatrixXd diff_cov;
      Eigen::MatrixXd covariance;
      {
        Eigen::VectorXd covariance_vector;
        jacobian_fwd(f, value_of(phi), covariance_vector, diff_cov);
        covariance = to_matrix(covariance_vector, theta_size_, theta_size_);
      }

      // Now compute the full gradient (using algorithm 5.1 of R & W)
      // Check: is there an efficient way to solve / divide a diagonal matrix?
      Eigen::MatrixXd Z;
      {
        Eigen::MatrixXd W_root_diag = W_root.asDiagonal();
        Z = W_root.asDiagonal() *
              L.transpose().triangularView<Eigen::Upper>()
               .solve(L.triangularView<Eigen::Lower>()
               .solve(W_root_diag));
      }

      Eigen::MatrixXd
        C = mdivide_left_tri<Eigen::Lower>(L, diag_pre_multiply(W_root, covariance));

      Eigen::VectorXd s2 = - 0.5 * (covariance.diagonal()
                 - (C * C.transpose()).diagonal())
                 .cwiseProduct(diff_likelihood.third_diff(theta));

      phi_adj_ = Eigen::VectorXd(phi_size_);
      for (int j = 0; j < phi_size_; j++) {
        Eigen::VectorXd j_col = diff_cov.col(j);
        C = to_matrix(j_col, theta_size_, theta_size_);

        double s1 = 0.5 * quad_form(C, a) - 0.5 * sum((Z * C).diagonal());
        Eigen::VectorXd b = multiply(C, l_grad);
        Eigen::VectorXd s3 = b - multiply(covariance, multiply(Z, b));
        phi_adj_[j] = (s1 + dot_product(s2, s3));
      }
    }

    void chain() { 
      for (int j = 0; j < phi_size_; j++)
        phi_[j]->adj_ += marginal_density_[0]->adj_ * phi_adj_[j];
    }
  };

  // A wrapper for the case where phi is passed as a var.
  // Note: the intial guess is also allowed to be passed as a var.
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
