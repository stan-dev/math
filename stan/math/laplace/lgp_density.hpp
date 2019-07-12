#ifndef STAN_MATH_LAPLACE_LGP_DENSITY_HPP
#define STAN_MATH_LAPLACE_LGP_DENSITY_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/quad_form_diag.hpp>
#include <stan/math/prim/mat/fun/diag_pre_multiply.hpp>
#include <stan/math/prim/mat/fun/diag_post_multiply.hpp>
#include <stan/math/prim/mat/fun/cholesky_decompose.hpp>
#include <stan/math/prim/mat/fun/sqrt.hpp>
#include <stan/math/rev/mat/fun/cholesky_decompose.hpp>


// Reference: Rasmussen and Williams,
// "Gaussian Processes for Machine Learning",
// The MIT Press, 2006.

namespace stan {
namespace math {

  /**
   * In R & W's notation, this is term a.
   * 
   * Input: initial guess, covaiance matrix, likelihood function,
   *        sum_log_diag_L (reference)
   */
  template <typename T0, typename T1, typename F>
  Eigen::Matrix<typename stan::return_type<T0, T1>::type, Eigen::Dynamic, 1>
  inline
  newton_update (const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
                 const Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>& K,
                 const F& likelihood_function,
                 typename stan::return_type<T0, T1>::type& sum_log_diag_L,
                 bool computing_density = 0) {
    typedef typename stan::return_type<T0, T1>::type scalar;

    Eigen::Matrix<scalar, Eigen::Dynamic, 1> gradient, hessian;
    likelihood_function.diff(theta, gradient, hessian);
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> W = - hessian;

    // WARNING: ill-defined derivatives if W has a zero element.
    // Seriously, how do we handle this??
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> W_root = sqrt(W);

    int M = theta.size();
    Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> L;
    {
      Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>
        B = Eigen::MatrixXd::Identity(M, M) + quad_form_diag(K, W_root);
        L = cholesky_decompose(B);
    }

    if (computing_density)
      sum_log_diag_L = sum(log(diagonal(L)));

    Eigen::Matrix<scalar, Eigen::Dynamic, 1>
      b = elt_multiply(W, theta) + gradient; 

    return b - diag_pre_multiply(W_root,
              mdivide_left(transpose(L),
                mdivide_left_tri_low(L,
                  diag_pre_multiply(W_root, multiply(K, b)))));
  }

  /**
   * Specialized Newton solver for latent Gaussian process, following
   * algorithm 3.1 in Gaussian Processes for Machine Learning, by
   * Rasmussen and Williams.
   * 
   * theta_0: the initial guess
   * phi: the hyper-parameters
   * x: data to be passed to the covariance structure.
   * D: a structure that contains the data and methods to compute
   *   the first and second-order derivatives of the log likelihood.
   * covariance: a structure that defines a method to construct the
   *             the covariance matrix, using phi and x.
   * . . .
   * When terms have not yet been defined, use the same notation
   * as Rasmussen and Williams.
   */
  template <typename T, typename D, typename K>
  Eigen::VectorXd
  gp_newton_solver (const Eigen::VectorXd& theta_0,
                    const Eigen::VectorXd& phi,
                    std::vector<Eigen::Matrix<T, 
                      Eigen::Dynamic, 1>>& x,
                    const D& diff_likelihood,
                    const K& covariance,
                    double function_tolerance = 1e-6,
                    long int max_num_steps = 100) {
    int M = theta_0.size();
    Eigen::MatrixXd Cov = covariance(phi, x, M);

    Eigen::VectorXd theta = theta_0;
    double sum_log_diag_L;  // place holder for Cholesky decomposition.
    double objective_old = - 1e+10;
    double objective_new;
    
    for (int i = 0; i <= max_num_steps; i++) {
      if (i == max_num_steps) {
        std::ostringstream message;
        message << "gp_newton_solver: max number of iterations:"
                << max_num_steps << " exceeded.";
        throw boost::math::evaluation_error(message.str());
      }

      Eigen::VectorXd a = newton_update(theta, Cov, diff_likelihood,
                                        sum_log_diag_L);
      theta = multiply(Cov, a);

      if (i != 0) objective_old = objective_new;
      objective_new = -0.5 * dot_product(a, theta)
                       + diff_likelihood.log_likelihood(theta);
      double objective_diff = abs(objective_new - objective_old);
      if (objective_diff < function_tolerance)
        std::cout << "iterations: " << i << std::endl;
      if (objective_diff < function_tolerance) break;

    }

    return theta;
  }

}  // namespace math
}  // namespace stan

#endif
