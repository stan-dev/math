#ifndef STAN_MATH_LAPLACE_LGP_DENSITY_HPP
#define STAN_MATH_LAPLACE_LGP_DENSITY_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/quad_form_diag.hpp>
#include <stan/math/prim/mat/fun/diag_pre_multiply.hpp>
#include <stan/math/prim/mat/fun/diag_post_multiply.hpp>
#include <stan/math/prim/mat/fun/cholesky_decompose.hpp>
#include <stan/math/prim/mat/fun/sqrt.hpp>
#include <stan/math/prim/mat/fun/exp.hpp>
#include <stan/math/rev/mat/fun/cholesky_decompose.hpp>


namespace stan {
namespace math {

  // TO DO: construct a general structure that takes in the
  // data and index for each data point, and returns the
  // gradient and hessian for a specified density.

  struct diff_poisson_log {
    Eigen::VectorXd n_samples_;
    Eigen::VectorXd sums_;

    diff_poisson_log() {}

    diff_poisson_log(const Eigen::VectorXd& n_samples,
                     const Eigen::VectorXd& sums)
      : n_samples_(n_samples), sums_(sums) { }

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
      gradient (const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta) const {
      return sums_ - elt_multiply(n_samples_, exp(theta));
    }

    /**
     * Since the Hessian is diagonal, it is stored inside a vector.
     * CHECK -- use above computation to avoid redundant calculation.
     */
    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    hessian (const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta) const {
      return - elt_multiply(n_samples_, exp(theta));
    }
  };

  struct spatial_covariance {
    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    operator() (const Eigen::Matrix<T, Eigen::Dynamic, 1>& phi,
                int M) const {
      int space_matters = true;
      using std::pow;
      T sigma = phi[0];
      T rho = phi[1];
      double exponent;

      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sigma(M, M);

      for (int i = 0; i < M; i++) {
        for (int j = 0; j < i; j++) {
          if (space_matters) {exponent = i - j;} else {exponent = 1;}
          Sigma(i, j) = pow(rho, exponent) * sigma;
          Sigma(j, i) = Sigma(i, j);
        }
        Sigma(i, i) = sigma;
      }

      return Sigma; 
    }
  };

  /**
   * In R & W's notation, this is term a.
   */
  template <typename T0, typename T1, typename D>
  Eigen::Matrix<typename stan::return_type<T0, T1>::type, Eigen::Dynamic, 1>
  newton_update (const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
                 const Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>& K,
                 const D& diff_likelihood,
                 typename stan::return_type<T0, T1>::type& sum_log_diag_L) {
    typedef typename stan::return_type<T0, T1>::type scalar;

    Eigen::Matrix<scalar, Eigen::Dynamic, 1>
      W = - diff_likelihood.hessian(theta);

    // WARNING: ill-defined derivatives if W has a zero element.
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> W_root = sqrt(W);

    int M = theta.size();
    Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> L;
    {
      Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>
        B = Eigen::MatrixXd::Identity(M, M) + quad_form_diag(K, W_root);
        L = cholesky_decompose(B);
    }

    sum_log_diag_L = sum(log(diagonal(L)));

    Eigen::Matrix<scalar, Eigen::Dynamic, 1>
      b = diag_pre_multiply(W, theta) + diff_likelihood.gradient(theta);

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
   * L: a structure that contains the data and methods to compute
   *   the first and second-order derivatives of the log likelihood.
   * phi: the hyper-parameters
   * covariance: a structure that defines a method to construct the
   *             the covariance matrix, using phi as an entry.
   * . . .
   * When terms have not yet been defined, use the same notation
   * as Rasmussen and Williams.
   */
  template <typename D, typename K>
  Eigen::VectorXd
  gp_newton_solver (const Eigen::VectorXd& theta_0,
                    const Eigen::VectorXd& phi,
                    const D& diff_likelihood,
                    const K& covariance,
                    double function_tolerance = 1e-6,
                    long int max_num_steps = 100) {
    int M = theta_0.size();
    Eigen::MatrixXd Cov = covariance(phi, M);
    Eigen::VectorXd theta = theta_0;
    double sum_log_diag_L;  // place holder for Cholesky decomp.

    for (int i = 0; i <= max_num_steps; i++) {
      if (i == max_num_steps) {
        std::ostringstream message;
        message << "gp_newton_solver: max number of iterations:"
                << max_num_steps << " exceeded.";
        throw boost::math::evaluation_error(message.str());
      }

      theta = multiply(Cov, newton_update(theta, Cov,
                                          diff_likelihood,
                                          sum_log_diag_L));

      // convergence criterion
      Eigen::VectorXd
        objective_grad = diff_likelihood.gradient(theta)
          - mdivide_left(Cov, theta);

      if (objective_grad.norm() < function_tolerance)
        std::cout << "iterations: " << i << std::endl;
      if (objective_grad.norm() < function_tolerance) break;
    }

    return theta;
  }

}  // namespace math
}  // namespace stan

#endif
