#ifndef STAN_MATH_LAPLACE_LGP_DENSE_SYSTEM_HPP
#define STAN_MATH_LAPLACE_LGP_DENSE_SYSTEM_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
// #include <stan/math/prim/mat/fun/mdivide_left.hpp>
#include <stan/math/rev/mat/fun/mdivide_left.hpp>
// #include <stan/math/rev/math/functor/jacobian.hpp>
#include <stan/math/rev/mat.hpp>
#include <iostream>
#include <string>
#include <vector>
// FIX ME - include specific files instead of rev/mat.hpp to
// improve compilation speed of unit tests.

namespace stan {
namespace math {

  // EXPERIMENT.
  // Functors for finding the mode of the conditional density
  // when doing a Poisson model with a latent gaussian
  // parameter. See Dan's experiment.

  /**
   * A function to constructs the covariance matrix,
   * based on the global paramter, phi.
   * Note the covariance matrix is assumed to have a specific structure,
   * namely same variance for all terms, and same correlation.
   */
  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
  covariance (Eigen::Matrix<T, Eigen::Dynamic, 1> phi, int M) {
    T sigma = phi[0];
    T rho = phi[1];

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sigma(M, M);
    
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < i; j++) {
        Sigma(i, j) = rho * sigma;
        Sigma(j, i) = rho * sigma;
      }
      Sigma(i, i) = sigma;
    }

    return Sigma;
  }

  /**
   * A structure for the parameters and data in a latent
   * Gaussian Poisson (lgp) model. Returns the gradient of 
   * the log density of the local parameter (theta) conditioned
   * on the observed data (y) and the global parameter (phi),
   * and also returns the conditional hessian.
   *
   * Note theta is not a member of the structure. This is because
   * in one HMC iteration, theta is computed using a Newton 
   * solver (i.e. unlike the other parameters, it is not fixed
   * during that iteration).
   */
  template<typename T0>
  struct lgp_dense_system {
    Eigen::Matrix<T0, Eigen::Dynamic, 1> phi_;
    Eigen::VectorXd n_samples_;  // number of samples for local parameter
    Eigen::VectorXd sums_;  // sums of observation for local parameter
    Eigen::Matrix<T0, Eigen::Dynamic, Eigen::Dynamic> Sigma_;  // covariance matrix

    // Constructors
    lgp_dense_system() {}

    lgp_dense_system(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& phi,
                     const Eigen::VectorXd n_samples,
                     const Eigen::VectorXd& sums)
      : phi_(phi), n_samples_(n_samples), sums_(sums) {
      int M = n_samples.size();
      Sigma_ = covariance(phi, M);
    }

    /**
     * Functions which retun class members
     */
    Eigen::Matrix<T0, Eigen::Dynamic, 1> get_phi() const { return phi_; }
    Eigen::VectorXd get_n_samples() const { return n_samples_; }
    Eigen::VectorXd get_sums() const { return sums_; }
    Eigen::MatrixXd get_Sigma() const { return Sigma_; }

    /**
     * An operator that returns the log conditional density, up to a
     * constant term. The returned object works as an
     * objective function we can optimize to find the mode. 
     */
    template <typename T1>
    typename stan::return_type<T0, T1>::type
    log_density (const Eigen::Matrix<T1, Eigen::Dynamic, 1> theta) const {
      Eigen::Matrix<T1, Eigen::Dynamic, 1> 
        poisson_term = elt_multiply(sums_, theta) - exp(theta);
      return sum(poisson_term) -
        0.5 * theta.transpose() * mdivide_left(Sigma_, theta);
    }

    /**
     * An operator that returns the gradient of the conditional density.
     */
    template <typename T1>
    Eigen::Matrix<typename stan::return_type<T0, T1>::type, 
                  Eigen::Dynamic, Eigen::Dynamic>
    cond_gradient(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta) const {
      return sums_ - elt_multiply(n_samples_, exp(theta)) -
        mdivide_left(Sigma_, theta);
    }

    /**
     * An operator that returns the hessian of the conditional density.
     * Required for Newton solver.
     */
    template <typename T1>
    Eigen::Matrix<typename stan::return_type<T0, T1>::type, 
                         Eigen::Dynamic, Eigen::Dynamic>
    cond_hessian(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta) const {
      Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> first_term = 
        elt_multiply(n_samples_, exp(theta)).asDiagonal();

      return - (first_term + inverse(Sigma_));
    }

    /**
     * A functor on which the Jacobian function can be called,
     * to compute the derivative of the objective function with
     * respect to phi.
     * NOTE: because we only want the derivatives with respect to
     * phi, we only compute the (second) term, which depends on phi.
     */
    struct deriv_objective {
      Eigen::VectorXd theta_;
      int M_;

      deriv_objective (const Eigen::VectorXd theta, int M) : 
        theta_(theta), M_(M) { }

      template <typename T>
      inline Eigen::Matrix<T, Eigen::Dynamic, 1>
      operator ()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& phi) const {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
          Sigma = covariance(phi, M_);

        return - mdivide_left(Sigma, theta_);
      }
    };

    /**
     * An operator that returns the Jacobian of theta with respect
     * to phi. This is done using the implicit function theorem.
     * The calculations are greatly simplified by the fact phi
     * is one dimensional and the hessian is a diagonal.
     */
    template<typename T1>
    Eigen::Matrix<typename stan::return_type<T0, T1>::type,
                  Eigen::Dynamic, Eigen::Dynamic>
    solver_gradient(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta) const {
      Eigen::VectorXd dummy;
      Eigen::MatrixXd phi_sensitivities;
      deriv_objective f(value_of(theta), theta.size());
      jacobian(f, phi_, dummy, phi_sensitivities);
      
      return - mdivide_left(cond_hessian(theta), phi_sensitivities);
    }
  };

}  // namespace math
}  // namespace stan

#endif
