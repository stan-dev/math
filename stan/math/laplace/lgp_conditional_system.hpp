#ifndef STAN_MATH_LAPLACE_LGP_CONDITIONAL_SYSTEM_HPP
#define STAN_MATH_LAPLACE_LGP_CONDITIONAL_SYSTEM_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/mat.hpp>
// #include <stan/math/prim/scal/meta/return_type.hpp>
// #include <stan/math/prim/mat/fun/exp.hpp>
// #include <stan/math/rev/mat/fun/multiply.hpp>
// #include <stan/math/rev/mat/fun/mdivide_left.hpp>
#include <iostream>
#include <string>
#include <vector>
// FIX ME - include specific files instead of rev/mat.hpp to
// improve compilation speed of unit tests.

namespace stan {
namespace math {

  // EXPERIMENT.
  // Functors for finding the mode of the conditional density
  // when doing a (simple) Poisson model with a latent gaussian
  // local parameter. The global parameter is one dimensional.
  // See Dan's experiment.

  /**
   * A structure for the parameters and data in a latent
   * Gaussian Poission (lgp) model. Returns the gradient of 
   * the density of the local parameter (theta) conditioned on 
   * the observed data (y) and the global parameter (phi), and
   * also returns the conditional hessian.
   * 
   * Note theta is not a member of the structure. This is because
   * in one HMC iteration, theta is computed using a Newton 
   * solver (i.e. unlike the other parameters, it is not fixed).
   */
  template<typename T0>
  struct lgp_conditional_system {
    // Eigen::Matrix<T0, Eigen::Dynamic, 1> theta_;  // local parameter
    T0 phi_;  // global parameter
    Eigen::VectorXd n_samples_;  // number of samples for local parameter
    Eigen::VectorXd sums_;  // sums of observation for local parameter

    lgp_conditional_system() {}

    lgp_conditional_system(const T0& phi,
                           const Eigen::VectorXd n_samples,
                           const Eigen::VectorXd& sums)
      : phi_(phi), n_samples_(n_samples), sums_(sums) {}

    /**
     * Functions which retun class members
     */
    T0 get_phi() { return phi_; }
    Eigen::VectorXd get_n_samples() { return n_samples; }
    Eigen::VectorXd get_sums() { return sums_; }
    
    /**
     * An operator that returns the conditional density.
     */

    /**
     * An operator that returns the gradient of the conditional density.
     */
    template <typename T1>
    Eigen::Matrix<typename stan::return_type<T0, T1>::type, 
                         Eigen::Dynamic, 1>
    cond_gradient(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta) {
      return sums_ - elt_multiply(n_samples_, exp(theta))
        - theta / (phi_ * phi_);
    }

    /**
     * An operator that returns the hessian of the conditional density.
     * Required for Newton solver.
     */
    template <typename T1>
    Eigen::Matrix<typename stan::return_type<T0, T1>::type, 
                         Eigen::Dynamic, 1>
    cond_hessian(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta) {
      // addition to vector is element-wise.
      return - add(elt_multiply(n_samples_, exp(theta)), 1/(phi_ * phi_));
    }
    
    /**
     * An operator that returns the Jacobian of theta with respect
     * to phi. This is done using the implicit function theorem.
     * The calculations are greatly simplified by the fact phi
     * is one dimensional and the hessian is a diagonal.
     */
    template<typename T1>
    Eigen::Matrix<typename stan::return_type<T0, T1>::type,
                  Eigen::Dynamic, 1>
    solver_gradient(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta) {
      return - 2 / (phi_ * phi_ * phi_) 
        * elt_divide(theta, add(n_samples_ - 1 / (phi * phi));
    }
  };

}  // namespace math
}  // namespace stan

#endif
