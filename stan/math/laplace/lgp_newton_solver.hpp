#ifndef STAN_MATH_LAPLACE_LGP_NEWTON_SOLVER_HPP
#define STAN_MATH_LAPLACE_LGP_NEWTON_SOLVER_HPP

#include <stan/math/laplace/lgp_conditional_system.hpp>


namespace stan {
namespace math {

  /**
   * The vari class for the latent gaussian poisson (lgp)
   * Newton solver.
   */
  template <typename T>
  struct lgp_newton_solver_vari : public vari {
    /** global parameter */
    vari** phi_;
    /** vector of solution */
    vari** theta_;
    /** size of solution */
    int theta_size_;
    /** Jacobian of the solution with respect to the global parameter */
    double* J_;

    lgp_newton_solver_vari(const T& phi,
                           const lgp_conditional_system<T>& system,
                           const Eigen::VectorXd& theta_dbl)
      : vari(theta_dbl(0)),
        phi_(ChainableStack::instance_.memalloc_.alloc_array<vari*>(
            theta.size())),
        theta_size_(theta.size()),
        theta_(ChainableStack::instance_.memalloc_.alloc_array<vari*>(
          theta_size_)) {
      using Eigen::Map;
      using Eigen::MatrixXd;

      phi_ = phi.vi_;

      theta_[0] = this;
      for (int i = 0; i < theta.size; i++)
        theta_[i] = new vari(theta_dbl(i), false);

      // Compute the Jacobian and store in array, using
      // the operator in the lgp_conditional_system structure.
      Map<VectorXd>(&J_[0], theta_size_) = system.solver_gradient(theta_);
    }

    void chain() {
      for (int i = 0; i < theta_size_; i++)
        phi_[i]->adj_ += theta_[i]->adj_ * J_[i];
    }
  };
  
  /**
   * Newton solver for lgp Newton solver.
   * First definition of the function considers the case where
   * the global parameter, phi, has type double.
   */
  template<typename T>  // template for variables
  Eigen::Matrix<T, Eigen::Dynamic, 1> lgp_newton_solver(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta_0,  // initial guess
    const lgp_conditional_system<double>& system,
    double tol,
    long int max_num_steps) {  // NOLINT(runtime/int)

    Eigen::VectorXd theta_dbl = value_of(theta);
    Eigen::VectorXd gradient; 

    for (int i = 0; i <= max_num_steps; i++) {

      // check if the max number of steps has been reached
      if (i == max_num_steps) {
        std::ostringstream message;
        message << "lgp_newton_solver: max number of iterations:"
                << max_num_steps << " exceeded.";
        throw boost::math::evaluation_error(message.str());
      }

      gradient = system.cond_gradient(theta));
      theta -= mdivide_left(gradient, system.cond_hessian(theta));

      // Check solution is a root of the gradient
      if (gradient.norm() <= tol) break;
    }

    return theta_dbl;
  }

  /**
   * lgp Newton solver.
   * Case where the global parameter is a var.
   */
  template <typename T1, typename T2>
  Eigen::Matrix<T2, Eigen::Dynamic, 1> lgp_newton_solver(
    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta_0,  // initial guess
    const lgp_conditional_system<T2>& system,
    double tol,
    long int max_num_steps) {  // NOLINT(runtime/int)

    lgp_conditional_system<double> 
      system_dbl(value_of(system.get_phi()),
                 system.get_n_samples(), 
                 system.get_sums());

    Eigen::VectorXd theta_dbl = lgp_newton_solver(theta_0, system_dbl,
                                                  tol, max_num_steps);

    // construct vari
    lgp_newton_solver_vari<T>* vi0
      = new lgp_solver_vari<T>(system.get_phi(), system, theta_dbl);
    
    Eigen::Matrix<var, Eigen::Dynamic, 1> theta(theta_dbl.size());
    theta(0) = var(vi0->theta_[0]);
    for (int i = 1; i < theta_dbl.size(); i++)
      theta(i) = var(vi0->theta_[i]);

    return theta;
  }

}  // namespace math
}  // namespace stan

#endif
