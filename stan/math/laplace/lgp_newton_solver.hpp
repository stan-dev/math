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
    /** size of solution */
    int theta_size_;
    /** vector of solution */
    vari** theta_;
    /** Jacobian of the solution with respect to the global parameter */
    Eigen::VectorXd J_;
    // double* J_;

    lgp_newton_solver_vari(const T& phi,
                           const lgp_conditional_system<double>& system,
                           const Eigen::VectorXd& theta_dbl)
      : vari(theta_dbl(0)),
        phi_(ChainableStack::instance_.memalloc_.alloc_array<vari*>(1)),
        theta_size_(theta_dbl.size()),
        theta_(ChainableStack::instance_.memalloc_.alloc_array<vari*>(
          theta_size_)) {
      using Eigen::Map;
      using Eigen::VectorXd;

      *phi_ = phi.vi_;  // CHECK - should be phi_ = phi.vi_?

      theta_[0] = this;
      for (int i = 0; i < theta_size_; i++)
        theta_[i] = new vari(theta_dbl(i), false);

      // Compute the Jacobian and store in array, using
      // the operator in the lgp_conditional_system structure.
      J_ = system.solver_gradient(theta_dbl);

      // CHECK - more memory efficient?
      // std::cout << system.solver_gradient(theta_dbl) << std::endl;
      // Map<VectorXd>(&J_[0], theta_size_) = system.solver_gradient(theta_dbl);
      // std::cout << "marker d" << std::endl;
    }

    void chain() {
      for (int i = 0; i < theta_size_; i++)
        phi_[0]->adj_ += theta_[i]->adj_ * J_[i];
    }
  };
  
  /**
   * Newton solver for lgp model.
   * In this instantiation, the global parameter phi has type double.
   * The initial guess can be passed as a parameter or fixed data.
   */
  template<typename T>  // template for variables
  Eigen::Matrix<T, Eigen::Dynamic, 1> lgp_newton_solver(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta_0,  // initial guess
    const lgp_conditional_system<double>& system,
    double tol = 1e-3,
    long int max_num_steps = 100) {  // NOLINT(runtime/int)

    Eigen::VectorXd theta_dbl = value_of(theta_0);
    Eigen::MatrixXd gradient; 

    for (int i = 0; i <= max_num_steps; i++) {

      // std::cout << "Iteration: " << i << " yields theta = " << theta_dbl << std::endl;

      // check if the max number of steps has been reached
      if (i == max_num_steps) {
        std::ostringstream message;
        message << "lgp_newton_solver: max number of iterations:"
                << max_num_steps << " exceeded.";
        throw boost::math::evaluation_error(message.str());
      }

      gradient = system.cond_gradient(theta_dbl);
      theta_dbl -= elt_divide(gradient, system.cond_hessian(theta_dbl));

      // std::cout << "gradient: " << gradient << std::endl;
      // std::cout << "cond_hessian: " << system.cond_hessian(theta_dbl) << std::endl;

      // Should really be an mdivide. However, the hessian is diagonal and 
      // stored as a vector, so need to use elt_divide.
      // theta_dbl -= mdivide_left(gradient, system.cond_hessian(theta_dbl));

      // Check solution is a root of the gradient
      if (gradient.norm() <= tol) break;
    }

    return theta_dbl;
  }

  /**
   * lgp Newton solver.
   * In this instantiation, phi is of type var.
   * The initial guess can be passed as a parameter or fixed data.
   */
  template <typename T1, typename T2>
  Eigen::Matrix<T2, Eigen::Dynamic, 1> lgp_newton_solver(
    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta_0,  // initial guess
    const lgp_conditional_system<T2>& system,
    double tol = 1e-6,
    long int max_num_steps = 100) {  // NOLINT(runtime/int)

    lgp_conditional_system<double> 
      system_dbl(value_of(system.get_phi()),
                 system.get_n_samples(), 
                 system.get_sums());

    Eigen::VectorXd theta_dbl 
      = lgp_newton_solver(value_of(theta_0), system_dbl, tol, max_num_steps);

    // construct vari
    lgp_newton_solver_vari<T2>* vi0
      = new lgp_newton_solver_vari<T2>(system.get_phi(), system_dbl, theta_dbl);

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta(theta_dbl.size());
    theta(0) = var(vi0->theta_[0]);
    for (int i = 1; i < theta_dbl.size(); i++)
      theta(i) = var(vi0->theta_[i]);

    return theta;
  }

  /**
   * Wrapper of the lgp solver for use in the Stan language. This is
   * because we cannot pass an object of type lgp_conditional_system.
   * Note the wrapper automatically handles cases where phi is double
   * or var.
   */
  template <typename T1, typename T2>
  Eigen::Matrix<T2, Eigen::Dynamic, 1> lgp_newton_solver(
    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta_0,  // initial guess
    const T2& phi,
    const std::vector<int>& n_samples,
    const std::vector<int>& sums) {

    return lgp_newton_solver(theta_0, 
                             lgp_conditional_system<T2>(phi, 
                                                        to_vector(n_samples), 
                                                        to_vector(sums)));
  }

}  // namespace math
}  // namespace stan

#endif
