#include <stan/math/rev/core.hpp>
#include <stan/math/laplace/lgp_conditional_system.hpp>
#include <stan/math/laplace/lgp_newton_solver.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>

// functor for algebraic solver
struct lgp_functor {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename stan::return_type<T0, T1>::type,
                       Eigen::Dynamic, 1>
  operator ()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
              const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
              const std::vector<double>& dat,
              const std::vector<int>& dat_int,
              std::ostream* pstream__) const {
    typedef typename stan::return_type<T0, T1>::type scalar;
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> fgrad;
    int dim_theta = 2;

    Eigen::VectorXd n_samples(dim_theta);
    n_samples(0) = dat[0];
    n_samples(1) = dat[1];

    Eigen::VectorXd sums(dim_theta);
    sums(0) = dat[2];
    sums(1) = dat[3];

    return sums - stan::math::elt_multiply(n_samples,
                                           stan::math::exp(theta))
      - theta / (phi(0) * phi(0));
  }
};

TEST(Laplace, lgp_conditional_system) {
  // R code to generate data and results are in
  // make_data.R with seed 1954.
  using stan::math::lgp_conditional_system;

  int dim_theta = 2;
  Eigen::VectorXd theta(dim_theta);
  theta << -0.7268203, 1.3347728;
  double phi = 2;
  Eigen::VectorXd n_samples(dim_theta);
  for (int i = 0; i < dim_theta; i++) n_samples(i) = 5;
  Eigen::VectorXd sums(dim_theta);
  sums << 3, 10;

  lgp_conditional_system<double> system(phi, n_samples, sums);

  // Test evaluation of the density
  EXPECT_FLOAT_EQ(6.595955, system.log_density(theta));
  
  // Test evaluation of the gradient
  Eigen::VectorXd cond_grad = system.cond_gradient(theta);
  EXPECT_FLOAT_EQ(0.7644863, cond_grad(0));
  EXPECT_FLOAT_EQ(-9.3293567, cond_grad(1));
  
  // Test evaluation of the hessian
  Eigen::VectorXd cond_hessian = system.cond_hessian(theta);
  EXPECT_FLOAT_EQ(-2.667219, cond_hessian(0));
  EXPECT_FLOAT_EQ(-19.245664, cond_hessian(1));
}

TEST(laplace, lgp_newton_solver) {
  // Compute mode using Powell's solver and use it as
  // a benchmark.
  using stan::math::var;
  using stan::math::algebra_solver;
  using stan::math::value_of;
  using stan::math::lgp_conditional_system;
  
  int dim_theta = 2;
  Eigen::VectorXd theta_0(dim_theta);  // initial guess
  theta_0 << -1, 1;
  Eigen::VectorXd n_samples(dim_theta);
  n_samples << 5, 5;
  Eigen::VectorXd sums(dim_theta);
  sums << 3, 10;

  std::vector<double> dat(2 * dim_theta);
  dat[0] = n_samples(0);
  dat[1] = n_samples(1);
  dat[2] = sums(0);
  dat[3] = sums(1);

  std::vector<int> dummy_int;
  std::vector<double> powell_solution(2);
  std::vector<double> solver_gradient(2);

  for (int k = 0; k < dim_theta; k++) {
    var phi = 2;
    Eigen::Matrix<var, Eigen::Dynamic, 1> phi_v(1);
    phi_v << phi;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta 
      = algebra_solver(lgp_functor(), theta_0, phi_v, dat, dummy_int);
    // returns -0.472228 0.6761  // CHECK - should this return true theta?

    AVEC parameters = createAVEC(phi);
    VEC g;
    theta(k).grad(parameters, g);
    solver_gradient[k] = g[0];
    powell_solution[k] = value_of(theta(k));  // a little redundant...
  }
  // solver gradient is -0.035052 0.0167667
  
  double phi_dbl = 2;
  
  // Evaluate gradients using finite differentiation
  double diff = 1e-6;
  double phi_min = phi_dbl - diff;
  double phi_max = phi_dbl + diff;
  lgp_conditional_system<double> system_min(phi_min, n_samples, sums);
  lgp_conditional_system<double> system_max(phi_max, n_samples, sums);
  Eigen::VectorXd finite_diff(2);
  finite_diff = (lgp_newton_solver(theta_0, system_max)
    - lgp_newton_solver(theta_0, system_min)) / (2 * diff);
  // result is -0.035052 0.0167667, which agrees with Powell's result.
  

  // Test newton solver with only double argument.
  lgp_conditional_system<double> system(phi_dbl, n_samples, sums);
  Eigen::VectorXd theta_dbl = lgp_newton_solver(theta_0, system);

  EXPECT_FLOAT_EQ(powell_solution[0], theta_dbl(0));
  EXPECT_FLOAT_EQ(powell_solution[1], theta_dbl(1));

  // Test newton solver with only double argument, and line search method.
  double tol = 1e-3;
  int max_num_steps = 100;
  bool line_search = true;
  theta_dbl = lgp_newton_solver(theta_0, system,
                                tol, max_num_steps, line_search);
  EXPECT_FLOAT_EQ(powell_solution[0], theta_dbl(0));
  EXPECT_FLOAT_EQ(powell_solution[1], theta_dbl(1));

  // Test lgp_conditional_system computes the correct gradient
  Eigen::VectorXd solution(2);
  solution << -0.472228, 0.6761;
  EXPECT_FLOAT_EQ(solver_gradient[0], system.solver_gradient(solution)(0));
  EXPECT_FLOAT_EQ(solver_gradient[1], system.solver_gradient(solution)(1));

  // Test newton solver with var argument.
  for (int k = 0; k < dim_theta; k++) {
    var phi = 2;
    lgp_conditional_system<var> system_v(phi, n_samples, sums);

    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> theta
      = lgp_newton_solver(theta_0, system_v);

    AVEC parameters = createAVEC(phi);
    VEC g;
    theta(k).grad(parameters, g);
    EXPECT_FLOAT_EQ(solver_gradient[k], g[0]);
    
    // check solution (redundant)
    EXPECT_FLOAT_EQ(powell_solution[0], value_of(theta(0)));
    EXPECT_FLOAT_EQ(powell_solution[1], value_of(theta(1)));
  }

  std::vector<int> n_samples_array(2);
  n_samples_array[0] = 5;
  n_samples_array[1] = 5;
  
  std::vector<int> sums_array(2);
  sums_array[0] = 3;
  sums_array[1] = 10;
  
  // Test newton solver wrapper
  for (int k = 0; k < dim_theta; k++) {
    var phi = 2;
    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> theta
      = lgp_newton_solver(theta_0, phi, n_samples_array, sums_array);
    
    AVEC parameters = createAVEC(phi);
    VEC g;
    theta(k).grad(parameters, g);
    EXPECT_FLOAT_EQ(solver_gradient[k], g[0]);
    
    // check solution (redundant)
    EXPECT_FLOAT_EQ(powell_solution[0], value_of(theta(0)));
    EXPECT_FLOAT_EQ(powell_solution[1], value_of(theta(1)));
  }
  
  // Test newton solver wrapper with tuning parameters
  int is_line_search = 0;
  for (int k = 0; k < dim_theta; k++) {
    var phi = 2;
    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> theta
    = lgp_newton_solver(theta_0, phi, n_samples_array, sums_array,
                        tol, max_num_steps, is_line_search);

    AVEC parameters = createAVEC(phi);
    VEC g;
    theta(k).grad(parameters, g);
    EXPECT_FLOAT_EQ(solver_gradient[k], g[0]);

    // check solution (redundant)
    EXPECT_FLOAT_EQ(powell_solution[0], value_of(theta(0)));
    EXPECT_FLOAT_EQ(powell_solution[1], value_of(theta(1)));
  }
  
  // Repeat test, this time using line_search
  is_line_search = 1;
  for (int k = 0; k < dim_theta; k++) {
    var phi = 2;
    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> theta
      = lgp_newton_solver(theta_0, phi, n_samples_array, sums_array,
                          tol, max_num_steps, is_line_search);

    AVEC parameters = createAVEC(phi);
    VEC g;
    theta(k).grad(parameters, g);
    EXPECT_FLOAT_EQ(solver_gradient[k], g[0]);

    // check solution (redundant)
    EXPECT_FLOAT_EQ(powell_solution[0], value_of(theta(0)));
    EXPECT_FLOAT_EQ(powell_solution[1], value_of(theta(1)));
  }
}

