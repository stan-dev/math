#include <stan/math/rev/core.hpp>
#include <stan/math/laplace/lgp_dense_system.hpp>
#include <stan/math/laplace/lgp_dense_newton_solver.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>



TEST(Laplace, lgp_dense_system) {
  // Since the dense function generalizes the diagonal case (semi-misleadingly
  // labeled conditional), we reproduce the result from the old function on
  // this new case.
  // One caveat: the global parameter in the diagonal case was the standard
  // deviation, but here it is easier to think in terms of variance.
  // Hence the first element of phi is not 2, but 4.
  // R code to generate results are in make_data.R with seed 1954.
  using stan::math::lgp_dense_system;
  using std::cout;
  using std::endl;

  int dim_theta = 2;
  Eigen::VectorXd theta(dim_theta);
  theta << -0.7268203, 1.3347728;

  Eigen::VectorXd phi(2);
  phi << 4, 0;  // global parameters: sigma and rho

  Eigen::VectorXd n_samples(dim_theta);
  for (int i = 0; i < dim_theta; i++) n_samples(i) = 5;
  Eigen::VectorXd sums(dim_theta);
  sums << 3, 10;

  lgp_dense_system<double> system(phi, theta, n_samples, sums);

  // Test evaulation of log densiy (up to a const)
  EXPECT_FLOAT_EQ(6.595955, system.log_density(theta));

  // Test evaluation of the gradient
  Eigen::VectorXd cond_grad = system.cond_gradient(theta);
  EXPECT_FLOAT_EQ(0.7644863, cond_grad(0));
  EXPECT_FLOAT_EQ(-9.3293567, cond_grad(1));

  Eigen::MatrixXd cond_hessian = system.cond_hessian(theta);
  EXPECT_FLOAT_EQ(-2.667219, cond_hessian(0, 0));
  EXPECT_FLOAT_EQ(0, cond_hessian(0, 1));
  EXPECT_FLOAT_EQ(0, cond_hessian(1, 0));
  EXPECT_FLOAT_EQ(-19.245664, cond_hessian(1, 1));
  
  // cout << system.cond_hessian(theta) << endl << endl;
}
