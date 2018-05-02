#include <stan/math/rev/core.hpp>
#include <stan/math/laplace/lgp_conditional_system.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>

TEST(Laplace, conditional_system) {
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

  // Test evaluation of the gradient
  Eigen::VectorXd cond_grad = system.cond_gradient(theta);

  EXPECT_FLOAT_EQ(0.7644863, cond_grad(0));
  EXPECT_FLOAT_EQ(-9.3293567, cond_grad(1));
  
  // Test evaluation of the hessian
  Eigen::VectorXd cond_hessian = system.cond_hessian(theta);

  EXPECT_FLOAT_EQ(-2.667219, cond_hessian(0));
  EXPECT_FLOAT_EQ(-19.245664, cond_hessian(1));
}
