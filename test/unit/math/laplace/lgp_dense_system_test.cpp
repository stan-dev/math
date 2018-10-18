#include <stan/math/rev/core.hpp>
#include <stan/math/laplace/lgp_conditional_system.hpp>
#include <stan/math/laplace/lgp_newton_solver.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>

TEST(Laplace, lgp_conditional_system) {
  // // R code to generate data and results are in
  // // draft_make_data.R with seed 1954.  
  // using Eigen::VectorXd 
  // using stan::math::lgp_conditional_system;
  // 
  // int dim_theta = 2;
  // VectorXd theta(dim_theta);
  // theta << -0.7268203, 1.3347728;
  // VectorXd phi(2);  // variance and correlation
  // phi << 1, 0;
  // 
  // VectorXd n_samples(dim_theta);
  // for (int i = 0; i < dim_theta; i++) n_samples(i) = 5;
  // Eigen::VectorXd sums(dim_theta);
  // sums << 3, 10;
  // 
  // lgp_conditional_system<double> system(phi, n_samples, sums);
  // 
  // // Test evaluation of the density
  // EXPECT_FLOAT_EQ(6.595955, system.log_density(theta));
}
