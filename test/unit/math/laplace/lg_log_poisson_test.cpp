#include <stan/math.hpp>
#include <stan/math/laplace/laplace_marginal_poisson.hpp>

#include <test/unit/math/laplace/lgp_utility.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>


TEST(laplace, likelihood_differentiation) {
  using stan::math::diff_poisson_log;
  using stan::math::var;
  using stan::math::to_vector;

  Eigen::VectorXd theta(2);
  theta << 1, 1;
  std::vector<int> n_samples = {1, 1};
  std::vector<int> sums = {1, 0};

  diff_poisson_log diff_functor(to_vector(n_samples),
                                to_vector(sums));
  double log_density = diff_functor.log_likelihood(theta);
  Eigen::VectorXd gradient, hessian;
  diff_functor.diff(theta, gradient, hessian);
  Eigen::VectorXd third_tensor = diff_functor.third_diff(theta);

  EXPECT_FLOAT_EQ(-4.436564, log_density);
  EXPECT_FLOAT_EQ(-1.718282, gradient(0));
  EXPECT_FLOAT_EQ(-2.718282, gradient(1));
  EXPECT_FLOAT_EQ(-2.718282, hessian(0));
  EXPECT_FLOAT_EQ(-2.718282, hessian(1));
  EXPECT_FLOAT_EQ(-2.718282, third_tensor(0));
  EXPECT_FLOAT_EQ(-2.718282, third_tensor(1));
}

TEST(laplace, lg_poisson_dim_2) {
  using stan::math::laplace_marginal_poisson;
  using stan::math::var;
  using stan::math::to_vector;
  using stan::math::value_of;

  int dim_phi = 2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi(dim_phi);
  phi << 1.6, 0.45;

  int dim_theta = 2;
  Eigen::VectorXd theta_0(dim_theta);
  theta_0 << 0, 0;

  int dim_x = 2;
  std::vector<Eigen::VectorXd> x(dim_theta);
  Eigen::VectorXd x_0(2);
  x_0 <<  0.05100797, 0.16086164;
  Eigen::VectorXd x_1(2);
  x_1 << -0.59823393, 0.98701425;
  x[0] = x_0;
  x[1] = x_1;

  std::vector<int> n_samples = {1, 1};
  std::vector<int> sums = {1, 0};

  squared_kernel_functor K;
  var target = laplace_marginal_poisson(theta_0, phi, x, n_samples, sums);

  // How to test this? The best way would be to generate a few
  // benchmarks using gpstuff.
  VEC g;
  AVEC parm_vec = createAVEC(phi(0), phi(1));
  target.grad(parm_vec, g);

  // finite diff test
  double diff = 1e-7;
  Eigen::VectorXd phi_dbl = value_of(phi);
  Eigen::VectorXd phi_1l = phi_dbl, phi_1u = phi_dbl,
    phi_2l = phi_dbl, phi_2u = phi_dbl;
  phi_1l(0) -= diff;
  phi_1u(0) += diff;
  phi_2l(1) -= diff;
  phi_2u(1) += diff;

  double target_1u = laplace_marginal_poisson(theta_0, phi_1u, x,
                                              n_samples, sums),
         target_1l = laplace_marginal_poisson(theta_0, phi_1l, x,
                                              n_samples, sums),
         target_2u = laplace_marginal_poisson(theta_0, phi_2u, x,
                                              n_samples, sums),
         target_2l = laplace_marginal_poisson(theta_0, phi_2l, x,
                                              n_samples, sums);

  VEC g_finite(dim_phi);
  g_finite[0] = (target_1u - target_1l) / (2 * diff);
  g_finite[1] = (target_2u - target_2l) / (2 * diff);

  double tol = 1.1e-4;
  EXPECT_NEAR(g_finite[0], g[0], tol);
  EXPECT_NEAR(g_finite[1], g[1], tol);
}
