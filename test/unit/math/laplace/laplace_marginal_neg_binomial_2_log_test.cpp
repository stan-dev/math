#include <stan/math.hpp>
#include <stan/math/laplace/laplace_likelihood.hpp>
#include <stan/math/laplace/laplace_marginal.hpp>
#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

TEST(laplace, likelihood_differentiation) {
  using stan::math::diff_neg_binomial_2_log;
  using stan::math::var;

  Eigen::VectorXd theta(2);
  theta << 1, 1;
  int n_theta = theta.size();
  Eigen::VectorXd eta(1);
  eta << 1.2;

  Eigen::VectorXd y(2);
  y << 0, 1;
  std::vector<int> y_index(2);
  y_index[0] = 0;
  y_index[1] = 1;

  // Eigen::Matrix<var, Eigen::Dynamic, 1> theta_v = theta;
  diff_neg_binomial_2_log diff_functor(y, y_index, n_theta);

  double log_density = diff_functor.log_likelihood(theta, eta);

  // benchmark against R
  EXPECT_FLOAT_EQ(-3.023328, log_density);

  Eigen::VectorXd gradient, hessian;
  diff_functor.diff(theta, eta, gradient, hessian);

  Eigen::VectorXd third_diff = diff_functor.third_diff(theta, eta);

  // Benchmark against finite diff
  double epsilon = 1e-6;
  Eigen::VectorXd theta_l0 = theta, theta_u0 = theta,
                  theta_l1 = theta, theta_u1 = theta;
  theta_u0(0) += epsilon;
  theta_l0(0) -= epsilon;
  theta_u1(1) += epsilon;
  theta_l1(1) -= epsilon;

  Eigen::VectorXd finite_gradient(2);
  finite_gradient(0) =
    (diff_functor.log_likelihood(theta_u0, eta)
      - diff_functor.log_likelihood(theta_l0, eta)) / (2 * epsilon);

  finite_gradient(1) =
    (diff_functor.log_likelihood(theta_u1, eta)
      - diff_functor.log_likelihood(theta_l1, eta)) / (2 * epsilon);

  Eigen::VectorXd gradient_l0, gradient_u0, gradient_l1, gradient_u1;
  Eigen::VectorXd hessian_l0, hessian_u0, hessian_l1, hessian_u1;
  Eigen::VectorXd hessian_dummy;
  diff_functor.diff(theta_l0, eta, gradient_l0, hessian_l0);
  diff_functor.diff(theta_u0, eta, gradient_u0, hessian_u0);
  diff_functor.diff(theta_l1, eta, gradient_l1, hessian_l1);
  diff_functor.diff(theta_u1, eta, gradient_u1, hessian_u1);

  Eigen::VectorXd finite_hessian(2);
  finite_hessian(0) = (gradient_u0 - gradient_l0)(0) / (2 * epsilon);
  finite_hessian(1) = (gradient_u1 - gradient_l1)(1) / (2 * epsilon);

  Eigen::VectorXd finite_third_diff(2);
  finite_third_diff(0) = (hessian_u0 - hessian_l0)(0) / (2 * epsilon);
  finite_third_diff(1) = (hessian_u1 - hessian_l1)(1) / (2 * epsilon);


  EXPECT_FLOAT_EQ(finite_gradient(0), gradient(0));
  EXPECT_FLOAT_EQ(finite_gradient(1), gradient(1));
  EXPECT_FLOAT_EQ(finite_hessian(0), hessian(0));
  EXPECT_FLOAT_EQ(finite_hessian(1), hessian(1));
  EXPECT_FLOAT_EQ(finite_third_diff(0), third_diff(0));
  EXPECT_FLOAT_EQ(finite_third_diff(1), third_diff(1));

  // derivatives wrt eta
  Eigen::VectorXd diff_eta = diff_functor.diff_eta(theta, eta);

  Eigen::VectorXd eta_l(1), eta_u(1);
  eta_l(0) = eta(0) - epsilon;
  eta_u(0) = eta(0) + epsilon;
  double finite_gradient_eta =
    (diff_functor.log_likelihood(theta, eta_u)
      - diff_functor.log_likelihood(theta, eta_l)) / (2 * epsilon);

  EXPECT_FLOAT_EQ(finite_gradient_eta,  diff_eta(0));

  Eigen::MatrixXd diff_theta_eta = diff_functor.diff_theta_eta(theta, eta);

  Eigen::VectorXd gradient_theta_l,
                  gradient_theta_u,
                  hessian_theta_u,
                  hessian_theta_l;

  diff_functor.diff(theta, eta_l, gradient_theta_l, hessian_theta_l);
  diff_functor.diff(theta, eta_u, gradient_theta_u, hessian_theta_u);
  Eigen::VectorXd finite_gradient_theta_eta
    = (gradient_theta_u - gradient_theta_l) / (2 * epsilon);

  EXPECT_FLOAT_EQ(finite_gradient_theta_eta(0), diff_theta_eta(0, 0));
  EXPECT_FLOAT_EQ(finite_gradient_theta_eta(1), diff_theta_eta(1, 0));

  Eigen::VectorXd W_root = (-hessian).cwiseSqrt();
  Eigen::MatrixXd diff2_theta_eta
    = diff_functor.diff2_theta_eta(theta, eta, W_root);

  Eigen::VectorXd finite_hessian_theta_eta
   = ((-hessian_theta_u).cwiseSqrt() - (-hessian_theta_l).cwiseSqrt())
       / (2 * epsilon);

  EXPECT_FLOAT_EQ(finite_hessian_theta_eta(0), diff2_theta_eta(0, 0));
  EXPECT_FLOAT_EQ(finite_hessian_theta_eta(1), diff2_theta_eta(1, 0));
}

TEST(laplace, neg_binomial_2_log_dbl) {
  using stan::math::to_vector;
  using stan::math::diff_neg_binomial_2_log;
  using stan::math::sqr_exp_kernel_functor;
  using stan::math::laplace_marginal_density;

  int dim_phi = 2, dim_eta = 1, dim_theta = 2;
  Eigen::VectorXd phi(dim_phi), eta(dim_eta), theta_0(dim_theta);
  phi << 1.6, 0.45;
  eta << 1;
  theta_0 << 0, 0;

  std::vector<Eigen::VectorXd> x(dim_theta);
  Eigen::VectorXd x_0(2), x_1(2);
  x_0 <<  0.05100797, 0.16086164;
  x_1 << -0.59823393, 0.98701425;
  x[0] = x_0;
  x[1] = x_1;

  std::vector<double> delta;
  std::vector<int> delta_int;
  std::vector<int> y_index = {1, 1};
  Eigen::VectorXd y = to_vector({1, 0});

  diff_neg_binomial_2_log diff_functor(y, y_index, dim_theta);
  stan::math::sqr_exp_kernel_functor K;

  double log_p = laplace_marginal_density(diff_functor, K, phi, eta, x, delta,
                                          delta_int, theta_0);

  // TODO: add test.


}
