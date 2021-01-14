#include <stan/math.hpp>
#include <stan/math/laplace/laplace_likelihood.hpp>
#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

TEST(laplace, likelihood_differentiation) {
  using stan::math::diff_student_t;
  using stan::math::var;

  double test_tolerance = 2e-4;

  Eigen::VectorXd theta(2);
  theta << 0, 0;  // -2.45809, -3.6127;
  Eigen::VectorXd eta(2);
  eta << 1.2, 1;  // nu, sigma

  Eigen::VectorXd y(2);
  y << -2.655953, -4.2044;
  std::vector<int> y_index(2);
  y_index[0] = 0;
  y_index[1] = 1;

  // Eigen::Matrix<var, Eigen::Dynamic, 1> theta_v = theta;
  diff_student_t diff_functor(y, y_index);

  double log_density = diff_functor.log_likelihood(theta, eta);

  // benchmark against R
  EXPECT_NEAR(-7.375673, log_density, test_tolerance);




  // diff_logistic_log diff_functor(n_samples, y);
  // double log_density = diff_functor.log_likelihood(theta);
  // Eigen::VectorXd gradient, hessian;
  // diff_functor.diff(theta, gradient, hessian);
  // Eigen::VectorXd third_tensor = diff_functor.third_diff(theta);
  //
  // EXPECT_NEAR(-2.566843, log_density, test_tolerance);

  // finite diff calculations for first-order derivatives
  // double diff = 1e-12;
  // Eigen::VectorXd theta_1u = theta;
  // Eigen::VectorXd theta_1l = theta;
  // Eigen::VectorXd theta_2u = theta;
  // Eigen::VectorXd theta_2l = theta;
  // theta_1u(0) = theta(0) + diff;
  // theta_1l(0) = theta(0) - diff;
  // theta_2u(1) = theta(1) + diff;
  // theta_2l(1) = theta(1) - diff;
  // double diff_1 = (diff_functor.log_likelihood(theta_1u)
  //                    - diff_functor.log_likelihood(theta_1l)) / (2 * diff);
  // double diff_2 = (diff_functor.log_likelihood(theta_2u)
  //                    - diff_functor.log_likelihood(theta_2l)) / (2 * diff);

  // EXPECT_NEAR(diff_1, gradient(0), test_tolerance);
  // EXPECT_NEAR(diff_2, gradient(1), test_tolerance);
  //
  // // finite diff calculation for second-order derivatives
  // Eigen::VectorXd gradient_1u, gradient_1l, hessian_1u, hessian_1l,
  // gradient_2u, gradient_2l, hessian_2u, hessian_2l;
  // diff_functor.diff(theta_1u, gradient_1u, hessian_1u);
  // diff_functor.diff(theta_1l, gradient_1l, hessian_1l);
  // diff_functor.diff(theta_2u, gradient_2u, hessian_2u);
  // diff_functor.diff(theta_2l, gradient_2l, hessian_2l);
  //
  // double diff_grad_1 = (gradient_1u(0) - gradient_1l(0)) / (2 * diff);
  // double diff_grad_2 = (gradient_2u(1) - gradient_2l(1)) / (2 * diff);
  //
  // EXPECT_NEAR(diff_grad_1, hessian(0), test_tolerance);
  // EXPECT_NEAR(diff_grad_2, hessian(1), test_tolerance);
  //
  // // finite diff calculation for third-order derivatives
  // double diff_hess_1 = (hessian_1u(0) - hessian_1l(0)) / (2 * diff);
  // double diff_hess_2 = (hessian_2u(1) - hessian_2l(1)) / (2 * diff);
  //
  // EXPECT_NEAR(diff_hess_1, third_tensor(0), test_tolerance);
  // EXPECT_NEAR(diff_hess_2, third_tensor(1), test_tolerance);
}
