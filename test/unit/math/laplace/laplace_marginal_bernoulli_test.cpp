#include <stan/math.hpp>
#include <stan/math/laplace/laplace_likelihood.hpp>
#include <stan/math/laplace/laplace_marginal_bernoulli.hpp>

#include <test/unit/math/laplace/lgp_utility.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

TEST(laplace, likelihood_differentiation) {
  using stan::math::diff_logistic_log;
  using stan::math::var;

  double test_tolerance = 2e-4;

  Eigen::VectorXd theta(2);
  theta << -2.45809, -3.6127;
  Eigen::VectorXd y(2), n_samples(2);
  y << 1, 0;
  n_samples << 1, 1;
  Eigen::Matrix<var, Eigen::Dynamic, 1> theta_v = theta;

  diff_logistic_log diff_functor(n_samples, y);
  double log_density = diff_functor.log_likelihood(theta);
  Eigen::VectorXd gradient, hessian;
  diff_functor.diff(theta, gradient, hessian);
  Eigen::VectorXd third_tensor = diff_functor.third_diff(theta);

  EXPECT_NEAR(-2.566843, log_density, test_tolerance);

  // finite diff calculations for first-order derivatives
  double diff = 1e-12;
  Eigen::VectorXd theta_1u = theta;
  Eigen::VectorXd theta_1l = theta;
  Eigen::VectorXd theta_2u = theta;
  Eigen::VectorXd theta_2l = theta;
  theta_1u(0) = theta(0) + diff;
  theta_1l(0) = theta(0) - diff;
  theta_2u(1) = theta(1) + diff;
  theta_2l(1) = theta(1) - diff;
  double diff_1 = (diff_functor.log_likelihood(theta_1u)
                     - diff_functor.log_likelihood(theta_1l)) / (2 * diff);
  double diff_2 = (diff_functor.log_likelihood(theta_2u)
                     - diff_functor.log_likelihood(theta_2l)) / (2 * diff);

  EXPECT_NEAR(diff_1, gradient(0), test_tolerance);
  EXPECT_NEAR(diff_2, gradient(1), test_tolerance);

  // finite diff calculation for second-order derivatives
  Eigen::VectorXd gradient_1u, gradient_1l, hessian_1u, hessian_1l,
  gradient_2u, gradient_2l, hessian_2u, hessian_2l;
  diff_functor.diff(theta_1u, gradient_1u, hessian_1u);
  diff_functor.diff(theta_1l, gradient_1l, hessian_1l);
  diff_functor.diff(theta_2u, gradient_2u, hessian_2u);
  diff_functor.diff(theta_2l, gradient_2l, hessian_2l);

  double diff_grad_1 = (gradient_1u(0) - gradient_1l(0)) / (2 * diff);
  double diff_grad_2 = (gradient_2u(1) - gradient_2l(1)) / (2 * diff);

  EXPECT_NEAR(diff_grad_1, hessian(0), test_tolerance);
  EXPECT_NEAR(diff_grad_2, hessian(1), test_tolerance);

  // finite diff calculation for third-order derivatives
  double diff_hess_1 = (hessian_1u(0) - hessian_1l(0)) / (2 * diff);
  double diff_hess_2 = (hessian_2u(1) - hessian_2l(1)) / (2 * diff);

  EXPECT_NEAR(diff_hess_1, third_tensor(0), test_tolerance);
  EXPECT_NEAR(diff_hess_2, third_tensor(1), test_tolerance);
}

TEST(laplace, logistic_lgm_dim500) {
  using stan::math::var;
  using stan::math::to_vector;
  using stan::math::diff_logistic_log;
  using stan::math::sqr_exp_kernel_functor;

  int dim_theta = 500;
  int n_observations = 500;
  std::string data_directory = "test/unit/math/laplace/aki_synth_data/";
  std::vector<double> x1(dim_theta), x2(dim_theta);
  std::vector<int> y(n_observations);
  read_in_data(dim_theta, n_observations, data_directory, x1, x2, y);

  // Look a some of the data.
  // std::cout << "x_1: " << x1[0] << " " << x2[0] << std::endl
  //           << "x_2: " << x1[1] << " " << x2[1] << std::endl
  //           << "y_1: " << y[0] << " y_2: " << y[1] << std::endl; 

  int dim_x = 2;
  std::vector<Eigen::VectorXd> x(dim_theta);
  for (int i = 0; i < dim_theta; i++) {
    Eigen::VectorXd coordinate(dim_x);
    coordinate << x1[i], x2[i];
    x[i] = coordinate;
  }
  std::vector<int> n_samples = stan::math::rep_array(1, dim_theta);

  Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);

  Eigen::VectorXd theta_laplace, W_root, a, l_grad;
  Eigen::MatrixXd L, covariance;
  std::vector<double> delta;
  std::vector<int> delta_int;

  // CASE 1: phi is passed as a double.
  Eigen::VectorXd phi(2);
  phi << 1.6, 1;  // standard deviation, length scale

  auto start_optimization = std::chrono::system_clock::now();

  double marginal_density
    = laplace_marginal_density(
      diff_logistic_log(to_vector(n_samples), to_vector(y)),
      sqr_exp_kernel_functor(),
      phi, x, delta, delta_int,
      covariance, theta_laplace, W_root, L, a, l_grad,
      theta_0, 1e-3, 100);

  auto end_optimization = std::chrono::system_clock::now();
  std::chrono::duration<double>
    elapsed_time_optimization = end_optimization - start_optimization;

  std::cout << "LAPLACE MARGINAL FOR DOUBLE: " << std::endl
            << "density: " << marginal_density << std::endl
            << "time: " << elapsed_time_optimization.count()
            << std::endl << std::endl;

  // Expected output
  // density: -195.368
  // time: 0.059645

  // CASE 2: phi is passed as a var
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi_v2 = phi;

  start_optimization = std::chrono::system_clock::now();
  var marginal_density_v
    = laplace_marginal_density(
        diff_logistic_log(to_vector(n_samples), to_vector(y)),
        sqr_exp_kernel_functor(),
        phi_v2, x, delta, delta_int,
        theta_0, 1e-3, 100);

  VEC g2;
  AVEC parm_vec2 = createAVEC(phi_v2(0), phi_v2(1));
  marginal_density_v.grad(parm_vec2, g2);

  end_optimization = std::chrono::system_clock::now();
  elapsed_time_optimization = end_optimization - start_optimization;

  std::cout << "LAPLACE MARGINAL AND VARI CLASS" << std::endl
            << "density: " << value_of(marginal_density_v) << std::endl
            << "autodiff grad: " << g2[0] << " " << g2[1]
            << std::endl
            << "total time: " << elapsed_time_optimization.count()
            << std::endl << std::endl;

  // EXPECTED
  // density: -195.368
  // autodiff grad: 21.9495 -32.5123
  // total time: 0.147897

  // TO DO -- get total time from GPStuff and do more comparisons.

  // CASE 3: use wrapper function and compare result.
  using stan::math::laplace_marginal_bernoulli;
  using stan::math::value_of;

  double marginal_density_v2
    = laplace_marginal_bernoulli(y, n_samples,
                                 phi, x, delta, delta_int,
                                 theta_0, 1e-3, 100);

  EXPECT_FLOAT_EQ(marginal_density, marginal_density_v2);

  marginal_density_v2
    = laplace_marginal_bernoulli(y, n_samples,
                                 sqr_exp_kernel_functor(),
                                 phi, x, delta, delta_int,
                                 theta_0, 1e-3, 100);

  EXPECT_FLOAT_EQ(marginal_density, marginal_density_v2);
}
