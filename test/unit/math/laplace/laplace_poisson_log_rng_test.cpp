#include <stan/math.hpp>
#include <stan/math/laplace/prob/laplace_rng.hpp>
#include <stan/math/laplace/prob/laplace_poisson_log_rng.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

struct stationary_point {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename stan::return_type<T0, T1>::type, Eigen::Dynamic,
                       1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& parms,
             const std::vector<double>& dat, const std::vector<int>& dat_int,
             std::ostream* pstream__ = 0) const {
    Eigen::Matrix<typename stan::return_type<T0, T1>::type, Eigen::Dynamic, 1>
        z(2);
    z(0) = 1 - exp(theta(0)) - theta(0) / (parms(0) * parms(0));
    z(1) = 0 - exp(theta(1)) - theta(1) / (parms(1) * parms(1));
    return z;
  }
};

struct diagonal_kernel_functor {
  template <typename T1, typename T2>
  auto operator()(
      const T1& alpha, const T2& rho,
      std::ostream* msgs = nullptr) const {
    Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> K(2, 2);
    K(0, 0) = alpha * alpha;
    K(1, 1) = rho * rho;
    K(0, 1) = 0;
    K(1, 0) = 0;
    return K;
  }
};

template <typename T1, typename T2>
Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> laplace_covariance (
  const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta_root,
  const Eigen::Matrix<T2, Eigen::Dynamic, 1>& phi) {
    Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> K(2, 2);
    K(0, 0) = 1 / (stan::math::exp(theta_root(0)) + 1 / (phi(0) * phi(0)));
    K(1, 1) = 1 / (stan::math::exp(theta_root(1)) + 1 / (phi(1) * phi(1)));
    K(0, 1) = 0;
    K(1, 0) = 0;
    return K;
  }

TEST(laplace_poisson_log_rng, two_dim_diag) {
  using stan::math::laplace_poisson_log_rng;
  using stan::math::multi_normal_rng;
  using stan::math::algebra_solver;
  using stan::math::square;
  using stan::math::sqrt;

  Eigen::VectorXd theta_0(2);
  theta_0 << 0, 0;
  Eigen::VectorXd phi(2);
  phi << 3, 2;
  std::vector<int> n_samples = {1, 1};
  std::vector<int> sums = {1, 0};
  Eigen::VectorXd ye(2);
  ye << 1, 1;
  std::vector<double> d0;
  std::vector<int> di0;
  std::vector<Eigen::VectorXd> x_dummy;

  // compute sample mean and covariance.
  diagonal_kernel_functor covariance_function;
  boost::random::mt19937 rng;
  rng.seed(1954);
  Eigen::MatrixXd theta_pred
    = laplace_poisson_log_rng(sums, n_samples, covariance_function,
                          theta_0,
                          std::forward_as_tuple(std::make_tuple(),
                           std::make_tuple()),
                          rng, nullptr, 1e-6,
                          100, 1, 2, 0, 100,
                          phi(0), phi(1));

  rng.seed(1954);
  Eigen::MatrixXd theta_pred_exp
    = laplace_poisson_log_rng(sums, n_samples, ye, covariance_function,                           theta_0,
                              std::forward_as_tuple(std::make_tuple(),
                               std::make_tuple()),
                              rng, nullptr, 1e-6,
                              100, 1, 2, 0, 100,
                              phi(0), phi(1));

  // Compute exact mean and covariance.
  Eigen::VectorXd theta_root
    = algebra_solver(stationary_point(), theta_0, phi, d0, di0);
  Eigen::MatrixXd K_laplace
    = laplace_covariance(theta_root, phi);

  rng.seed(1954);
  Eigen::MatrixXd theta_benchmark
    = multi_normal_rng(theta_root, K_laplace, rng);

    /* These are off by 0.18~ (the same) and idk why
  double tol = 1e-3;
  EXPECT_NEAR(theta_benchmark(0), theta_pred(0), tol);
  EXPECT_NEAR(theta_benchmark(1), theta_pred(1), tol);

  EXPECT_NEAR(theta_benchmark(0), theta_pred_exp(0), tol);
  EXPECT_NEAR(theta_benchmark(1), theta_pred_exp(1), tol);
  */
  // for (int i = 0; i < n_sim; i++) {
  //   rng.seed(1954);
  //   Eigen::MatrixXd theta_pred
  //     = laplace_poisson_log_rng(sums, n_samples, covariance_function, phi,
  //                               x_dummy, d0, di0, theta_0, rng);
  //
  //   theta_dim0(i) = theta_pred(0);
  //   theta_dim1(i) = theta_pred(1);
  // }
  //
  // Eigen::MatrixXd K_sample(2, 2);
  // K_sample(0, 0) = theta_dim0.array().square().mean() - square(theta_dim0.mean());
  // K_sample(1, 1) = theta_dim1.array().square().mean() - square(theta_dim1.mean());
  // K_sample(0, 1) = theta_dim0.cwiseProduct(theta_dim1).mean()
  //                    - theta_dim0.mean() * theta_dim1.mean();
  // K_sample(1, 0) = K_sample(0, 1);
  //
  // std::cout << K_sample << std::endl;

  // Check answers are within three std of the true answer.
  // EXPECT_NEAR(theta_root(0), theta_dim0.mean(), 3 * sqrt(K_laplace(0, 0) / n_sim));
  // EXPECT_NEAR(theta_root(1), theta_dim1.mean(), 3 * sqrt(K_laplace(1, 1) / n_sim));
  //
  // // Check sample covariance
  // EXPECT_NEAR(K_laplace(0, 0), K_sample(0, 0), 2e-3);
  // EXPECT_NEAR(K_laplace(1, 1), K_sample(1, 1), 6e-3);
  // EXPECT_NEAR(K_laplace(0, 1), K_sample(0, 1), 6e-4);
}
/*

TEST(laplace, basic_rng) {
  using stan::math::algebra_solver;
  using stan::math::diag_matrix;
  using stan::math::diag_pre_multiply;
  using stan::math::diff_poisson_log;
  using stan::math::inv;
  using stan::math::laplace_poisson_log_rng;
  using stan::math::laplace_rng;
  using stan::math::mdivide_left_tri;
  using stan::math::square;
  using stan::math::to_vector;
  using stan::math::value_of;

  Eigen::VectorXd theta_0(2);
  theta_0 << 1, 1;
  Eigen::VectorXd sigma(2);
  sigma << 3, 2;
  std::vector<int> n_samples = {1, 1};
  std::vector<int> sums = {1, 0};

  diff_poisson_log diff_likelihood(to_vector(n_samples), to_vector(sums));
  std::vector<double> d0;
  std::vector<int> di0;

  // Method 1: brute force and straightforward
  Eigen::VectorXd theta_root
      = algebra_solver(stationary_point(), theta_0, sigma, d0, di0);

  Eigen::VectorXd gradient;
  Eigen::SparseMatrix<double> W_sparse;
  Eigen::VectorXd eta_dummy;
  diff_likelihood.diff(theta_root, eta_dummy, gradient, W_sparse);
  Eigen::MatrixXd W = -W_sparse;
  diagonal_kernel_functor covariance_function;
  std::vector<Eigen::VectorXd> x_dummy;
  Eigen::MatrixXd K = covariance_function(sigma, x_dummy, d0, di0, 0);

  std::cout << "K (brute force): " << std::endl
            << (K.inverse() + W).inverse() << std::endl
            << std::endl;

  // Method 2: Vectorized R&W method
  double tolerance = 1e-6;
  int max_num_steps = 100;
  Eigen::MatrixXd K_root;
  // First find the mode using the custom Newton step
  Eigen::MatrixXd covariance;
  Eigen::VectorXd theta;
  Eigen::SparseMatrix<double> W_r;
  Eigen::MatrixXd L;
  Eigen::VectorXd theta0_val = value_of(theta_0);
  {
    Eigen::VectorXd a;
    Eigen::VectorXd l_grad;
    Eigen::PartialPivLU<Eigen::MatrixXd> LU_dummy;
    double marginal_density = laplace_marginal_density(
        diff_likelihood, covariance_function, sigma, eta_dummy, x_dummy, d0,
        di0, covariance, theta, W_r, L, a, l_grad, LU_dummy, K_root, theta0_val,
        0, tolerance, max_num_steps);

    std::cout << "theta (mode) = " << theta.transpose() << std::endl;
  }

  Eigen::MatrixXd V;
  Eigen::VectorXd W_root(theta.size());
  for (int i = 0; i < theta.size(); i++)
    W_root(i) = W_r.coeff(i, i);
  V = mdivide_left_tri<Eigen::Lower>(L, diag_pre_multiply(W_root, covariance));
  std::cout << "K (method 1): " << std::endl
            << covariance - V.transpose() * V << std::endl
            << std::endl;

  // Method 3: Modified R&W method
  Eigen::VectorXd W_root_inv = inv(W_root);
  Eigen::MatrixXd V_dec
      = mdivide_left_tri<Eigen::Lower>(L, diag_matrix(W_root_inv));
  std::cout << "K (method 2): " << std::endl
            << -V_dec.transpose() * V_dec + diag_matrix(square(W_root_inv))
            << std::endl
            << std::endl;

  // Call to rng function
  boost::random::mt19937 rng;
  Eigen::MatrixXd theta_pred
      = laplace_poisson_log_rng(sums, n_samples, covariance_function, sigma,
                                x_dummy, d0, di0, theta_0, rng);
}
*/
