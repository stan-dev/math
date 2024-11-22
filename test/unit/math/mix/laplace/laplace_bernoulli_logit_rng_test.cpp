#include <stan/math.hpp>
#include <stan/math/mix.hpp>

#include <Eigen/Sparse>

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
    z(0) = 1 / (1 + exp(theta(0))) - theta(0) / (parms(0) * parms(0));
    z(1) = -1 / (1 + exp(-theta(1))) - theta(1) / (parms(1) * parms(1));
    return z;
  }
};

struct diagonal_kernel_functor {
  template <typename T1, typename T2>
  auto operator()(const T1& arg1, const T2& arg2,
                  std::ostream* msgs = nullptr) const {
    Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, Eigen::Dynamic>
        K(2, 2);
    K(0, 0) = arg1 * arg1;
    K(1, 1) = arg2 * arg2;
    K(0, 1) = 0;
    K(1, 0) = 0;
    return K;
  }
};

template <typename T1, typename T2>
Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> laplace_covariance(
    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
    const Eigen::Matrix<T2, Eigen::Dynamic, 1>& phi) {
  using stan::math::exp;
  using stan::math::square;
  Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> K(2, 2);
  K(0, 0)
      = -1
        / (-1 / (phi(0) * phi(0)) - exp(theta(0)) / square(1 + exp(theta(0))));
  K(1, 1) = -1
            / (-1 / (phi(1) * phi(1))
               - exp(-theta(1)) / square(1 + exp(-theta(1))));
  K(0, 1) = 0;
  K(1, 0) = 0;
  return K;
}

TEST(laplace_bernoulli_logit_rng, two_dim_diag) {
  using stan::math::algebra_solver;
  using stan::math::laplace_marginal_bernoulli_logit_rng;
  using stan::math::multi_normal_rng;
  using stan::math::sqrt;
  using stan::math::square;

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
  Eigen::MatrixXd theta_pred = laplace_marginal_bernoulli_logit_rng(
      sums, n_samples, theta_0, covariance_function, std::make_tuple(),
      std::make_tuple(), rng, nullptr, std::forward_as_tuple(phi(0), phi(1)));

  // Compute exact mean and covariance
  Eigen::VectorXd theta_root
      = algebra_solver(stationary_point(), theta_0, phi, d0, di0);
  Eigen::MatrixXd K_laplace = laplace_covariance(theta_root, phi);

  rng.seed(1954);
  Eigen::MatrixXd theta_benchmark
      = multi_normal_rng(theta_root, K_laplace, rng);

  double tol = 1e-3;
  EXPECT_NEAR(theta_benchmark(0), theta_pred(0), tol);
  EXPECT_NEAR(theta_benchmark(1), theta_pred(1), tol);
}

// TEST(laplace, bernoulli_logit_basic_rng) {
//   // make sure the right covariance function is computed
//   // and compare results.
//   using stan::math::diff_poisson_log;
//   using stan::math::laplace_bernoulli_logit_rng;
//   using stan::math::laplace_poisson_log_rng;
//   using stan::math::laplace_rng;
//
//   using stan::math::algebra_solver;
//   using stan::math::diag_matrix;
//   using stan::math::diag_pre_multiply;
//   using stan::math::inv;
//   using stan::math::mdivide_left_tri;
//   using stan::math::square;
//   using stan::math::to_vector;
//   using stan::math::value_of;
//
//   Eigen::VectorXd theta_0(2);
//   theta_0 << 1, 1;
//   Eigen::VectorXd sigma(2);
//   sigma << 3, 2;
//   std::vector<int> n_samples = {1, 1};
//   std::vector<int> sums = {1, 0};
//
//   diff_poisson_log laplace_likelihood(to_vector(n_samples), to_vector(sums));
//   std::vector<double> d0;
//   std::vector<int> di0;
//
//   // Method 1: brute force and straightforward
//   Eigen::VectorXd theta_root
//       = algebra_solver(stationary_point(), theta_0, sigma, d0, di0);
//
//   Eigen::VectorXd gradient, eta_dummy;
//   Eigen::SparseMatrix<double> W_sparse;
//   laplace_likelihood.diff(theta_root, eta_dummy, gradient, W_sparse);
//   Eigen::MatrixXd W = -W_sparse;
//   diagonal_kernel_functor covariance_function;
//   std::vector<Eigen::VectorXd> x_dummy;
//   Eigen::MatrixXd x_dummay_mat;
//   Eigen::MatrixXd K = covariance_function(sigma, x_dummy, d0, di0, 0);
//
//   std::cout << "K (brute force): " << std::endl
//             << (K.inverse() + W).inverse() << std::endl
//             << std::endl;
//
//   // Method 2: Vectorized R&W method
//   double tolerance = 1e-6;
//   int max_num_steps = 100;
//
//   // First find the mode using the custom Newton step
//   Eigen::MatrixXd covariance;
//   Eigen::VectorXd theta;
//   // Eigen::VectorXd W_root;
//   Eigen::SparseMatrix<double> W_r;
//   Eigen::MatrixXd L;
//   Eigen::MatrixXd K_root;
//   Eigen::VectorXd theta0_val = value_of(theta_0);
//   {
//     Eigen::VectorXd a;
//     Eigen::VectorXd l_grad;
//     Eigen::PartialPivLU<Eigen::MatrixXd> LU_dummy;
//     double marginal_density = laplace_marginal_density(
//         laplace_likelihood, covariance_function, sigma, eta_dummy, x_dummy,
//         d0, di0, covariance, theta, W_r, L, a, l_grad, LU_dummy, K_root,
//         theta0_val, 0, tolerance, max_num_steps);
//   }
//
//   Eigen::VectorXd W_root(theta.size());
//   for (int i = 0; i < theta.size(); i++)
//     W_root(i) = W_r.coeff(i, i);
//   Eigen::MatrixXd V;
//   V = mdivide_left_tri<Eigen::Lower>(L, diag_pre_multiply(W_root,
//   covariance)); std::cout << "K (method 1): " << std::endl
//             << covariance - V.transpose() * V << std::endl
//             << std::endl;
//
//   // Method 3: Modified R&W method
//   Eigen::VectorXd W_root_inv = inv(W_root);
//   Eigen::MatrixXd V_dec
//       = mdivide_left_tri<Eigen::Lower>(L, diag_matrix(W_root_inv));
//   std::cout << "K (method 2): " << std::endl
//             << -V_dec.transpose() * V_dec + diag_matrix(square(W_root_inv))
//             << std::endl
//             << std::endl;
//
//   // Check calls to rng functions compile
//   boost::random::mt19937 rng;
//   Eigen::MatrixXd theta_pred
//       = laplace_base_rng(laplace_likelihood, covariance_function, sigma,
//       eta_dummy,
//                          x_dummy, x_dummy, d0, di0, theta_0, rng);
//
//   theta_pred
//       = laplace_bernoulli_logit_rng(sums, n_samples, covariance_function,
//       sigma,
//                                     x_dummay_mat, d0, di0, theta_0, rng);
//
//   // Bonus: make the distribution with a poisson rng also runs.
//   theta_pred = laplace_poisson_log_rng(sums, n_samples, covariance_function,
//                                        sigma, x_dummy, d0, di0, theta_0,
//                                        rng);
// }
