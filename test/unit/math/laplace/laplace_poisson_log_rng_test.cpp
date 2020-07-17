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
  template<typename T0, typename T1>
  inline Eigen::Matrix<typename stan::return_type<T0, T1>::type,
                       Eigen::Dynamic, 1>
  operator() (const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
           const Eigen::Matrix<T1, Eigen::Dynamic, 1>& parms,
           const std::vector<double>& dat,
           const std::vector<int>& dat_int,
           std::ostream* pstream__ = 0) const {
    Eigen::Matrix<typename stan::return_type<T0, T1>::type,
                  Eigen::Dynamic, 1> z(2);
    z(0) = 1 - exp(theta(0)) - theta(0) / (parms(0) * parms(0));
    z(1) = - exp(theta(1)) - theta(1) / (parms(1) * parms(1));
    return z;
  }
};

struct diagonal_kernel_functor {
  template <typename T1, typename T2>
  Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>
  operator() (const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
              const T2& x,
              const std::vector<double>& delta,
              const std::vector<int>& delta_int,
              std::ostream* msgs = nullptr) const {
    Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> K(2, 2);
    K(0, 0) = phi(0) * phi(0);
    K(1, 1) = phi(1) * phi(1);
    K(0, 1) = 0;
    K(1, 0) = 0;
    return K;
  }
};

TEST(laplace, basic_rng) {
  using stan::math::algebra_solver;
  using stan::math::diff_poisson_log;
  using stan::math::to_vector;
  using stan::math::diag_matrix;
  using stan::math::laplace_rng;
  using stan::math::laplace_poisson_log_rng;
  using stan::math::value_of;
  using stan::math::mdivide_left_tri;
  using stan::math::diag_pre_multiply;
  using stan::math::inv;
  using stan::math::square;


  Eigen::VectorXd theta_0(2);
  theta_0 << 1, 1;
  Eigen::VectorXd sigma(2);
  sigma << 3, 2;
  std::vector<int> n_samples = {1, 1};
  std::vector<int> sums = {1, 0};

  diff_poisson_log diff_likelihood(to_vector(n_samples),
                                   to_vector(sums));
  std::vector<double> d0;
  std::vector<int> di0;


  // Method 1: brute force and straightforward
  Eigen::VectorXd theta_root
    = algebra_solver(stationary_point(),
                     theta_0, sigma, d0, di0);

  Eigen::VectorXd gradient, W;
  diff_likelihood.diff(theta_root, gradient, W);
  W = -W;
  diagonal_kernel_functor covariance_function;
  std::vector<Eigen::VectorXd> x_dummy;
  Eigen::MatrixXd K = covariance_function(sigma, x_dummy, d0, di0, 0);

  std::cout << "K (brute force): "
            << std::endl
            << (K.inverse() + diag_matrix(W)).inverse()
            << std::endl << std::endl;

  // Method 2: Vectorized R&W method
  double tolerance = 1e-6;
  int max_num_steps = 100;

  // First find the mode using the custom Newton step
  Eigen::MatrixXd covariance;
  Eigen::VectorXd theta;
  Eigen::VectorXd W_root;
  Eigen::MatrixXd L;
  {
    Eigen::VectorXd a;
    Eigen::VectorXd l_grad;
    double marginal_density
      = laplace_marginal_density(diff_likelihood,
                                 covariance_function,
                                 sigma, x_dummy, d0, di0,
                                 covariance, theta, W_root, L, a, l_grad,
                                 value_of(theta_0), 0,
                                 tolerance, max_num_steps);
  }

  Eigen::MatrixXd V;
  V = mdivide_left_tri<Eigen::Lower>(L,
        diag_pre_multiply(W_root, covariance));
  std::cout << "K (method 1): " << std::endl
            << covariance - V.transpose() * V << std::endl
            << std::endl;

  // Method 3: Modified R&W method
  Eigen::VectorXd W_root_inv = inv(W_root);
  Eigen::MatrixXd V_dec = mdivide_left_tri<Eigen::Lower>(L,
                            diag_matrix(W_root_inv));
  std::cout << "K (method 2): " << std::endl
            << - V_dec.transpose() * V_dec + diag_matrix(square(W_root_inv))
            << std::endl << std::endl;


  // Call to rng function
  boost::random::mt19937 rng;
  Eigen::MatrixXd theta_pred
   = laplace_rng(diff_likelihood, covariance_function,
                 sigma, x_dummy, d0, di0, theta_0, rng);

  theta_pred = laplace_poisson_log_rng(sums, n_samples, covariance_function,
                                       sigma, x_dummy, d0, di0, theta_0, rng);
}
