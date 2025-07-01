#include <stan/math.hpp>
#include <stan/math/mix.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <vector>

namespace {
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
  using stan::math::laplace_latent_bernoulli_logit_rng;
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
  Eigen::VectorXd mean(2);
  mean << 0, 0;
  std::vector<double> d0;
  std::vector<int> di0;
  std::vector<Eigen::VectorXd> x_dummy;
  boost::random::mt19937 rng;
  rng.seed(1954);
  Eigen::MatrixXd theta_pred = laplace_latent_bernoulli_logit_rng(
      sums, n_samples, mean, diagonal_kernel_functor{},
      std::forward_as_tuple(phi(0), phi(1)), rng, nullptr);

  // Compute exact mean and covariance
  Eigen::VectorXd theta_root
      = algebra_solver(stationary_point{}, theta_0, phi, d0, di0);
  Eigen::MatrixXd K_laplace = laplace_covariance(theta_root, phi);

  rng.seed(1954);
  Eigen::MatrixXd theta_benchmark
      = multi_normal_rng(theta_root, K_laplace, rng);

  double tol = 1e-3;
  EXPECT_NEAR(theta_benchmark(0), theta_pred(0), tol);
  EXPECT_NEAR(theta_benchmark(1), theta_pred(1), tol);
}
}  // namespace
