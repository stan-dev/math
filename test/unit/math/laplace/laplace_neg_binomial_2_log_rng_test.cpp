#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <vector>

namespace {
struct stationary_point_nb {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename stan::return_type<T0, T1>::type, Eigen::Dynamic,
                       1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& parms,
             const std::vector<double>& dat, const std::vector<int>& dat_int,
             std::ostream* pstream__ = 0) const {
    using stan::math::exp;
    Eigen::Matrix<typename stan::return_type<T0, T1>::type, Eigen::Dynamic, 1>
        z(2);
    double eta = dat[0];
    Eigen::Matrix<T0, Eigen::Dynamic, 1> exp_theta = exp(theta);
    std::vector<int> y = {1, 0};

    z(0) = -exp_theta(0) * (y[0] + eta) / (eta + exp_theta(0)) + y[0]
           - theta(0) / parms(0);
    z(1) = -exp_theta(1) * (y[1] + eta) / (eta + exp_theta(1)) + y[1]
           - theta(1) / parms(1);
    // z(0) = 1 - (1 - eta) / (1 + eta * exp(theta(0))) - theta(0) / parms(0);
    // z(1) = 0 - (0 - eta) / (1 + eta * exp(theta(1))) - theta(1) / parms(1);
    return z;
  }
};

struct diagonal_kernel_nb_functor {
  template <typename T1, typename T2>
  auto operator()(const T1& alpha, const T2& rho,
                  std::ostream* msgs = nullptr) const {
    Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> K(2, 2);
    K(0, 0) = alpha;
    K(1, 1) = rho;
    K(0, 1) = 0;
    K(1, 0) = 0;
    return K;
  }
};

template <typename T1, typename T2>
Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> laplace_covariance_nb(
    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
    const Eigen::Matrix<T2, Eigen::Dynamic, 1>& phi, const double& eta) {
  using stan::math::exp;
  using stan::math::square;
  Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> K(2, 2);
  Eigen::Matrix<T1, Eigen::Dynamic, 1> exp_theta = exp(theta);
  std::vector<int> y = {1, 0};
  K(0, 0) = 1
            / ((eta * exp_theta(0) * (y[0] + eta) / square(eta + exp_theta(0)))
               + 1 / phi(0));
  K(1, 1) = 1
            / ((eta * exp_theta(1) * (y[1] + eta) / square(eta + exp_theta(1)))
               + 1 / phi(1));

  // K(0, 0) = 1 / (1 / phi(0) + (1 - eta) * exp(theta(0))
  //                                          / square(1 + eta *
  //                                          exp(theta(0))));
  // K(1, 1) = 1 / (1 / phi(1) + (0 - eta) * exp(theta(1))
  //                                          / square(1 + eta *
  //                                          exp(theta(1))));
  K(0, 1) = 0;
  K(1, 0) = 0;
  return K;
}

TEST(laplace_latent_neg_binomial_2_log_rng, count_two_dim_diag) {
  using stan::math::algebra_solver;
  using stan::math::laplace_latent_neg_binomial_2_log_rng;
  using stan::math::laplace_latent_tol_neg_binomial_2_log_rng;
  using stan::math::multi_normal_rng;
  using stan::math::sqrt;
  using stan::math::square;

  std::vector<int> y = {1, 0};
  std::vector<int> y_index = {1, 2};
  Eigen::VectorXd theta_0{{1, 1}};
  Eigen::VectorXd phi{{3, 2}};

  double eta = 2;
  std::vector<double> d0 = {eta};
  std::vector<int> di0;

  Eigen::VectorXd theta_root
      = algebra_solver(stationary_point_nb(), theta_0, phi, d0, di0);
  Eigen::MatrixXd K_laplace = laplace_covariance_nb(theta_root, phi, eta);

  boost::random::mt19937 rng;
  rng.seed(1954);
  Eigen::MatrixXd theta_benchmark
      = stan::math::multi_normal_rng(theta_root, K_laplace, rng);

  rng.seed(1954);
  Eigen::MatrixXd theta_pred = laplace_latent_neg_binomial_2_log_rng(
      y, y_index, eta, 0, diagonal_kernel_nb_functor{},
      std::forward_as_tuple(phi(0), phi(1)), rng, nullptr);

  double tol = 1e-3;
  EXPECT_NEAR(theta_benchmark(0), theta_pred(0), tol);
  EXPECT_NEAR(theta_benchmark(1), theta_pred(1), tol);

  int n_sim = 5e5;
  Eigen::VectorXd theta_dim0(n_sim);
  Eigen::VectorXd theta_dim1(n_sim);
  for (int i = 0; i < n_sim; i++) {
    rng.seed(2025 + i);
    Eigen::MatrixXd theta_pred = laplace_latent_neg_binomial_2_log_rng(
        y, y_index, eta, 0, diagonal_kernel_nb_functor{},
        std::forward_as_tuple(phi(0), phi(1)), rng, nullptr);

    theta_dim0(i) = theta_pred(0);
    theta_dim1(i) = theta_pred(1);
  }

  Eigen::MatrixXd K_sample(2, 2);
  K_sample(0, 0)
      = theta_dim0.array().square().mean() - square(theta_dim0.mean());
  K_sample(1, 1)
      = theta_dim1.array().square().mean() - square(theta_dim1.mean());
  K_sample(0, 1) = theta_dim0.cwiseProduct(theta_dim1).mean()
                   - theta_dim0.mean() * theta_dim1.mean();
  K_sample(1, 0) = K_sample(0, 1);

  // Check answers are within three std of the true answer.
  EXPECT_NEAR(theta_root(0), theta_dim0.mean(),
              3 * sqrt(K_laplace(0, 0) / n_sim));
  EXPECT_NEAR(theta_root(1), theta_dim1.mean(),
              3 * sqrt(K_laplace(1, 1) / n_sim));

  // Check sample covariance
  EXPECT_NEAR(K_laplace(0, 0), K_sample(0, 0), 6e-3);
  EXPECT_NEAR(K_laplace(1, 1), K_sample(1, 1), 6e-3);
  EXPECT_NEAR(K_laplace(0, 1), K_sample(0, 1), 1e-3);
}
}  // namespace
