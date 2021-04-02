#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsSkewDoubleExponential, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << -100, INFINITY, 11.3;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.1, 0.5;
  Eigen::VectorXd y_value(N);
  y_value << NAN, 0.5, 0.99;

  Eigen::VectorXd mu(N);
  mu << -10.3, 0.8, 2.1;
  Eigen::VectorXd mu_size(N - 1);
  mu_size << 0.3, 0.8;
  Eigen::VectorXd mu_value(N);
  mu_value << 0.3, INFINITY, 0.5;

  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 11.0;
  Eigen::VectorXd sigma_size(N - 1);
  sigma_size << 0.3, 0.8;
  Eigen::VectorXd sigma_value1(N);
  sigma_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd sigma_value2(N);
  sigma_value2 << 0.3, INFINITY, 0.5;

  Eigen::VectorXd tau(N);
  tau << 0.3, 0.8, 0.9;
  Eigen::VectorXd tau_size(N - 1);
  tau_size << 0.3, 0.8;
  Eigen::VectorXd tau_value1(N);
  tau_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd tau_value2(N);
  tau_value2 << 0.3, 1.1, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);

  stan::math::matrix_cl<double> mu_cl(mu);
  stan::math::matrix_cl<double> mu_size_cl(mu_size);
  stan::math::matrix_cl<double> mu_value_cl(mu_value);

  stan::math::matrix_cl<double> sigma_cl(sigma);
  stan::math::matrix_cl<double> sigma_size_cl(sigma_size);
  stan::math::matrix_cl<double> sigma_value1_cl(sigma_value1);
  stan::math::matrix_cl<double> sigma_value2_cl(sigma_value2);

  stan::math::matrix_cl<double> tau_cl(tau);
  stan::math::matrix_cl<double> tau_size_cl(tau_size);
  stan::math::matrix_cl<double> tau_value1_cl(tau_value1);
  stan::math::matrix_cl<double> tau_value2_cl(tau_value2);

  EXPECT_NO_THROW(
      stan::math::skew_double_exponential_lpdf(y_cl, mu_cl, sigma_cl, tau_cl));

  EXPECT_THROW(stan::math::skew_double_exponential_lpdf(y_size_cl, mu_cl,
                                                        sigma_cl, tau_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::skew_double_exponential_lpdf(y_cl, mu_size_cl,
                                                        sigma_cl, tau_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::skew_double_exponential_lpdf(y_cl, mu_cl,
                                                        sigma_size_cl, tau_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::skew_double_exponential_lpdf(y_cl, mu_cl, sigma_cl,
                                                        tau_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::skew_double_exponential_lpdf(y_value_cl, mu_cl,
                                                        sigma_cl, tau_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::skew_double_exponential_lpdf(y_cl, mu_value_cl,
                                                        sigma_cl, tau_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::skew_double_exponential_lpdf(
                   y_cl, mu_cl, sigma_value1_cl, tau_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::skew_double_exponential_lpdf(y_cl, mu_cl, sigma_cl,
                                                        tau_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::skew_double_exponential_lpdf(
                   y_cl, mu_cl, sigma_value2_cl, tau_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::skew_double_exponential_lpdf(y_cl, mu_cl, sigma_cl,
                                                        tau_value2_cl),
               std::domain_error);
}

auto skew_double_exponential_lpdf_functor
    = [](const auto& y, const auto& mu, const auto& sigma, const auto& tau) {
        return stan::math::skew_double_exponential_lpdf(y, mu, sigma, tau);
      };
auto skew_double_exponential_lpdf_functor_propto =
    [](const auto& y, const auto& mu, const auto& sigma, const auto& tau) {
      return stan::math::skew_double_exponential_lpdf<true>(y, mu, sigma, tau);
    };

TEST(ProbDistributionsSkewDoubleExponential, opencl_matches_cpu_small) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 1;
  Eigen::VectorXd mu(N);
  mu << -10.3, 0.8, 2.1;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 11.0;
  Eigen::VectorXd tau(N);
  tau << 0.3, 0.8, 0.9;

  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_double_exponential_lpdf_functor, y, mu, sigma, tau);
  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_double_exponential_lpdf_functor_propto, y, mu, sigma, tau);
  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_double_exponential_lpdf_functor, y.transpose().eval(),
      mu.transpose().eval(), sigma.transpose().eval(), tau.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_double_exponential_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), sigma.transpose().eval(), tau.transpose().eval());
}

TEST(ProbDistributionsSkewDoubleExponential, opencl_broadcast_y) {
  int N = 3;

  double y = 0.3;
  Eigen::VectorXd mu(N);
  mu << -10.3, 0.8, 2.1;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 11.0;
  Eigen::VectorXd tau(N);
  tau << 0.3, 0.8, 0.9;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      skew_double_exponential_lpdf_functor, y, mu, sigma, tau);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      skew_double_exponential_lpdf_functor_propto, y, mu, sigma, tau);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      skew_double_exponential_lpdf_functor, y, mu.transpose().eval(),
      sigma.transpose().eval(), tau);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      skew_double_exponential_lpdf_functor_propto, y, mu,
      sigma.transpose().eval(), tau.transpose().eval());
}

TEST(ProbDistributionsSkewDoubleExponential, opencl_broadcast_mu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 1;
  double mu = 3.2;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 11.0;
  Eigen::VectorXd tau(N);
  tau << 0.3, 0.8, 0.9;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      skew_double_exponential_lpdf_functor, y, mu, sigma, tau);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      skew_double_exponential_lpdf_functor_propto, y, mu, sigma, tau);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      skew_double_exponential_lpdf_functor, y, mu, sigma.transpose().eval(),
      tau.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      skew_double_exponential_lpdf_functor_propto, y.transpose().eval(), mu,
      sigma, tau.transpose().eval());
}

TEST(ProbDistributionsSkewDoubleExponential, opencl_broadcast_sigma) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 1;
  Eigen::VectorXd mu(N);
  mu << -10.3, 0.8, 2.1;
  double sigma = 0.8;
  Eigen::VectorXd tau(N);
  tau << 0.3, 0.8, 0.9;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      skew_double_exponential_lpdf_functor, y, mu, sigma, tau);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      skew_double_exponential_lpdf_functor_propto, y, mu, sigma, tau);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      skew_double_exponential_lpdf_functor, y.transpose().eval(), mu, sigma,
      tau.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      skew_double_exponential_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), sigma, tau);
}

TEST(ProbDistributionsSkewDoubleExponential, opencl_broadcast_tau) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 1;
  Eigen::VectorXd mu(N);
  mu << -10.3, 0.8, 2.1;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 11.0;
  double tau = 0.9;

  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      skew_double_exponential_lpdf_functor, y, mu, sigma, tau);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      skew_double_exponential_lpdf_functor_propto, y, mu, sigma, tau);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      skew_double_exponential_lpdf_functor, y.transpose().eval(),
      mu.transpose().eval(), sigma, tau);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      skew_double_exponential_lpdf_functor_propto, y, mu.transpose().eval(),
      sigma.transpose().eval(), tau);
}

TEST(ProbDistributionsSkewDoubleExponential, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> sigma
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs().array()
        + 0.1;
  Eigen::Matrix<double, Eigen::Dynamic, 1> tau
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_double_exponential_lpdf_functor, y, mu, sigma, tau);
  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_double_exponential_lpdf_functor_propto, y, mu, sigma, tau);
  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_double_exponential_lpdf_functor, y.transpose().eval(),
      mu.transpose().eval(), sigma.transpose().eval(), tau.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_double_exponential_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), sigma.transpose().eval(), tau.transpose().eval());
}

#endif
