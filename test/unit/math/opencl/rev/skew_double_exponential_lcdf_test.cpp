#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsSkewDoubleExponentialLcdf, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.3, 0.8;
  Eigen::VectorXd y_value(N);
  y_value << 0.3, NAN, 0.5;

  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu_size(N - 1);
  mu_size << 0.3, 0.8;
  Eigen::VectorXd mu_value(N);
  mu_value << 0.3, -INFINITY, 0.5;

  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma_size(N - 1);
  sigma_size << 0.3, 0.8;
  Eigen::VectorXd sigma_value(N);
  sigma_value << 0.3, 0, 0.5;

  Eigen::VectorXd tau(N);
  tau << 0.4, 0.4, 1.0;
  Eigen::VectorXd tau_size(N - 1);
  tau_size << 0.3, 0.8;
  Eigen::VectorXd tau_value(N);
  tau_value << 0.3, 0.1, 1.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);
  stan::math::matrix_cl<double> mu_cl(mu);
  stan::math::matrix_cl<double> mu_size_cl(mu_size);
  stan::math::matrix_cl<double> mu_value_cl(mu_value);
  stan::math::matrix_cl<double> sigma_cl(sigma);
  stan::math::matrix_cl<double> sigma_size_cl(sigma_size);
  stan::math::matrix_cl<double> sigma_value_cl(sigma_value);
  stan::math::matrix_cl<double> tau_cl(tau);
  stan::math::matrix_cl<double> tau_size_cl(tau_size);
  stan::math::matrix_cl<double> tau_value_cl(tau_value);

  EXPECT_NO_THROW(
      stan::math::skew_double_exponential_lcdf(y_cl, mu_cl, sigma_cl, tau_cl));

  EXPECT_THROW(stan::math::skew_double_exponential_lcdf(y_size_cl, mu_cl,
                                                        sigma_cl, tau_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::skew_double_exponential_lcdf(y_cl, mu_size_cl,
                                                        sigma_cl, tau_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::skew_double_exponential_lcdf(y_cl, mu_cl,
                                                        sigma_size_cl, tau_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::skew_double_exponential_lcdf(y_cl, mu_cl, sigma_cl,
                                                        tau_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::skew_double_exponential_lcdf(y_value_cl, mu_cl,
                                                        sigma_cl, tau_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::skew_double_exponential_lcdf(y_cl, mu_value_cl,
                                                        sigma_cl, tau_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::skew_double_exponential_lcdf(y_cl, mu_cl,
                                                        sigma_value_cl, tau_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::skew_double_exponential_lcdf(y_cl, mu_cl, sigma_cl,
                                                        tau_value_cl),
               std::domain_error);
}

auto skew_double_exponential_lcdf_functor
    = [](const auto& y, const auto& mu, const auto& sigma, const auto& tau) {
        return stan::math::skew_double_exponential_lcdf(y, mu, sigma, tau);
      };

TEST(ProbDistributionsSkewDoubleExponentialLcdf, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << -0.3, 1.8, 1.4;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd tau(N);
  tau << 0.3, 0.4, 0.9;

  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_double_exponential_lcdf_functor, y, mu, sigma, tau);
  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_double_exponential_lcdf_functor, y.transpose().eval(),
      mu.transpose().eval(), sigma.transpose().eval(), tau.transpose().eval());
}

TEST(ProbDistributionsSkewDoubleExponentialLcdf,
     opencl_matches_cpu_small_y_neg_inf) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << -INFINITY, 1.8, 1.4;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd tau(N);
  tau << 0.3, 0.4, 0.9;

  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_double_exponential_lcdf_functor, y, mu, sigma, tau);
  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_double_exponential_lcdf_functor, y.transpose().eval(),
      mu.transpose().eval(), sigma.transpose().eval(), tau.transpose().eval());
}

TEST(ProbDistributionsSkewDoubleExponentialLcdf, opencl_broadcast_y) {
  int N = 3;

  double y_scal = 12.3;
  Eigen::VectorXd mu(N);
  mu << 0.5, 1.2, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd tau(N);
  tau << 0.3, 0.4, 0.9;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      skew_double_exponential_lcdf_functor, y_scal, mu, sigma, tau);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      skew_double_exponential_lcdf_functor, y_scal, mu.transpose().eval(),
      sigma, tau.transpose().eval());
}

TEST(ProbDistributionsSkewDoubleExponentialLcdf, opencl_broadcast_mu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double mu_scal = 12.3;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd tau(N);
  tau << 0.3, 0.4, 0.9;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      skew_double_exponential_lcdf_functor, y, mu_scal, sigma, tau);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      skew_double_exponential_lcdf_functor, y.transpose().eval(), mu_scal,
      sigma, tau.transpose().eval());
}

TEST(ProbDistributionsSkewDoubleExponentialLcdf, opencl_matches_cpu_big) {
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
      skew_double_exponential_lcdf_functor, y, mu, sigma, tau);
  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_double_exponential_lcdf_functor, y.transpose().eval(),
      mu.transpose().eval(), sigma.transpose().eval(), tau.transpose().eval());
}

#endif
