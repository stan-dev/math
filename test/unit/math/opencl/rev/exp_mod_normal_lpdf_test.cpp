#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsExpModNormal, error_checking) {
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

  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.8, 1.3;
  Eigen::VectorXd lambda_size(N - 1);
  lambda_size << 0.3, 0.8;
  Eigen::VectorXd lambda_value1(N);
  lambda_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd lambda_value2(N);
  lambda_value2 << 0.3, INFINITY, 0.5;

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

  stan::math::matrix_cl<double> lambda_cl(lambda);
  stan::math::matrix_cl<double> lambda_size_cl(lambda_size);
  stan::math::matrix_cl<double> lambda_value1_cl(lambda_value1);
  stan::math::matrix_cl<double> lambda_value2_cl(lambda_value2);

  EXPECT_NO_THROW(
      stan::math::exp_mod_normal_lpdf(y_cl, mu_cl, sigma_cl, lambda_cl));

  EXPECT_THROW(
      stan::math::exp_mod_normal_lpdf(y_size_cl, mu_cl, sigma_cl, lambda_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::exp_mod_normal_lpdf(y_cl, mu_size_cl, sigma_cl, lambda_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::exp_mod_normal_lpdf(y_cl, mu_cl, sigma_size_cl, lambda_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::exp_mod_normal_lpdf(y_cl, mu_cl, sigma_cl, lambda_size_cl),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::exp_mod_normal_lpdf(y_value_cl, mu_cl, sigma_cl, lambda_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::exp_mod_normal_lpdf(y_cl, mu_value_cl, sigma_cl, lambda_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::exp_mod_normal_lpdf(y_cl, mu_cl, sigma_value1_cl, lambda_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::exp_mod_normal_lpdf(y_cl, mu_cl, sigma_cl, lambda_value1_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::exp_mod_normal_lpdf(y_cl, mu_cl, sigma_value2_cl, lambda_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::exp_mod_normal_lpdf(y_cl, mu_cl, sigma_cl, lambda_value2_cl),
      std::domain_error);
}

auto exp_mod_normal_lpdf_functor
    = [](const auto& y, const auto& mu, const auto& sigma, const auto& lambda) {
        return stan::math::exp_mod_normal_lpdf(y, mu, sigma, lambda);
      };
auto exp_mod_normal_lpdf_functor_propto
    = [](const auto& y, const auto& mu, const auto& sigma, const auto& lambda) {
        return stan::math::exp_mod_normal_lpdf<true>(y, mu, sigma, lambda);
      };

TEST(ProbDistributionsExpModNormal, opencl_matches_cpu_small) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 1;
  Eigen::VectorXd mu(N);
  mu << -10.3, 0.8, 2.1;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 11.0;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.8, 1.30;

  stan::math::test::compare_cpu_opencl_prim_rev(exp_mod_normal_lpdf_functor, y,
                                                mu, sigma, lambda);
  stan::math::test::compare_cpu_opencl_prim_rev(
      exp_mod_normal_lpdf_functor_propto, y, mu, sigma, lambda);
  stan::math::test::compare_cpu_opencl_prim_rev(
      exp_mod_normal_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval(), lambda.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      exp_mod_normal_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), sigma.transpose().eval(),
      lambda.transpose().eval());
}

TEST(ProbDistributionsExpModNormal, opencl_broadcast_y) {
  int N = 3;

  double y = 0.3;
  Eigen::VectorXd mu(N);
  mu << -10.3, 0.8, 2.1;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 11.0;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.8, 1.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      exp_mod_normal_lpdf_functor, y, mu, sigma, lambda);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      exp_mod_normal_lpdf_functor_propto, y, mu, sigma, lambda);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      exp_mod_normal_lpdf_functor, y, mu.transpose().eval(),
      sigma.transpose().eval(), lambda);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      exp_mod_normal_lpdf_functor_propto, y, mu, sigma.transpose().eval(),
      lambda.transpose().eval());
}

TEST(ProbDistributionsExpModNormal, opencl_broadcast_mu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 1;
  double mu = 3.2;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 11.0;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.8, 1.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      exp_mod_normal_lpdf_functor, y, mu, sigma, lambda);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      exp_mod_normal_lpdf_functor_propto, y, mu, sigma, lambda);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      exp_mod_normal_lpdf_functor, y, mu, sigma.transpose().eval(),
      lambda.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      exp_mod_normal_lpdf_functor_propto, y.transpose().eval(), mu, sigma,
      lambda.transpose().eval());
}

TEST(ProbDistributionsExpModNormal, opencl_broadcast_sigma) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 1;
  Eigen::VectorXd mu(N);
  mu << -10.3, 0.8, 2.1;
  double sigma = 0.8;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.8, 1.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      exp_mod_normal_lpdf_functor, y, mu, sigma, lambda);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      exp_mod_normal_lpdf_functor_propto, y, mu, sigma, lambda);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      exp_mod_normal_lpdf_functor, y.transpose().eval(), mu, sigma,
      lambda.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      exp_mod_normal_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), sigma, lambda);
}

TEST(ProbDistributionsExpModNormal, opencl_broadcast_lambda) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 1;
  Eigen::VectorXd mu(N);
  mu << -10.3, 0.8, 2.1;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 11.0;
  double lambda = 1.2;

  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      exp_mod_normal_lpdf_functor, y, mu, sigma, lambda);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      exp_mod_normal_lpdf_functor_propto, y, mu, sigma, lambda);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      exp_mod_normal_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma, lambda);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      exp_mod_normal_lpdf_functor_propto, y, mu.transpose().eval(),
      sigma.transpose().eval(), lambda);
}

TEST(ProbDistributionsExpModNormal, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> sigma
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs().array()
        + 0.1;
  Eigen::Matrix<double, Eigen::Dynamic, 1> lambda
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(exp_mod_normal_lpdf_functor, y,
                                                mu, sigma, lambda);
  stan::math::test::compare_cpu_opencl_prim_rev(
      exp_mod_normal_lpdf_functor_propto, y, mu, sigma, lambda);
  stan::math::test::compare_cpu_opencl_prim_rev(
      exp_mod_normal_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval(), lambda.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      exp_mod_normal_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), sigma.transpose().eval(),
      lambda.transpose().eval());
}

#endif
