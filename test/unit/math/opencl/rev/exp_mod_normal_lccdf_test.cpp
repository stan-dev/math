#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsDoubleExpModNormalLccdf, error_checking) {
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

  Eigen::VectorXd lambda(N);
  lambda << 0.4, 0.4, 1.4;
  Eigen::VectorXd lambda_size(N - 1);
  lambda_size << 0.3, 0.8;
  Eigen::VectorXd lambda_value(N);
  lambda_value << 0.3, 0, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);
  stan::math::matrix_cl<double> mu_cl(mu);
  stan::math::matrix_cl<double> mu_size_cl(mu_size);
  stan::math::matrix_cl<double> mu_value_cl(mu_value);
  stan::math::matrix_cl<double> sigma_cl(sigma);
  stan::math::matrix_cl<double> sigma_size_cl(sigma_size);
  stan::math::matrix_cl<double> sigma_value_cl(sigma_value);
  stan::math::matrix_cl<double> lambda_cl(lambda);
  stan::math::matrix_cl<double> lambda_size_cl(lambda_size);
  stan::math::matrix_cl<double> lambda_value_cl(lambda_value);

  EXPECT_NO_THROW(
      stan::math::exp_mod_normal_lccdf(y_cl, mu_cl, sigma_cl, lambda_cl));

  EXPECT_THROW(
      stan::math::exp_mod_normal_lccdf(y_size_cl, mu_cl, sigma_cl, lambda_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::exp_mod_normal_lccdf(y_cl, mu_size_cl, sigma_cl, lambda_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::exp_mod_normal_lccdf(y_cl, mu_cl, sigma_size_cl, lambda_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::exp_mod_normal_lccdf(y_cl, mu_cl, sigma_cl, lambda_size_cl),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::exp_mod_normal_lccdf(y_value_cl, mu_cl, sigma_cl, lambda_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::exp_mod_normal_lccdf(y_cl, mu_value_cl, sigma_cl, lambda_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::exp_mod_normal_lccdf(y_cl, mu_cl, sigma_value_cl, lambda_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::exp_mod_normal_lccdf(y_cl, mu_cl, sigma_cl, lambda_value_cl),
      std::domain_error);
}

auto exp_mod_normal_lccdf_functor
    = [](const auto& y, const auto& mu, const auto& sigma, const auto& lambda) {
        return stan::math::exp_mod_normal_lccdf(y, mu, sigma, lambda);
      };

TEST(ProbDistributionsDoubleExpModNormalLccdf, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << -0.3, 1.8, 1.4;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.4, 1.1;

  stan::math::test::compare_cpu_opencl_prim_rev(exp_mod_normal_lccdf_functor, y,
                                                mu, sigma, lambda);
  stan::math::test::compare_cpu_opencl_prim_rev(
      exp_mod_normal_lccdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval(), lambda.transpose().eval());
}

TEST(ProbDistributionsDoubleExpModNormalLccdf,
     opencl_matches_cpu_small_y_pos_inf) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << -0.3, 1.8, INFINITY;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.4, 1.1;

  stan::math::test::compare_cpu_opencl_prim_rev(exp_mod_normal_lccdf_functor, y,
                                                mu, sigma, lambda);
  stan::math::test::compare_cpu_opencl_prim_rev(
      exp_mod_normal_lccdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval(), lambda.transpose().eval());
}

TEST(ProbDistributionsDoubleExpModNormalLccdf,
     opencl_matches_cpu_small_y_neg_inf) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << -0.3, 1.8, -INFINITY;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.4, 1.1;

  stan::math::test::compare_cpu_opencl_prim_rev(exp_mod_normal_lccdf_functor, y,
                                                mu, sigma, lambda);
  stan::math::test::compare_cpu_opencl_prim_rev(
      exp_mod_normal_lccdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval(), lambda.transpose().eval());
}

TEST(ProbDistributionsDoubleExpModNormalLccdf, opencl_broadcast_y) {
  int N = 3;

  double y_scal = 12.3;
  Eigen::VectorXd mu(N);
  mu << 0.5, 1.2, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.4, 1.1;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      exp_mod_normal_lccdf_functor, y_scal, mu, sigma, lambda);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      exp_mod_normal_lccdf_functor, y_scal, mu.transpose().eval(), sigma,
      lambda.transpose().eval());
}

TEST(ProbDistributionsDoubleExpModNormalLccdf, opencl_matches_cpu_big) {
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

  stan::math::test::compare_cpu_opencl_prim_rev(exp_mod_normal_lccdf_functor, y,
                                                mu, sigma, lambda);
  stan::math::test::compare_cpu_opencl_prim_rev(
      exp_mod_normal_lccdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval(), lambda.transpose().eval());
}

#endif
