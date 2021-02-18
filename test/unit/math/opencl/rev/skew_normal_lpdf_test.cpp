#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsSkewNormal, error_checking) {
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

  Eigen::VectorXd alpha(N);
  alpha << 0.3, -0.8, 1.0;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value(N);
  alpha_value << 0.3, INFINITY, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);
  stan::math::matrix_cl<double> mu_cl(mu);
  stan::math::matrix_cl<double> mu_size_cl(mu_size);
  stan::math::matrix_cl<double> mu_value_cl(mu_value);
  stan::math::matrix_cl<double> sigma_cl(sigma);
  stan::math::matrix_cl<double> sigma_size_cl(sigma_size);
  stan::math::matrix_cl<double> sigma_value_cl(sigma_value);
  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value_cl(alpha_value);

  EXPECT_NO_THROW(
      stan::math::skew_normal_lpdf(y_cl, mu_cl, sigma_cl, alpha_cl));

  EXPECT_THROW(
      stan::math::skew_normal_lpdf(y_size_cl, mu_cl, sigma_cl, alpha_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::skew_normal_lpdf(y_cl, mu_size_cl, sigma_cl, alpha_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::skew_normal_lpdf(y_cl, mu_cl, sigma_size_cl, alpha_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::skew_normal_lpdf(y_cl, mu_cl, sigma_cl, alpha_size_cl),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::skew_normal_lpdf(y_value_cl, mu_cl, sigma_cl, alpha_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::skew_normal_lpdf(y_cl, mu_value_cl, sigma_cl, alpha_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::skew_normal_lpdf(y_cl, mu_cl, sigma_value_cl, alpha_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::skew_normal_lpdf(y_cl, mu_cl, sigma_cl, alpha_value_cl),
      std::domain_error);
}

auto skew_normal_lpdf_functor
    = [](const auto& y, const auto& mu, const auto& sigma, const auto& alpha) {
        return stan::math::skew_normal_lpdf(y, mu, sigma, alpha);
      };
auto skew_normal_lpdf_functor_propto
    = [](const auto& y, const auto& mu, const auto& sigma, const auto& alpha) {
        return stan::math::skew_normal_lpdf<true>(y, mu, sigma, alpha);
      };

TEST(ProbDistributionsSkewNormal, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.3, 1.5;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, -0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(skew_normal_lpdf_functor, y, mu,
                                                sigma, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(skew_normal_lpdf_functor_propto,
                                                y, mu, sigma, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_normal_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval(), alpha.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_normal_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), sigma.transpose().eval(),
      alpha.transpose().eval());
}

TEST(ProbDistributionsSkewNormal, opencl_broadcast_y) {
  int N = 3;

  double y_scal = 12.3;
  Eigen::VectorXd mu(N);
  mu << 0.5, 1.2, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, -0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      skew_normal_lpdf_functor, y_scal, mu, sigma, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      skew_normal_lpdf_functor_propto, y_scal, mu, sigma, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      skew_normal_lpdf_functor, y_scal, mu.transpose().eval(),
      sigma.transpose().eval(), alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      skew_normal_lpdf_functor_propto, y_scal, mu, sigma.transpose().eval(),
      alpha.transpose().eval());
}

TEST(ProbDistributionsSkewNormal, opencl_broadcast_mu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double mu_scal = 12.3;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, -0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      skew_normal_lpdf_functor, y, mu_scal, sigma, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      skew_normal_lpdf_functor_propto, y, mu_scal, sigma, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      skew_normal_lpdf_functor, y, mu_scal, sigma.transpose().eval(),
      alpha.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      skew_normal_lpdf_functor_propto, y.transpose().eval(), mu_scal, sigma,
      alpha.transpose().eval());
}

TEST(ProbDistributionsSkewNormal, opencl_broadcast_sigma) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.3, 1.5;
  double sigma_scal = 12.3;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, -0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      skew_normal_lpdf_functor, y, mu, sigma_scal, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      skew_normal_lpdf_functor_propto, y, mu, sigma_scal, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      skew_normal_lpdf_functor, y.transpose().eval(), mu, sigma_scal,
      alpha.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      skew_normal_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), sigma_scal, alpha);
}

TEST(ProbDistributionsSkewNormal, opencl_broadcast_alpha) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.3, 1.5;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  double alpha_scal = 12.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      skew_normal_lpdf_functor, y, mu, sigma, alpha_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      skew_normal_lpdf_functor_propto, y, mu, sigma, alpha_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      skew_normal_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma, alpha_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      skew_normal_lpdf_functor_propto, y, mu.transpose().eval(),
      sigma.transpose().eval(), alpha_scal);
}

TEST(ProbDistributionsSkewNormal, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> sigma
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(skew_normal_lpdf_functor, y, mu,
                                                sigma, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(skew_normal_lpdf_functor_propto,
                                                y, mu, sigma, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_normal_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval(), alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(
      skew_normal_lpdf_functor_propto, y, mu.transpose().eval(),
      sigma.transpose().eval(), alpha.transpose().eval());
}

#endif
