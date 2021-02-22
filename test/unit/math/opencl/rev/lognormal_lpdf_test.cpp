#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsLognormal, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 10;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.1, 0.5;
  Eigen::VectorXd y_value1(N);
  y_value1 << -0.1, 0.5, 0.99;
  Eigen::VectorXd y_value2(N);
  y_value2 << 0.1, 1.1, NAN;

  Eigen::VectorXd mu(N);
  mu << -1.3, 0.8, 1.5;
  Eigen::VectorXd mu_size(N - 1);
  mu_size << 0.3, 0.8;
  Eigen::VectorXd mu_value1(N);
  mu_value1 << 0.3, NAN, 0.5;
  Eigen::VectorXd mu_value2(N);
  mu_value2 << 0.3, INFINITY, 0.5;

  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma_size(N - 1);
  sigma_size << 0.3, 0.8;
  Eigen::VectorXd sigma_value1(N);
  sigma_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd sigma_value2(N);
  sigma_value2 << 0.3, INFINITY, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value1_cl(y_value1);
  stan::math::matrix_cl<double> y_value2_cl(y_value2);

  stan::math::matrix_cl<double> mu_cl(mu);
  stan::math::matrix_cl<double> mu_size_cl(mu_size);
  stan::math::matrix_cl<double> mu_value1_cl(mu_value1);
  stan::math::matrix_cl<double> mu_value2_cl(mu_value2);

  stan::math::matrix_cl<double> sigma_cl(sigma);
  stan::math::matrix_cl<double> sigma_size_cl(sigma_size);
  stan::math::matrix_cl<double> sigma_value1_cl(sigma_value1);
  stan::math::matrix_cl<double> sigma_value2_cl(sigma_value2);

  EXPECT_NO_THROW(stan::math::lognormal_lpdf(y_cl, mu_cl, sigma_cl));

  EXPECT_THROW(stan::math::lognormal_lpdf(y_size_cl, mu_cl, sigma_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::lognormal_lpdf(y_cl, mu_size_cl, sigma_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::lognormal_lpdf(y_cl, mu_cl, sigma_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::lognormal_lpdf(y_value1_cl, mu_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::lognormal_lpdf(y_cl, mu_value1_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::lognormal_lpdf(y_cl, mu_cl, sigma_value1_cl),
               std::domain_error);

  EXPECT_THROW(stan::math::lognormal_lpdf(y_value2_cl, mu_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::lognormal_lpdf(y_cl, mu_value2_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::lognormal_lpdf(y_cl, mu_cl, sigma_value2_cl),
               std::domain_error);
}

auto lognormal_lpdf_functor
    = [](const auto& y, const auto& mu, const auto& sigma) {
        return stan::math::lognormal_lpdf(y, mu, sigma);
      };
auto lognormal_lpdf_functor_propto
    = [](const auto& y, const auto& mu, const auto& sigma) {
        return stan::math::lognormal_lpdf<true>(y, mu, sigma);
      };

TEST(ProbDistributionsLognormal, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.1, 0.5, 1;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.3;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.6, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(lognormal_lpdf_functor, y, mu,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(lognormal_lpdf_functor_propto,
                                                y, mu, sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(
      lognormal_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      lognormal_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), sigma.transpose().eval());
}
TEST(ProbDistributionsLognormal, opencl_matches_cpu_small_zero_y) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 1;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.6, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.3;

  stan::math::test::compare_cpu_opencl_prim_rev(lognormal_lpdf_functor, y, mu,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(lognormal_lpdf_functor_propto,
                                                y, mu, sigma);
}

TEST(ProbDistributionsLognormal, opencl_broadcast_y) {
  int N = 3;

  double y = 0.3;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(lognormal_lpdf_functor,
                                                         y, mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      lognormal_lpdf_functor_propto, y, mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      lognormal_lpdf_functor, y, mu.transpose().eval(), sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      lognormal_lpdf_functor_propto, y, mu, sigma.transpose().eval());
}

TEST(ProbDistributionsLognormal, opencl_broadcast_mu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double mu = 0.3;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(lognormal_lpdf_functor,
                                                         y, mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      lognormal_lpdf_functor_propto, y, mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      lognormal_lpdf_functor, y.transpose().eval(), mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      lognormal_lpdf_functor_propto, y, mu, sigma.transpose().eval());
}

TEST(ProbDistributionsLognormal, opencl_broadcast_sigma) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  double sigma = 0.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(lognormal_lpdf_functor,
                                                         y, mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      lognormal_lpdf_functor_propto, y, mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      lognormal_lpdf_functor, y.transpose().eval(), mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      lognormal_lpdf_functor_propto, y, mu.transpose().eval(), sigma);
}

TEST(ProbDistributionsLognormal, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> sigma
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(lognormal_lpdf_functor, y, mu,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(lognormal_lpdf_functor_propto,
                                                y, mu, sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(
      lognormal_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      lognormal_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), sigma.transpose().eval());
}

#endif
