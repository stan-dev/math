#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsGumbel, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << -2.0, 0.5, 1;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.1, 0.5;
  Eigen::VectorXd y_value(N);
  y_value << NAN, 0.5, 0.99;

  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu_size(N - 1);
  mu_size << 0.3, 0.8;
  Eigen::VectorXd mu_value1(N);
  mu_value1 << 0.3, -INFINITY, 0.5;
  Eigen::VectorXd mu_value2(N);
  mu_value2 << 0.3, NAN, 0.5;

  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, INFINITY;
  Eigen::VectorXd beta_size(N - 1);
  beta_size << 0.3, 0.8;
  Eigen::VectorXd beta_value1(N);
  beta_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd beta_value2(N);
  beta_value2 << 0.3, NAN, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);

  stan::math::matrix_cl<double> mu_cl(mu);
  stan::math::matrix_cl<double> mu_size_cl(mu_size);
  stan::math::matrix_cl<double> mu_value1_cl(mu_value1);
  stan::math::matrix_cl<double> mu_value2_cl(mu_value2);

  stan::math::matrix_cl<double> beta_cl(beta);
  stan::math::matrix_cl<double> beta_size_cl(beta_size);
  stan::math::matrix_cl<double> beta_value1_cl(beta_value1);
  stan::math::matrix_cl<double> beta_value2_cl(beta_value2);

  EXPECT_NO_THROW(stan::math::gumbel_lpdf(y_cl, mu_cl, beta_cl));

  EXPECT_THROW(stan::math::gumbel_lpdf(y_size_cl, mu_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::gumbel_lpdf(y_cl, mu_size_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::gumbel_lpdf(y_cl, mu_cl, beta_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::gumbel_lpdf(y_value_cl, mu_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::gumbel_lpdf(y_cl, mu_value1_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::gumbel_lpdf(y_cl, mu_cl, beta_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::gumbel_lpdf(y_cl, mu_value2_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::gumbel_lpdf(y_cl, mu_cl, beta_value2_cl),
               std::domain_error);
}

auto gumbel_lpdf_functor = [](const auto& y, const auto& mu, const auto& beta) {
  return stan::math::gumbel_lpdf(y, mu, beta);
};
auto gumbel_lpdf_functor_propto
    = [](const auto& y, const auto& mu, const auto& beta) {
        return stan::math::gumbel_lpdf<true>(y, mu, beta);
      };

TEST(ProbDistributionsGumbel, opencl_matches_cpu_small) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 1;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(gumbel_lpdf_functor, y, mu,
                                                beta);
  stan::math::test::compare_cpu_opencl_prim_rev(gumbel_lpdf_functor_propto, y,
                                                mu, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      gumbel_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      gumbel_lpdf_functor_propto, y.transpose().eval(), mu.transpose().eval(),
      beta.transpose().eval());
}

TEST(ProbDistributionsGumbel, opencl_broadcast_y) {
  int N = 3;

  double y = 0.3;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(gumbel_lpdf_functor, y,
                                                         mu, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      gumbel_lpdf_functor_propto, y, mu, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      gumbel_lpdf_functor, y, mu.transpose().eval(), beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      gumbel_lpdf_functor_propto, y, mu, beta.transpose().eval());
}

TEST(ProbDistributionsGumbel, opencl_broadcast_mu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double mu = 0.3;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(gumbel_lpdf_functor, y,
                                                         mu, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      gumbel_lpdf_functor_propto, y, mu, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      gumbel_lpdf_functor, y.transpose().eval(), mu, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      gumbel_lpdf_functor_propto, y, mu, beta.transpose().eval());
}

TEST(ProbDistributionsGumbel, opencl_broadcast_beta) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  double beta = 0.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(gumbel_lpdf_functor, y,
                                                         mu, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      gumbel_lpdf_functor_propto, y, mu, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      gumbel_lpdf_functor, y.transpose().eval(), mu, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      gumbel_lpdf_functor_propto, y, mu.transpose().eval(), beta);
}

TEST(ProbDistributionsGumbel, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> beta
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(gumbel_lpdf_functor, y, mu,
                                                beta);
  stan::math::test::compare_cpu_opencl_prim_rev(gumbel_lpdf_functor_propto, y,
                                                mu, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      gumbel_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      gumbel_lpdf_functor_propto, y.transpose().eval(), mu.transpose().eval(),
      beta.transpose().eval());
}

#endif
