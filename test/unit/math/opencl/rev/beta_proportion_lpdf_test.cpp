#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsBetaProportion, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.3, 0.8;
  Eigen::VectorXd y_value(N);
  y_value << 0.3, -0.8, 0.5;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 0.4;
  Eigen::VectorXd mu_size(N - 1);
  mu_size << 0.3, 0.8;
  Eigen::VectorXd mu_value(N);
  mu_value << 0.3, 1, 0.5;
  Eigen::VectorXd kappa(N);
  kappa << 0.3, 0.8, 3.0;
  Eigen::VectorXd kappa_size(N - 1);
  kappa_size << 0.3, 0.8;
  Eigen::VectorXd kappa_value(N);
  kappa_value << 0.3, -0.8, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);
  stan::math::matrix_cl<double> mu_cl(mu);
  stan::math::matrix_cl<double> mu_size_cl(mu_size);
  stan::math::matrix_cl<double> mu_value_cl(mu_value);
  stan::math::matrix_cl<double> kappa_cl(kappa);
  stan::math::matrix_cl<double> kappa_size_cl(kappa_size);
  stan::math::matrix_cl<double> kappa_value_cl(kappa_value);

  EXPECT_NO_THROW(stan::math::beta_proportion_lpdf(y_cl, mu_cl, kappa_cl));

  EXPECT_THROW(stan::math::beta_proportion_lpdf(y_size_cl, mu_cl, kappa_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::beta_proportion_lpdf(y_cl, mu_size_cl, kappa_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::beta_proportion_lpdf(y_cl, mu_cl, kappa_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::beta_proportion_lpdf(y_value_cl, mu_cl, kappa_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_proportion_lpdf(y_cl, mu_value_cl, kappa_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_proportion_lpdf(y_cl, mu_cl, kappa_value_cl),
               std::domain_error);
}

auto beta_proportion_lpdf_functor
    = [](const auto& y, const auto& mu, const auto& kappa) {
        return stan::math::beta_proportion_lpdf(y, mu, kappa);
      };
auto beta_proportion_lpdf_functor_propto
    = [](const auto& y, const auto& mu, const auto& kappa) {
        return stan::math::beta_proportion_lpdf<true>(y, mu, kappa);
      };

TEST(ProbDistributionsBetaProportion, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.1, 0.3, 0.4;
  Eigen::VectorXd kappa(N);
  kappa << 0.3, 1.8, 3.0;

  stan::math::test::compare_cpu_opencl_prim_rev(
      beta_proportion_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      kappa.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      beta_proportion_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), kappa.transpose().eval());
}

TEST(ProbDistributionsBetaProportion, opencl_broadcast_y) {
  int N = 3;

  double y_scal = 0.7;
  Eigen::VectorXd mu(N);
  mu << 0.1, 0.3, 0.4;
  Eigen::VectorXd kappa(N);
  kappa << 0.3, 1.8, 3.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      beta_proportion_lpdf_functor, y_scal, mu, kappa);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      beta_proportion_lpdf_functor_propto, y_scal, mu, kappa);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      beta_proportion_lpdf_functor, y_scal, mu.transpose().eval(), kappa);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      beta_proportion_lpdf_functor_propto, y_scal, mu,
      kappa.transpose().eval());
}

TEST(ProbDistributionsBetaProportion, opencl_broadcast_mu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double mu_scal = 0.9;
  Eigen::VectorXd kappa(N);
  kappa << 0.3, 1.8, 3.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      beta_proportion_lpdf_functor, y, mu_scal, kappa);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      beta_proportion_lpdf_functor_propto, y, mu_scal, kappa);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      beta_proportion_lpdf_functor, y.transpose().eval(), mu_scal, kappa);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      beta_proportion_lpdf_functor_propto, y, mu_scal,
      kappa.transpose().eval());
}

TEST(ProbDistributionsBetaProportion, opencl_broadcast_kappa) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.1, 0.3, 0.4;
  double kappa_scal = 3.2;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      beta_proportion_lpdf_functor, y, mu, kappa_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      beta_proportion_lpdf_functor_propto, y, mu, kappa_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      beta_proportion_lpdf_functor, y.transpose().eval(), mu, kappa_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      beta_proportion_lpdf_functor_propto, y, mu.transpose().eval(),
      kappa_scal);
}

TEST(ProbDistributionsBetaProportion, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> kappa
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(beta_proportion_lpdf_functor, y,
                                                mu, kappa);
  stan::math::test::compare_cpu_opencl_prim_rev(
      beta_proportion_lpdf_functor_propto, y, mu, kappa);
  stan::math::test::compare_cpu_opencl_prim_rev(
      beta_proportion_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      kappa.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      beta_proportion_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), kappa.transpose().eval());
}

#endif
