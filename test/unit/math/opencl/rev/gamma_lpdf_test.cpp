#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsGamma, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 2.0, 0.5, 1;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.1, 0.5;
  Eigen::VectorXd y_value(N);
  y_value << NAN, 0.5, 0.99;

  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value1(N);
  alpha_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd alpha_value2(N);
  alpha_value2 << 0.3, INFINITY, 0.5;

  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;
  Eigen::VectorXd beta_size(N - 1);
  beta_size << 0.3, 0.8;
  Eigen::VectorXd beta_value1(N);
  beta_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd beta_value2(N);
  beta_value2 << 0.3, INFINITY, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);

  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value1_cl(alpha_value1);
  stan::math::matrix_cl<double> alpha_value2_cl(alpha_value2);

  stan::math::matrix_cl<double> beta_cl(beta);
  stan::math::matrix_cl<double> beta_size_cl(beta_size);
  stan::math::matrix_cl<double> beta_value1_cl(beta_value1);
  stan::math::matrix_cl<double> beta_value2_cl(beta_value2);

  EXPECT_NO_THROW(stan::math::gamma_lpdf(y_cl, alpha_cl, beta_cl));

  EXPECT_THROW(stan::math::gamma_lpdf(y_size_cl, alpha_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::gamma_lpdf(y_cl, alpha_size_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::gamma_lpdf(y_cl, alpha_cl, beta_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::gamma_lpdf(y_value_cl, alpha_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::gamma_lpdf(y_cl, alpha_value1_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::gamma_lpdf(y_cl, alpha_cl, beta_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::gamma_lpdf(y_cl, alpha_value2_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::gamma_lpdf(y_cl, alpha_cl, beta_value2_cl),
               std::domain_error);
}

auto gamma_lpdf_functor
    = [](const auto& y, const auto& alpha, const auto& beta) {
        return stan::math::gamma_lpdf(y, alpha, beta);
      };
auto gamma_lpdf_functor_propto
    = [](const auto& y, const auto& alpha, const auto& beta) {
        return stan::math::gamma_lpdf<true>(y, alpha, beta);
      };

TEST(ProbDistributionsGamma, opencl_matches_cpu_small) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.1, 0.5, 1;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(gamma_lpdf_functor, y, alpha,
                                                beta);
  stan::math::test::compare_cpu_opencl_prim_rev(gamma_lpdf_functor_propto, y,
                                                alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      gamma_lpdf_functor, y.transpose().eval(), alpha.transpose().eval(),
      beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      gamma_lpdf_functor_propto, y.transpose().eval(), alpha.transpose().eval(),
      beta.transpose().eval());
}

TEST(ProbDistributionsGamma, opencl_broadcast_y) {
  int N = 3;

  double y = 0.3;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(gamma_lpdf_functor, y,
                                                         alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      gamma_lpdf_functor_propto, y, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      gamma_lpdf_functor, y, alpha.transpose().eval(), beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      gamma_lpdf_functor_propto, y, alpha, beta.transpose().eval());
}

TEST(ProbDistributionsGamma, opencl_broadcast_alpha) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double alpha = 0.3;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(gamma_lpdf_functor, y,
                                                         alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      gamma_lpdf_functor_propto, y, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      gamma_lpdf_functor, y.transpose().eval(), alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      gamma_lpdf_functor_propto, y, alpha, beta.transpose().eval());
}

TEST(ProbDistributionsGamma, opencl_broadcast_beta) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;
  double beta = 0.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(gamma_lpdf_functor, y,
                                                         alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      gamma_lpdf_functor_propto, y, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      gamma_lpdf_functor, y.transpose().eval(), alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      gamma_lpdf_functor_propto, y, alpha.transpose().eval(), beta);
}

TEST(ProbDistributionsGamma, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> beta
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(gamma_lpdf_functor, y, alpha,
                                                beta);
  stan::math::test::compare_cpu_opencl_prim_rev(gamma_lpdf_functor_propto, y,
                                                alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      gamma_lpdf_functor, y.transpose().eval(), alpha.transpose().eval(),
      beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      gamma_lpdf_functor_propto, y.transpose().eval(), alpha.transpose().eval(),
      beta.transpose().eval());
}

#endif
