#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsExponential, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 12.3;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.1, 0.5;
  Eigen::VectorXd y_value(N);
  y_value << -0.1, 0.5, 0.99;

  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 12.4;
  Eigen::VectorXd beta_size(N - 1);
  beta_size << 0.3, 0.8;
  Eigen::VectorXd beta_value1(N);
  beta_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd beta_value2(N);
  beta_value2 << 0.3, INFINITY, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);

  stan::math::matrix_cl<double> beta_cl(beta);
  stan::math::matrix_cl<double> beta_size_cl(beta_size);
  stan::math::matrix_cl<double> beta_value1_cl(beta_value1);
  stan::math::matrix_cl<double> beta_value2_cl(beta_value2);

  EXPECT_NO_THROW(stan::math::exponential_lpdf(y_cl, beta_cl));

  EXPECT_THROW(stan::math::exponential_lpdf(y_size_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::exponential_lpdf(y_cl, beta_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::exponential_lpdf(y_value_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::exponential_lpdf(y_cl, beta_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::exponential_lpdf(y_cl, beta_value2_cl),
               std::domain_error);
}

auto exponential_lpdf_functor = [](const auto& y, const auto& beta) {
  return stan::math::exponential_lpdf(y, beta);
};
auto exponential_lpdf_functor_propto = [](const auto& y, const auto& beta) {
  return stan::math::exponential_lpdf<true>(y, beta);
};

TEST(ProbDistributionsExponential, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 12.3;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 21.0;

  stan::math::test::compare_cpu_opencl_prim_rev(exponential_lpdf_functor, y,
                                                beta);
  stan::math::test::compare_cpu_opencl_prim_rev(exponential_lpdf_functor_propto,
                                                y, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(exponential_lpdf_functor,
                                                y.transpose().eval(), beta);
  stan::math::test::compare_cpu_opencl_prim_rev(exponential_lpdf_functor_propto,
                                                y, beta.transpose().eval());
}

TEST(ProbDistributionsExponential, opencl_broadcast_y) {
  int N = 3;

  double y = 0.5;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      exponential_lpdf_functor, y, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      exponential_lpdf_functor_propto, y, beta);
}

TEST(ProbDistributionsExponential, opencl_broadcast_beta) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double beta = 9.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      exponential_lpdf_functor, y, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      exponential_lpdf_functor_propto, y, beta);
}

TEST(ProbDistributionsExponential, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> beta
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(exponential_lpdf_functor, y,
                                                beta);
  stan::math::test::compare_cpu_opencl_prim_rev(exponential_lpdf_functor_propto,
                                                y, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      exponential_lpdf_functor, y.transpose().eval(), beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(exponential_lpdf_functor_propto,
                                                y.transpose().eval(),
                                                beta.transpose().eval());
}

#endif
