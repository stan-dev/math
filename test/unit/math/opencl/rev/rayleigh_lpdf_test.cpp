#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsRayleigh, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.1, 0.5, 12.3;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.1, 0.5;
  Eigen::VectorXd y_value(N);
  y_value << -0.1, 0.5, 0.99;

  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 12.4;
  Eigen::VectorXd sigma_size(N - 1);
  sigma_size << 0.3, 0.8;
  Eigen::VectorXd sigma_value(N);
  sigma_value << 0.3, -0.8, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);

  stan::math::matrix_cl<double> sigma_cl(sigma);
  stan::math::matrix_cl<double> sigma_size_cl(sigma_size);
  stan::math::matrix_cl<double> sigma_value_cl(sigma_value);

  EXPECT_NO_THROW(stan::math::rayleigh_lpdf(y_cl, sigma_cl));

  EXPECT_THROW(stan::math::rayleigh_lpdf(y_size_cl, sigma_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::rayleigh_lpdf(y_cl, sigma_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::rayleigh_lpdf(y_value_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::rayleigh_lpdf(y_cl, sigma_value_cl),
               std::domain_error);
}

auto rayleigh_lpdf_functor = [](const auto& y, const auto& sigma) {
  return stan::math::rayleigh_lpdf(y, sigma);
};
auto rayleigh_lpdf_functor_propto = [](const auto& y, const auto& sigma) {
  return stan::math::rayleigh_lpdf<true>(y, sigma);
};

TEST(ProbDistributionsRayleigh, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.1, 0.5, 12.3;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 21.0;

  stan::math::test::compare_cpu_opencl_prim_rev(rayleigh_lpdf_functor, y,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(rayleigh_lpdf_functor_propto, y,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(rayleigh_lpdf_functor,
                                                y.transpose().eval(), sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(rayleigh_lpdf_functor_propto, y,
                                                sigma.transpose().eval());
}

TEST(ProbDistributionsRayleigh, opencl_broadcast_y) {
  int N = 3;

  double y = 0.5;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(rayleigh_lpdf_functor,
                                                         y, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      rayleigh_lpdf_functor_propto, y, sigma);
}

TEST(ProbDistributionsRayleigh, opencl_broadcast_sigma) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double sigma = 9.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(rayleigh_lpdf_functor,
                                                         y, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      rayleigh_lpdf_functor_propto, y, sigma);
}

TEST(ProbDistributionsRayleigh, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> sigma
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(rayleigh_lpdf_functor, y,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(rayleigh_lpdf_functor_propto, y,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(
      rayleigh_lpdf_functor, y.transpose().eval(), sigma.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(rayleigh_lpdf_functor_propto,
                                                y.transpose().eval(),
                                                sigma.transpose().eval());
}

#endif
