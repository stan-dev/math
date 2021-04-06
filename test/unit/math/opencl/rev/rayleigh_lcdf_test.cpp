#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsRayleighLcdf, error_checking) {
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

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);
  stan::math::matrix_cl<double> mu_cl(mu);
  stan::math::matrix_cl<double> mu_size_cl(mu_size);
  stan::math::matrix_cl<double> mu_value_cl(mu_value);

  EXPECT_NO_THROW(stan::math::rayleigh_lcdf(y_cl, mu_cl));

  EXPECT_THROW(stan::math::rayleigh_lcdf(y_size_cl, mu_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::rayleigh_lcdf(y_cl, mu_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::rayleigh_lcdf(y_value_cl, mu_cl), std::domain_error);
  EXPECT_THROW(stan::math::rayleigh_lcdf(y_cl, mu_value_cl), std::domain_error);
}

auto rayleigh_lcdf_functor = [](const auto& y, const auto& mu) {
  return stan::math::rayleigh_lcdf(y, mu);
};

TEST(ProbDistributionsRayleighLcdf, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(rayleigh_lcdf_functor, y, mu);
  stan::math::test::compare_cpu_opencl_prim_rev(
      rayleigh_lcdf_functor, y.transpose().eval(), mu.transpose().eval());
}

TEST(ProbDistributionsRayleighLcdf, opencl_broadcast_y) {
  int N = 3;

  double y_scal = 12.3;
  Eigen::VectorXd mu(N);
  mu << 0.5, 1.2, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(rayleigh_lcdf_functor,
                                                         y_scal, mu);
}

TEST(ProbDistributionsRayleighLcdf, opencl_broadcast_mu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double mu_scal = 12.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(rayleigh_lcdf_functor,
                                                         y, mu_scal);
}

TEST(ProbDistributionsRayleighLcdf, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(rayleigh_lcdf_functor, y, mu);
  stan::math::test::compare_cpu_opencl_prim_rev(
      rayleigh_lcdf_functor, y.transpose().eval(), mu.transpose().eval());
}
#endif
