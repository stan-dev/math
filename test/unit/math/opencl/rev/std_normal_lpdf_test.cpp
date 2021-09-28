#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsStdNormal, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd y_value(N);
  y_value << 0.3, NAN, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_value_cl(y_value);

  EXPECT_NO_THROW(stan::math::std_normal_lpdf(y_cl));
  EXPECT_THROW(stan::math::std_normal_lpdf(y_value_cl), std::domain_error);
}

auto std_normal_lpdf_functor
    = [](const auto& y) { return stan::math::std_normal_lpdf(y); };
auto std_normal_lpdf_functor_propto
    = [](const auto& y) { return stan::math::std_normal_lpdf<true>(y); };

TEST(ProbDistributionsStdNormal, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(std_normal_lpdf_functor, y);
  stan::math::test::compare_cpu_opencl_prim_rev(std_normal_lpdf_functor_propto,
                                                y);
  stan::math::test::compare_cpu_opencl_prim_rev(std_normal_lpdf_functor,
                                                y.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(std_normal_lpdf_functor_propto,
                                                y.transpose().eval());
}

TEST(ProbDistributionsStdNormal, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(std_normal_lpdf_functor, y);
  stan::math::test::compare_cpu_opencl_prim_rev(std_normal_lpdf_functor_propto,
                                                y);
  stan::math::test::compare_cpu_opencl_prim_rev(std_normal_lpdf_functor,
                                                y.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(std_normal_lpdf_functor_propto,
                                                y.transpose().eval());
}

#endif
