#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsStdNormalLcdf, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd y_value(N);
  y_value << 0.3, NAN, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_value_cl(y_value);

  EXPECT_NO_THROW(stan::math::std_normal_lcdf(y_cl));

  EXPECT_THROW(stan::math::std_normal_lcdf(y_value_cl), std::domain_error);
}

auto std_normal_lcdf_functor
    = [](const auto& y) { return stan::math::std_normal_lcdf(y); };

TEST(ProbDistributionsStdNormalLcdf, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(std_normal_lcdf_functor, y);
  stan::math::test::compare_cpu_opencl_prim_rev(std_normal_lcdf_functor,
                                                y.transpose().eval());
}

TEST(ProbDistributionsStdNormalLcdf, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1);

  stan::math::test::compare_cpu_opencl_prim_rev(std_normal_lcdf_functor, y);
  stan::math::test::compare_cpu_opencl_prim_rev(std_normal_lcdf_functor,
                                                y.transpose().eval());
}

TEST(ProbDistributionsStdNormalLcdf, opencl_matches_cpu_tail_branches) {
  Eigen::VectorXd y(12);
  y << -1e100, -40, -6, -4 * stan::math::SQRT_TWO, -1, 0, 0.3, 1, 4, 8, 40,
      1e100;
  stan::math::test::compare_cpu_opencl_prim_rev(std_normal_lcdf_functor, y);
  // Check each value separately so the extreme tail cannot dominate the sum.
  for (Eigen::Index i = 0; i < y.size(); ++i) {
    SCOPED_TRACE(y[i]);
    stan::math::test::compare_cpu_opencl_prim_rev(std_normal_lcdf_functor,
                                                  y.segment(i, 1).eval());
  }
}

TEST(ProbDistributionsStdNormalLcdf, empty_and_extreme_gradient) {
  using namespace stan::math;
  matrix_cl<double> empty(Eigen::VectorXd(0));
  EXPECT_EQ(0, std_normal_lcdf(empty));
  EXPECT_EQ(0, std_normal_lccdf(empty));
  nested_rev_autodiff nested;
  var_value<matrix_cl<double>> y(
      to_matrix_cl(Eigen::VectorXd::Constant(1, -1.5e308)));
  auto lp = std_normal_lcdf(y);
  lp.grad();
  EXPECT_NEAR(1.5e308, from_matrix_cl(y.adj())(0, 0), 1.5e296);
}
#endif
