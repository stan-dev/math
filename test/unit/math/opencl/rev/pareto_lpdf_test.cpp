#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsPareto, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << -0.3, 0.8, 11.0;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.3, 0.8;
  Eigen::VectorXd y_value(N);
  y_value << 0.3, NAN, 0.5;

  Eigen::VectorXd y_min(N);
  y_min << 0.3, 0.8, 1.0;
  Eigen::VectorXd y_min_size(N - 1);
  y_min_size << 0.3, 0.8;
  Eigen::VectorXd y_min_value1(N);
  y_min_value1 << 0.3, INFINITY, 0.5;
  Eigen::VectorXd y_min_value2(N);
  y_min_value2 << 0.3, -0.6, 0.5;

  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value1(N);
  alpha_value1 << 0.3, 0, 0.5;
  Eigen::VectorXd alpha_value2(N);
  alpha_value2 << 0.3, INFINITY, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);
  stan::math::matrix_cl<double> y_min_cl(y_min);
  stan::math::matrix_cl<double> y_min_size_cl(y_min_size);
  stan::math::matrix_cl<double> y_min_value1_cl(y_min_value1);
  stan::math::matrix_cl<double> y_min_value2_cl(y_min_value2);
  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value1_cl(alpha_value1);
  stan::math::matrix_cl<double> alpha_value2_cl(alpha_value2);

  EXPECT_NO_THROW(stan::math::pareto_lpdf(y_cl, y_min_cl, alpha_cl));

  EXPECT_THROW(stan::math::pareto_lpdf(y_size_cl, y_min_cl, alpha_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::pareto_lpdf(y_cl, y_min_size_cl, alpha_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::pareto_lpdf(y_cl, y_min_cl, alpha_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::pareto_lpdf(y_value_cl, y_min_cl, alpha_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::pareto_lpdf(y_cl, y_min_value1_cl, alpha_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::pareto_lpdf(y_cl, y_min_cl, alpha_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::pareto_lpdf(y_cl, y_min_value2_cl, alpha_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::pareto_lpdf(y_cl, y_min_cl, alpha_value2_cl),
               std::domain_error);
}

auto pareto_lpdf_functor
    = [](const auto& y, const auto& y_min, const auto& alpha) {
        return stan::math::pareto_lpdf(y, y_min, alpha);
      };
auto pareto_lpdf_functor_propto
    = [](const auto& y, const auto& y_min, const auto& alpha) {
        return stan::math::pareto_lpdf<true>(y, y_min, alpha);
      };

TEST(ProbDistributionsPareto, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd y_min(N);
  y_min << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(pareto_lpdf_functor, y, y_min,
                                                alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(pareto_lpdf_functor_propto, y,
                                                y_min, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(
      pareto_lpdf_functor, y.transpose().eval(), y_min.transpose().eval(),
      alpha.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      pareto_lpdf_functor_propto, y.transpose().eval(),
      y_min.transpose().eval(), alpha.transpose().eval());
}

TEST(ProbDistributionsPareto, opencl_broadcast_y) {
  int N = 3;

  double y_scal = 12.3;
  Eigen::VectorXd y_min(N);
  y_min << 0.5, 1.2, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(pareto_lpdf_functor,
                                                         y_scal, y_min, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      pareto_lpdf_functor_propto, y_scal, y_min, alpha);

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      pareto_lpdf_functor, y_scal, y_min.transpose().eval(), alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      pareto_lpdf_functor_propto, y_scal, y_min, alpha.transpose().eval());
}

TEST(ProbDistributionsPareto, opencl_broadcast_y_min) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double y_min_scal = 12.3;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(pareto_lpdf_functor, y,
                                                         y_min_scal, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      pareto_lpdf_functor_propto, y, y_min_scal, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      pareto_lpdf_functor, y.transpose().eval(), y_min_scal, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      pareto_lpdf_functor_propto, y, y_min_scal, alpha.transpose().eval());
}

TEST(ProbDistributionsPareto, opencl_broadcast_alpha) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd y_min(N);
  y_min << 0.3, 0.8, 1.0;
  double alpha_scal = 12.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(pareto_lpdf_functor, y,
                                                         y_min, alpha_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      pareto_lpdf_functor_propto, y, y_min, alpha_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      pareto_lpdf_functor, y.transpose().eval(), y_min, alpha_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      pareto_lpdf_functor_propto, y, y_min.transpose().eval(), alpha_scal);
}

TEST(ProbDistributionsPareto, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_min
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(pareto_lpdf_functor, y, y_min,
                                                alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(pareto_lpdf_functor_propto, y,
                                                y_min, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(
      pareto_lpdf_functor, y.transpose().eval(), y_min.transpose().eval(),
      alpha.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      pareto_lpdf_functor_propto, y.transpose().eval(),
      y_min.transpose().eval(), alpha.transpose().eval());
}

#endif
