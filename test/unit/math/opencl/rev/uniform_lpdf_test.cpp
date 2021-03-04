#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsUniform, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.3, 0.8;
  Eigen::VectorXd y_value(N);
  y_value << 0.3, NAN, 0.5;

  Eigen::VectorXd alpha(N);
  alpha << 0.2, -0.8, -2.0;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value(N);
  alpha_value << 0.2, INFINITY, -1.5;

  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, -1.0;
  Eigen::VectorXd beta_size(N - 1);
  beta_size << 0.3, 0.8;
  Eigen::VectorXd beta_value1(N);
  beta_value1 << 0.3, INFINITY, 0.5;
  Eigen::VectorXd beta_value2(N);
  beta_value2 << 0.3, -0.8, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);
  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value_cl(alpha_value);
  stan::math::matrix_cl<double> beta_cl(beta);
  stan::math::matrix_cl<double> beta_size_cl(beta_size);
  stan::math::matrix_cl<double> beta_value1_cl(beta_value1);
  stan::math::matrix_cl<double> beta_value2_cl(beta_value2);

  EXPECT_NO_THROW(stan::math::uniform_lpdf(y_cl, alpha_cl, beta_cl));

  EXPECT_THROW(stan::math::uniform_lpdf(y_size_cl, alpha_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::uniform_lpdf(y_cl, alpha_size_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::uniform_lpdf(y_cl, alpha_cl, beta_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::uniform_lpdf(y_value_cl, alpha_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::uniform_lpdf(y_cl, alpha_value_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::uniform_lpdf(y_cl, alpha_cl, beta_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::uniform_lpdf(y_cl, alpha_cl, beta_value2_cl),
               std::domain_error);
}

auto uniform_lpdf_functor
    = [](const auto& y, const auto& alpha, const auto& beta) {
        return stan::math::uniform_lpdf(y, alpha, beta);
      };
auto uniform_lpdf_functor_propto
    = [](const auto& y, const auto& alpha, const auto& beta) {
        return stan::math::uniform_lpdf<true>(y, alpha, beta);
      };

TEST(ProbDistributionsUniform, opencl_matches_cpu_small) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.2, 0.4, -1.0;
  Eigen::VectorXd beta(N);
  beta << 0.4, 1.5, 5.0;

  stan::math::test::compare_cpu_opencl_prim_rev(uniform_lpdf_functor, y, alpha,
                                                beta);
  stan::math::test::compare_cpu_opencl_prim_rev(uniform_lpdf_functor_propto, y,
                                                alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      uniform_lpdf_functor, y.transpose().eval(), alpha.transpose().eval(),
      beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      uniform_lpdf_functor_propto, y.transpose().eval(),
      alpha.transpose().eval(), beta.transpose().eval());
}

TEST(ProbDistributionsUniform, opencl_matches_cpu_small_y_out_of_bounds) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 14.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.2, 0.4, -1.0;
  Eigen::VectorXd beta(N);
  beta << 0.4, 1.5, 5.0;

  stan::math::test::compare_cpu_opencl_prim_rev(uniform_lpdf_functor, y, alpha,
                                                beta);
  stan::math::test::compare_cpu_opencl_prim_rev(uniform_lpdf_functor_propto, y,
                                                alpha, beta);
}

TEST(ProbDistributionsUniform, opencl_broadcast_y) {
  int N = 3;

  double y_scal = 0.7;
  Eigen::VectorXd alpha(N);
  alpha << 0.5, -1.2, -1.0;
  Eigen::VectorXd beta(N);
  beta << 1.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(uniform_lpdf_functor,
                                                         y_scal, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      uniform_lpdf_functor_propto, y_scal, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      uniform_lpdf_functor, y_scal, alpha.transpose().eval(), beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      uniform_lpdf_functor_propto, y_scal, alpha, beta.transpose().eval());
}

TEST(ProbDistributionsUniform, opencl_broadcast_alpha) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double alpha_scal = -12.3;
  Eigen::VectorXd beta(N);
  beta << 1.3, 2.8, 3.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(uniform_lpdf_functor,
                                                         y, alpha_scal, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      uniform_lpdf_functor_propto, y, alpha_scal, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      uniform_lpdf_functor, y.transpose().eval(), alpha_scal, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      uniform_lpdf_functor_propto, y, alpha_scal, beta.transpose().eval());
}

TEST(ProbDistributionsUniform, opencl_broadcast_beta) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << -0.3, 0.5, -1.0;
  double beta_scal = 12.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(uniform_lpdf_functor,
                                                         y, alpha, beta_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      uniform_lpdf_functor_propto, y, alpha, beta_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      uniform_lpdf_functor, y.transpose().eval(), alpha, beta_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      uniform_lpdf_functor_propto, y, alpha.transpose().eval(), beta_scal);
}

TEST(ProbDistributionsUniform, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> beta
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs()
        + alpha.array();
  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs()
            * (beta.array() - alpha.array())
        + beta.array();

  stan::math::test::compare_cpu_opencl_prim_rev(uniform_lpdf_functor, y, alpha,
                                                beta);
  stan::math::test::compare_cpu_opencl_prim_rev(uniform_lpdf_functor_propto, y,
                                                alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      uniform_lpdf_functor, y.transpose().eval(), alpha.transpose().eval(),
      beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      uniform_lpdf_functor_propto, y.transpose().eval(),
      alpha.transpose().eval(), beta.transpose().eval());
}

#endif
