#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsUniformLccdf, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.3, 0.8;
  Eigen::VectorXd y_value(N);
  y_value << 0.3, NAN, 0.5;

  Eigen::VectorXd alpha(N);
  alpha << 0.1, 0.8, -1.0;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value(N);
  alpha_value << 0.1, -INFINITY, 0.5;

  Eigen::VectorXd beta(N);
  beta << 0.3, 1.8, 1.0;
  Eigen::VectorXd beta_size(N - 1);
  beta_size << 0.3, 0.8;
  Eigen::VectorXd beta_value(N);
  beta_value << 0.3, 0.8, INFINITY;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);
  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value_cl(alpha_value);
  stan::math::matrix_cl<double> beta_cl(beta);
  stan::math::matrix_cl<double> beta_size_cl(beta_size);
  stan::math::matrix_cl<double> beta_value_cl(beta_value);

  EXPECT_NO_THROW(stan::math::uniform_lccdf(y_cl, alpha_cl, beta_cl));

  EXPECT_THROW(stan::math::uniform_lccdf(y_size_cl, alpha_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::uniform_lccdf(y_cl, alpha_size_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::uniform_lccdf(y_cl, alpha_cl, beta_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::uniform_lccdf(y_value_cl, alpha_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::uniform_lccdf(y_cl, alpha_value_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::uniform_lccdf(y_cl, alpha_cl, beta_value_cl),
               std::domain_error);
}

auto uniform_lccdf_functor
    = [](const auto& y, const auto& alpha, const auto& beta) {
        return stan::math::uniform_lccdf(y, alpha, beta);
      };

TEST(ProbDistributionsUniformLccdf, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.1, 0.5, -1.0;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(uniform_lccdf_functor, y, alpha,
                                                beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      uniform_lccdf_functor, y.transpose().eval(), alpha.transpose().eval(),
      beta.transpose().eval());
}

TEST(ProbDistributionsUniformLccdf, opencl_matches_cpu_small_y_neg_inf) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.3, -INFINITY, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.1, 0.5, -1.0;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(uniform_lccdf_functor, y, alpha,
                                                beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      uniform_lccdf_functor, y.transpose().eval(), alpha.transpose().eval(),
      beta.transpose().eval());
}

TEST(ProbDistributionsUniformLccdf, opencl_broadcast_y) {
  int N = 3;

  double y_scal = 12.3;
  Eigen::VectorXd alpha(N);
  alpha << 0.1, 0.5, -1.0;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(uniform_lccdf_functor,
                                                         y_scal, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      uniform_lccdf_functor, y_scal, alpha.transpose().eval(), beta);
}

TEST(ProbDistributionsUniformLccdf, opencl_broadcast_alpha) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double alpha_scal = -0.1;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(uniform_lccdf_functor,
                                                         y, alpha_scal, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      uniform_lccdf_functor, y.transpose().eval(), alpha_scal, beta);
}

TEST(ProbDistributionsUniformLccdf, opencl_broadcast_beta) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.1, 0.5, -1.0;
  double beta_scal = 12.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(uniform_lccdf_functor,
                                                         y, alpha, beta_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      uniform_lccdf_functor, y.transpose().eval(), alpha, beta_scal);
}

TEST(ProbDistributionsUniformLccdf, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> beta
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs()
        + alpha.array();

  stan::math::test::compare_cpu_opencl_prim_rev(uniform_lccdf_functor, y, alpha,
                                                beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      uniform_lccdf_functor, y.transpose().eval(), alpha.transpose().eval(),
      beta.transpose().eval());
}
#endif
