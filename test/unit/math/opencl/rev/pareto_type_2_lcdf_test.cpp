#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsParetoType2Lcdf, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.3, 0.8;
  Eigen::VectorXd y_value(N);
  y_value << 0.3, -3, 0.5;

  Eigen::VectorXd mu(N);
  mu << 0.2, 0.5, -1.0;
  Eigen::VectorXd mu_size(N - 1);
  mu_size << 0.3, 0.8;
  Eigen::VectorXd mu_value(N);
  mu_value << 0.2, 1.3, 0.5;

  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.8, 1.0;
  Eigen::VectorXd lambda_size(N - 1);
  lambda_size << 0.3, 0.8;
  Eigen::VectorXd lambda_value(N);
  lambda_value << 0.3, 0, 0.5;

  Eigen::VectorXd alpha(N);
  alpha << 0.4, 0.4, 1.4;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value(N);
  alpha_value << 0.3, 0, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);
  stan::math::matrix_cl<double> mu_cl(mu);
  stan::math::matrix_cl<double> mu_size_cl(mu_size);
  stan::math::matrix_cl<double> mu_value_cl(mu_value);
  stan::math::matrix_cl<double> lambda_cl(lambda);
  stan::math::matrix_cl<double> lambda_size_cl(lambda_size);
  stan::math::matrix_cl<double> lambda_value_cl(lambda_value);
  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value_cl(alpha_value);

  EXPECT_NO_THROW(
      stan::math::pareto_type_2_lcdf(y_cl, mu_cl, lambda_cl, alpha_cl));

  EXPECT_THROW(
      stan::math::pareto_type_2_lcdf(y_size_cl, mu_cl, lambda_cl, alpha_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::pareto_type_2_lcdf(y_cl, mu_size_cl, lambda_cl, alpha_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::pareto_type_2_lcdf(y_cl, mu_cl, lambda_size_cl, alpha_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::pareto_type_2_lcdf(y_cl, mu_cl, lambda_cl, alpha_size_cl),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::pareto_type_2_lcdf(y_value_cl, mu_cl, lambda_cl, alpha_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::pareto_type_2_lcdf(y_cl, mu_value_cl, lambda_cl, alpha_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::pareto_type_2_lcdf(y_cl, mu_cl, lambda_value_cl, alpha_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::pareto_type_2_lcdf(y_cl, mu_cl, lambda_cl, alpha_value_cl),
      std::domain_error);
}

auto pareto_type_2_lcdf_functor
    = [](const auto& y, const auto& mu, const auto& lambda, const auto& alpha) {
        return stan::math::pareto_type_2_lcdf(y, mu, lambda, alpha);
      };

TEST(ProbDistributionsParetoType2Lcdf, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.5, 1.8, 1.4;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.4, 1.1;

  stan::math::test::compare_cpu_opencl_prim_rev(pareto_type_2_lcdf_functor, y,
                                                mu, lambda, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(
      pareto_type_2_lcdf_functor, y.transpose().eval(), mu.transpose().eval(),
      lambda.transpose().eval(), alpha.transpose().eval());
}

TEST(ProbDistributionsParetoType2Lcdf, opencl_broadcast_y) {
  int N = 3;

  double y_scal = 12.3;
  Eigen::VectorXd mu(N);
  mu << 0.5, 1.2, 1.0;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.4, 1.1;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      pareto_type_2_lcdf_functor, y_scal, mu, lambda, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      pareto_type_2_lcdf_functor, y_scal, mu.transpose().eval(), lambda,
      alpha.transpose().eval());
}

TEST(ProbDistributionsParetoType2Lcdf, opencl_broadcast_mu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double mu_scal = 0.1;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.4, 1.1;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      pareto_type_2_lcdf_functor, y, mu_scal, lambda, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      pareto_type_2_lcdf_functor, y.transpose().eval(), mu_scal, lambda,
      alpha.transpose().eval());
}

TEST(ProbDistributionsParetoType2Lcdf, opencl_broadcast_lambda) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.5, 1.8, 1.1;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  double lambda_scal = 12.3;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.4, 1.1;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      pareto_type_2_lcdf_functor, y, mu, lambda_scal, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      pareto_type_2_lcdf_functor, y.transpose().eval(), mu, lambda_scal,
      alpha.transpose().eval());
}

TEST(ProbDistributionsParetoType2Lcdf, opencl_broadcast_alpha) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.5, 1.8, 1.1;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.4, 1.1;
  double alpha_scal = 12.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      pareto_type_2_lcdf_functor, y, mu, lambda, alpha_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      pareto_type_2_lcdf_functor, y.transpose().eval(), mu,
      lambda.transpose().eval(), alpha_scal);
}

TEST(ProbDistributionsParetoType2Lcdf, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> mu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs()
        + mu.array();
  Eigen::Matrix<double, Eigen::Dynamic, 1> lambda
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs().array()
        + 0.1;
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(pareto_type_2_lcdf_functor, y,
                                                mu, lambda, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(
      pareto_type_2_lcdf_functor, y.transpose().eval(), mu.transpose().eval(),
      lambda.transpose().eval(), alpha.transpose().eval());
}

#endif
