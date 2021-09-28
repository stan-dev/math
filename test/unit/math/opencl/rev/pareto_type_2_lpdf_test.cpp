#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsParetoType2, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << -10, INFINITY, 11.3;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.1, 0.5;
  Eigen::VectorXd y_value(N);
  y_value << 10, 0.5, 9.99;

  Eigen::VectorXd mu(N);
  mu << -10.3, 0.8, 2.1;
  Eigen::VectorXd mu_size(N - 1);
  mu_size << 0.3, 0.8;
  Eigen::VectorXd mu_value(N);
  mu_value << -20.3, NAN, 0.5;

  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.8, 11.0;
  Eigen::VectorXd lambda_size(N - 1);
  lambda_size << 0.3, 0.8;
  Eigen::VectorXd lambda_value1(N);
  lambda_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd lambda_value2(N);
  lambda_value2 << 0.3, INFINITY, 0.5;

  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 13.0;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value1(N);
  alpha_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd alpha_value2(N);
  alpha_value2 << 0.3, INFINITY, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);

  stan::math::matrix_cl<double> mu_cl(mu);
  stan::math::matrix_cl<double> mu_size_cl(mu_size);
  stan::math::matrix_cl<double> mu_value_cl(mu_value);

  stan::math::matrix_cl<double> lambda_cl(lambda);
  stan::math::matrix_cl<double> lambda_size_cl(lambda_size);
  stan::math::matrix_cl<double> lambda_value1_cl(lambda_value1);
  stan::math::matrix_cl<double> lambda_value2_cl(lambda_value2);

  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value1_cl(alpha_value1);
  stan::math::matrix_cl<double> alpha_value2_cl(alpha_value2);

  EXPECT_NO_THROW(
      stan::math::pareto_type_2_lpdf(y_cl, mu_cl, lambda_cl, alpha_cl));

  EXPECT_THROW(
      stan::math::pareto_type_2_lpdf(y_size_cl, mu_cl, lambda_cl, alpha_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::pareto_type_2_lpdf(y_cl, mu_size_cl, lambda_cl, alpha_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::pareto_type_2_lpdf(y_cl, mu_cl, lambda_size_cl, alpha_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::pareto_type_2_lpdf(y_cl, mu_cl, lambda_cl, alpha_size_cl),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::pareto_type_2_lpdf(y_value_cl, mu_cl, lambda_cl, alpha_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::pareto_type_2_lpdf(y_cl, mu_value_cl, lambda_cl, alpha_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::pareto_type_2_lpdf(y_cl, mu_cl, lambda_value1_cl, alpha_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::pareto_type_2_lpdf(y_cl, mu_cl, lambda_cl, alpha_value1_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::pareto_type_2_lpdf(y_cl, mu_cl, lambda_value2_cl, alpha_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::pareto_type_2_lpdf(y_cl, mu_cl, lambda_cl, alpha_value2_cl),
      std::domain_error);
}

auto pareto_type_2_lpdf_functor
    = [](const auto& y, const auto& mu, const auto& lambda, const auto& alpha) {
        return stan::math::pareto_type_2_lpdf(y, mu, lambda, alpha);
      };
auto pareto_type_2_lpdf_functor_propto
    = [](const auto& y, const auto& mu, const auto& lambda, const auto& alpha) {
        return stan::math::pareto_type_2_lpdf<true>(y, mu, lambda, alpha);
      };

TEST(ProbDistributionsParetoType2, opencl_matches_cpu_small) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.9, 31;
  Eigen::VectorXd mu(N);
  mu << -10.3, 0.8, 21.0;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.8, 11.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 13.0;

  stan::math::test::compare_cpu_opencl_prim_rev(pareto_type_2_lpdf_functor, y,
                                                mu, lambda, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(
      pareto_type_2_lpdf_functor_propto, y, mu, lambda, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(
      pareto_type_2_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      lambda.transpose().eval(), alpha.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      pareto_type_2_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), lambda.transpose().eval(),
      alpha.transpose().eval());
}

TEST(ProbDistributionsParetoType2, opencl_broadcast_y) {
  int N = 3;

  double y = 30.3;
  Eigen::VectorXd mu(N);
  mu << -10.3, 0.8, 21.0;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.8, 11.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 13.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      pareto_type_2_lpdf_functor, y, mu, lambda, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      pareto_type_2_lpdf_functor_propto, y, mu, lambda, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      pareto_type_2_lpdf_functor, y, mu.transpose().eval(),
      lambda.transpose().eval(), alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      pareto_type_2_lpdf_functor_propto, y, mu, lambda.transpose().eval(),
      alpha.transpose().eval());
}

TEST(ProbDistributionsParetoType2, opencl_broadcast_mu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 1;
  double mu = -1.1;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.8, 11.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 13.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      pareto_type_2_lpdf_functor, y, mu, lambda, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      pareto_type_2_lpdf_functor_propto, y, mu, lambda, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      pareto_type_2_lpdf_functor, y, mu, lambda.transpose().eval(),
      alpha.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      pareto_type_2_lpdf_functor_propto, y.transpose().eval(), mu, lambda,
      alpha.transpose().eval());
}

TEST(ProbDistributionsParetoType2, opencl_broadcast_lambda) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.9, 31;
  Eigen::VectorXd mu(N);
  mu << -10.3, 0.8, 21.0;
  double lambda = 0.8;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 13.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      pareto_type_2_lpdf_functor, y, mu, lambda, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      pareto_type_2_lpdf_functor_propto, y, mu, lambda, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      pareto_type_2_lpdf_functor, y.transpose().eval(), mu, lambda,
      alpha.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      pareto_type_2_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), lambda, alpha);
}

TEST(ProbDistributionsParetoType2, opencl_broadcast_alpha) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.9, 31;
  Eigen::VectorXd mu(N);
  mu << -10.3, 0.8, 21.0;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.8, 11.0;
  double alpha = 1.2;

  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      pareto_type_2_lpdf_functor, y, mu, lambda, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      pareto_type_2_lpdf_functor_propto, y, mu, lambda, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      pareto_type_2_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      lambda, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      pareto_type_2_lpdf_functor_propto, y, mu.transpose().eval(),
      lambda.transpose().eval(), alpha);
}

TEST(ProbDistributionsParetoType2, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu
      = -Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs()
        + y.array();
  Eigen::Matrix<double, Eigen::Dynamic, 1> lambda
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(pareto_type_2_lpdf_functor, y,
                                                mu, lambda, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(
      pareto_type_2_lpdf_functor_propto, y, mu, lambda, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(
      pareto_type_2_lpdf_functor, y, mu.transpose().eval(),
      lambda.transpose().eval(), alpha.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      pareto_type_2_lpdf_functor_propto, y.transpose().eval(),
      mu.transpose().eval(), lambda.transpose().eval(),
      alpha.transpose().eval());
}

#endif
