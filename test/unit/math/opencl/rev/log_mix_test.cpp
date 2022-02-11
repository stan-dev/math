#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(OpenCLLogMix, error_checking) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd theta(N);
  theta << 0.1, 0.4, 0.5;
  Eigen::VectorXd theta_size1(N - 1);
  theta_size1 << 0.4, 0.5;
  Eigen::VectorXd theta_value1(N);
  theta_value1 << -0.1, 0.5, 0.6;
  Eigen::VectorXd theta_value2(N);
  theta_value2 << 0.1, 0.4, 1.4;

  Eigen::MatrixXd lambda1(N, M);
  lambda1 << 0.3, 0.8, -1.0, 0.12, 0.5, 12.3;
  Eigen::VectorXd lambda2(N);
  lambda2 << 0.3, 0.8, 1.0;
  Eigen::MatrixXd lambda_size1(N - 1, M);
  lambda_size1 << 0.3, 0.8, 0.3, 0.8;
  Eigen::VectorXd lambda_value(N);
  lambda_value << 0.3, 0.0, NAN;

  stan::math::matrix_cl<double> theta_cl(theta);
  stan::math::matrix_cl<double> theta_size1_cl(theta_size1);
  stan::math::matrix_cl<double> theta_value1_cl(theta_value1);
  stan::math::matrix_cl<double> theta_value2_cl(theta_value2);

  stan::math::matrix_cl<double> lambda1_cl(lambda1);
  stan::math::matrix_cl<double> lambda2_cl(lambda2);
  stan::math::matrix_cl<double> lambda_size1_cl(lambda_size1);
  stan::math::matrix_cl<double> lambda_value_cl(lambda_value);

  EXPECT_NO_THROW(stan::math::log_mix(theta_cl, lambda1_cl));
  EXPECT_NO_THROW(stan::math::log_mix(theta_cl, lambda2_cl));

  EXPECT_THROW(stan::math::log_mix(theta_size1_cl, lambda1_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::log_mix(theta_cl, lambda_size1_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::log_mix(theta_value1_cl, lambda1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::log_mix(theta_value2_cl, lambda1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::log_mix(theta_cl, lambda_value_cl),
               std::domain_error);
}

auto log_mix_functor = [](const auto& theta, const auto& lambda) {
  return stan::math::log_mix(theta, lambda);
};

TEST(OpenCLLogMix, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd theta1(N);
  theta1 << 0.5, 0.4, 0.1;
  Eigen::RowVectorXd theta2(N);
  theta2 << 0.6, 0.2, 0.2;
  Eigen::VectorXd lambda1(N);
  lambda1 << 0.5, 0.1, 12.3;
  Eigen::RowVectorXd lambda2(N);
  lambda2 << 2.1, 3.4, 2.3;
  std::vector<Eigen::VectorXd> lambda3{lambda1, lambda2};
  std::vector<Eigen::RowVectorXd> lambda4{lambda2, lambda1};

  stan::math::test::compare_cpu_opencl_prim_rev(log_mix_functor, theta1,
                                                lambda1);
  stan::math::test::compare_cpu_opencl_prim_rev(log_mix_functor, theta1,
                                                lambda1);
  stan::math::test::compare_cpu_opencl_prim_rev(log_mix_functor, theta1,
                                                lambda3);
  stan::math::test::compare_cpu_opencl_prim_rev(log_mix_functor, theta1,
                                                lambda3);
  stan::math::test::compare_cpu_opencl_prim_rev(log_mix_functor, theta1,
                                                lambda4);
  stan::math::test::compare_cpu_opencl_prim_rev(log_mix_functor, theta1,
                                                lambda4);
}

TEST(OpenCLLogMix, opencl_matches_cpu_big) {
  int N = 153;
  int M = 11;

  Eigen::VectorXd theta1
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::VectorXd lambda1;
  std::vector<Eigen::VectorXd> lambda3;
  std::vector<Eigen::RowVectorXd> lambda4;

  for (int i = 0; i < M; i++) {
    lambda1 = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
    lambda3.push_back(lambda1);
    lambda4.push_back(lambda1);
  }
  Eigen::RowVectorXd theta2 = theta1;
  Eigen::RowVectorXd lambda2 = lambda1;

  stan::math::test::compare_cpu_opencl_prim_rev(log_mix_functor, theta1,
                                                lambda1);
  stan::math::test::compare_cpu_opencl_prim_rev(log_mix_functor, theta1,
                                                lambda1);
  stan::math::test::compare_cpu_opencl_prim_rev(log_mix_functor, theta1,
                                                lambda3);
  stan::math::test::compare_cpu_opencl_prim_rev(log_mix_functor, theta1,
                                                lambda3);
  stan::math::test::compare_cpu_opencl_prim_rev(log_mix_functor, theta1,
                                                lambda4);
  stan::math::test::compare_cpu_opencl_prim_rev(log_mix_functor, theta1,
                                                lambda4);
}

#endif
