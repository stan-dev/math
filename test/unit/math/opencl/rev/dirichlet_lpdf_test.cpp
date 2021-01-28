#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsDirichlet, error_checking) {
  int N = 3;
  int M = 2;

  Eigen::MatrixXd theta1(N, M);
  theta1 << 0.1, 0.3, 0.4, 0.7, 0.5, 0;
  Eigen::VectorXd theta2(N);
  theta2 << 0.1, 0.4, 0.5;
  Eigen::MatrixXd theta_size1(N - 1, M);
  theta_size1 << 0.4, 0.5, 0.6, 0.5;
  Eigen::MatrixXd theta_size2(N, M + 1);
  theta_size2 << 0.4, 0.2, 0.4, 0.2, 0.4, 0.2, 0.4, 0.4, 0.4;
  Eigen::VectorXd theta_value1(N);
  theta_value1 << -0.1, 0.5, 0.6;
  Eigen::VectorXd theta_value2(N);
  theta_value2 << 0.1, 0.4, 0.4;

  Eigen::MatrixXd alpha1(N, M);
  alpha1 << 0.3, 0.8, 1.0, 0.12, 0.5, 12.3;
  Eigen::VectorXd alpha2(N);
  alpha2 << 0.3, 0.8, 1.0;
  Eigen::MatrixXd alpha_size1(N - 1, M);
  alpha_size1 << 0.3, 0.8, 0.3, 0.8;
  Eigen::MatrixXd alpha_size2(N, M + 1);
  alpha_size2 << 0.3, 0.8, 1.0, 0.12, 0.5, 12.3, 1, 2, 3;
  Eigen::VectorXd alpha_value(N);
  alpha_value << 0.3, 0.0, 0.5;

  stan::math::matrix_cl<double> theta1_cl(theta1);
  stan::math::matrix_cl<double> theta2_cl(theta2);
  stan::math::matrix_cl<double> theta_size1_cl(theta_size1);
  stan::math::matrix_cl<double> theta_size2_cl(theta_size2);
  stan::math::matrix_cl<double> theta_value1_cl(theta_value1);
  stan::math::matrix_cl<double> theta_value2_cl(theta_value2);

  stan::math::matrix_cl<double> alpha1_cl(alpha1);
  stan::math::matrix_cl<double> alpha2_cl(alpha2);
  stan::math::matrix_cl<double> alpha_size1_cl(alpha_size1);
  stan::math::matrix_cl<double> alpha_size2_cl(alpha_size2);
  stan::math::matrix_cl<double> alpha_value_cl(alpha_value);

  EXPECT_NO_THROW(stan::math::dirichlet_lpdf(theta1_cl, alpha1_cl));
  EXPECT_NO_THROW(stan::math::dirichlet_lpdf(theta1_cl, alpha2_cl));
  EXPECT_NO_THROW(stan::math::dirichlet_lpdf(theta2_cl, alpha1_cl));

  EXPECT_THROW(stan::math::dirichlet_lpdf(theta_size1_cl, alpha1_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::dirichlet_lpdf(theta_size2_cl, alpha1_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::dirichlet_lpdf(theta1_cl, alpha_size1_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::dirichlet_lpdf(theta1_cl, alpha_size2_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::dirichlet_lpdf(theta_value1_cl, alpha1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::dirichlet_lpdf(theta_value2_cl, alpha1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::dirichlet_lpdf(theta1_cl, alpha_value_cl),
               std::domain_error);
}

auto dirichlet_lpdf_functor = [](const auto& theta, const auto& alpha) {
  return stan::math::dirichlet_lpdf(theta, alpha);
};
auto dirichlet_lpdf_functor_propto = [](const auto& theta, const auto& alpha) {
  return stan::math::dirichlet_lpdf<true>(theta, alpha);
};

TEST(ProbDistributionsDirichlet, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd theta1(N);
  theta1 << 0.5, 0.4, 0.1;
  Eigen::VectorXd theta2(N);
  theta2 << 0.6, 0.2, 0.2;
  std::vector<Eigen::VectorXd> theta{theta1, theta2};
  Eigen::VectorXd alpha1(N);
  alpha1 << 0.5, 0.1, 12.3;
  Eigen::VectorXd alpha2(N);
  alpha2 << 2.1, 3.4, 2.3;
  std::vector<Eigen::VectorXd> alpha{alpha1, alpha2};

  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor, theta,
                                                alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor_propto,
                                                theta, alpha);

  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor, theta1,
                                                alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor_propto,
                                                theta1, alpha);

  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor, theta,
                                                alpha1);
  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor_propto,
                                                theta, alpha1);

  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor, theta1,
                                                alpha1);
  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor_propto,
                                                theta1, alpha1);
}

TEST(ProbDistributionsDirichlet, opencl_matches_cpu_big) {
  int N = 153;
  int M = 11;
  Eigen::VectorXd theta1;
  Eigen::VectorXd alpha1;
  std::vector<Eigen::VectorXd> theta;
  std::vector<Eigen::VectorXd> alpha;

  for (int i = 0; i < M; i++) {
    theta1 = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
    theta1 /= theta1.sum();
    alpha1 = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
    theta.push_back(theta1);
    alpha.push_back(alpha1);
  }

  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor, theta,
                                                alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor_propto,
                                                theta, alpha);

  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor, theta1,
                                                alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor_propto,
                                                theta1, alpha);

  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor, theta,
                                                alpha1);
  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor_propto,
                                                theta, alpha1);

  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor, theta1,
                                                alpha1);
  stan::math::test::compare_cpu_opencl_prim_rev(dirichlet_lpdf_functor_propto,
                                                theta1, alpha1);
}

#endif
