#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(muProbDistributionsNegBinomial, error_checking) {
  int N = 3;

  std::vector<int> n{1, 0, 12};
  std::vector<int> n_size{1, 0, 1, 0};
  std::vector<int> n_value{0, 1, -3};

  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.3;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value1(N);
  alpha_value1 << 0.3, -0.3, 0.5;
  Eigen::VectorXd alpha_value2(N);
  alpha_value2 << 0.3, INFINITY, 0.5;

  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.3;
  Eigen::VectorXd beta_size(N - 1);
  beta_size << 0.3, 0.8;
  Eigen::VectorXd beta_value1(N);
  beta_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd beta_value2(N);
  beta_value2 << 0.3, INFINITY, 0.5;

  stan::math::matrix_cl<int> n_cl(n);
  stan::math::matrix_cl<int> n_size_cl(n_size);
  stan::math::matrix_cl<int> n_value_cl(n_value);
  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value1_cl(alpha_value1);
  stan::math::matrix_cl<double> alpha_value2_cl(alpha_value2);
  stan::math::matrix_cl<double> beta_cl(beta);
  stan::math::matrix_cl<double> beta_size_cl(beta_size);
  stan::math::matrix_cl<double> beta_value1_cl(beta_value1);
  stan::math::matrix_cl<double> beta_value2_cl(beta_value2);

  EXPECT_NO_THROW(stan::math::neg_binomial_lpmf(n_cl, alpha_cl, beta_cl));

  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_size_cl, alpha_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_cl, alpha_size_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_cl, alpha_cl, beta_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_value_cl, alpha_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_cl, alpha_value1_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_cl, alpha_value2_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_cl, alpha_cl, beta_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_cl, alpha_cl, beta_value2_cl),
               std::domain_error);
}

auto neg_binomial_lpmf_functor
    = [](const auto& n, const auto& alpha, const auto& beta) {
        return stan::math::neg_binomial_lpmf(n, alpha, beta);
      };
auto neg_binomial_lpmf_functor_propto
    = [](const auto& n, const auto& alpha, const auto& beta) {
        return stan::math::neg_binomial_lpmf<true>(n, alpha, beta);
      };

TEST(muProbDistributionsNegBinomial, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  std::vector<int> n{1, 0, 12};
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.5, 1.8;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.3;

  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_lpmf_functor, n,
                                                alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_lpmf_functor_propto, n, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_lpmf_functor, n,
                                                alpha.transpose().eval(),
                                                beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_lpmf_functor_propto, n, alpha.transpose().eval(),
      beta.transpose().eval());
}

TEST(muProbDistributionsNegBinomial, opencl_broadcast_n) {
  int N = 3;

  int n = 2;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.5, 1.8;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_lpmf_functor, n, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_lpmf_functor_propto, n, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_lpmf_functor, n, alpha.transpose().eval(), beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_lpmf_functor_propto, n, alpha, beta.transpose().eval());
}

TEST(muProbDistributionsNegBinomial, opencl_broadcast_alpha) {
  int N = 3;

  std::vector<int> n{1, 0, 12};
  double alpha = 0.4;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_lpmf_functor, n, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_lpmf_functor_propto, n, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_lpmf_functor, n, alpha, beta.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_lpmf_functor_propto, n, alpha, beta.transpose().eval());
}

TEST(muProbDistributionsNegBinomial, opencl_broadcast_beta) {
  int N = 3;

  std::vector<int> n{1, 0, 12};
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.5, 1.8;
  double beta = 0.4;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_lpmf_functor, n, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_lpmf_functor_propto, n, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_lpmf_functor, n, alpha.transpose().eval(), beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_lpmf_functor_propto, n, alpha.transpose().eval(), beta);
}

TEST(muProbDistributionsNegBinomial, opencl_matches_cpu_big) {
  int N = 153;

  std::vector<int> n(N);
  for (int i = 0; i < N; i++) {
    n[i] = Eigen::Array<int, Eigen::Dynamic, 1>::Random(1, 1).abs()(0) % 1000;
  }
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> beta
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_lpmf_functor, n,
                                                alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_lpmf_functor_propto, n, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_lpmf_functor, n,
                                                alpha.transpose().eval(),
                                                beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_lpmf_functor_propto, n, alpha.transpose().eval(),
      beta.transpose().eval());
}

#endif
