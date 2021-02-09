#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsBetaBinomial, error_checking) {
  int N = 3;

  std::vector<int> n{2, 0, 12};
  std::vector<int> n_size{1, 0, 1, 0};
  std::vector<int> N_{2, 0, 123};
  std::vector<int> N_size{100, 100, 100, 100};
  std::vector<int> N_value{2, -1, 23};
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 1.8, 1.3;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value1(N);
  alpha_value1 << 0, 0.4, 0.5;
  Eigen::VectorXd alpha_value2(N);
  alpha_value2 << 0.3, INFINITY, 0.5;
  Eigen::VectorXd beta(N);
  beta << 0.3, 1.8, 1.2;
  Eigen::VectorXd beta_size(N - 1);
  beta_size << 0.3, 0.8;
  Eigen::VectorXd beta_value1(N);
  beta_value1 << 0, 0.4, 0.5;
  Eigen::VectorXd beta_value2(N);
  beta_value2 << 0.3, INFINITY, 0.5;

  stan::math::matrix_cl<int> n_cl(n);
  stan::math::matrix_cl<int> n_size_cl(n_size);
  stan::math::matrix_cl<int> N_cl(N_);
  stan::math::matrix_cl<int> N_size_cl(N_size);
  stan::math::matrix_cl<int> N_value_cl(N_value);
  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value1_cl(alpha_value1);
  stan::math::matrix_cl<double> alpha_value2_cl(alpha_value2);
  stan::math::matrix_cl<double> beta_cl(beta);
  stan::math::matrix_cl<double> beta_size_cl(beta_size);
  stan::math::matrix_cl<double> beta_value1_cl(beta_value1);
  stan::math::matrix_cl<double> beta_value2_cl(beta_value2);

  EXPECT_NO_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_cl, beta_cl));

  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_size_cl, N_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_size_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_size_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_cl, beta_size_cl),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_value_cl, alpha_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_value1_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_value2_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_cl, beta_value1_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_cl, beta_value2_cl),
      std::domain_error);
}

auto beta_binomial_lpmf_functor
    = [](const auto& n, const auto& N, const auto& alpha, const auto& beta) {
        return stan::math::beta_binomial_lpmf(n, N, alpha, beta);
      };
auto beta_binomial_lpmf_functor_propto
    = [](const auto& n, const auto& N, const auto& alpha, const auto& beta) {
        return stan::math::beta_binomial_lpmf<true>(n, N, alpha, beta);
      };

TEST(ProbDistributionsBetaBinomial, opencl_matches_cpu_small) {
  int N_ = 3;

  std::vector<int> n{2, 0, 12};
  std::vector<int> N{2, 0, 123};
  Eigen::VectorXd alpha(N_);
  alpha << 0.3, 1.8, 1.3;
  Eigen::VectorXd beta(N_);
  beta << 0.3, 1.8, 1.2;

  stan::math::test::compare_cpu_opencl_prim_rev(beta_binomial_lpmf_functor, n,
                                                N, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      beta_binomial_lpmf_functor_propto, n, N, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(beta_binomial_lpmf_functor, n,
                                                N, alpha.transpose().eval(),
                                                beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      beta_binomial_lpmf_functor_propto, n, N, alpha.transpose().eval(),
      beta.transpose().eval());
}

TEST(ProbDistributionsBetaBinomial, opencl_broadcast_n) {
  int N_ = 3;

  int n = 1;
  std::vector<int> N{2, 0, 123};
  Eigen::VectorXd alpha(N_);
  alpha << 0.3, 1.8, 1.3;
  Eigen::VectorXd beta(N_);
  beta << 0.3, 1.8, 1.2;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      beta_binomial_lpmf_functor, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      beta_binomial_lpmf_functor_propto, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      beta_binomial_lpmf_functor, n, N, alpha.transpose().eval(),
      beta.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      beta_binomial_lpmf_functor_propto, n, N, alpha.transpose().eval(),
      beta.transpose().eval());
}

TEST(ProbDistributionsBetaBinomial, opencl_broadcast_N) {
  int N_ = 3;

  std::vector<int> n{2, 0, 12};
  int N = 15;
  Eigen::VectorXd alpha(N_);
  alpha << 0.3, 1.8, 1.3;
  Eigen::VectorXd beta(N_);
  beta << 0.3, 1.8, 1.2;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      beta_binomial_lpmf_functor, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      beta_binomial_lpmf_functor_propto, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      beta_binomial_lpmf_functor, n, N, alpha.transpose().eval(),
      beta.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      beta_binomial_lpmf_functor_propto, n, N, alpha.transpose().eval(),
      beta.transpose().eval());
}

TEST(ProbDistributionsBetaBinomial, opencl_broadcast_alpha) {
  int N_ = 3;

  std::vector<int> n{2, 0, 12};
  std::vector<int> N{2, 0, 123};
  double alpha = 1.1;
  Eigen::VectorXd beta(N_);
  beta << 0.3, 1.8, 1.2;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      beta_binomial_lpmf_functor, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      beta_binomial_lpmf_functor_propto, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      beta_binomial_lpmf_functor, n, N, alpha, beta.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      beta_binomial_lpmf_functor_propto, n, N, alpha, beta.transpose().eval());
}

TEST(ProbDistributionsBetaBinomial, opencl_broadcast_beta) {
  int N_ = 3;

  std::vector<int> n{2, 0, 12};
  std::vector<int> N{2, 0, 123};
  Eigen::VectorXd alpha(N_);
  alpha << 0.3, 1.8, 1.3;
  double beta = 1.2;
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      beta_binomial_lpmf_functor, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      beta_binomial_lpmf_functor_propto, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      beta_binomial_lpmf_functor, n, N, alpha.transpose().eval(), beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      beta_binomial_lpmf_functor_propto, n, N, alpha.transpose().eval(), beta);
}

TEST(ProbDistributionsBetaBinomial, opencl_matches_cpu_big) {
  int N_ = 153;

  std::vector<int> n(N_);
  std::vector<int> N(N_);
  for (int i = 0; i < N_; i++) {
    n[i] = Eigen::Array<int, Eigen::Dynamic, 1>::Random(1, 1).abs()(0) % 123;
    N[i] = Eigen::Array<int, Eigen::Dynamic, 1>::Random(1, 1).abs()(0) % 123
           + 123;
  }
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N_, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> beta
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N_, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(beta_binomial_lpmf_functor, n,
                                                N, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      beta_binomial_lpmf_functor_propto, n, N, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(beta_binomial_lpmf_functor, n,
                                                N, alpha.transpose().eval(),
                                                beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      beta_binomial_lpmf_functor_propto, n, N, alpha.transpose().eval(),
      beta.transpose().eval());
}

#endif
