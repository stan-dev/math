#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <stan/math/opencl/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsBinomialLogitGLM, error_checking) {
  using stan::math::binomial_logit_glm_lpmf;
  using stan::math::matrix_cl;

  int N = 3;
  int M = 2;

  std::vector<int> n{1, 0, 1};
  std::vector<int> n_size{1, 0, 1, 0};
  std::vector<int> n_value{0, 1, -23};

  std::vector<int> trials{10, 5, 2};
  std::vector<int> trials_size{1, 0, 1, 0};
  std::vector<int> trials_value{5, 1, -1};

  Eigen::MatrixXd x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Eigen::MatrixXd x_size1(N - 1, M);
  x_size1 << -12, 46, -42, 24;
  Eigen::MatrixXd x_size2(N, M - 1);
  x_size2 << -12, 46, -42;
  Eigen::MatrixXd x_value(N, M);
  x_value << -12, 46, -42, 24, 25, -INFINITY;
  Eigen::VectorXd beta(M);
  beta << 0.3, 2;
  Eigen::VectorXd beta_size(M + 1);
  beta_size << 0.3, 2, 0.4;
  Eigen::VectorXd beta_value(M);
  beta_value << 0.3, INFINITY;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, -0.8, 1.8;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, -0.8;
  Eigen::VectorXd alpha_value(N);
  alpha_value << 0.3, -0.8, NAN;

  matrix_cl<double> x_cl(x);
  matrix_cl<double> x_size1_cl(x_size1);
  matrix_cl<double> x_size2_cl(x_size2);
  matrix_cl<double> x_value_cl(x_value);
  matrix_cl<int> n_cl(n);
  matrix_cl<int> n_size_cl(n_size);
  matrix_cl<int> n_value_cl(n_value);
  matrix_cl<int> trials_cl(trials);
  matrix_cl<int> trials_size_cl(trials_size);
  matrix_cl<int> trials_value_cl(trials_value);
  matrix_cl<double> beta_cl(beta);
  matrix_cl<double> beta_size_cl(beta_size);
  matrix_cl<double> beta_value_cl(beta_value);
  matrix_cl<double> alpha_cl(alpha);
  matrix_cl<double> alpha_size_cl(alpha_size);
  matrix_cl<double> alpha_value_cl(alpha_value);

  EXPECT_NO_THROW(
      binomial_logit_glm_lpmf(n_cl, trials_cl, x_cl, alpha_cl, beta_cl));

  EXPECT_THROW(
      binomial_logit_glm_lpmf(n_size_cl, trials_cl, x_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      binomial_logit_glm_lpmf(n_cl, trials_size_cl, x_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      binomial_logit_glm_lpmf(n_cl, trials_cl, x_size1_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      binomial_logit_glm_lpmf(n_cl, trials_cl, x_size2_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      binomial_logit_glm_lpmf(n_cl, trials_cl, x_cl, alpha_size_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      binomial_logit_glm_lpmf(n_cl, trials_cl, x_cl, alpha_cl, beta_size_cl),
      std::invalid_argument);

  EXPECT_THROW(
      binomial_logit_glm_lpmf(n_value_cl, trials_cl, x_cl, alpha_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      binomial_logit_glm_lpmf(n_cl, trials_value_cl, x_cl, alpha_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      binomial_logit_glm_lpmf(n_cl, trials_cl, x_value_cl, alpha_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      binomial_logit_glm_lpmf(n_cl, trials_cl, x_cl, alpha_value_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      binomial_logit_glm_lpmf(n_cl, trials_cl, x_cl, alpha_cl, beta_value_cl),
      std::domain_error);
}

auto binomial_logit_glm_lpmf_functor
    = [](const auto& n, const auto& trials, const auto& x, const auto& alpha,
         const auto& beta) {
        return stan::math::binomial_logit_glm_lpmf(n, trials, x, alpha, beta);
      };
auto binomial_logit_glm_lpmf_functor_propto
    = [](const auto& n, const auto& trials, const auto& x, const auto& alpha,
         const auto& beta) {
        return stan::math::binomial_logit_glm_lpmf<true>(n, trials, x, alpha,
                                                         beta);
      };

TEST(ProbDistributionsBinomialLogitGLM, opencl_matches_cpu_small_simple) {
  int N = 3;
  int M = 2;

  std::vector<int> n{0, 1, 0};
  std::vector<int> trials{10, 4, 15};
  Eigen::MatrixXd x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Eigen::VectorXd beta(M);
  beta << 0.3, 2;
  double alpha = 0.3;

  stan::math::test::compare_cpu_opencl_prim_rev(binomial_logit_glm_lpmf_functor,
                                                n, trials, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      binomial_logit_glm_lpmf_functor_propto, n, trials, x, alpha, beta);
}

TEST(ProbDistributionsBinomialLogitGLM, opencl_broadcast_n) {
  int N = 3;
  int M = 2;

  int n_scal = 1;
  std::vector<int> trials{10, 4, 15};
  Eigen::MatrixXd x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Eigen::VectorXd beta(M);
  beta << 0.3, 2;
  double alpha = 0.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      binomial_logit_glm_lpmf_functor, n_scal, trials, x, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      binomial_logit_glm_lpmf_functor_propto, n_scal, trials, x, alpha, beta);
}

TEST(ProbDistributionsBinomialLogitGLM, opencl_matches_cpu_zero_instances) {
  int N = 0;
  int M = 2;

  std::vector<int> n{};
  std::vector<int> trials{};
  Eigen::MatrixXd x(N, M);
  Eigen::VectorXd beta(M);
  beta << 0.3, 2;
  double alpha = 0.3;

  stan::math::test::compare_cpu_opencl_prim_rev(binomial_logit_glm_lpmf_functor,
                                                n, trials, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      binomial_logit_glm_lpmf_functor_propto, n, trials, x, alpha, beta);
}

TEST(ProbDistributionsBinomialLogitGLM, opencl_matches_cpu_zero_attributes) {
  int N = 3;
  int M = 0;

  std::vector<int> n{0, 1, 0};
  std::vector<int> trials{10, 5, 4};
  Eigen::MatrixXd x(N, M);
  Eigen::VectorXd beta(M);
  double alpha = 0.3;

  stan::math::test::compare_cpu_opencl_prim_rev(binomial_logit_glm_lpmf_functor,
                                                n, trials, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      binomial_logit_glm_lpmf_functor_propto, n, trials, x, alpha, beta);
}

TEST(ProbDistributionsBinomialLogitGLM, opencl_matches_cpu_small_vector_alpha) {
  int N = 3;
  int M = 2;

  std::vector<int> n{0, 1, 0};
  std::vector<int> trials{0, 1, 0};
  Eigen::MatrixXd x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Eigen::VectorXd beta(M);
  beta << 0.3, 2;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, -0.8, 1.8;

  stan::math::test::compare_cpu_opencl_prim_rev(binomial_logit_glm_lpmf_functor,
                                                n, trials, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      binomial_logit_glm_lpmf_functor_propto, n, trials, x, alpha, beta);
}

TEST(ProbDistributionsBinomialLogitGLM, opencl_matches_cpu_big) {
  int N = 153;
  int M = 71;

  std::vector<int> n(N);
  std::vector<int> trials(N);
  for (int i = 0; i < N; i++) {
    n[i] = Eigen::ArrayXi::Random(1, 1).abs()(0);
    trials[i] = n[i] + Eigen::ArrayXi::Random(1, 1).abs()(0);
  }
  Eigen::MatrixXd x = Eigen::MatrixXd::Random(N, M);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(M);
  Eigen::VectorXd alpha = Eigen::VectorXd::Random(N);

  stan::math::test::compare_cpu_opencl_prim_rev(binomial_logit_glm_lpmf_functor,
                                                n, trials, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      binomial_logit_glm_lpmf_functor_propto, n, trials, x, alpha, beta);
}

#endif
