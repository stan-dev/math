#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <stan/math/opencl/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

using Eigen::Array;
using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::matrix_cl;
using stan::math::var;
using stan::test::expect_near_rel;
using std::vector;

TEST(ProbDistributionsBernoulliLogitGLM, error_checking) {
  int N = 3;
  int M = 2;

  vector<int> y{1, 0, 1};
  vector<int> y_size{1, 0, 1, 0};
  vector<int> y_value{0, 1, 23};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, Dynamic> x_size1(N - 1, M);
  x_size1 << -12, 46, -42, 24;
  Matrix<double, Dynamic, Dynamic> x_size2(N, M - 1);
  x_size2 << -12, 46, -42;
  Matrix<double, Dynamic, Dynamic> x_value(N, M);
  x_value << -12, 46, -42, 24, 25, -INFINITY;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  Matrix<double, Dynamic, 1> beta_size(M + 1, 1);
  beta_size << 0.3, 2, 0.4;
  Matrix<double, Dynamic, 1> beta_value(M, 1);
  beta_value << 0.3, INFINITY;
  Matrix<double, Dynamic, 1> alpha(N, 1);
  alpha << 0.3, -0.8, 1.8;
  Matrix<double, Dynamic, 1> alpha_size(N - 1, 1);
  alpha_size << 0.3, -0.8;
  Matrix<double, Dynamic, 1> alpha_value(N, 1);
  alpha_value << 0.3, -0.8, NAN;

  matrix_cl<double> x_cl(x);
  matrix_cl<double> x_size1_cl(x_size1);
  matrix_cl<double> x_size2_cl(x_size2);
  matrix_cl<double> x_value_cl(x_value);
  matrix_cl<int> y_cl(y);
  matrix_cl<int> y_size_cl(y_size);
  matrix_cl<int> y_value_cl(y_value);
  matrix_cl<double> beta_cl(beta);
  matrix_cl<double> beta_size_cl(beta_size);
  matrix_cl<double> beta_value_cl(beta_value);
  matrix_cl<double> alpha_cl(alpha);
  matrix_cl<double> alpha_size_cl(alpha_size);
  matrix_cl<double> alpha_value_cl(alpha_value);

  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_lpmf(y_cl, x_cl, alpha_cl, beta_cl));

  EXPECT_THROW(
      stan::math::bernoulli_logit_glm_lpmf(y_size_cl, x_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::bernoulli_logit_glm_lpmf(y_cl, x_size1_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::bernoulli_logit_glm_lpmf(y_cl, x_size2_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::bernoulli_logit_glm_lpmf(y_cl, x_cl, alpha_size_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::bernoulli_logit_glm_lpmf(y_cl, x_cl, alpha_cl, beta_size_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::bernoulli_logit_glm_lpmf(y_cl, x_cl, alpha_cl, beta_size_cl),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::bernoulli_logit_glm_lpmf(y_value_cl, x_cl, alpha_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::bernoulli_logit_glm_lpmf(y_cl, x_value_cl, alpha_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::bernoulli_logit_glm_lpmf(y_cl, x_cl, alpha_value_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::bernoulli_logit_glm_lpmf(y_cl, x_cl, alpha_cl, beta_value_cl),
      std::domain_error);
}

auto bernoulli_logit_glm_lpmf_functor
    = [](const auto& y, const auto& x, const auto& alpha, const auto& beta) {
        return stan::math::bernoulli_logit_glm_lpmf(y, x, alpha, beta);
      };
auto bernoulli_logit_glm_lpmf_functor_propto
    = [](const auto& y, const auto& x, const auto& alpha, const auto& beta) {
        return stan::math::bernoulli_logit_glm_lpmf<true>(y, x, alpha, beta);
      };

TEST(ProbDistributionsBernoulliLogitGLM, opencl_matches_cpu_small_simple) {
  int N = 3;
  int M = 2;

  vector<int> y{0, 1, 0};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  double alpha = 0.3;

  stan::math::test::compare_cpu_opencl_prim_rev(
      bernoulli_logit_glm_lpmf_functor, y, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      bernoulli_logit_glm_lpmf_functor_propto, y, x, alpha, beta);
}

TEST(ProbDistributionsBernoulliLogitGLM, opencl_broadcast_y) {
  int N = 3;
  int M = 2;

  int y_scal = 1;
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  double alpha = 0.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      bernoulli_logit_glm_lpmf_functor, y_scal, x, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      bernoulli_logit_glm_lpmf_functor_propto, y_scal, x, alpha, beta);
}

TEST(ProbDistributionsBernoulliLogitGLM, opencl_matches_cpu_zero_instances) {
  int N = 0;
  int M = 2;

  vector<int> y{};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  double alpha = 0.3;

  stan::math::test::compare_cpu_opencl_prim_rev(
      bernoulli_logit_glm_lpmf_functor, y, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      bernoulli_logit_glm_lpmf_functor_propto, y, x, alpha, beta);
}

TEST(ProbDistributionsBernoulliLogitGLM, opencl_matches_cpu_zero_attributes) {
  int N = 3;
  int M = 0;

  vector<int> y{0, 1, 0};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  Matrix<double, Dynamic, 1> beta(M, 1);
  double alpha = 0.3;

  stan::math::test::compare_cpu_opencl_prim_rev(
      bernoulli_logit_glm_lpmf_functor, y, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      bernoulli_logit_glm_lpmf_functor_propto, y, x, alpha, beta);
}

TEST(ProbDistributionsBernoulliLogitGLM,
     opencl_matches_cpu_small_vector_alpha) {
  int N = 3;
  int M = 2;

  vector<int> y{0, 1, 0};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  Matrix<double, Dynamic, 1> alpha(N, 1);
  alpha << 0.3, -0.8, 1.8;

  stan::math::test::compare_cpu_opencl_prim_rev(
      bernoulli_logit_glm_lpmf_functor, y, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      bernoulli_logit_glm_lpmf_functor_propto, y, x, alpha, beta);
}

TEST(ProbDistributionsBernoulliLogitGLM, opencl_matches_cpu_big) {
  int N = 153;
  int M = 71;

  vector<int> y(N);
  for (int i = 0; i < N; i++) {
    y[i] = Array<int, Dynamic, 1>::Random(1, 1).abs()(0) % 2;
  }
  Matrix<double, Dynamic, Dynamic> x
      = Matrix<double, Dynamic, Dynamic>::Random(N, M);
  Matrix<double, Dynamic, 1> beta = Matrix<double, Dynamic, 1>::Random(M, 1);
  Matrix<double, Dynamic, 1> alpha = Matrix<double, Dynamic, 1>::Random(N, 1);

  stan::math::test::compare_cpu_opencl_prim_rev(
      bernoulli_logit_glm_lpmf_functor, y, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      bernoulli_logit_glm_lpmf_functor_propto, y, x, alpha, beta);
}

#endif
