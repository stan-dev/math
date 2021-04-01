#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <stan/math/opencl/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

using Eigen::Array;
using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::matrix_cl;
using stan::math::var;
using std::vector;

TEST(ProbDistributionsCategoricalLogitGLM, error_checking) {
  int N = 3;
  int M = 2;
  int C = 3;

  vector<int> y{1, 3, 2};
  vector<int> y_size{1, 3, 2, 3};
  vector<int> y_value{1, 2, -23};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, Dynamic> x_size1(N - 1, M);
  x_size1 << -12, 46, -42, 24;
  Matrix<double, Dynamic, Dynamic> x_size2(N, M - 1);
  x_size2 << -12, 46, -42;
  Matrix<double, Dynamic, Dynamic> x_value(N, M);
  x_value << -12, 46, -42, 24, 25, -INFINITY;
  Matrix<double, Dynamic, Dynamic> beta(M, C);
  beta << 0.3, 2, 0.4, -0.1, -1.3, 1;
  Matrix<double, Dynamic, Dynamic> beta_size1(M + 1, C);
  beta_size1 << 0.3, 2, 0.4, -0.1, -1.3, 1, -0.1, -1.3, 1;
  Matrix<double, Dynamic, Dynamic> beta_size2(M, C + 1);
  beta_size2 << 0.3, 2, 0.4, -0.1, -1.3, 1, -1.3, 1;
  Matrix<double, Dynamic, Dynamic> beta_value(M, C);
  beta_value << 0.3, 2, 0.4, -0.1, -1.3, NAN;
  Matrix<double, Dynamic, 1> alpha(C, 1);
  alpha << 0.3, -0.8, 1.8;
  Matrix<double, Dynamic, 1> alpha_size(C - 1, 1);
  alpha_size << 0.3, -0.8;
  Matrix<double, Dynamic, 1> alpha_value(C, 1);
  alpha_value << 0.3, -0.8, NAN;

  matrix_cl<double> x_cl(x);
  matrix_cl<double> x_size1_cl(x_size1);
  matrix_cl<double> x_size2_cl(x_size2);
  matrix_cl<double> x_value_cl(x_value);
  matrix_cl<int> y_cl(y);
  matrix_cl<int> y_size_cl(y_size);
  matrix_cl<int> y_value_cl(y_value);
  matrix_cl<double> beta_cl(beta);
  matrix_cl<double> beta_size1_cl(beta_size1);
  matrix_cl<double> beta_size2_cl(beta_size2);
  matrix_cl<double> beta_value_cl(beta_value);
  matrix_cl<double> alpha_cl(alpha);
  matrix_cl<double> alpha_size_cl(alpha_size);
  matrix_cl<double> alpha_value_cl(alpha_value);

  EXPECT_NO_THROW(
      stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha_cl, beta_cl));

  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y_size_cl, x_cl, alpha_cl,
                                                      beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y_cl, x_size1_cl,
                                                      alpha_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y_cl, x_size2_cl,
                                                      alpha_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha_size_cl,
                                                      beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha_cl,
                                                      beta_size1_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha_cl,
                                                      beta_size2_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y_value_cl, x_cl,
                                                      alpha_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y_cl, x_value_cl,
                                                      alpha_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y_cl, x_cl,
                                                      alpha_value_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha_cl,
                                                      beta_value_cl),
               std::domain_error);
}

auto categorical_logit_glm_lpmf_functor
    = [](const auto& y, const auto& x, const auto& alpha, const auto& beta) {
        return stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta);
      };
auto categorical_logit_glm_lpmf_functor_propto
    = [](const auto& y, const auto& x, const auto& alpha, const auto& beta) {
        return stan::math::categorical_logit_glm_lpmf<true>(y, x, alpha, beta);
      };

TEST(ProbDistributionsCategoricalLogitGLM, opencl_matches_cpu_small_simple) {
  int N = 3;
  int M = 2;
  int C = 3;

  vector<int> y{1, 3, 2};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, Dynamic> beta(M, C);
  beta << 0.3, 2, 0.4, -0.1, -1.3, 1;
  Matrix<double, Dynamic, 1> alpha(C);
  alpha << 0.3, -2, 0.8;

  stan::math::test::compare_cpu_opencl_prim_rev(
      categorical_logit_glm_lpmf_functor, y, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      categorical_logit_glm_lpmf_functor_propto, y, x, alpha, beta);
}

TEST(ProbDistributionsCategoricalLogitGLM, opencl_broadcast_y) {
  int N = 3;
  int M = 2;
  int C = 3;

  int y_scal = 1;
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, Dynamic> beta(M, C);
  beta << 0.3, 2, 0.4, -0.1, -1.3, 1;
  Matrix<double, Dynamic, 1> alpha(C);
  alpha << 0.3, -2, 0.8;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      categorical_logit_glm_lpmf_functor, y_scal, x, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      categorical_logit_glm_lpmf_functor_propto, y_scal, x, alpha, beta);
}

TEST(ProbDistributionsCategoricalLogitGLM, opencl_matches_cpu_zero_instances) {
  int N = 0;
  int M = 2;
  int C = 3;

  vector<int> y{};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  Matrix<double, Dynamic, Dynamic> beta(M, C);
  beta << 0.3, 2, 0.4, -0.1, -1.3, 1;
  Matrix<double, Dynamic, 1> alpha(C);
  alpha << 0.3, -2, 0.8;

  stan::math::test::compare_cpu_opencl_prim_rev(
      categorical_logit_glm_lpmf_functor, y, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      categorical_logit_glm_lpmf_functor_propto, y, x, alpha, beta);
}

TEST(ProbDistributionsCategoricalLogitGLM, opencl_matches_cpu_zero_attributes) {
  int N = 3;
  int M = 0;
  int C = 3;

  vector<int> y{1, 3, 2};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  Matrix<double, Dynamic, Dynamic> beta(M, C);
  Matrix<double, Dynamic, 1> alpha(C);
  alpha << 0.3, -2, 0.8;

  stan::math::test::compare_cpu_opencl_prim_rev(
      categorical_logit_glm_lpmf_functor, y, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      categorical_logit_glm_lpmf_functor_propto, y, x, alpha, beta);
}

TEST(ProbDistributionsCategoricalLogitGLM, opencl_matches_cpu_single_class) {
  int N = 3;
  int M = 2;
  int C = 1;

  vector<int> y{1, 1, 1};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, Dynamic> beta(M, C);
  beta << 100000.3, 200000.5;
  Matrix<double, Dynamic, 1> alpha(C);
  alpha << 100000.3;

  stan::math::test::compare_cpu_opencl_prim_rev(
      categorical_logit_glm_lpmf_functor, y, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      categorical_logit_glm_lpmf_functor_propto, y, x, alpha, beta);
}

TEST(ProbDistributionsCategoricalLogitGLM, opencl_matches_cpu_all_vars) {
  int N = 5;
  int M = 3;
  int C = 2;

  vector<int> y(N);
  for (int i = 0; i < N; i++) {
    y[i] = Array<int, Dynamic, 1>::Random(1, 1).abs()(0) % C + 1;
  }
  Matrix<double, Dynamic, Dynamic> x
      = Matrix<double, Dynamic, Dynamic>::Random(N, M);
  Matrix<double, Dynamic, Dynamic> beta
      = Matrix<double, Dynamic, Dynamic>::Random(M, C);
  Matrix<double, Dynamic, 1> alpha = Matrix<double, Dynamic, 1>::Random(C, 1);

  stan::math::test::compare_cpu_opencl_prim_rev(
      categorical_logit_glm_lpmf_functor, y, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      categorical_logit_glm_lpmf_functor_propto, y, x, alpha, beta);
}

TEST(ProbDistributionsCategoricalLogitGLM, opencl_matches_cpu_big) {
  int N = 153;
  int M = 71;
  int C = 43;

  vector<int> y(N);
  for (int i = 0; i < N; i++) {
    y[i] = Array<int, Dynamic, 1>::Random(1, 1).abs()(0) % C + 1;
  }
  Matrix<double, Dynamic, Dynamic> x
      = Matrix<double, Dynamic, Dynamic>::Random(N, M);
  Matrix<double, Dynamic, Dynamic> beta
      = Matrix<double, Dynamic, Dynamic>::Random(M, C);
  Matrix<double, Dynamic, 1> alpha = Matrix<double, Dynamic, 1>::Random(C, 1);

  stan::math::test::compare_cpu_opencl_prim_rev(
      categorical_logit_glm_lpmf_functor, y, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      categorical_logit_glm_lpmf_functor_propto, y, x, alpha, beta);
}

#endif
