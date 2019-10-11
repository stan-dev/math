#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <stan/math/opencl/opencl.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/expect_near_rel.hpp>
#include <vector>

using Eigen::Array;
using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::matrix_cl;
using stan::math::var;
using stan::test::expect_near_rel;
using std::vector;

TEST(ProbDistributionsCategoricalLogitGLM, error_checking) {
  double eps = 1e-9;
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
  matrix_cl<int> y_cl(y, N, 1);
  matrix_cl<int> y_size_cl(y_size, N + 1, 1);
  matrix_cl<int> y_value_cl(y_value, N, 1);

  EXPECT_NO_THROW(
      stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha, beta));

  EXPECT_THROW(
      stan::math::categorical_logit_glm_lpmf(y_size_cl, x_cl, alpha, beta),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::categorical_logit_glm_lpmf(y_cl, x_size1_cl, alpha, beta),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::categorical_logit_glm_lpmf(y_cl, x_size2_cl, alpha, beta),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha_size, beta),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha, beta_size1),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha, beta_size2),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::categorical_logit_glm_lpmf(y_value_cl, x_cl, alpha, beta),
      std::domain_error);
  EXPECT_THROW(
      stan::math::categorical_logit_glm_lpmf(y_cl, x_value_cl, alpha, beta),
      std::domain_error);
  EXPECT_THROW(
      stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha_value, beta),
      std::domain_error);
  EXPECT_THROW(
      stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha, beta_value),
      std::domain_error);
}

TEST(ProbDistributionsCategoricalLogitGLM, gpu_matches_cpu_small_simple) {
  double eps = 1e-9;
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

  matrix_cl<double> x_cl(x);
  matrix_cl<int> y_cl(y, N, 1);

  expect_near_rel(
      "categorical_logit_glm_lpmf (OpenCL)",
      stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha, beta),
      stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta));
  expect_near_rel(
      "categorical_logit_glm_lpmf (OpenCL)",
      stan::math::categorical_logit_glm_lpmf<true>(y_cl, x_cl, alpha, beta),
      stan::math::categorical_logit_glm_lpmf<true>(y, x, alpha, beta));

  Matrix<var, Dynamic, Dynamic> beta_var1 = beta;
  Matrix<var, Dynamic, Dynamic> beta_var2 = beta;
  Matrix<var, Dynamic, 1> alpha_var1 = alpha;
  Matrix<var, Dynamic, 1> alpha_var2 = alpha;

  var res1 = stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha_var1,
                                                    beta_var1);
  var res2
      = stan::math::categorical_logit_glm_lpmf(y, x, alpha_var2, beta_var2);

  (res1 + res2).grad();

  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)", res1.val(),
                  res2.val());

  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)",
                  alpha_var1.adj().eval(), alpha_var2.adj().eval());
  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)", beta_var1.adj().eval(),
                  beta_var2.adj().eval());
}

TEST(ProbDistributionsCategoricalLogitGLM, gpu_matches_cpu_zero_instances) {
  double eps = 1e-9;
  int N = 0;
  int M = 2;
  int C = 3;

  vector<int> y{};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  Matrix<double, Dynamic, Dynamic> beta(M, C);
  beta << 0.3, 2, 0.4, -0.1, -1.3, 1;
  Matrix<double, Dynamic, 1> alpha(C);
  alpha << 0.3, -2, 0.8;

  matrix_cl<double> x_cl(x);
  matrix_cl<int> y_cl(y, N, 1);

  expect_near_rel(
      "categorical_logit_glm_lpmf (OpenCL)",
      stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha, beta),
      stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta));
  expect_near_rel(
      "categorical_logit_glm_lpmf (OpenCL)",
      stan::math::categorical_logit_glm_lpmf<true>(y_cl, x_cl, alpha, beta),
      stan::math::categorical_logit_glm_lpmf<true>(y, x, alpha, beta));

  Matrix<var, Dynamic, Dynamic> beta_var1 = beta;
  Matrix<var, Dynamic, Dynamic> beta_var2 = beta;
  Matrix<var, Dynamic, 1> alpha_var1 = alpha;
  Matrix<var, Dynamic, 1> alpha_var2 = alpha;

  var res1 = stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha_var1,
                                                    beta_var1);
  var res2
      = stan::math::categorical_logit_glm_lpmf(y, x, alpha_var2, beta_var2);

  (res1 + res2).grad();

  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)", res1.val(),
                  res2.val());

  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)",
                  alpha_var1.adj().eval(), alpha_var2.adj().eval());
  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)", beta_var1.adj().eval(),
                  beta_var2.adj().eval());
}

TEST(ProbDistributionsCategoricalLogitGLM, gpu_matches_cpu_zero_attributes) {
  double eps = 1e-9;
  int N = 3;
  int M = 0;
  int C = 3;

  vector<int> y{1, 3, 2};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  Matrix<double, Dynamic, Dynamic> beta(M, C);
  Matrix<double, Dynamic, 1> alpha(C);
  alpha << 0.3, -2, 0.8;

  matrix_cl<double> x_cl(x);
  matrix_cl<int> y_cl(y, N, 1);

  expect_near_rel(
      "categorical_logit_glm_lpmf (OpenCL)",
      stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha, beta),
      stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta));
  expect_near_rel(
      "categorical_logit_glm_lpmf (OpenCL)",
      stan::math::categorical_logit_glm_lpmf<true>(y_cl, x_cl, alpha, beta),
      stan::math::categorical_logit_glm_lpmf<true>(y, x, alpha, beta));

  Matrix<var, Dynamic, Dynamic> beta_var1 = beta;
  Matrix<var, Dynamic, Dynamic> beta_var2 = beta;
  Matrix<var, Dynamic, 1> alpha_var1 = alpha;
  Matrix<var, Dynamic, 1> alpha_var2 = alpha;

  var res1 = stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha_var1,
                                                    beta_var1);
  var res2
      = stan::math::categorical_logit_glm_lpmf(y, x, alpha_var2, beta_var2);

  (res1 + res2).grad();

  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)", res1.val(),
                  res2.val());

  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)",
                  alpha_var1.adj().eval(), alpha_var2.adj().eval());
  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)", beta_var1.adj().eval(),
                  beta_var2.adj().eval());
}

TEST(ProbDistributionsCategoricalLogitGLM, gpu_matches_cpu_single_class) {
  double eps = 1e-9;
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

  matrix_cl<double> x_cl(x);
  matrix_cl<int> y_cl(y, N, 1);

  expect_near_rel(
      "categorical_logit_glm_lpmf (OpenCL)",
      stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha, beta),
      stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta));
  expect_near_rel(
      "categorical_logit_glm_lpmf (OpenCL)",
      stan::math::categorical_logit_glm_lpmf<true>(y_cl, x_cl, alpha, beta),
      stan::math::categorical_logit_glm_lpmf<true>(y, x, alpha, beta));

  Matrix<var, Dynamic, Dynamic> beta_var1 = beta;
  Matrix<var, Dynamic, Dynamic> beta_var2 = beta;
  Matrix<var, Dynamic, 1> alpha_var1 = alpha;
  Matrix<var, Dynamic, 1> alpha_var2 = alpha;

  var res1 = stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha_var1,
                                                    beta_var1);
  var res2
      = stan::math::categorical_logit_glm_lpmf(y, x, alpha_var2, beta_var2);

  (res1 + res2).grad();

  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)", res1.val(),
                  res2.val());

  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)",
                  alpha_var1.adj().eval(), alpha_var2.adj().eval());
  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)", beta_var1.adj().eval(),
                  beta_var2.adj().eval());
}

TEST(ProbDistributionsCategoricalLogitGLM, gpu_matches_cpu_big) {
  double eps = 1e-9;
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

  matrix_cl<double> x_cl(x);
  matrix_cl<int> y_cl(y, N, 1);

  expect_near_rel(
      "categorical_logit_glm_lpmf (OpenCL)",
      stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha, beta),
      stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta));
  expect_near_rel(
      "categorical_logit_glm_lpmf (OpenCL)",
      stan::math::categorical_logit_glm_lpmf(y_cl, x_cl, alpha, beta),
      stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta));

  Matrix<var, Dynamic, Dynamic> beta_var1 = beta;
  Matrix<var, Dynamic, Dynamic> beta_var2 = beta;
  Matrix<var, Dynamic, 1> alpha_var1 = alpha;
  Matrix<var, Dynamic, 1> alpha_var2 = alpha;

  var res1 = stan::math::categorical_logit_glm_lpmf<true>(
      y_cl, x_cl, alpha_var1, beta_var1);
  var res2 = stan::math::categorical_logit_glm_lpmf<true>(y, x, alpha_var2,
                                                          beta_var2);

  (res1 + res2).grad();

  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)", res1.val(),
                  res2.val());

  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)", beta_var1.adj().eval(),
                  beta_var2.adj().eval());
  expect_near_rel("categorical_logit_glm_lpmf (OpenCL)",
                  alpha_var1.adj().eval(), alpha_var2.adj().eval());
}

#endif
