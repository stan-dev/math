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

TEST(ProbDistributionsOrderedLogisitcGLM, error_checking) {
  double eps = 1e-9;
  int N = 3;
  int M = 2;
  int C = 5;

  vector<int> y{1, 3, 2};
  vector<int> y_size{1, 3, 1, 2};
  vector<int> y_value1{0, 1, 2};
  vector<int> y_value2{1, 2, 6};
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
  Matrix<double, Dynamic, 1> cuts(C - 1, 1);
  cuts << -0.3, 0.8, 1.8, 3.0;
  Matrix<double, Dynamic, 1> cuts_value1(C - 1, 1);
  cuts_value1 << 0.3, -0.8, 3.0, NAN;
  Matrix<double, Dynamic, 1> cuts_value2(C - 1, 1);
  cuts_value2 << 0.3, -0.8, 0, 4.5;

  matrix_cl<double> x_cl(x);
  matrix_cl<double> x_size1_cl(x_size1);
  matrix_cl<double> x_size2_cl(x_size2);
  matrix_cl<double> x_value_cl(x_value);
  matrix_cl<int> y_cl(y, N, 1);
  matrix_cl<int> y_size_cl(y_size, N + 1, 1);
  matrix_cl<int> y_value1_cl(y_value1, N, 1);
  matrix_cl<int> y_value2_cl(y_value2, N, 1);

  EXPECT_NO_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta, cuts));

  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_size_cl, x_cl, beta, cuts),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_cl, x_size1_cl, beta, cuts),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_cl, x_size2_cl, beta, cuts),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta_size, cuts),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_value1_cl, x_cl, beta, cuts),
      std::domain_error);
  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_value2_cl, x_cl, beta, cuts),
      std::domain_error);
  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_cl, x_value_cl, beta, cuts),
      std::domain_error);
  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta_value, cuts),
      std::domain_error);
  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta, cuts_value1),
      std::domain_error);
  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta, cuts_value2),
      std::domain_error);
}

TEST(ProbDistributionsOrderedLogisitcGLM, gpu_matches_cpu_small_simple) {
  double eps = 1e-9;
  int N = 3;
  int M = 2;
  int C = 5;

  vector<int> y{2, 1, 5};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  Matrix<double, Dynamic, 1> cuts(C - 1, 1);
  cuts << -0.4, 0.1, 0.3, 4.5;

  matrix_cl<double> x_cl(x);
  matrix_cl<int> y_cl(y, N, 1);

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)",
                  stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta, cuts),
                  stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts));
  expect_near_rel(
      "ordered_logistic_glm_lpmf (OpenCL)",
      stan::math::ordered_logistic_glm_lpmf<true>(y_cl, x_cl, beta, cuts),
      stan::math::ordered_logistic_glm_lpmf<true>(y, x, beta, cuts));

  Matrix<var, Dynamic, 1> beta_var1 = beta;
  Matrix<var, Dynamic, 1> beta_var2 = beta;
  Matrix<var, Dynamic, 1> cuts_var1 = cuts;
  Matrix<var, Dynamic, 1> cuts_var2 = cuts;

  var res1
      = stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta_var1, cuts_var1);
  var res2 = stan::math::ordered_logistic_glm_lpmf(y, x, beta_var2, cuts_var2);

  (res1 + res2).grad();

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", res1.val(), res2.val());

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", beta_var1.adj().eval(),
                  beta_var2.adj().eval());
  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", cuts_var1.adj().eval(),
                  cuts_var2.adj().eval());
}

TEST(ProbDistributionsOrderedLogisitcGLM, gpu_matches_cpu_broadcast_y) {
  double eps = 1e-9;
  int N = 3;
  int M = 2;
  int C = 5;

  int y = 1;
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  Matrix<double, Dynamic, 1> cuts(C - 1, 1);
  cuts << -0.4, 0.1, 0.3, 4.5;

  matrix_cl<double> x_cl(x);
  matrix_cl<int> y_cl(y);

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)",
                  stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta, cuts),
                  stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts));
  expect_near_rel(
      "ordered_logistic_glm_lpmf (OpenCL)",
      stan::math::ordered_logistic_glm_lpmf<true>(y_cl, x_cl, beta, cuts),
      stan::math::ordered_logistic_glm_lpmf<true>(y, x, beta, cuts));

  Matrix<var, Dynamic, 1> beta_var1 = beta;
  Matrix<var, Dynamic, 1> beta_var2 = beta;
  Matrix<var, Dynamic, 1> cuts_var1 = cuts;
  Matrix<var, Dynamic, 1> cuts_var2 = cuts;

  var res1
      = stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta_var1, cuts_var1);
  var res2 = stan::math::ordered_logistic_glm_lpmf(y, x, beta_var2, cuts_var2);

  (res1 + res2).grad();

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", res1.val(), res2.val());

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", beta_var1.adj().eval(),
                  beta_var2.adj().eval());
  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", cuts_var1.adj().eval(),
                  cuts_var2.adj().eval());
}

TEST(ProbDistributionsOrderedLogisitcGLM, gpu_matches_cpu_zero_instances) {
  double eps = 1e-9;
  int N = 0;
  int M = 2;
  int C = 5;

  vector<int> y{};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  Matrix<double, Dynamic, 1> cuts(C - 1, 1);
  cuts << -0.4, 0.1, 0.3, 4.5;

  matrix_cl<double> x_cl(x);
  matrix_cl<int> y_cl(y, N, 1);

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)",
                  stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta, cuts),
                  stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts));
  expect_near_rel(
      "ordered_logistic_glm_lpmf (OpenCL)",
      stan::math::ordered_logistic_glm_lpmf<true>(y_cl, x_cl, beta, cuts),
      stan::math::ordered_logistic_glm_lpmf<true>(y, x, beta, cuts));

  Matrix<var, Dynamic, 1> beta_var1 = beta;
  Matrix<var, Dynamic, 1> beta_var2 = beta;
  Matrix<var, Dynamic, 1> cuts_var1 = cuts;
  Matrix<var, Dynamic, 1> cuts_var2 = cuts;

  var res1
      = stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta_var1, cuts_var1);
  var res2 = stan::math::ordered_logistic_glm_lpmf(y, x, beta_var2, cuts_var2);

  (res1 + res2).grad();

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", res1.val(), res2.val());

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", beta_var1.adj().eval(),
                  beta_var2.adj().eval());
  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", cuts_var1.adj().eval(),
                  cuts_var2.adj().eval());
}

TEST(ProbDistributionsOrderedLogisitcGLM, gpu_matches_cpu_zero_attributes) {
  double eps = 1e-9;
  int N = 3;
  int M = 0;
  int C = 5;

  vector<int> y{2, 1, 5};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  Matrix<double, Dynamic, 1> beta(M, 1);
  Matrix<double, Dynamic, 1> cuts(C - 1, 1);
  cuts << -0.4, 0.1, 0.3, 4.5;

  matrix_cl<double> x_cl(x);
  matrix_cl<int> y_cl(y, N, 1);

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)",
                  stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta, cuts),
                  stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts));
  expect_near_rel(
      "ordered_logistic_glm_lpmf (OpenCL)",
      stan::math::ordered_logistic_glm_lpmf<true>(y_cl, x_cl, beta, cuts),
      stan::math::ordered_logistic_glm_lpmf<true>(y, x, beta, cuts));

  Matrix<var, Dynamic, 1> beta_var1 = beta;
  Matrix<var, Dynamic, 1> beta_var2 = beta;
  Matrix<var, Dynamic, 1> cuts_var1 = cuts;
  Matrix<var, Dynamic, 1> cuts_var2 = cuts;

  var res1
      = stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta_var1, cuts_var1);
  var res2 = stan::math::ordered_logistic_glm_lpmf(y, x, beta_var2, cuts_var2);

  (res1 + res2).grad();

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", res1.val(), res2.val());

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", beta_var1.adj().eval(),
                  beta_var2.adj().eval());
  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", cuts_var1.adj().eval(),
                  cuts_var2.adj().eval());
}

TEST(ProbDistributionsOrderedLogisitcGLM, gpu_matches_cpu_single_class) {
  double eps = 1e-9;
  int N = 3;
  int M = 2;
  int C = 1;

  vector<int> y{1, 1, 1};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  Matrix<double, Dynamic, 1> cuts(C - 1, 1);

  matrix_cl<double> x_cl(x);
  matrix_cl<int> y_cl(y, N, 1);

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)",
                  stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta, cuts),
                  stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts));
  expect_near_rel(
      "ordered_logistic_glm_lpmf (OpenCL)",
      stan::math::ordered_logistic_glm_lpmf<true>(y_cl, x_cl, beta, cuts),
      stan::math::ordered_logistic_glm_lpmf<true>(y, x, beta, cuts));

  Matrix<var, Dynamic, 1> beta_var1 = beta;
  Matrix<var, Dynamic, 1> beta_var2 = beta;
  Matrix<var, Dynamic, 1> cuts_var1 = cuts;
  Matrix<var, Dynamic, 1> cuts_var2 = cuts;

  var res1
      = stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta_var1, cuts_var1);
  var res2 = stan::math::ordered_logistic_glm_lpmf(y, x, beta_var2, cuts_var2);

  (res1 + res2).grad();

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", res1.val(), res2.val());

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", beta_var1.adj().eval(),
                  beta_var2.adj().eval());
  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", cuts_var1.adj().eval(),
                  cuts_var2.adj().eval());
}

TEST(ProbDistributionsOrderedLogisitcGLM, gpu_matches_cpu_big) {
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
  Matrix<double, Dynamic, 1> beta = Matrix<double, Dynamic, 1>::Random(M, 1);
  Matrix<double, Dynamic, 1> cuts
      = Array<double, Dynamic, 1>::Random(C - 1, 1).abs() + 0.00001;
  for (int i = 1; i < C - 1; i++) {
    cuts[i] += cuts[i - 1];
  }

  matrix_cl<double> x_cl(x);
  matrix_cl<int> y_cl(y, N, 1);

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)",
                  stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta, cuts),
                  stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts));
  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)",
                  stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta, cuts),
                  stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts));

  Matrix<var, Dynamic, 1> beta_var1 = beta;
  Matrix<var, Dynamic, 1> beta_var2 = beta;
  Matrix<var, Dynamic, 1> cuts_var1 = cuts;
  Matrix<var, Dynamic, 1> cuts_var2 = cuts;

  var res1 = stan::math::ordered_logistic_glm_lpmf<true>(y_cl, x_cl, beta_var1,
                                                         cuts_var1);
  var res2
      = stan::math::ordered_logistic_glm_lpmf<true>(y, x, beta_var2, cuts_var2);

  (res1 + res2).grad();

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", res1.val(), res2.val());

  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", beta_var1.adj().eval(),
                  beta_var2.adj().eval());
  expect_near_rel("ordered_logistic_glm_lpmf (OpenCL)", cuts_var1.adj().eval(),
                  cuts_var2.adj().eval());
}

#endif
