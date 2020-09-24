#ifdef STAN_OPENCL
#include <stan/math/opencl/rev/opencl.hpp>
#include <stan/math.hpp>
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

TEST(ProbDistributionsPoissonLogGLM, error_checking) {
  int N = 3;
  int M = 2;

  vector<int> y{0, 1, 5};
  vector<int> y_size{0, 1, 5, 0};
  vector<int> y_value{1, 4, -23};
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
      stan::math::poisson_log_glm_lpmf(y_cl, x_cl, alpha_cl, beta_cl));

  EXPECT_THROW(
      stan::math::poisson_log_glm_lpmf(y_size_cl, x_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::poisson_log_glm_lpmf(y_cl, x_size1_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::poisson_log_glm_lpmf(y_cl, x_size2_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::poisson_log_glm_lpmf(y_cl, x_cl, alpha_size_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::poisson_log_glm_lpmf(y_cl, x_cl, alpha_cl, beta_size_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::poisson_log_glm_lpmf(y_cl, x_cl, alpha_cl, beta_size_cl),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::poisson_log_glm_lpmf(y_value_cl, x_cl, alpha_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::poisson_log_glm_lpmf(y_cl, x_value_cl, alpha_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::poisson_log_glm_lpmf(y_cl, x_cl, alpha_value_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::poisson_log_glm_lpmf(y_cl, x_cl, alpha_cl, beta_value_cl),
      std::domain_error);
}

auto poisson_log_glm_lpmf_functor
    = [](const auto& y, const auto& x, const auto& alpha, const auto& beta) {
        return stan::math::poisson_log_glm_lpmf(y, x, alpha, beta);
      };
auto poisson_log_glm_lpmf_functor_propto
    = [](const auto& y, const auto& x, const auto& alpha, const auto& beta) {
        return stan::math::poisson_log_glm_lpmf<true>(y, x, alpha, beta);
      };

TEST(ProbDistributionsPoissonLogGLM, opencl_matches_cpu_small_simple) {
  int N = 3;
  int M = 2;

  vector<int> y{0, 1, 5};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  double alpha = 0.3;

  stan::math::test::compare_cpu_opencl_prim_rev(poisson_log_glm_lpmf_functor, y,
                                                x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      poisson_log_glm_lpmf_functor_propto, y, x, alpha, beta);
}

TEST(ProbDistributionsPoissonLogGLM, opencl_broadcast_y) {
  int N = 3;
  int M = 2;

  int y = 4;
  vector<int> y_vec{y, y, y};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  double alpha = 0.3;

  matrix_cl<double> x_cl(x);
  matrix_cl<int> y_vec_cl(y_vec);
  matrix_cl<double> beta_cl(beta);

  expect_near_rel(
      "poisson_log_glm_lpmf (OpenCL)",
      stan::math::poisson_log_glm_lpmf(y, x_cl, alpha, beta_cl),
      stan::math::poisson_log_glm_lpmf(y_vec_cl, x_cl, alpha, beta_cl));
  expect_near_rel(
      "poisson_log_glm_lpmf (OpenCL)",
      stan::math::poisson_log_glm_lpmf<true>(y, x_cl, alpha, beta_cl),
      stan::math::poisson_log_glm_lpmf<true>(y_vec_cl, x_cl, alpha, beta_cl));

  Matrix<var, Dynamic, Dynamic> x_var1 = x;
  Matrix<var, Dynamic, Dynamic> x_var2 = x;
  Matrix<var, Dynamic, 1> beta_var1 = beta;
  Matrix<var, Dynamic, 1> beta_var2 = beta;
  auto x_var1_cl = to_matrix_cl(x_var1);
  auto x_var2_cl = to_matrix_cl(x_var2);
  auto beta_var1_cl = stan::math::to_matrix_cl(beta_var1);
  auto beta_var2_cl = stan::math::to_matrix_cl(beta_var2);
  var alpha_var1 = alpha;
  var alpha_var2 = alpha;

  var res1 = stan::math::poisson_log_glm_lpmf(y, x_var1_cl, alpha_var1,
                                              beta_var1_cl);
  var res2 = stan::math::poisson_log_glm_lpmf(y_vec_cl, x_var2_cl, alpha_var2,
                                              beta_var2_cl);

  (res1 + res2).grad();

  expect_near_rel("poisson_log_glm_lpmf (OpenCL)", res1.val(), res2.val());

  expect_near_rel("bernoulli_logit_glm_lpmf (OpenCL)", x_var1.adj().eval(),
                  x_var2.adj().eval());
  expect_near_rel("poisson_log_glm_lpmf (OpenCL)", alpha_var1.adj(),
                  alpha_var2.adj());
  expect_near_rel("poisson_log_glm_lpmf (OpenCL)", beta_var1.adj().eval(),
                  beta_var2.adj().eval());
}

TEST(ProbDistributionsPoissonLogGLM, opencl_matches_cpu_zero_instances) {
  int N = 0;
  int M = 2;

  vector<int> y{};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  double alpha = 0.3;

  stan::math::test::compare_cpu_opencl_prim_rev(poisson_log_glm_lpmf_functor, y,
                                                x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      poisson_log_glm_lpmf_functor_propto, y, x, alpha, beta);
}

TEST(ProbDistributionsPoissonLogGLM, opencl_matches_cpu_zero_attributes) {
  int N = 3;
  int M = 0;

  vector<int> y{0, 1, 5};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  Matrix<double, Dynamic, 1> beta(M, 1);
  double alpha = 0.3;

  stan::math::test::compare_cpu_opencl_prim_rev(poisson_log_glm_lpmf_functor, y,
                                                x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      poisson_log_glm_lpmf_functor_propto, y, x, alpha, beta);
}

TEST(ProbDistributionsPoissonLogGLM, opencl_matches_cpu_small_vector_alpha) {
  int N = 3;
  int M = 2;

  vector<int> y{0, 1, 5};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  Matrix<double, Dynamic, 1> alpha(N, 1);
  alpha << 0.3, -0.8, 1.8;

  stan::math::test::compare_cpu_opencl_prim_rev(poisson_log_glm_lpmf_functor, y,
                                                x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      poisson_log_glm_lpmf_functor_propto, y, x, alpha, beta);
}

TEST(ProbDistributionsPoissonLogGLM, opencl_matches_cpu_big) {
  int N = 153;
  int M = 71;

  vector<int> y(N);
  for (int i = 0; i < N; i++) {
    y[i] = Array<int, Dynamic, 1>::Random(1, 1).abs()(0);
  }
  Matrix<double, Dynamic, Dynamic> x
      = Matrix<double, Dynamic, Dynamic>::Random(N, M);
  Matrix<double, Dynamic, 1> beta = Matrix<double, Dynamic, 1>::Random(M, 1);
  Matrix<double, Dynamic, 1> alpha = Matrix<double, Dynamic, 1>::Random(N, 1);

  stan::math::test::compare_cpu_opencl_prim_rev(poisson_log_glm_lpmf_functor, y,
                                                x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      poisson_log_glm_lpmf_functor_propto, y, x, alpha, beta);
}

#endif
