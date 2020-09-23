#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <stan/math/opencl/opencl.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

using Eigen::Array;
using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::matrix_cl;
using stan::math::var;
using std::vector;

TEST(ProbDistributionsOrderedLogisitcGLM, error_checking) {
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
  matrix_cl<int> y_cl(y);
  matrix_cl<int> y_size_cl(y_size);
  matrix_cl<int> y_value1_cl(y_value1);
  matrix_cl<int> y_value2_cl(y_value2);
  matrix_cl<double> beta_cl(beta);
  matrix_cl<double> beta_size_cl(beta_size);
  matrix_cl<double> beta_value_cl(beta_value);
  matrix_cl<double> cuts_cl(cuts);
  matrix_cl<double> cuts_value1_cl(cuts_value1);
  matrix_cl<double> cuts_value2_cl(cuts_value2);

  EXPECT_NO_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta_cl, cuts_cl));

  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_size_cl, x_cl, beta_cl, cuts_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_cl, x_size1_cl, beta_cl, cuts_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_cl, x_size2_cl, beta_cl, cuts_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta_size_cl, cuts_cl),
      std::invalid_argument);

  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y_value1_cl, x_cl, beta_cl,
                                                     cuts_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y_value2_cl, x_cl, beta_cl,
                                                     cuts_cl),
               std::domain_error);
  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_cl, x_value_cl, beta_cl, cuts_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta_value_cl, cuts_cl),
      std::domain_error);
  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta_cl,
                                                     cuts_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y_cl, x_cl, beta_cl,
                                                     cuts_value2_cl),
               std::domain_error);
}

auto ordered_logistic_glm_lpmf_functor
    = [](const auto& y, const auto& x, const auto& beta, const auto& cuts) {
        return stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts);
      };
auto ordered_logistic_glm_lpmf_functor_propto
    = [](const auto& y, const auto& x, const auto& beta, const auto& cuts) {
        return stan::math::ordered_logistic_glm_lpmf<true>(y, x, beta, cuts);
      };

TEST(ProbDistributionsOrderedLogisitcGLM, opencl_matches_cpu_small_simple) {
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

  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_glm_lpmf_functor, y, x, beta, cuts);
  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_glm_lpmf_functor_propto, y, x, beta, cuts);
}

TEST(ProbDistributionsOrderedLogisitcGLM, opencl_matches_cpu_broadcast_y) {
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
  matrix_cl<double> beta_cl(beta);
  matrix_cl<double> cuts_cl(cuts);

  EXPECT_NEAR_REL(
      stan::math::ordered_logistic_glm_lpmf(y, x_cl, beta_cl, cuts_cl),
      stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts));
  EXPECT_NEAR_REL(
      stan::math::ordered_logistic_glm_lpmf<true>(y, x_cl, beta_cl, cuts_cl),
      stan::math::ordered_logistic_glm_lpmf<true>(y, x, beta, cuts));

  Matrix<var, Dynamic, Dynamic> x_var1 = x;
  Matrix<var, Dynamic, Dynamic> x_var2 = x;
  Matrix<var, Dynamic, 1> beta_var1 = beta;
  Matrix<var, Dynamic, 1> beta_var2 = beta;
  Matrix<var, Dynamic, 1> cuts_var1 = cuts;
  Matrix<var, Dynamic, 1> cuts_var2 = cuts;
  auto x_var1_cl = stan::math::to_matrix_cl(x_var1);
  auto beta_var1_cl = stan::math::to_matrix_cl(beta_var1);
  auto cuts_var1_cl = stan::math::to_matrix_cl(cuts_var1);

  var res1 = stan::math::ordered_logistic_glm_lpmf(y, x_var1_cl, beta_var1_cl,
                                                   cuts_var1_cl);
  var res2
      = stan::math::ordered_logistic_glm_lpmf(y, x_var2, beta_var2, cuts_var2);

  (res1 + res2).grad();

  EXPECT_NEAR_REL(res1.val(), res2.val());

  EXPECT_NEAR_REL(x_var1.adj().eval(), x_var2.adj().eval());
  EXPECT_NEAR_REL(beta_var1.adj().eval(), beta_var2.adj().eval());
  EXPECT_NEAR_REL(cuts_var1.adj().eval(), cuts_var2.adj().eval());
}

TEST(ProbDistributionsOrderedLogisitcGLM, opencl_matches_cpu_zero_instances) {
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
  matrix_cl<int> y_cl(y);
  matrix_cl<double> beta_cl(beta);
  matrix_cl<double> cuts_cl(cuts);

  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_glm_lpmf_functor, y, x, beta, cuts);
  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_glm_lpmf_functor_propto, y, x, beta, cuts);
}

TEST(ProbDistributionsOrderedLogisitcGLM, opencl_matches_cpu_zero_attributes) {
  int N = 3;
  int M = 0;
  int C = 5;

  vector<int> y{2, 1, 5};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  Matrix<double, Dynamic, 1> beta(M, 1);
  Matrix<double, Dynamic, 1> cuts(C - 1, 1);
  cuts << -0.4, 0.1, 0.3, 4.5;

  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_glm_lpmf_functor, y, x, beta, cuts);
  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_glm_lpmf_functor_propto, y, x, beta, cuts);
}

TEST(ProbDistributionsOrderedLogisitcGLM, opencl_matches_cpu_single_class) {
  int N = 3;
  int M = 2;
  int C = 1;

  vector<int> y{1, 1, 1};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << -12, 46, -42, 24, 25, 27;
  Matrix<double, Dynamic, 1> beta(M, 1);
  beta << 0.3, 2;
  Matrix<double, Dynamic, 1> cuts(C - 1, 1);

  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_glm_lpmf_functor, y, x, beta, cuts);
  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_glm_lpmf_functor_propto, y, x, beta, cuts);
}

TEST(ProbDistributionsOrderedLogisitcGLM, opencl_matches_cpu_big) {
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

  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_glm_lpmf_functor, y, x, beta, cuts);
  stan::math::test::compare_cpu_opencl_prim_rev(
      ordered_logistic_glm_lpmf_functor_propto, y, x, beta, cuts);
}

#endif
