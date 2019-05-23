#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <Eigen/Core>

using stan::math::var;
using stan::math::ordered_logistic_lpmf;
using stan::math::ordered_logistic_glm_lpmf;
using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::Array;
using Eigen::Dynamic;
using Eigen::VectorXd;
using Eigen::RowVectorXd;

template<bool propto, typename T_x, typename T_beta, typename T_cuts>
typename stan::return_type<T_x, T_beta, T_cuts>::type
ordered_logistic_glm_simple_lpmf(
        const Matrix<int, Dynamic, 1>& y, const Matrix<T_x, Dynamic, Dynamic>& x,
        const T_beta& beta, const T_cuts& cuts) {
  typedef typename stan::return_type<T_x, T_beta>::type T_x_beta;
  using stan::math::as_column_vector_or_scalar;

  auto& beta_col = as_column_vector_or_scalar(beta);

  Eigen::Matrix<T_x_beta, Dynamic, 1> location = x.template cast<T_x_beta>() * beta_col.template cast<T_x_beta>();

  return ordered_logistic_lpmf<propto>(y, location, cuts);
}

TEST(ProbDistributionsOrderedLogisticGLM, glm_matches_ordered_logistic_doubles) {
  double eps = 1e-13;
  int N = 5;
  int M = 2;
  int C = 4;
  Matrix<int, Dynamic, 1> y(N);
  y << 1, 4, 3, 3, 2;
  VectorXd cuts(C - 1);
  cuts << 0.9, 1.1, 7;
  VectorXd beta(M);
  beta << 1.1, 0.4;
  MatrixXd x(N, M);
  x << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
  EXPECT_FLOAT_EQ(ordered_logistic_glm_lpmf(y, x, beta, cuts), ordered_logistic_glm_simple_lpmf<false>(y, x, beta, cuts));
}

TEST(ProbDistributionsOrderedLogisticGLM, glm_matches_ordered_logistic_vars) {
  double eps = 1e-13;
  int N = 5;
  int M = 2;
  int C = 3;
  Matrix<int, Dynamic, 1> y(N);
  y << 1, 1, 2, 4, 4;
  Matrix<var, Dynamic, 1> cuts1(C), cuts2(C);
  cuts1 << 0.9, 1.1, 7;
  cuts2 << 0.9, 1.1, 7;
  Matrix<var, Dynamic, 1> beta1(M), beta2(M);
  beta1 << 1.1, 0.4;
  beta2 << 1.1, 0.4;
  Matrix<var, Dynamic, Dynamic> x1(N, M), x2(N, M);
  x1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
  x2 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
  var res1 = ordered_logistic_glm_lpmf(y, x1, beta1, cuts1);
  var res2 = ordered_logistic_glm_simple_lpmf<false>(y, x2, beta2, cuts2);
  (res1 + res2).grad();

  Matrix<double, Dynamic, Dynamic> x_adj(N, M);

  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      x_adj(j, i) = x2(j, i).adj();
      EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
    }
  }
  for (int i = 0; i < M; i++) {
    EXPECT_NEAR(beta1[i].adj(), beta2[i].adj(), eps);
  }
  for (int i = 0; i < C; i++) {
    EXPECT_NEAR(cuts1[i].adj(), cuts2[i].adj(), eps);
  }
}

TEST(ProbDistributionsOrderedLogisticGLM, glm_matches_ordered_logistic_vars_big) {
  double eps = 1e-7;
  int N = 155;
  int M = 15;
  int C = 10;
  Matrix<int, Dynamic, 1> y(N);
  for (int i = 0; i < N; i++) {
    y[i] = Matrix<unsigned int, Dynamic, 1>::Random(1)[0] % (C + 1) + 1;
  }
  VectorXd cuts_double = (VectorXd::Random(C).array() + 1) / 2 / C;
  for (int i = 1; i < C; i++) {
    cuts_double[i] += cuts_double[i - 1];
  }
  Matrix<var, Dynamic, 1> cuts1 = cuts_double, cuts2 = cuts_double;
  VectorXd beta_double = VectorXd::Random(M);
  Matrix<var, Dynamic, 1> beta1 = beta_double, beta2 = beta_double;
  MatrixXd x_double = MatrixXd::Random(N, M);
  Matrix<var, Dynamic, Dynamic> x1 = x_double, x2 = x_double;
  var res1 = ordered_logistic_glm_lpmf(y, x1, beta1, cuts1);
  var res2 = ordered_logistic_glm_simple_lpmf<false>(y, x2, beta2, cuts2);
  (res1 + res2).grad();

  Matrix<double, Dynamic, Dynamic> x_adj(N, M);

  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      x_adj(j, i) = x2(j, i).adj();
      EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
    }
  }
  for (int i = 0; i < M; i++) {
    EXPECT_NEAR(beta1[i].adj(), beta2[i].adj(), eps);
  }
  for (int i = 0; i < C; i++) {
    EXPECT_NEAR(cuts1[i].adj(), cuts2[i].adj(), eps);
  }
}

TEST(ProbDistributionsOrderedLogisticGLM, glm_interfaces) {
  int N = 5;
  int M = 2;
  int C = 3;
  Matrix<int, Dynamic, 1> y(N);
  y << 1, 1, 2, 4, 4;
  std::vector<int> y_vec = {1, 1, 2, 4, 4};
  Matrix<double, Dynamic, 1> cuts_double(C);
  cuts_double << 0.9, 1.1, 7;
  Matrix<var, Dynamic, 1> cuts_var = cuts_double;
  std::vector<double> cuts_vec_double = {0.5, 1, 4};
  std::vector<var> cuts_vec_var;
  for (int i = 0; i < C; i++) {
    cuts_vec_var.emplace_back(cuts_vec_double[i]);
  }
  Matrix<double, Dynamic, 1> beta_double(M);
  beta_double << 1.1, 0.4;
  Matrix<var, Dynamic, 1> beta_var = beta_double;
  std::vector<double> beta_vec_double = {1.1, 0.4};
  std::vector<var> beta_vec_var;
  for (int i = 0; i < M; i++) {
    beta_vec_var.emplace_back(beta_vec_double[i]);
  }
  Matrix<double, Dynamic, Dynamic> x_double(N, M);
  x_double << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
  Matrix<var, Dynamic, Dynamic> x_var = x_double;

  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_double, cuts_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_var, beta_double, cuts_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_var, cuts_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_double, cuts_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_var, cuts_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_var, beta_double, cuts_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_var, beta_var, cuts_double));

  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_double, cuts_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_var, beta_double, cuts_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_var, cuts_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_double, cuts_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_var, cuts_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_var, beta_double, cuts_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_var, beta_var, cuts_double));

  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_vec_double, cuts_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_var, beta_vec_double, cuts_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_vec_var, cuts_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_vec_double, cuts_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_vec_var, cuts_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_var, beta_vec_double, cuts_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_var, beta_vec_var, cuts_double));

  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_vec_double, cuts_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_var, beta_vec_double, cuts_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_vec_var, cuts_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_vec_double, cuts_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_vec_var, cuts_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_var, beta_vec_double, cuts_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_var, beta_vec_var, cuts_double));


  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_double, cuts_vec_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_var, beta_double, cuts_vec_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_var, cuts_vec_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_double, cuts_vec_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_var, cuts_vec_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_var, beta_double, cuts_vec_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_var, beta_var, cuts_vec_double));

  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_double, cuts_vec_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_var, beta_double, cuts_vec_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_var, cuts_vec_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_double, cuts_vec_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_var, cuts_vec_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_var, beta_double, cuts_vec_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_var, beta_var, cuts_vec_double));

  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_vec_double, cuts_vec_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_var, beta_vec_double, cuts_vec_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_vec_var, cuts_vec_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_vec_double, cuts_vec_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_double, beta_vec_var, cuts_vec_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_var, beta_vec_double, cuts_vec_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x_var, beta_vec_var, cuts_vec_double));

  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_vec_double, cuts_vec_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_var, beta_vec_double, cuts_vec_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_vec_var, cuts_vec_double));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_vec_double, cuts_vec_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_double, beta_vec_var, cuts_vec_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_var, beta_vec_double, cuts_vec_var));
  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y_vec, x_var, beta_vec_var, cuts_vec_double));
}

TEST(ProbDistributionsOrderedLogisticGLM, glm_errors) {
  int N = 5;
  int M = 2;
  int C = 3;
  Matrix<int, Dynamic, 1> y(N);
  y << 1, 1, 2, 4, 4;
  Matrix<int, Dynamic, 1> y_size(N + 1);
  y_size << 1, 1, 2, 4, 4, 2;
  Matrix<int, Dynamic, 1> y_val1(N);
  y_val1 << 1, 1, 2, 4, 5;
  Matrix<int, Dynamic, 1> y_val2(N);
  y_val2 << 1, 1, 2, 4, 0;
  Matrix<double, Dynamic, 1> cuts(C);
  cuts << 0.9, 1.1, 7;
  Matrix<double, Dynamic, 1> cuts_val1(C);
  cuts_val1 << 0.9, 1.1, INFINITY;
  Matrix<double, Dynamic, 1> cuts_val2(C);
  cuts_val2 << -INFINITY, 1.1, 4;
  Matrix<double, Dynamic, 1> cuts_val3(C);
  cuts_val3 << 0, 1.1, 0.5;
  Matrix<double, Dynamic, 1> cuts_val4(C);
  cuts_val4 << 0, NAN, 0.5;
  Matrix<double, Dynamic, 1> beta(M);
  beta << 1.1, 0.4;
  Matrix<double, Dynamic, 1> beta_size(M + 1);
  beta_size << 1.1, 0.4, 0;
  Matrix<double, Dynamic, 1> beta_val(M);
  beta_val << 1.1, INFINITY;
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
  Matrix<double, Dynamic, Dynamic> x_size1(N + 1, M);
  x_size1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 9, 8;
  Matrix<double, Dynamic, Dynamic> x_size2(N, M + 1);
  x_size2 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 4, 5, 6, 7, 9;
  Matrix<double, Dynamic, Dynamic> x_val(N, M);
  x_val << 1, 2, 3, 4, 5, 6, 7, INFINITY, 9, 0;

  EXPECT_NO_THROW(ordered_logistic_glm_lpmf(y, x, beta, cuts));

  EXPECT_THROW(ordered_logistic_glm_lpmf(y_size, x, beta, cuts), std::invalid_argument);
  EXPECT_THROW(ordered_logistic_glm_lpmf(y, x_size1, beta, cuts), std::invalid_argument);
  EXPECT_THROW(ordered_logistic_glm_lpmf(y, x_size2, beta, cuts), std::invalid_argument);
  EXPECT_THROW(ordered_logistic_glm_lpmf(y, x, beta_size, cuts), std::invalid_argument);

  EXPECT_THROW(ordered_logistic_glm_lpmf(y_val1, x, beta, cuts), std::domain_error);
  EXPECT_THROW(ordered_logistic_glm_lpmf(y_val2, x, beta, cuts), std::domain_error);
  EXPECT_THROW(ordered_logistic_glm_lpmf(y, x_val, beta, cuts), std::domain_error);
  EXPECT_THROW(ordered_logistic_glm_lpmf(y, x, beta_val, cuts), std::domain_error);
  EXPECT_THROW(ordered_logistic_glm_lpmf(y, x, beta, cuts_val1), std::domain_error);
  EXPECT_THROW(ordered_logistic_glm_lpmf(y, x, beta, cuts_val2), std::domain_error);
  EXPECT_THROW(ordered_logistic_glm_lpmf(y, x, beta, cuts_val3), std::domain_error);
  EXPECT_THROW(ordered_logistic_glm_lpmf(y, x, beta, cuts_val4), std::domain_error);
}
