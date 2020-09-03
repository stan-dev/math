#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <Eigen/Core>
#include <vector>

template <bool propto, typename T_x, typename T_beta, typename T_cuts>
stan::return_type_t<T_x, T_beta, T_cuts> ordered_logistic_glm_simple_lpmf(
    const std::vector<int>& y,
    const Eigen::Matrix<T_x, Eigen::Dynamic, Eigen::Dynamic>& x,
    const T_beta& beta, const T_cuts& cuts) {
  using T_x_beta = stan::return_type_t<T_x, T_beta>;
  using stan::math::as_column_vector_or_scalar;

  auto& beta_col = as_column_vector_or_scalar(beta);

  Eigen::Matrix<T_x_beta, Eigen::Dynamic, 1> location
      = x.template cast<T_x_beta>() * beta_col.template cast<T_x_beta>();

  return stan::math::ordered_logistic_lpmf<propto>(y, location, cuts);
}

TEST(ProbDistributionsOrderedLogisticGLM,
     glm_matches_ordered_logistic_doubles) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using std::vector;
  double eps = 1e-13;
  int N = 5;
  int M = 2;
  int C = 4;
  vector<int> y{1, 4, 3, 3, 2};
  VectorXd cuts(C - 1);
  cuts << 0.9, 1.1, 7;
  VectorXd beta(M);
  beta << 1.1, 0.4;
  MatrixXd x(N, M);
  x << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
  EXPECT_FLOAT_EQ(stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts),
                  ordered_logistic_glm_simple_lpmf<false>(y, x, beta, cuts));
}

TEST(ProbDistributionsOrderedLogisticGLM,
     glm_matches_ordered_logistic_doubles_broadcast_y) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using std::vector;
  double eps = 1e-13;
  int N = 5;
  int M = 2;
  int C = 4;
  VectorXd cuts(C - 1);
  cuts << 0.9, 1.1, 7;
  VectorXd beta(M);
  beta << 1.1, 0.4;
  MatrixXd x(N, M);
  x << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
  for (int y_scal = 1; y_scal <= C; y_scal++) {
    vector<int> y(N, y_scal);
    EXPECT_FLOAT_EQ(
        stan::math::ordered_logistic_glm_lpmf(y_scal, x, beta, cuts),
        stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts));
  }
}

TEST(ProbDistributionsOrderedLogisticGLM, glm_matches_ordered_logistic_vars) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::var;
  using std::vector;
  double eps = 1e-13;
  int N = 5;
  int M = 2;
  int C = 3;
  vector<int> y{1, 1, 2, 4, 4};
  Matrix<var, Dynamic, 1> cuts1(C), cuts2(C);
  cuts1 << 0.9, 1.1, 7;
  cuts2 << 0.9, 1.1, 7;
  Matrix<var, Dynamic, 1> beta1(M), beta2(M);
  beta1 << 1.1, 0.4;
  beta2 << 1.1, 0.4;
  Matrix<var, Dynamic, Dynamic> x1(N, M), x2(N, M);
  x1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
  x2 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
  var res1 = stan::math::ordered_logistic_glm_lpmf(y, x1, beta1, cuts1);
  var res2 = ordered_logistic_glm_simple_lpmf<false>(y, x2, beta2, cuts2);
  (res1 + res2).grad();

  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
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

TEST(ProbDistributionsOrderedLogisticGLM,
     glm_matches_ordered_logistic_vars_broadcast_y) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::var;
  using std::vector;
  double eps = 1e-13;
  int N = 5;
  int M = 2;
  int C = 3;
  Matrix<var, Dynamic, 1> cuts1(C), cuts2(C);
  cuts1 << 0.9, 1.1, 7;
  cuts2 << 0.9, 1.1, 7;
  Matrix<var, Dynamic, 1> beta1(M), beta2(M);
  beta1 << 1.1, 0.4;
  beta2 << 1.1, 0.4;
  Matrix<var, Dynamic, Dynamic> x1(N, M), x2(N, M);
  x1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
  x2 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;

  for (int y_scal = 1; y_scal <= C; y_scal++) {
    vector<int> y(N, y_scal);
    var res1 = ordered_logistic_glm_lpmf(y, x1, beta1, cuts1);
    var res2 = ordered_logistic_glm_lpmf(y_scal, x2, beta2, cuts2);
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
}

TEST(ProbDistributionsOrderedLogisticGLM,
     glm_matches_ordered_logistic_vars_broadcast_x) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::var;
  using std::vector;
  double eps = 1e-13;
  int N = 5;
  int M = 2;
  int C = 3;
  vector<int> y{1, 1, 2, 4, 4};
  Matrix<var, Dynamic, 1> cuts1(C), cuts2(C);
  cuts1 << 0.9, 1.1, 7;
  cuts2 << 0.9, 1.1, 7;
  Matrix<var, Dynamic, 1> beta1(M), beta2(M);
  beta1 << 1.1, 0.4;
  beta2 << 1.1, 0.4;
  RowVectorXd x_double(M);
  x_double << 1, 2;
  Matrix<var, 1, Dynamic> x_row = x_double;
  Matrix<var, Dynamic, Dynamic> x = x_double.replicate(N, 1);

  var res1 = ordered_logistic_glm_lpmf(y, x_row, beta1, cuts1);
  var res2 = ordered_logistic_glm_lpmf(y, x, beta2, cuts2);
  (res1 + res2).grad();

  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < M; i++) {
    double x_sum = 0;
    for (int j = 0; j < N; j++) {
      x_sum += x(j, i).adj();
    }
    EXPECT_NEAR(x_row[i].adj(), x_sum, eps);
  }
  for (int i = 0; i < M; i++) {
    EXPECT_NEAR(beta1[i].adj(), beta2[i].adj(), eps);
  }
  for (int i = 0; i < C; i++) {
    EXPECT_NEAR(cuts1[i].adj(), cuts2[i].adj(), eps);
  }
}

TEST(ProbDistributionsOrderedLogisticGLM,
     glm_matches_ordered_logistic_single_instance) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::var;
  using std::vector;
  double eps = 1e-13;
  int N = 1;
  int M = 2;
  int C = 3;
  vector<int> y{1};
  Matrix<var, Dynamic, 1> cuts1(C), cuts2(C);
  cuts1 << 0.9, 1.1, 7;
  cuts2 << 0.9, 1.1, 7;
  Matrix<var, Dynamic, 1> beta1(M), beta2(M);
  beta1 << 1.1, 0.4;
  beta2 << 1.1, 0.4;
  Matrix<var, Dynamic, Dynamic> x1(N, M), x2(N, M);
  x1 << 1, 2;
  x2 << 1, 2;
  var res1 = ordered_logistic_glm_lpmf(y, x1, beta1, cuts1);
  var res2 = ordered_logistic_glm_simple_lpmf<false>(y, x2, beta2, cuts2);
  (res1 + res2).grad();

  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < M; i++) {
    EXPECT_NEAR(x1(0, i).adj(), x2(0, i).adj(), eps);
  }
  for (int i = 0; i < M; i++) {
    EXPECT_NEAR(beta1[i].adj(), beta2[i].adj(), eps);
  }
  for (int i = 0; i < C; i++) {
    EXPECT_NEAR(cuts1[i].adj(), cuts2[i].adj(), eps);
  }
}

TEST(ProbDistributionsOrderedLogisticGLM,
     glm_matches_ordered_logistic_zero_instances) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::var;
  using std::vector;
  double eps = 1e-13;
  int N = 0;
  int M = 2;
  int C = 3;
  vector<int> y{};
  Matrix<var, Dynamic, 1> cuts1(C), cuts2(C);
  cuts1 << 0.9, 1.1, 7;
  cuts2 << 0.9, 1.1, 7;
  Matrix<var, Dynamic, 1> beta1(M), beta2(M);
  beta1 << 1.1, 0.4;
  beta2 << 1.1, 0.4;
  Matrix<var, Dynamic, Dynamic> x1(N, M), x2(N, M);
  var res1 = ordered_logistic_glm_lpmf(y, x1, beta1, cuts1);
  var res2 = ordered_logistic_glm_simple_lpmf<false>(y, x2, beta2, cuts2);
  (res1 + res2).grad();

  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < M; i++) {
    EXPECT_NEAR(beta1[i].adj(), beta2[i].adj(), eps);
  }
  for (int i = 0; i < C; i++) {
    EXPECT_NEAR(cuts1[i].adj(), cuts2[i].adj(), eps);
  }
}

TEST(ProbDistributionsOrderedLogisticGLM,
     glm_matches_ordered_logistic_zero_attributes) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::var;
  using std::vector;
  double eps = 1e-13;
  int N = 5;
  int M = 0;
  int C = 3;
  vector<int> y{1, 1, 2, 4, 4};
  Matrix<var, Dynamic, 1> cuts1(C), cuts2(C);
  cuts1 << 0.9, 1.1, 7;
  cuts2 << 0.9, 1.1, 7;
  Matrix<var, Dynamic, 1> beta1(M), beta2(M);
  Matrix<var, Dynamic, Dynamic> x1(N, M), x2(N, M);
  var res1 = ordered_logistic_glm_lpmf(y, x1, beta1, cuts1);
  var res2 = ordered_logistic_glm_simple_lpmf<false>(y, x2, beta2, cuts2);
  (res1 + res2).grad();

  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < C; i++) {
    EXPECT_NEAR(cuts1[i].adj(), cuts2[i].adj(), eps);
  }
}

TEST(ProbDistributionsOrderedLogisticGLM,
     glm_matches_ordered_logistic_single_class) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::var;
  using std::vector;
  double eps = 1e-13;
  int N = 5;
  int M = 2;
  int C = 1;
  vector<int> y{1, 1, 1, 1, 1};
  Matrix<var, Dynamic, 1> cuts1(C), cuts2(C);
  cuts1 << 0.9;
  cuts2 << 0.9;
  Matrix<var, Dynamic, 1> beta1(M), beta2(M);
  beta1 << 1.1, 0.4;
  beta2 << 1.1, 0.4;
  Matrix<var, Dynamic, Dynamic> x1(N, M), x2(N, M);
  x1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
  x2 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
  var res1 = ordered_logistic_glm_lpmf(y, x1, beta1, cuts1);
  var res2 = ordered_logistic_glm_simple_lpmf<false>(y, x2, beta2, cuts2);
  (res1 + res2).grad();

  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
    }
  }
  for (int i = 0; i < M; i++) {
    EXPECT_NEAR(beta1[i].adj(), beta2[i].adj(), eps);
  }
  EXPECT_NEAR(cuts1[0].adj(), cuts2[0].adj(), eps);
}

TEST(ProbDistributionsOrderedLogisticGLM,
     glm_matches_ordered_logistic_vars_big) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::var;
  using std::vector;
  double eps = 1e-7;
  int N = 155;
  int M = 15;
  int C = 10;
  vector<int> y(N);
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
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::var;
  using std::vector;
  int N = 5;
  int M = 2;
  int C = 3;
  vector<int> y = {1, 1, 2, 4, 4};
  int y_scal = 1;
  VectorXd cuts_double(C);
  cuts_double << 0.9, 1.1, 7;
  Matrix<var, Dynamic, 1> cuts_var = cuts_double;
  VectorXd beta_double(M);
  beta_double << 1.1, 0.4;
  Matrix<var, Dynamic, 1> beta_var = beta_double;
  Matrix<double, Dynamic, Dynamic> x_double(N, M);
  x_double << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
  Matrix<var, Dynamic, Dynamic> x_var = x_double;
  Matrix<double, 1, Dynamic> x_row_double(M);
  x_row_double << 1, 2;
  Matrix<var, 1, Dynamic> x_row_var = x_row_double;

  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(
      y, x_double, beta_double, cuts_double));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y, x_var, beta_double,
                                                        cuts_double));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y, x_double, beta_var,
                                                        cuts_double));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y, x_double,
                                                        beta_double, cuts_var));
  EXPECT_NO_THROW(
      stan::math::ordered_logistic_glm_lpmf(y, x_double, beta_var, cuts_var));
  EXPECT_NO_THROW(
      stan::math::ordered_logistic_glm_lpmf(y, x_var, beta_double, cuts_var));
  EXPECT_NO_THROW(
      stan::math::ordered_logistic_glm_lpmf(y, x_var, beta_var, cuts_double));

  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(
      y_scal, x_double, beta_double, cuts_double));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(
      y_scal, x_var, beta_double, cuts_double));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y_scal, x_double,
                                                        beta_var, cuts_double));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y_scal, x_double,
                                                        beta_double, cuts_var));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y_scal, x_double,
                                                        beta_var, cuts_var));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y_scal, x_var,
                                                        beta_double, cuts_var));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y_scal, x_var, beta_var,
                                                        cuts_double));

  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(
      y, x_row_double, beta_double, cuts_double));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(
      y, x_row_var, beta_double, cuts_double));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y, x_row_double,
                                                        beta_var, cuts_double));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y, x_row_double,
                                                        beta_double, cuts_var));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y, x_row_double,
                                                        beta_var, cuts_var));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y, x_row_var,
                                                        beta_double, cuts_var));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y, x_row_var, beta_var,
                                                        cuts_double));

  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(
      y_scal, x_row_double, beta_double, cuts_double));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(
      y_scal, x_row_var, beta_double, cuts_double));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y_scal, x_row_double,
                                                        beta_var, cuts_double));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y_scal, x_row_double,
                                                        beta_double, cuts_var));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y_scal, x_row_double,
                                                        beta_var, cuts_var));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y_scal, x_row_var,
                                                        beta_double, cuts_var));
  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y_scal, x_row_var,
                                                        beta_var, cuts_double));
}

TEST(ProbDistributionsOrderedLogisticGLM, glm_errors) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::var;
  using std::vector;
  int N = 5;
  int M = 2;
  int C = 3;
  vector<int> y{1, 1, 2, 4, 4};
  vector<int> y_size{1, 1, 2, 4, 4, 2};
  vector<int> y_val1{1, 1, 2, 4, 5};
  vector<int> y_val2{1, 1, 2, 4, 0};
  VectorXd cuts(C);
  cuts << 0.9, 1.1, 7;
  VectorXd cuts_val1(C);
  cuts_val1 << 0.9, 1.1, INFINITY;
  VectorXd cuts_val2(C);
  cuts_val2 << -INFINITY, 1.1, 4;
  VectorXd cuts_val3(C);
  cuts_val3 << 0, 1.1, 0.5;
  VectorXd cuts_val4(C);
  cuts_val4 << 0, NAN, 0.5;
  VectorXd beta(M);
  beta << 1.1, 0.4;
  VectorXd beta_size(M + 1);
  beta_size << 1.1, 0.4, 0;
  VectorXd beta_val(M);
  beta_val << 1.1, INFINITY;
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
  Matrix<double, Dynamic, Dynamic> x_size1(N + 1, M);
  x_size1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 9, 8;
  Matrix<double, Dynamic, Dynamic> x_size2(N, M + 1);
  x_size2 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 4, 5, 6, 7, 9;
  Matrix<double, 1, Dynamic> x_size3(M + 1);
  x_size3 << 1, 2, 3;
  Matrix<double, Dynamic, Dynamic> x_val(N, M);
  x_val << 1, 2, 3, 4, 5, 6, 7, INFINITY, 9, 0;

  EXPECT_NO_THROW(stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts));

  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y_size, x, beta, cuts),
               std::invalid_argument);
  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y, x_size1, beta, cuts),
               std::invalid_argument);
  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y, x_size2, beta, cuts),
               std::invalid_argument);
  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y, x_size3, beta, cuts),
               std::invalid_argument);
  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y, x, beta_size, cuts),
               std::invalid_argument);

  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y_val1, x, beta, cuts),
               std::domain_error);
  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y_val2, x, beta, cuts),
               std::domain_error);
  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y, x_val, beta, cuts),
               std::domain_error);
  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y, x, beta_val, cuts),
               std::domain_error);
  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts_val1),
               std::domain_error);
  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts_val2),
               std::domain_error);
  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts_val3),
               std::domain_error);
  EXPECT_THROW(stan::math::ordered_logistic_glm_lpmf(y, x, beta, cuts_val4),
               std::domain_error);
}
