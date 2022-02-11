#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

template <bool propto, typename T_x, typename T_alpha, typename T_beta>
stan::return_type_t<T_x, T_alpha, T_beta> categorical_logit_glm_simple_lpmf(
    const std::vector<int>& y, const T_x& x, const T_alpha& alpha,
    const T_beta& beta) {
  using T_x_beta = stan::return_type_t<T_x, T_beta>;
  using T_return = stan::return_type_t<T_x, T_beta, T_alpha>;

  const size_t N_instances = x.rows();

  const auto& alpha_row
      = stan::math::as_column_vector_or_scalar(alpha).transpose();

  auto tmp = stan::math::to_ref(
      stan::math::multiply(x, beta)
      + stan::math::rep_matrix<std::decay_t<decltype(alpha_row)>>(alpha_row,
                                                                  x.rows()));

  T_return lpmf = 0;
  // iterate overt instances
  for (int i = 0; i < N_instances; i++) {
    lpmf += stan::math::categorical_logit_lpmf<propto>(
        y[i], tmp.row(i).transpose().eval());
  }
  return lpmf;
}

TEST(ProbDistributionsCategoricalLogitGLM,
     glm_matches_categorical_logit_doubles) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::categorical_logit_glm_lpmf;
  using stan::math::var;
  using std::vector;
  double eps = 1e-13;
  const size_t N_instances = 5;
  const size_t N_attributes = 2;
  const size_t N_classes = 3;
  vector<int> y{1, 3, 1, 2, 2};
  MatrixXd x(N_instances, N_attributes);
  x << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  MatrixXd beta(N_attributes, N_classes);
  beta << 0.3, 2, 0.4, -0.1, -1.3, 1;
  VectorXd alpha(N_classes);
  alpha << 0.5, -2, 4;

  EXPECT_NEAR((categorical_logit_glm_simple_lpmf<false>(y, x, alpha, beta)),
              (categorical_logit_glm_lpmf(y, x, alpha, beta)), eps);
  EXPECT_NEAR((categorical_logit_glm_simple_lpmf<true>(y, x, alpha, beta)),
              (categorical_logit_glm_lpmf<true>(y, x, alpha, beta)), eps);
}

template <class T>
class ProbDistributionsCategoricalLogitGLM
    : public stan::math::test::VarMatrixTypedTests<T> {};

TYPED_TEST_SUITE(ProbDistributionsCategoricalLogitGLM,
                 stan::math::test::VarMatImpls);

TYPED_TEST(ProbDistributionsCategoricalLogitGLM,
           glm_matches_categorical_logit_vars) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::categorical_logit_glm_lpmf;
  using stan::math::var;
  using std::vector;
  using matrix_v = typename TypeParam::matrix_v;
  using vector_v = typename TypeParam::vector_v;
  using row_vector_v = typename TypeParam::row_vector_v;
  double eps = 1e-13;
  const size_t N_instances = 5;
  const size_t N_attributes = 2;
  const size_t N_classes = 3;
  vector<int> y{1, 3, 1, 2, 2};
  Matrix<double, Dynamic, Dynamic> x1_val(N_instances, N_attributes);
  x1_val << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  matrix_v x1 = x1_val;
  matrix_v x2 = x1.val();
  Matrix<double, Dynamic, Dynamic> beta1_val(N_attributes, N_classes);
  beta1_val << 0.3, 2, 0.4, -0.1, -1.3, 1;
  matrix_v beta1 = beta1_val;
  matrix_v beta2 = beta1.val();
  Matrix<double, Dynamic, 1> alpha1_val(N_classes);
  alpha1_val << 0.5, -2, 4;
  vector_v alpha1 = alpha1_val;
  vector_v alpha2 = alpha1.val();

  var res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1.adj()(j, i), x2.adj()(j, i), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
  }

  stan::math::set_zero_all_adjoints();

  res1 = categorical_logit_glm_simple_lpmf<true>(y, x1, alpha1, beta1);
  res2 = categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1.adj()(j, i), x2.adj()(j, i), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
  }
}

TYPED_TEST(ProbDistributionsCategoricalLogitGLM, single_instance) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::categorical_logit_glm_lpmf;
  using stan::math::var;
  using std::vector;
  using matrix_v = typename TypeParam::matrix_v;
  using vector_v = typename TypeParam::vector_v;
  using row_vector_v = typename TypeParam::row_vector_v;

  double eps = 1e-13;
  const size_t N_instances = 1;
  const size_t N_attributes = 2;
  const size_t N_classes = 3;
  vector<int> y{1};
  Matrix<double, Dynamic, Dynamic> x1_val(N_instances, N_attributes);
  x1_val << -12, 46;
  matrix_v x1 = x1_val;
  matrix_v x2 = x1.val();
  Matrix<double, Dynamic, Dynamic> beta1_val(N_attributes, N_classes);
  beta1_val << 0.3, 2, 0.4, -0.1, -1.3, 1;
  matrix_v beta1 = beta1_val;
  matrix_v beta2 = beta1_val;
  Matrix<double, Dynamic, 1> alpha1_val(N_classes);
  alpha1_val << 0.5, -2, 4;
  vector_v alpha1 = alpha1_val;
  vector_v alpha2 = alpha1.val();
  var res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1.adj()(j, i), x2.adj()(j, i), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
  }

  stan::math::set_zero_all_adjoints();

  res1 = categorical_logit_glm_simple_lpmf<true>(y, x1, alpha1, beta1);
  res2 = categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1.adj()(j, i), x2.adj()(j, i), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
  }
}

TYPED_TEST(ProbDistributionsCategoricalLogitGLM, zero_instances) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::categorical_logit_glm_lpmf;
  using stan::math::var;
  using std::vector;
  using matrix_v = typename TypeParam::matrix_v;
  using vector_v = typename TypeParam::vector_v;
  using row_vector_v = typename TypeParam::row_vector_v;
  double eps = 1e-13;
  const size_t N_instances = 0;
  const size_t N_attributes = 2;
  const size_t N_classes = 3;
  std::vector<int> y{};
  matrix_v x1 = Matrix<double, Dynamic, Dynamic>(N_instances, N_attributes);
  matrix_v x2 = x1.val();
  Matrix<double, Dynamic, Dynamic> beta1_val(N_attributes, N_classes);
  beta1_val << 0.3, 2, 0.4, -0.1, -1.3, 1;
  matrix_v beta1 = beta1_val;
  matrix_v beta2 = beta1_val;
  Matrix<double, Dynamic, 1> alpha1_val(N_classes);
  alpha1_val << 0.5, -2, 4;
  vector_v alpha1 = alpha1_val;
  vector_v alpha2 = alpha1.val();

  var res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
  }

  stan::math::set_zero_all_adjoints();

  res1 = categorical_logit_glm_simple_lpmf<true>(y, x1, alpha1, beta1);
  res2 = categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
  }
}

TYPED_TEST(ProbDistributionsCategoricalLogitGLM, single_class) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::categorical_logit_glm_lpmf;
  using stan::math::var;
  using std::vector;
  using matrix_v = typename TypeParam::matrix_v;
  using vector_v = typename TypeParam::vector_v;
  using row_vector_v = typename TypeParam::row_vector_v;
  double eps = 1e-13;
  const size_t N_instances = 5;
  const size_t N_attributes = 2;
  const size_t N_classes = 1;
  vector<int> y{1, 1, 1, 1, 1};
  Matrix<double, Dynamic, Dynamic> x1_val(N_instances, N_attributes);
  x1_val << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  matrix_v x1 = x1_val;
  matrix_v x2 = x1.val();
  Matrix<double, Dynamic, Dynamic> beta1_val(N_attributes, N_classes);
  beta1_val << 0.3, 2;
  matrix_v beta1 = beta1_val;
  matrix_v beta2 = beta1_val;
  Matrix<double, Dynamic, 1> alpha1_val(N_classes);
  alpha1_val << 0.5;
  vector_v alpha1 = alpha1_val;
  vector_v alpha2 = alpha1.val();

  var res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1.adj()(j, i), x2.adj()(j, i), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
  }

  stan::math::set_zero_all_adjoints();

  res1 = categorical_logit_glm_simple_lpmf<true>(y, x1, alpha1, beta1);
  res2 = categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1.adj()(j, i), x2.adj()(j, i), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
  }
}

TYPED_TEST(ProbDistributionsCategoricalLogitGLM, zero_attributes) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::categorical_logit_glm_lpmf;
  using stan::math::var;
  using std::vector;
  using matrix_v = typename TypeParam::matrix_v;
  using vector_v = typename TypeParam::vector_v;
  using row_vector_v = typename TypeParam::row_vector_v;
  double eps = 1e-13;
  const size_t N_instances = 5;
  const size_t N_attributes = 0;
  const size_t N_classes = 3;
  vector<int> y{1, 3, 1, 2, 2};
  Matrix<double, Dynamic, Dynamic> x_val(N_instances, N_attributes);
  matrix_v x1 = x_val;
  matrix_v x2 = x_val;
  Matrix<double, Dynamic, Dynamic> beta_val(N_attributes, N_classes);
  matrix_v beta1 = beta_val;
  matrix_v beta2 = beta_val;
  Matrix<double, Dynamic, 1> alpha1_val(N_classes);
  alpha1_val << 0.5, -2, 4;
  vector_v alpha1 = alpha1_val;
  vector_v alpha2 = alpha1.val();

  var res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
  }

  stan::math::set_zero_all_adjoints();

  res1 = categorical_logit_glm_simple_lpmf<true>(y, x1, alpha1, beta1);
  res2 = categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
  }
}

TYPED_TEST(ProbDistributionsCategoricalLogitGLM, x_broadcasting) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::categorical_logit_glm_lpmf;
  using stan::math::var;
  using std::vector;
  using matrix_v = typename TypeParam::matrix_v;
  using vector_v = typename TypeParam::vector_v;
  using row_vector_v = typename TypeParam::row_vector_v;
  double eps = 1e-13;
  const size_t N_instances = 5;
  const size_t N_attributes = 2;
  const size_t N_classes = 3;
  vector<int> y{1, 3, 1, 2, 2};
  Matrix<double, 1, Dynamic> x_double(N_attributes);
  x_double << -12, 46;
  row_vector_v x_row = x_double;
  matrix_v x = x_double.replicate(N_instances, 1);
  Matrix<double, Dynamic, Dynamic> beta_val(N_attributes, N_classes);
  beta_val << 0.3, 2, 0.4, -0.1, -1.3, 1;
  matrix_v beta1 = beta_val;
  matrix_v beta2 = beta_val;
  Matrix<double, Dynamic, 1> alpha_val(N_classes);
  alpha_val << 0.5, -2, 4;
  vector_v alpha1 = alpha_val;
  vector_v alpha2 = alpha_val;
  var res1 = categorical_logit_glm_lpmf(y, x_row, alpha1, beta1);
  var res2 = categorical_logit_glm_lpmf(y, x, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    double x_sum = 0;
    for (int j = 0; j < N_instances; j++) {
      x_sum += x.adj()(j, i);
    }
    EXPECT_NEAR(x_row.adj()[i], x_sum, eps);
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
  }

  stan::math::set_zero_all_adjoints();

  res1 = categorical_logit_glm_lpmf<true>(y, x_row, alpha1, beta1);
  res2 = categorical_logit_glm_lpmf<true>(y, x, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    double x_sum = 0;
    for (int j = 0; j < N_instances; j++) {
      x_sum += x.adj()(j, i);
    }
    EXPECT_NEAR(x_row.adj()[i], x_sum, eps);
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
  }
}

TYPED_TEST(ProbDistributionsCategoricalLogitGLM, y_broadcasting) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::categorical_logit_glm_lpmf;
  using stan::math::var;
  using std::vector;
  using matrix_v = typename TypeParam::matrix_v;
  using vector_v = typename TypeParam::vector_v;
  using row_vector_v = typename TypeParam::row_vector_v;
  double eps = 1e-13;
  const size_t N_instances = 5;
  const size_t N_attributes = 2;
  const size_t N_classes = 3;
  Matrix<double, Dynamic, Dynamic> x_val(N_instances, N_attributes);
  x_val << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  matrix_v x1 = x_val;
  matrix_v x2 = x_val;
  Matrix<double, Dynamic, Dynamic> beta_val(N_attributes, N_classes);
  beta_val << 0.3, 2, 0.4, -0.1, -1.3, 1;
  matrix_v beta1 = beta_val;
  matrix_v beta2 = beta_val;
  Matrix<double, Dynamic, 1> alpha_val(N_classes);
  alpha_val << 0.5, -2, 4;
  vector_v alpha1 = alpha_val;
  vector_v alpha2 = alpha_val;
  for (int y_scal = 1; y_scal <= N_classes; y_scal++) {
    vector<int> y(N_instances, y_scal);
    var res1 = categorical_logit_glm_lpmf(y_scal, x1, alpha1, beta1);
    var res2 = categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
    (res1 + res2).grad();
    EXPECT_NEAR(res1.val(), res2.val(), eps);
    for (int i = 0; i < N_attributes; i++) {
      for (int j = 0; j < N_instances; j++) {
        EXPECT_NEAR(x1.adj()(j, i), x2.adj()(j, i), eps);
      }
      for (int j = 0; j < N_classes; j++) {
        EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
      }
    }
    for (int i = 0; i < N_classes; i++) {
      EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
    }
    stan::math::set_zero_all_adjoints();
    res1 = categorical_logit_glm_lpmf<true>(y_scal, x1, alpha1, beta1);
    res2 = categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
    (res1 + res2).grad();
    EXPECT_NEAR(res1.val(), res2.val(), eps);
    for (int i = 0; i < N_attributes; i++) {
      for (int j = 0; j < N_instances; j++) {
        EXPECT_NEAR(x1.adj()(j, i), x2.adj()(j, i), eps);
      }
      for (int j = 0; j < N_classes; j++) {
        EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
      }
    }
    for (int i = 0; i < N_classes; i++) {
      EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
    }
  }
}

TYPED_TEST(ProbDistributionsCategoricalLogitGLM,
           glm_matches_categorical_logit_vars_big) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::categorical_logit_glm_lpmf;
  using stan::math::var;
  using std::vector;
  using matrix_v = typename TypeParam::matrix_v;
  using vector_v = typename TypeParam::vector_v;
  using row_vector_v = typename TypeParam::row_vector_v;
  double eps = 1e-11;
  const size_t N_instances = 89;
  const size_t N_attributes = 23;
  const size_t N_classes = 11;
  vector<int> y(N_instances);
  for (int i = 0; i < N_instances; i++) {
    y[i] = Matrix<int, Dynamic, 1>::Random(1)[0] % N_classes + 1;
  }
  MatrixXd x_double = MatrixXd::Random(N_instances, N_attributes);
  matrix_v x1 = x_double;
  matrix_v x2 = x_double;
  MatrixXd beta_double = MatrixXd::Random(N_attributes, N_classes);
  matrix_v beta1 = beta_double;
  matrix_v beta2 = beta_double;
  RowVectorXd alpha_double = RowVectorXd::Random(N_classes);
  vector_v alpha1 = alpha_double.transpose();
  vector_v alpha2 = alpha_double.transpose();

  var res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1.adj()(j, i), x2.adj()(j, i), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
  }

  stan::math::set_zero_all_adjoints();

  res1 = categorical_logit_glm_simple_lpmf<true>(y, x1, alpha1, beta1);
  res2 = categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1.adj()(j, i), x2.adj()(j, i), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
  }
}

TYPED_TEST(ProbDistributionsCategoricalLogitGLM, glm_interfaces) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::categorical_logit_glm_lpmf;
  using stan::math::var;
  using std::vector;
  using matrix_v = typename TypeParam::matrix_v;
  using vector_v = typename TypeParam::vector_v;
  using row_vector_v = typename TypeParam::row_vector_v;
  const size_t N_instances = 5;
  const size_t N_attributes = 2;
  const size_t N_classes = 3;
  int y_scal = 1;
  vector<int> y{1, 3, 1, 2, 2};
  MatrixXd x_double(N_instances, N_attributes);
  x_double << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  MatrixXd beta_double(N_attributes, N_classes);
  beta_double << 0.3, 2, 0.4, -0.1, -1.3, 1;
  VectorXd alpha_double(N_classes);
  alpha_double << 0.5, -2, 4;

  RowVectorXd x_double_row = x_double.row(0);
  matrix_v x_var = x_double;
  row_vector_v x_var_row = x_double_row;
  matrix_v beta_var = beta_double;
  vector_v alpha_var = alpha_double;

  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y, x_double, alpha_double, beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y, x_var, alpha_double, beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y, x_double, alpha_var, beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y, x_double, alpha_double, beta_var));
  EXPECT_NO_THROW(categorical_logit_glm_lpmf(y, x_double, alpha_var, beta_var));
  EXPECT_NO_THROW(categorical_logit_glm_lpmf(y, x_var, alpha_double, beta_var));
  EXPECT_NO_THROW(categorical_logit_glm_lpmf(y, x_var, alpha_var, beta_double));
  EXPECT_NO_THROW(categorical_logit_glm_lpmf(y, x_var, alpha_var, beta_var));

  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_double, alpha_double, beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_var, alpha_double, beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_double, alpha_var, beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_double, alpha_double, beta_var));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_double, alpha_var, beta_var));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_var, alpha_double, beta_var));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_var, alpha_var, beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_var, alpha_var, beta_var));

  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y, x_double_row, alpha_double, beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y, x_var_row, alpha_double, beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y, x_double_row, alpha_var, beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y, x_double_row, alpha_double, beta_var));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y, x_double_row, alpha_var, beta_var));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y, x_var_row, alpha_double, beta_var));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y, x_var_row, alpha_var, beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y, x_var_row, alpha_var, beta_var));

  EXPECT_NO_THROW(categorical_logit_glm_lpmf(y_scal, x_double_row, alpha_double,
                                             beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_var_row, alpha_double, beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_double_row, alpha_var, beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_double_row, alpha_double, beta_var));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_double_row, alpha_var, beta_var));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_var_row, alpha_double, beta_var));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_var_row, alpha_var, beta_double));
  EXPECT_NO_THROW(
      categorical_logit_glm_lpmf(y_scal, x_var_row, alpha_var, beta_var));
}

TEST(ProbDistributionsCategoricalLogitGLM, glm_errors) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::categorical_logit_glm_lpmf;
  using stan::math::var;
  using std::vector;
  const size_t N_instances = 5;
  const size_t N_attributes = 2;
  const size_t N_classes = 3;
  vector<int> y{1, 3, 1, 2, 2};
  vector<int> y_size{1, 3, 1, 2, 2, 1};
  vector<int> y_val1{1, 3, 1, 2, 4};
  vector<int> y_val2{1, 3, 1, 2, 0};
  MatrixXd x(N_instances, N_attributes);
  x << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  MatrixXd x_size1(N_instances + 1, N_attributes);
  x_size1 << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18, 1, 2;
  MatrixXd x_size2(N_instances, N_attributes + 1);
  x_size2 << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18, 27, -14, -11, 5, 18;
  MatrixXd x_val1(N_instances, N_attributes);
  x_val1 << -12, 46, -42, 24, 25, 27, -14, -11, 5, NAN;
  MatrixXd x_val2(N_instances, N_attributes);
  x_val2 << -12, 46, -42, 24, 25, 27, -14, -11, 5, INFINITY;
  MatrixXd beta(N_attributes, N_classes);
  beta << 0.3, 2, 0.4, -0.1, -1.3, 1;
  MatrixXd beta_size1(N_attributes + 1, N_classes);
  beta_size1 << 0.3, 2, 0.4, -0.1, -1.3, 1, -0.1, -1.3, 1;
  MatrixXd beta_size2(N_attributes, N_classes + 1);
  beta_size2 << 0.3, 2, 0.4, -0.1, -1.3, 1, -1.3, 1;
  MatrixXd beta_val1(N_attributes, N_classes);
  beta_val1 << 0.3, 2, 0.4, -0.1, -1.3, NAN;
  MatrixXd beta_val2(N_attributes, N_classes);
  beta_val2 << 0.3, 2, 0.4, -0.1, -1.3, INFINITY;
  VectorXd alpha(N_classes);
  alpha << 0.5, -2, 4;
  VectorXd alpha_size(N_classes + 1);
  alpha_size << 0.5, -2, 4, 3;
  VectorXd alpha_val1(N_classes);
  alpha_val1 << 0.5, -2, NAN;
  VectorXd alpha_val2(N_classes);
  alpha_val2 << 0.5, -2, INFINITY;

  EXPECT_NO_THROW(categorical_logit_glm_lpmf(y, x, alpha, beta));
  EXPECT_THROW(categorical_logit_glm_lpmf(y_size, x, alpha, beta),
               std::invalid_argument);
  EXPECT_THROW(categorical_logit_glm_lpmf(y_val1, x, alpha, beta),
               std::domain_error);
  EXPECT_THROW(categorical_logit_glm_lpmf(y_val2, x, alpha, beta),
               std::domain_error);
  EXPECT_THROW(categorical_logit_glm_lpmf(y, x_size1, alpha, beta),
               std::invalid_argument);
  EXPECT_THROW(categorical_logit_glm_lpmf(y, x_size2, alpha, beta),
               std::invalid_argument);
  EXPECT_THROW(categorical_logit_glm_lpmf(y, x_val1, alpha, beta),
               std::domain_error);
  EXPECT_THROW(categorical_logit_glm_lpmf(y, x_val2, alpha, beta),
               std::domain_error);
  EXPECT_THROW(categorical_logit_glm_lpmf(y, x, alpha_size, beta),
               std::invalid_argument);
  EXPECT_THROW(categorical_logit_glm_lpmf(y, x, alpha_val1, beta),
               std::domain_error);
  EXPECT_THROW(categorical_logit_glm_lpmf(y, x, alpha_val2, beta),
               std::domain_error);
  EXPECT_THROW(categorical_logit_glm_lpmf(y, x, alpha, beta_size1),
               std::invalid_argument);
  EXPECT_THROW(categorical_logit_glm_lpmf(y, x, alpha, beta_size2),
               std::invalid_argument);
  EXPECT_THROW(categorical_logit_glm_lpmf(y, x, alpha, beta_val1),
               std::domain_error);
  EXPECT_THROW(categorical_logit_glm_lpmf(y, x, alpha, beta_val2),
               std::domain_error);
}
