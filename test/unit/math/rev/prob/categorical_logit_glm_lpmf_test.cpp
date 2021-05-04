#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

template <bool propto, typename T_x, typename T_alpha, typename T_beta>
stan::return_type_t<T_x, T_alpha, T_beta> categorical_logit_glm_simple_lpmf(
    const std::vector<int>& y,
    const Eigen::Matrix<T_x, Eigen::Dynamic, Eigen::Dynamic>& x,
    const T_alpha& alpha,
    const Eigen::Matrix<T_beta, Eigen::Dynamic, Eigen::Dynamic>& beta) {
  using T_x_beta = stan::return_type_t<T_x, T_beta>;
  using T_return = stan::return_type_t<T_x, T_beta, T_alpha>;

  const size_t N_instances = x.rows();

  const auto& alpha_row
      = stan::math::as_column_vector_or_scalar(alpha).transpose();

  Eigen::Matrix<T_return, Eigen::Dynamic, Eigen::Dynamic> tmp
      = (x.template cast<T_x_beta>() * beta.template cast<T_x_beta>())
            .array()
            .rowwise()
        + alpha_row.array();

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

TEST(ProbDistributionsCategoricalLogitGLM, glm_matches_categorical_logit_vars) {
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
  Matrix<var, Dynamic, Dynamic> x1(N_instances, N_attributes),
      x2(N_instances, N_attributes);
  x1 << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  x2 << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  Matrix<var, Dynamic, Dynamic> beta1(N_attributes, N_classes),
      beta2(N_attributes, N_classes);
  beta1 << 0.3, 2, 0.4, -0.1, -1.3, 1;
  beta2 << 0.3, 2, 0.4, -0.1, -1.3, 1;
  Matrix<var, Dynamic, 1> alpha1(N_classes), alpha2(N_classes);
  alpha1 << 0.5, -2, 4;
  alpha2 << 0.5, -2, 4;

  var res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1(i, j).adj(), beta2(i, j).adj(), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
  }

  stan::math::set_zero_all_adjoints();

  res1 = categorical_logit_glm_simple_lpmf<true>(y, x1, alpha1, beta1);
  res2 = categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1(i, j).adj(), beta2(i, j).adj(), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
  }
}

TEST(ProbDistributionsCategoricalLogitGLM, single_instance) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::categorical_logit_glm_lpmf;
  using stan::math::var;
  using std::vector;
  double eps = 1e-13;
  const size_t N_instances = 1;
  const size_t N_attributes = 2;
  const size_t N_classes = 3;
  vector<int> y{1};
  Matrix<var, Dynamic, Dynamic> x1(N_instances, N_attributes),
      x2(N_instances, N_attributes);
  x1 << -12, 46;
  x2 << -12, 46;
  Matrix<var, Dynamic, Dynamic> beta1(N_attributes, N_classes),
      beta2(N_attributes, N_classes);
  beta1 << 0.3, 2, 0.4, -0.1, -1.3, 1;
  beta2 << 0.3, 2, 0.4, -0.1, -1.3, 1;
  Matrix<var, Dynamic, 1> alpha1(N_classes), alpha2(N_classes);
  alpha1 << 0.5, -2, 4;
  alpha2 << 0.5, -2, 4;

  var res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1(i, j).adj(), beta2(i, j).adj(), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
  }

  stan::math::set_zero_all_adjoints();

  res1 = categorical_logit_glm_simple_lpmf<true>(y, x1, alpha1, beta1);
  res2 = categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1(i, j).adj(), beta2(i, j).adj(), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
  }
}

TEST(ProbDistributionsCategoricalLogitGLM, zero_instances) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::categorical_logit_glm_lpmf;
  using stan::math::var;
  using std::vector;
  double eps = 1e-13;
  const size_t N_instances = 0;
  const size_t N_attributes = 2;
  const size_t N_classes = 3;
  std::vector<int> y{};
  Matrix<var, Dynamic, Dynamic> x1(N_instances, N_attributes),
      x2(N_instances, N_attributes);
  Matrix<var, Dynamic, Dynamic> beta1(N_attributes, N_classes),
      beta2(N_attributes, N_classes);
  beta1 << 0.3, 2, 0.4, -0.1, -1.3, 1;
  beta2 << 0.3, 2, 0.4, -0.1, -1.3, 1;
  Matrix<var, Dynamic, 1> alpha1(N_classes), alpha2(N_classes);
  alpha1 << 0.5, -2, 4;
  alpha2 << 0.5, -2, 4;

  var res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1(i, j).adj(), beta2(i, j).adj(), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
  }

  stan::math::set_zero_all_adjoints();

  res1 = categorical_logit_glm_simple_lpmf<true>(y, x1, alpha1, beta1);
  res2 = categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1(i, j).adj(), beta2(i, j).adj(), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
  }
}

TEST(ProbDistributionsCategoricalLogitGLM, single_class) {
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
  const size_t N_classes = 1;
  vector<int> y{1, 1, 1, 1, 1};
  Matrix<var, Dynamic, Dynamic> x1(N_instances, N_attributes),
      x2(N_instances, N_attributes);
  x1 << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  x2 << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  Matrix<var, Dynamic, Dynamic> beta1(N_attributes, N_classes),
      beta2(N_attributes, N_classes);
  beta1 << 0.3, 2;
  beta2 << 0.3, 2;
  Matrix<var, Dynamic, 1> alpha1(N_classes), alpha2(N_classes);
  alpha1 << 0.5;
  alpha2 << 0.5;

  var res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1(i, j).adj(), beta2(i, j).adj(), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
  }

  stan::math::set_zero_all_adjoints();

  res1 = categorical_logit_glm_simple_lpmf<true>(y, x1, alpha1, beta1);
  res2 = categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1(i, j).adj(), beta2(i, j).adj(), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
  }
}

TEST(ProbDistributionsCategoricalLogitGLM, zero_attributes) {
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
  const size_t N_attributes = 0;
  const size_t N_classes = 3;
  vector<int> y{1, 3, 1, 2, 2};
  Matrix<var, Dynamic, Dynamic> x(N_instances, N_attributes);
  Matrix<var, Dynamic, Dynamic> beta(N_attributes, N_classes);
  Matrix<var, Dynamic, 1> alpha1(N_classes), alpha2(N_classes);
  alpha1 << 0.5, -2, 4;
  alpha2 << 0.5, -2, 4;

  var res1 = categorical_logit_glm_simple_lpmf<false>(y, x, alpha1, beta);
  var res2 = categorical_logit_glm_lpmf(y, x, alpha2, beta);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
  }

  stan::math::set_zero_all_adjoints();

  res1 = categorical_logit_glm_simple_lpmf<true>(y, x, alpha1, beta);
  res2 = categorical_logit_glm_lpmf<true>(y, x, alpha2, beta);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
  }
}

TEST(ProbDistributionsCategoricalLogitGLM, x_broadcasting) {
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
  Matrix<double, 1, Dynamic> x_double(N_attributes);
  x_double << -12, 46;
  Matrix<var, 1, Dynamic> x_row = x_double;
  Matrix<var, Dynamic, Dynamic> x = x_double.replicate(N_instances, 1);
  Matrix<var, Dynamic, Dynamic> beta1(N_attributes, N_classes),
      beta2(N_attributes, N_classes);
  beta1 << 0.3, 2, 0.4, -0.1, -1.3, 1;
  beta2 << 0.3, 2, 0.4, -0.1, -1.3, 1;
  Matrix<var, Dynamic, 1> alpha1(N_classes), alpha2(N_classes);
  alpha1 << 0.5, -2, 4;
  alpha2 << 0.5, -2, 4;

  var res1 = categorical_logit_glm_lpmf(y, x_row, alpha1, beta1);
  var res2 = categorical_logit_glm_lpmf(y, x, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    double x_sum = 0;
    for (int j = 0; j < N_instances; j++) {
      x_sum += x(j, i).adj();
    }
    EXPECT_NEAR(x_row[i].adj(), x_sum, eps);
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1(i, j).adj(), beta2(i, j).adj(), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
  }

  stan::math::set_zero_all_adjoints();

  res1 = categorical_logit_glm_lpmf<true>(y, x_row, alpha1, beta1);
  res2 = categorical_logit_glm_lpmf<true>(y, x, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    double x_sum = 0;
    for (int j = 0; j < N_instances; j++) {
      x_sum += x(j, i).adj();
    }
    EXPECT_NEAR(x_row[i].adj(), x_sum, eps);
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1(i, j).adj(), beta2(i, j).adj(), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
  }
}

TEST(ProbDistributionsCategoricalLogitGLM, y_broadcasting) {
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
  Matrix<var, Dynamic, Dynamic> x1(N_instances, N_attributes),
      x2(N_instances, N_attributes);
  x1 << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  x2 << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  Matrix<var, Dynamic, Dynamic> beta1(N_attributes, N_classes),
      beta2(N_attributes, N_classes);
  beta1 << 0.3, 2, 0.4, -0.1, -1.3, 1;
  beta2 << 0.3, 2, 0.4, -0.1, -1.3, 1;
  Matrix<var, Dynamic, 1> alpha1(N_classes), alpha2(N_classes);
  alpha1 << 0.5, -2, 4;
  alpha2 << 0.5, -2, 4;

  for (int y_scal = 1; y_scal <= N_classes; y_scal++) {
    vector<int> y(N_instances, y_scal);
    var res1 = categorical_logit_glm_lpmf(y_scal, x1, alpha1, beta1);
    var res2 = categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
    (res1 + res2).grad();
    EXPECT_NEAR(res1.val(), res2.val(), eps);
    for (int i = 0; i < N_attributes; i++) {
      for (int j = 0; j < N_instances; j++) {
        EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
      }
      for (int j = 0; j < N_classes; j++) {
        EXPECT_NEAR(beta1(i, j).adj(), beta2(i, j).adj(), eps);
      }
    }
    for (int i = 0; i < N_classes; i++) {
      EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
    }

    stan::math::set_zero_all_adjoints();

    res1 = categorical_logit_glm_lpmf<true>(y_scal, x1, alpha1, beta1);
    res2 = categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
    (res1 + res2).grad();
    EXPECT_NEAR(res1.val(), res2.val(), eps);
    for (int i = 0; i < N_attributes; i++) {
      for (int j = 0; j < N_instances; j++) {
        EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
      }
      for (int j = 0; j < N_classes; j++) {
        EXPECT_NEAR(beta1(i, j).adj(), beta2(i, j).adj(), eps);
      }
    }
    for (int i = 0; i < N_classes; i++) {
      EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
    }
  }
}

TEST(ProbDistributionsCategoricalLogitGLM,
     glm_matches_categorical_logit_vars_big) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::categorical_logit_glm_lpmf;
  using stan::math::var;
  using std::vector;
  double eps = 1e-11;
  const size_t N_instances = 89;
  const size_t N_attributes = 23;
  const size_t N_classes = 11;
  vector<int> y(N_instances);
  for (int i = 0; i < N_instances; i++) {
    y[i] = Matrix<int, Dynamic, 1>::Random(1)[0] % N_classes + 1;
  }
  MatrixXd x_double = MatrixXd::Random(N_instances, N_attributes);
  Matrix<var, Dynamic, Dynamic> x1 = x_double, x2 = x_double;
  MatrixXd beta_double = MatrixXd::Random(N_attributes, N_classes);
  Matrix<var, Dynamic, Dynamic> beta1 = beta_double, beta2 = beta_double;
  RowVectorXd alpha_double = RowVectorXd::Random(N_classes);
  Matrix<var, Dynamic, 1> alpha1 = alpha_double, alpha2 = alpha_double;

  var res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1(i, j).adj(), beta2(i, j).adj(), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
  }

  stan::math::set_zero_all_adjoints();

  res1 = categorical_logit_glm_simple_lpmf<true>(y, x1, alpha1, beta1);
  res2 = categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < N_attributes; i++) {
    for (int j = 0; j < N_instances; j++) {
      EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
    }
    for (int j = 0; j < N_classes; j++) {
      EXPECT_NEAR(beta1(i, j).adj(), beta2(i, j).adj(), eps);
    }
  }
  for (int i = 0; i < N_classes; i++) {
    EXPECT_NEAR(alpha1[i].adj(), alpha2[i].adj(), eps);
  }
}

TEST(ProbDistributionsCategoricalLogitGLM, glm_interfaces) {
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
  int y_scal = 1;
  vector<int> y{1, 3, 1, 2, 2};
  MatrixXd x_double(N_instances, N_attributes);
  x_double << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  MatrixXd beta_double(N_attributes, N_classes);
  beta_double << 0.3, 2, 0.4, -0.1, -1.3, 1;
  VectorXd alpha_double(N_classes);
  alpha_double << 0.5, -2, 4;

  RowVectorXd x_double_row = x_double.row(0);
  Matrix<var, Dynamic, Dynamic> x_var = x_double;
  Matrix<var, 1, Dynamic> x_var_row = x_double_row;
  Matrix<var, Dynamic, Dynamic> beta_var = beta_double;
  Matrix<var, Dynamic, 1> alpha_var = alpha_double;

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
