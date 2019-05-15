#include <stan/math/rev/mat.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/err/check_consistent_size.hpp>
#include <stan/math/prim/mat/prob/categorical_logit_lpmf.hpp>
#include <vector>
#include <cmath>

using namespace Eigen;
using namespace stan::math;

namespace stan {
namespace math {
template <bool propto, typename T_x, typename T_alpha, typename T_beta>
typename return_type<T_x, T_alpha, T_beta>::type
categorical_logit_glm_simple_lpmf(
    const Matrix<int, Dynamic, 1>& y, const Matrix<T_x, Dynamic, Dynamic>& x,
    const T_alpha& alpha, const Matrix<T_beta, Dynamic, Dynamic>& beta) {
  typedef typename return_type<T_x, T_beta>::type T_x_beta;
  typedef typename return_type<T_x, T_beta, T_alpha>::type T_return;

  const size_t N_instances = x.rows();

  const auto& alpha_row = as_column_vector_or_scalar(alpha).transpose();

  Eigen::Matrix<T_return, Dynamic, Dynamic> tmp
      = (x.template cast<T_x_beta>() * beta.template cast<T_x_beta>())
            .array()
            .rowwise()
        + alpha_row.array();

  T_return lpmf = 0;
  // iterate overt instances
  for (int i = 0; i < N_instances; i++) {
    lpmf += categorical_logit_lpmf<propto>(y[i], tmp.row(i).transpose().eval());
  }
  return lpmf;
}
}  // namespace math
}  // namespace stan

TEST(ProbDistributionsCategoricalLogitGLM,
     glm_matches_categorical_logit_doubles) {
  double eps = 1e-13;
  const size_t N_instances = 5;
  const size_t N_attributes = 2;
  const size_t N_classes = 3;
  Matrix<int, Dynamic, 1> y(N_instances);
  y << 1, 3, 1, 2, 2;
  MatrixXd x(N_instances, N_attributes);
  x << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  MatrixXd beta(N_attributes, N_classes);
  beta << 0.3, 2, 0.4, -0.1, -1.3, 1;
  RowVectorXd alpha(N_classes);
  alpha << 0.5, -2, 4;
  std::vector<double> alpha_vec = {0.5, -2, 4};

  EXPECT_NEAR((categorical_logit_glm_simple_lpmf<false>(y, x, alpha, beta)),
              (stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta)), eps);
  EXPECT_NEAR((categorical_logit_glm_simple_lpmf<true>(y, x, alpha, beta)),
              (stan::math::categorical_logit_glm_lpmf<true>(y, x, alpha, beta)),
              eps);

  EXPECT_NEAR((categorical_logit_glm_simple_lpmf<false>(y, x, alpha_vec, beta)),
              (stan::math::categorical_logit_glm_lpmf(y, x, alpha_vec, beta)),
              eps);
  EXPECT_NEAR(
      (categorical_logit_glm_simple_lpmf<true>(y, x, alpha_vec, beta)),
      (stan::math::categorical_logit_glm_lpmf<true>(y, x, alpha_vec, beta)),
      eps);
}

TEST(ProbDistributionsCategoricalLogitGLM, glm_matches_categorical_logit_vars) {
  double eps = 1e-13;
  const size_t N_instances = 5;
  const size_t N_attributes = 2;
  const size_t N_classes = 3;
  Matrix<int, Dynamic, 1> y(N_instances);
  y << 1, 3, 1, 2, 2;
  Matrix<var, Dynamic, Dynamic> x1(N_instances, N_attributes),
      x2(N_instances, N_attributes);
  x1 << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  x2 << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  Matrix<var, Dynamic, Dynamic> beta1(N_attributes, N_classes),
      beta2(N_attributes, N_classes);
  beta1 << 0.3, 2, 0.4, -0.1, -1.3, 1;
  beta2 << 0.3, 2, 0.4, -0.1, -1.3, 1;
  Matrix<var, 1, Dynamic> alpha1(N_classes), alpha2(N_classes);
  alpha1 << 0.5, -2, 4;
  alpha2 << 0.5, -2, 4;
  std::vector<var> alpha_vec1 = {0.5, -2, 4};
  std::vector<var> alpha_vec2 = {0.5, -2, 4};

  var res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = stan::math::categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
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
  res2 = stan::math::categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
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

  res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha_vec1, beta1);
  res2 = stan::math::categorical_logit_glm_lpmf(y, x2, alpha_vec2, beta2);
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

  res1 = categorical_logit_glm_simple_lpmf<true>(y, x1, alpha_vec1, beta1);
  res2 = stan::math::categorical_logit_glm_lpmf<true>(y, x2, alpha_vec2, beta2);
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

TEST(ProbDistributionsCategoricalLogitGLM,
     glm_matches_categorical_logit_vars_big) {
  double eps = 1e-11;
  const size_t N_instances = 89;
  const size_t N_attributes = 23;
  const size_t N_classes = 11;
  Matrix<int, Dynamic, 1> y(N_instances);
  for (int i = 0; i < N_instances; i++) {
    y[i] = Matrix<int, Dynamic, 1>::Random(1)[0] % N_classes + 1;
  }
  MatrixXd x_double = MatrixXd::Random(N_instances, N_attributes);
  Matrix<var, Dynamic, Dynamic> x1 = x_double, x2 = x_double;
  MatrixXd beta_double = MatrixXd::Random(N_attributes, N_classes);
  Matrix<var, Dynamic, Dynamic> beta1 = beta_double, beta2 = beta_double;
  RowVectorXd alpha_double = RowVectorXd::Random(N_classes);
  Matrix<var, 1, Dynamic> alpha1 = alpha_double, alpha2 = alpha_double;
  std::vector<var> alpha_vec1, alpha_vec2;
  for (int i = 0; i < N_classes; i++) {
    alpha_vec1.emplace_back(alpha_double[i]);
    alpha_vec2.emplace_back(alpha_double[i]);
  }

  var res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = stan::math::categorical_logit_glm_lpmf(y, x2, alpha2, beta2);
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
  res2 = stan::math::categorical_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
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

  res1 = categorical_logit_glm_simple_lpmf<false>(y, x1, alpha_vec1, beta1);
  res2 = stan::math::categorical_logit_glm_lpmf(y, x2, alpha_vec2, beta2);
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

  res1 = categorical_logit_glm_simple_lpmf<true>(y, x1, alpha_vec1, beta1);
  res2 = stan::math::categorical_logit_glm_lpmf<true>(y, x2, alpha_vec2, beta2);
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
  const size_t N_instances = 5;
  const size_t N_attributes = 2;
  const size_t N_classes = 3;
  Matrix<int, Dynamic, 1> y(N_instances);
  y << 1, 3, 1, 2, 2;
  MatrixXd x_double(N_instances, N_attributes);
  x_double << -12, 46, -42, 24, 25, 27, -14, -11, 5, 18;
  MatrixXd beta_double(N_attributes, N_classes);
  beta_double << 0.3, 2, 0.4, -0.1, -1.3, 1;
  RowVectorXd alpha_double(N_classes);
  alpha_double << 0.5, -2, 4;
  std::vector<double> alpha_vec_double = {0.5, -2, 4};

  Matrix<var, Dynamic, Dynamic> x_var = x_double;
  Matrix<var, Dynamic, Dynamic> beta_var = beta_double;
  Matrix<var, 1, Dynamic> alpha_var = alpha_double;
  std::vector<var> alpha_vec_var;
  for (int i = 0; i < N_classes; i++) {
    alpha_vec_var.emplace_back(alpha_vec_double[i]);
  }

  EXPECT_NO_THROW(stan::math::categorical_logit_glm_lpmf(
      y, x_double, alpha_double, beta_double));
  EXPECT_NO_THROW(stan::math::categorical_logit_glm_lpmf(
      y, x_double, alpha_vec_double, beta_double));

  EXPECT_NO_THROW(stan::math::categorical_logit_glm_lpmf(y, x_var, alpha_double,
                                                         beta_double));
  EXPECT_NO_THROW(stan::math::categorical_logit_glm_lpmf(
      y, x_var, alpha_vec_double, beta_double));
  EXPECT_NO_THROW(stan::math::categorical_logit_glm_lpmf(y, x_double, alpha_var,
                                                         beta_double));
  EXPECT_NO_THROW(stan::math::categorical_logit_glm_lpmf(
      y, x_double, alpha_vec_var, beta_double));
  EXPECT_NO_THROW(stan::math::categorical_logit_glm_lpmf(
      y, x_double, alpha_double, beta_var));
  EXPECT_NO_THROW(stan::math::categorical_logit_glm_lpmf(
      y, x_double, alpha_vec_double, beta_var));

  EXPECT_NO_THROW(
      stan::math::categorical_logit_glm_lpmf(y, x_double, alpha_var, beta_var));
  EXPECT_NO_THROW(stan::math::categorical_logit_glm_lpmf(
      y, x_double, alpha_vec_var, beta_var));
  EXPECT_NO_THROW(
      stan::math::categorical_logit_glm_lpmf(y, x_var, alpha_double, beta_var));
  EXPECT_NO_THROW(stan::math::categorical_logit_glm_lpmf(
      y, x_var, alpha_vec_double, beta_var));
  EXPECT_NO_THROW(
      stan::math::categorical_logit_glm_lpmf(y, x_var, alpha_var, beta_double));
  EXPECT_NO_THROW(stan::math::categorical_logit_glm_lpmf(
      y, x_var, alpha_vec_var, beta_double));

  EXPECT_NO_THROW(
      stan::math::categorical_logit_glm_lpmf(y, x_var, alpha_var, beta_var));
  EXPECT_NO_THROW(stan::math::categorical_logit_glm_lpmf(
      y, x_var, alpha_vec_var, beta_var));
}

TEST(ProbDistributionsCategoricalLogitGLM, glm_errors) {
  const size_t N_instances = 5;
  const size_t N_attributes = 2;
  const size_t N_classes = 3;
  Matrix<int, Dynamic, 1> y(N_instances);
  y << 1, 3, 1, 2, 2;
  Matrix<int, Dynamic, 1> y_size(N_instances + 1);
  y_size << 1, 3, 1, 2, 2, 1;
  Matrix<int, Dynamic, 1> y_val1(N_instances);
  y_val1 << 1, 3, 1, 2, 4;
  Matrix<int, Dynamic, 1> y_val2(N_instances);
  y_val2 << 1, 3, 1, 2, 0;
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
  RowVectorXd alpha(N_classes);
  alpha << 0.5, -2, 4;
  RowVectorXd alpha_size(N_classes + 1);
  alpha_size << 0.5, -2, 4, 3;
  RowVectorXd alpha_val1(N_classes);
  alpha_val1 << 0.5, -2, NAN;
  RowVectorXd alpha_val2(N_classes);
  alpha_val2 << 0.5, -2, INFINITY;

  EXPECT_NO_THROW(stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta));
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y_size, x, alpha, beta),
               std::invalid_argument);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y_val1, x, alpha, beta),
               std::domain_error);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y_val2, x, alpha, beta),
               std::domain_error);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y, x_size1, alpha, beta),
               std::invalid_argument);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y, x_size2, alpha, beta),
               std::invalid_argument);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y, x_val1, alpha, beta),
               std::domain_error);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y, x_val2, alpha, beta),
               std::domain_error);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y, x, alpha_size, beta),
               std::invalid_argument);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y, x, alpha_val1, beta),
               std::domain_error);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y, x, alpha_val2, beta),
               std::domain_error);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta_size1),
               std::invalid_argument);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta_size2),
               std::invalid_argument);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta_val1),
               std::domain_error);
  EXPECT_THROW(stan::math::categorical_logit_glm_lpmf(y, x, alpha, beta_val2),
               std::domain_error);
}