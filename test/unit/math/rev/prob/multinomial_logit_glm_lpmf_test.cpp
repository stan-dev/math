#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <vector>

// Reference: sum multinomial_logit_lpmf over all N instances
template <bool propto, typename T_x, typename T_alpha, typename T_beta>
inline stan::return_type_t<T_x, T_alpha, T_beta>
multinomial_logit_glm_simple_lpmf(const std::vector<std::vector<int>>& y,
                                  const T_x& x, const T_alpha& alpha,
                                  const T_beta& beta) {
  using T_return = stan::return_type_t<T_x, T_alpha, T_beta>;
  const int N = x.rows();
  // alpha is 1xK; as_column_vector_or_scalar+transpose gives a row-vector
  // expression whose decay type lets rep_matrix deduce its Ret parameter.
  const auto& alpha_row
      = stan::math::as_column_vector_or_scalar(alpha).transpose();
  auto lin = stan::math::to_ref(
      stan::math::multiply(x, beta)
      + stan::math::rep_matrix<std::decay_t<decltype(alpha_row)>>(alpha_row,
                                                                  N));
  T_return lpmf = 0;
  for (int n = 0; n < N; ++n)
    lpmf += stan::math::multinomial_logit_lpmf<propto>(
        y[n], lin.row(n).transpose().eval());
  return lpmf;
}

TEST_F(AgradRev, MultinomialLogitGLM_matches_simple_doubles) {
  using stan::math::multinomial_logit_glm_lpmf;
  const size_t N = 5, M = 2, K = 3;
  std::vector<std::vector<int>> y{
      {2, 1, 3}, {0, 4, 1}, {3, 0, 2}, {1, 2, 0}, {0, 1, 4}};
  Eigen::MatrixXd x(N, M);
  x << -1.2, 0.5, 0.3, -0.7, 1.0, 0.2, -0.4, 0.8, 0.6, -1.1;
  Eigen::MatrixXd beta(M, K);
  beta << 0.3, -0.2, 0.1, -0.1, 0.4, -0.3;
  Eigen::RowVectorXd alpha(K);
  alpha << 0.2, -0.1, 0.5;
  const double eps = 1e-13;

  EXPECT_NEAR(multinomial_logit_glm_simple_lpmf<false>(y, x, alpha, beta),
              multinomial_logit_glm_lpmf(y, x, alpha, beta), eps);
  EXPECT_NEAR(multinomial_logit_glm_simple_lpmf<true>(y, x, alpha, beta),
              multinomial_logit_glm_lpmf<true>(y, x, alpha, beta), eps);
}

template <class T>
class ProbDistributionsMultinomialLogitGLM
    : public stan::math::test::VarMatrixTypedTests<T> {};

TYPED_TEST_SUITE(ProbDistributionsMultinomialLogitGLM,
                 stan::math::test::VarMatImpls);

TYPED_TEST(ProbDistributionsMultinomialLogitGLM, glm_matches_simple_vars) {
  using stan::math::multinomial_logit_glm_lpmf;
  using stan::math::var;
  using matrix_v = typename TypeParam::matrix_v;
  using row_vector_v = typename TypeParam::row_vector_v;
  const size_t N = 5, M = 2, K = 3;
  const double eps = 1e-13;
  std::vector<std::vector<int>> y{
      {2, 1, 3}, {0, 4, 1}, {3, 0, 2}, {1, 2, 0}, {0, 1, 4}};
  Eigen::MatrixXd x_val(N, M);
  x_val << -1.2, 0.5, 0.3, -0.7, 1.0, 0.2, -0.4, 0.8, 0.6, -1.1;
  Eigen::MatrixXd beta_val(M, K);
  beta_val << 0.3, -0.2, 0.1, -0.1, 0.4, -0.3;
  Eigen::RowVectorXd alpha_val(K);
  alpha_val << 0.2, -0.1, 0.5;

  // Reference always uses Matrix<var> so multinomial_logit_lpmf can extract
  // rows.
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> x1 = x_val,
                                                     beta1 = beta_val;
  Eigen::Matrix<var, 1, Eigen::Dynamic> alpha1 = alpha_val;
  matrix_v x2 = x_val;
  matrix_v beta2 = beta_val;
  row_vector_v alpha2 = alpha_val;

  var res1 = multinomial_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = multinomial_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j)
      EXPECT_NEAR(x1.adj()(j, i), x2.adj()(j, i), eps);
    for (size_t j = 0; j < K; ++j)
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
  }
  for (size_t i = 0; i < K; ++i)
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);

  stan::math::set_zero_all_adjoints();

  res1 = multinomial_logit_glm_simple_lpmf<true>(y, x1, alpha1, beta1);
  res2 = multinomial_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j)
      EXPECT_NEAR(x1.adj()(j, i), x2.adj()(j, i), eps);
    for (size_t j = 0; j < K; ++j)
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
  }
  for (size_t i = 0; i < K; ++i)
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
}

TYPED_TEST(ProbDistributionsMultinomialLogitGLM, glm_matches_simple_big) {
  using stan::math::multinomial_logit_glm_lpmf;
  using stan::math::var;
  using matrix_v = typename TypeParam::matrix_v;
  using row_vector_v = typename TypeParam::row_vector_v;
  const size_t N = 37, M = 11, K = 7;
  const double eps = 1e-10;

  std::vector<std::vector<int>> y(N, std::vector<int>(K));
  for (size_t n = 0; n < N; ++n)
    for (size_t k = 0; k < K; ++k)
      y[n][k] = (n + k) % 5;

  Eigen::MatrixXd x_val = Eigen::MatrixXd::Random(N, M);
  Eigen::MatrixXd beta_val = Eigen::MatrixXd::Random(M, K);
  Eigen::RowVectorXd alpha_val = Eigen::RowVectorXd::Random(K);

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> x1 = x_val,
                                                     beta1 = beta_val;
  Eigen::Matrix<var, 1, Eigen::Dynamic> alpha1 = alpha_val;
  matrix_v x2 = x_val;
  matrix_v beta2 = beta_val;
  row_vector_v alpha2 = alpha_val;

  var res1 = multinomial_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = multinomial_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j)
      EXPECT_NEAR(x1.adj()(j, i), x2.adj()(j, i), eps);
    for (size_t j = 0; j < K; ++j)
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
  }
  for (size_t i = 0; i < K; ++i)
    EXPECT_NEAR(alpha1.adj()[i], alpha2.adj()[i], eps);
}

TYPED_TEST(ProbDistributionsMultinomialLogitGLM, matrix_alpha_grads) {
  using stan::math::multinomial_logit_glm_lpmf;
  using stan::math::var;
  using matrix_v = typename TypeParam::matrix_v;
  const size_t N = 4, M = 2, K = 3;
  const double eps = 1e-13;
  std::vector<std::vector<int>> y{{2, 1, 3}, {0, 4, 1}, {3, 0, 2}, {1, 2, 0}};
  Eigen::MatrixXd x_val(N, M);
  x_val << -1.2, 0.5, 0.3, -0.7, 1.0, 0.2, -0.4, 0.8;
  Eigen::MatrixXd alpha_val(N, K);
  alpha_val << 0.2, -0.1, 0.5, -0.3, 0.4, 0.1, 0.1, 0.0, -0.2, 0.5, -0.3, 0.2;
  Eigen::MatrixXd beta_val(M, K);
  beta_val << 0.3, -0.2, 0.1, -0.1, 0.4, -0.3;
  // full x + matrix alpha: exercises the T_alpha_rows != 1 gradient path
  // Reference uses Matrix<var> so multinomial_logit_lpmf can extract Eigen
  // rows.
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> x1 = x_val,
                                                     alpha1 = alpha_val,
                                                     beta1 = beta_val;
  matrix_v x2 = x_val;
  matrix_v alpha2 = alpha_val;
  matrix_v beta2 = beta_val;

  var res1 = 0.0;
  {
    auto lin = stan::math::to_ref(stan::math::multiply(x1, beta1) + alpha1);
    for (size_t n = 0; n < N; ++n)
      res1 += stan::math::multinomial_logit_lpmf(y[n],
                                                 lin.row(n).transpose().eval());
  }
  var res2 = multinomial_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j)
      EXPECT_NEAR(x1.adj()(j, i), x2.adj()(j, i), eps);
    for (size_t j = 0; j < K; ++j)
      EXPECT_NEAR(beta1.adj()(i, j), beta2.adj()(i, j), eps);
  }
  for (size_t n = 0; n < N; ++n)
    for (size_t k = 0; k < K; ++k)
      EXPECT_NEAR(alpha1.adj()(n, k), alpha2.adj()(n, k), eps);
}

TYPED_TEST(ProbDistributionsMultinomialLogitGLM,
           extreme_broadcast_alpha_value_and_gradients) {
  using stan::math::multinomial_logit_glm_lpmf;
  using stan::math::var;
  using matrix_v = typename TypeParam::matrix_v;
  using row_vector_v = typename TypeParam::row_vector_v;
  std::vector<std::vector<int>> y{{0, 1}, {1, 0}};
  Eigen::MatrixXd x_val(2, 1);
  x_val << 0, 0;
  Eigen::RowVectorXd alpha_val(2);
  alpha_val << 0, -1000;
  Eigen::MatrixXd beta_val(1, 2);
  beta_val << 0, 0;

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> x1 = x_val,
                                                     beta1 = beta_val;
  Eigen::Matrix<var, 1, Eigen::Dynamic> alpha1 = alpha_val;
  matrix_v x2 = x_val;
  row_vector_v alpha2 = alpha_val;
  matrix_v beta2 = beta_val;

  var res1 = multinomial_logit_glm_simple_lpmf<false>(y, x1, alpha1, beta1);
  var res2 = multinomial_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), 1e-12);
  EXPECT_NEAR(-1000, res2.val(), 1e-12);
  EXPECT_TRUE(x1.adj().isApprox(x2.adj(), 1e-12));
  EXPECT_TRUE(alpha1.adj().isApprox(alpha2.adj(), 1e-12));
  EXPECT_TRUE(beta1.adj().isApprox(beta2.adj(), 1e-12));

  stan::math::set_zero_all_adjoints();
  res1 = multinomial_logit_glm_simple_lpmf<true>(y, x1, alpha1, beta1);
  res2 = multinomial_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), 1e-12);
  EXPECT_NEAR(-1000, res2.val(), 1e-12);
  EXPECT_TRUE(x1.adj().isApprox(x2.adj(), 1e-12));
  EXPECT_TRUE(alpha1.adj().isApprox(alpha2.adj(), 1e-12));
  EXPECT_TRUE(beta1.adj().isApprox(beta2.adj(), 1e-12));
}

TYPED_TEST(ProbDistributionsMultinomialLogitGLM,
           extreme_matrix_alpha_value_and_gradients) {
  using stan::math::multinomial_logit_glm_lpmf;
  using stan::math::var;
  using matrix_v = typename TypeParam::matrix_v;
  std::vector<std::vector<int>> y{{0, 1, 0}, {1, 0, 2}};
  Eigen::MatrixXd x_val(2, 1);
  x_val << 0, 0;
  Eigen::MatrixXd alpha_val(2, 3);
  alpha_val << 1000, 0, -1000, 0, -1000, -2000;
  Eigen::MatrixXd beta_val = Eigen::MatrixXd::Zero(1, 3);

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> x1 = x_val,
                                                     alpha1 = alpha_val,
                                                     beta1 = beta_val;
  matrix_v x2 = x_val;
  matrix_v alpha2 = alpha_val;
  matrix_v beta2 = beta_val;

  var res1 = 0;
  auto lin = stan::math::to_ref(stan::math::multiply(x1, beta1) + alpha1);
  for (int n = 0; n < x_val.rows(); ++n) {
    res1 += stan::math::multinomial_logit_lpmf(y[n],
                                               lin.row(n).transpose().eval());
  }
  var res2 = multinomial_logit_glm_lpmf(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), 1e-10);
  EXPECT_TRUE(x1.adj().isApprox(x2.adj(), 1e-10));
  EXPECT_TRUE(alpha1.adj().isApprox(alpha2.adj(), 1e-10));
  EXPECT_TRUE(beta1.adj().isApprox(beta2.adj(), 1e-10));

  stan::math::set_zero_all_adjoints();
  res1 = 0;
  for (int n = 0; n < x_val.rows(); ++n) {
    res1 += stan::math::multinomial_logit_lpmf<true>(
        y[n], lin.row(n).transpose().eval());
  }
  res2 = multinomial_logit_glm_lpmf<true>(y, x2, alpha2, beta2);
  (res1 + res2).grad();
  EXPECT_NEAR(res1.val(), res2.val(), 1e-10);
  EXPECT_TRUE(x1.adj().isApprox(x2.adj(), 1e-10));
  EXPECT_TRUE(alpha1.adj().isApprox(alpha2.adj(), 1e-10));
  EXPECT_TRUE(beta1.adj().isApprox(beta2.adj(), 1e-10));
}

TYPED_TEST(ProbDistributionsMultinomialLogitGLM,
           zero_total_all_neg_inf_and_zero_classes_have_zero_gradients) {
  using stan::math::multinomial_logit_glm_lpmf;
  using stan::math::var;
  using matrix_v = typename TypeParam::matrix_v;
  using row_vector_v = typename TypeParam::row_vector_v;

  Eigen::MatrixXd x_val(1, 1);
  x_val << 0.4;
  Eigen::MatrixXd alpha_val
      = Eigen::MatrixXd::Constant(1, 2, -stan::math::INFTY);
  Eigen::MatrixXd beta_val(1, 2);
  beta_val << 0.2, -0.3;
  matrix_v x = x_val;
  matrix_v alpha = alpha_val;
  matrix_v beta = beta_val;
  var res = multinomial_logit_glm_lpmf(std::vector<std::vector<int>>{{0, 0}}, x,
                                       alpha, beta);
  res.grad();
  EXPECT_EQ(0, res.val());
  EXPECT_TRUE(x.adj().isZero());
  EXPECT_TRUE(alpha.adj().isZero());
  EXPECT_TRUE(beta.adj().isZero());

  stan::math::set_zero_all_adjoints();
  Eigen::MatrixXd x0_val(2, 1);
  x0_val << 0.5, -0.2;
  matrix_v x0 = x0_val;
  row_vector_v alpha0 = Eigen::RowVectorXd(0);
  matrix_v beta0 = Eigen::MatrixXd(1, 0);
  var res0 = multinomial_logit_glm_lpmf(std::vector<std::vector<int>>(2), x0,
                                        alpha0, beta0);
  res0.grad();
  EXPECT_EQ(0, res0.val());
  EXPECT_TRUE(x0.adj().isZero());
}

TEST_F(AgradRev, MultinomialLogitGLM_eigen_expression_gradients) {
  using stan::math::multinomial_logit_glm_lpmf;
  using stan::math::var;
  std::vector<std::vector<int>> y{{2, 1, 3}, {0, 4, 1}};

  Eigen::MatrixXd x_storage_val(4, 2);
  x_storage_val << 0.1, 0.2, 1.0, -0.5, -0.3, 0.7, 0.4, -0.2;
  Eigen::MatrixXd x_add_val(2, 2);
  x_add_val << 0.2, -0.1, 0.3, 0.4;
  Eigen::VectorXd alpha_col_val(3), alpha_add_val(3);
  alpha_col_val << 0.1, -0.3, 0.2;
  alpha_add_val << 0.2, 0.1, -0.2;
  Eigen::MatrixXd beta_left_val(2, 2), beta_right_val(2, 3);
  beta_left_val << 1.0, 0.2, -0.1, 0.8;
  beta_right_val << 0.3, -0.2, 0.1, -0.1, 0.4, -0.3;

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> x_storage1 = x_storage_val;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> x_add1 = x_add_val;
  Eigen::Matrix<var, Eigen::Dynamic, 1> alpha_col1 = alpha_col_val;
  Eigen::Matrix<var, Eigen::Dynamic, 1> alpha_add1 = alpha_add_val;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> beta_left1 = beta_left_val;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> beta_right1
      = beta_right_val;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> x_storage2 = x_storage_val;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> x_add2 = x_add_val;
  Eigen::Matrix<var, Eigen::Dynamic, 1> alpha_col2 = alpha_col_val;
  Eigen::Matrix<var, Eigen::Dynamic, 1> alpha_add2 = alpha_add_val;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> beta_left2 = beta_left_val;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> beta_right2
      = beta_right_val;

  var expected = multinomial_logit_glm_lpmf(
      y, (x_storage1.middleRows(1, 2) + x_add1).eval(),
      (alpha_col1 + alpha_add1).transpose().eval(),
      (beta_left1 * beta_right1).eval());
  var actual = multinomial_logit_glm_lpmf(
      y, x_storage2.middleRows(1, 2) + x_add2,
      (alpha_col2 + alpha_add2).transpose(), beta_left2 * beta_right2);
  (expected + actual).grad();

  EXPECT_NEAR(expected.val(), actual.val(), 1e-12);
  EXPECT_TRUE(x_storage1.adj().isApprox(x_storage2.adj(), 1e-12));
  EXPECT_TRUE(x_add1.adj().isApprox(x_add2.adj(), 1e-12));
  EXPECT_TRUE(alpha_col1.adj().isApprox(alpha_col2.adj(), 1e-12));
  EXPECT_TRUE(alpha_add1.adj().isApprox(alpha_add2.adj(), 1e-12));
  EXPECT_TRUE(beta_left1.adj().isApprox(beta_left2.adj(), 1e-12));
  EXPECT_TRUE(beta_right1.adj().isApprox(beta_right2.adj(), 1e-12));
}

TYPED_TEST(ProbDistributionsMultinomialLogitGLM, interfaces) {
  using stan::math::multinomial_logit_glm_lpmf;
  using matrix_v = typename TypeParam::matrix_v;
  using row_vector_v = typename TypeParam::row_vector_v;
  const size_t N = 3, M = 2, K = 3;
  std::vector<std::vector<int>> y{{1, 2, 0}, {0, 1, 3}, {2, 0, 1}};
  Eigen::MatrixXd x_d(N, M);
  x_d << 1.0, -0.5, 0.3, 0.7, -0.2, 0.4;
  Eigen::MatrixXd beta_d(M, K);
  beta_d << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
  Eigen::RowVectorXd alpha_d(K);
  alpha_d << 0.1, 0.2, 0.3;

  matrix_v x_v = x_d;
  matrix_v beta_v = beta_d;
  row_vector_v alpha_v = alpha_d;

  EXPECT_NO_THROW(multinomial_logit_glm_lpmf(y, x_d, alpha_d, beta_d));
  EXPECT_NO_THROW(multinomial_logit_glm_lpmf(y, x_v, alpha_d, beta_d));
  EXPECT_NO_THROW(multinomial_logit_glm_lpmf(y, x_d, alpha_v, beta_d));
  EXPECT_NO_THROW(multinomial_logit_glm_lpmf(y, x_d, alpha_d, beta_v));
  EXPECT_NO_THROW(multinomial_logit_glm_lpmf(y, x_v, alpha_v, beta_d));
  EXPECT_NO_THROW(multinomial_logit_glm_lpmf(y, x_v, alpha_d, beta_v));
  EXPECT_NO_THROW(multinomial_logit_glm_lpmf(y, x_d, alpha_v, beta_v));
  EXPECT_NO_THROW(multinomial_logit_glm_lpmf(y, x_v, alpha_v, beta_v));
}
