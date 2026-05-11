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
