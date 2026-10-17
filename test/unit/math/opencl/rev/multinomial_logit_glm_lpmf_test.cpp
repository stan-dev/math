#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <stan/math/opencl/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <climits>
#include <limits>
#include <utility>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::matrix_cl;
using stan::math::var;
using std::vector;

namespace {
template <typename T_y>
using device_outcome_call_t = decltype(stan::math::multinomial_logit_glm_lpmf(
    std::declval<T_y>(), std::declval<matrix_cl<double>>(),
    std::declval<matrix_cl<double>>(), std::declval<matrix_cl<double>>()));

static_assert(stan::is_detected<matrix_cl<int>, device_outcome_call_t>::value,
              "matrix_cl<int> outcomes must be accepted");
static_assert(
    !stan::is_detected<matrix_cl<double>, device_outcome_call_t>::value,
    "non-int device outcome buffers must be rejected");
}  // namespace

TEST(ProbDistributionsMultinomialLogitGLM, error_checking) {
  int N = 3, M = 2, K = 3;

  vector<vector<int>> y{{1, 2, 0}, {0, 3, 1}, {2, 0, 2}};
  vector<vector<int>> y_bad_n{{1, 2, 0}, {0, 3, 1}};
  vector<vector<int>> y_bad_k{{1, 2}, {0, 3}, {2, 0}};
  vector<vector<int>> y_neg{{-1, 2, 0}, {0, 3, 1}, {2, 0, 2}};

  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << 1.0, -0.5, 0.3, 0.7, -0.2, 1.1;
  Matrix<double, Dynamic, Dynamic> x_bad_n(N + 1, M);
  Matrix<double, Dynamic, Dynamic> x_bad_m(N, M + 1);
  Matrix<double, Dynamic, Dynamic> x_inf = x;
  x_inf(0, 0) = stan::math::INFTY;

  Matrix<double, Dynamic, Dynamic> beta(M, K);
  beta << 0.3, -0.2, 0.1, -0.1, 0.4, -0.3;
  Matrix<double, Dynamic, Dynamic> beta_bad_m(M + 1, K);
  Matrix<double, Dynamic, Dynamic> beta_bad_k(M, K + 1);
  Matrix<double, Dynamic, Dynamic> beta_inf = beta;
  beta_inf(0, 0) = stan::math::INFTY;

  Matrix<double, 1, Dynamic> alpha(1, K);
  alpha << 0.1, -0.3, 0.2;
  Matrix<double, 1, Dynamic> alpha_bad_k(1, K + 1);
  Matrix<double, 1, Dynamic> alpha_inf = alpha;
  alpha_inf(0) = stan::math::INFTY;

  Matrix<double, Dynamic, Dynamic> alpha_mat(N, K);
  alpha_mat << 0.1, -0.3, 0.2, -0.2, 0.1, 0.3, 0.0, 0.2, -0.1;
  Matrix<double, Dynamic, Dynamic> alpha_mat_bad_n(N + 1, K);
  Matrix<double, Dynamic, Dynamic> alpha_mat_bad_k(N, K + 1);

  matrix_cl<double> x_cl(x), x_bad_n_cl(x_bad_n), x_bad_m_cl(x_bad_m),
      x_inf_cl(x_inf);
  matrix_cl<double> beta_cl(beta), beta_bad_m_cl(beta_bad_m),
      beta_bad_k_cl(beta_bad_k), beta_inf_cl(beta_inf);
  matrix_cl<double> alpha_cl(alpha), alpha_bad_k_cl(alpha_bad_k),
      alpha_inf_cl(alpha_inf);
  matrix_cl<double> alpha_mat_cl(alpha_mat),
      alpha_mat_bad_n_cl(alpha_mat_bad_n), alpha_mat_bad_k_cl(alpha_mat_bad_k);

  EXPECT_NO_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x_cl, alpha_cl, beta_cl));
  EXPECT_NO_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x_cl, alpha_mat_cl, beta_cl));

  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y_bad_n, x_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y_bad_k, x_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y_neg, x_cl, alpha_cl, beta_cl),
      std::domain_error);

  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x_bad_n_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x_bad_m_cl, alpha_cl, beta_cl),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x_cl, alpha_bad_k_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(stan::math::multinomial_logit_glm_lpmf(
                   y, x_cl, alpha_mat_bad_n_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multinomial_logit_glm_lpmf(
                   y, x_cl, alpha_mat_bad_k_cl, beta_cl),
               std::invalid_argument);

  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x_cl, alpha_cl, beta_bad_m_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x_cl, alpha_cl, beta_bad_k_cl),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x_inf_cl, alpha_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x_cl, alpha_inf_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x_cl, alpha_cl, beta_inf_cl),
      std::domain_error);
}

TEST(ProbDistributionsMultinomialLogitGLM,
     opencl_matches_cpu_small_broadcast_alpha) {
  int N = 3, M = 2, K = 3;
  vector<vector<int>> y{{1, 2, 0}, {0, 3, 1}, {2, 0, 2}};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << 1.0, -0.5, 0.3, 0.7, -0.2, 1.1;
  Matrix<double, 1, Dynamic> alpha(1, K);
  alpha << 0.1, -0.3, 0.2;
  Matrix<double, Dynamic, Dynamic> beta(M, K);
  beta << 0.3, -0.2, 0.1, -0.1, 0.4, -0.3;

  auto f = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf(y, x_, a_, b_);
  };
  auto f_propto = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf<true>(y, x_, a_, b_);
  };
  stan::math::test::compare_cpu_opencl_prim_rev(f, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(f_propto, x, alpha, beta);
}

TEST(ProbDistributionsMultinomialLogitGLM,
     opencl_matches_cpu_small_matrix_alpha) {
  int N = 3, M = 2, K = 3;
  vector<vector<int>> y{{1, 2, 0}, {0, 3, 1}, {2, 0, 2}};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << 1.0, -0.5, 0.3, 0.7, -0.2, 1.1;
  Matrix<double, Dynamic, Dynamic> alpha(N, K);
  alpha << 0.1, -0.3, 0.2, -0.2, 0.1, 0.3, 0.0, 0.2, -0.1;
  Matrix<double, Dynamic, Dynamic> beta(M, K);
  beta << 0.3, -0.2, 0.1, -0.1, 0.4, -0.3;

  auto f = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf(y, x_, a_, b_);
  };
  auto f_propto = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf<true>(y, x_, a_, b_);
  };
  stan::math::test::compare_cpu_opencl_prim_rev(f, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(f_propto, x, alpha, beta);
}

TEST(ProbDistributionsMultinomialLogitGLM, opencl_matches_cpu_zero_instances) {
  int M = 2, K = 3;
  vector<vector<int>> y{};
  Matrix<double, Dynamic, Dynamic> x(0, M);
  Matrix<double, 1, Dynamic> alpha(1, K);
  alpha << 0.1, -0.3, 0.2;
  Matrix<double, Dynamic, Dynamic> beta(M, K);
  beta << 0.3, -0.2, 0.1, -0.1, 0.4, -0.3;

  auto f = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf(y, x_, a_, b_);
  };
  auto f_propto = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf<true>(y, x_, a_, b_);
  };
  stan::math::test::compare_cpu_opencl_prim_rev(f, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(f_propto, x, alpha, beta);
}

TEST(ProbDistributionsMultinomialLogitGLM, opencl_matches_cpu_zero_attributes) {
  int N = 3, K = 3;
  vector<vector<int>> y{{1, 2, 0}, {0, 3, 1}, {2, 0, 2}};
  Matrix<double, Dynamic, Dynamic> x(N, 0);
  Matrix<double, 1, Dynamic> alpha(1, K);
  alpha << 0.1, -0.3, 0.2;
  Matrix<double, Dynamic, Dynamic> beta(0, K);

  auto f = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf(y, x_, a_, b_);
  };
  auto f_propto = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf<true>(y, x_, a_, b_);
  };
  stan::math::test::compare_cpu_opencl_prim_rev(f, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(f_propto, x, alpha, beta);
}

TEST(ProbDistributionsMultinomialLogitGLM, opencl_matches_cpu_single_instance) {
  int N = 1, M = 3, K = 4;
  vector<vector<int>> y{{2, 0, 3, 1}};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << 0.5, -1.2, 0.8;
  Matrix<double, 1, Dynamic> alpha(1, K);
  alpha << 0.1, -0.3, 0.5, -0.2;
  Matrix<double, Dynamic, Dynamic> beta(M, K);
  beta << 0.2, -0.1, 0.4, 0.0, -0.3, 0.5, -0.2, 0.1, 0.1, 0.0, 0.3, -0.4;

  auto f = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf(y, x_, a_, b_);
  };
  stan::math::test::compare_cpu_opencl_prim_rev(f, x, alpha, beta);
}

TEST(ProbDistributionsMultinomialLogitGLM, opencl_matches_cpu_all_zeros_y) {
  int N = 2, M = 2, K = 3;
  vector<vector<int>> y{{0, 0, 0}, {0, 0, 0}};
  Matrix<double, Dynamic, Dynamic> x
      = Matrix<double, Dynamic, Dynamic>::Random(N, M);
  Matrix<double, 1, Dynamic> alpha = Matrix<double, 1, Dynamic>::Random(1, K);
  Matrix<double, Dynamic, Dynamic> beta
      = Matrix<double, Dynamic, Dynamic>::Random(M, K);

  auto f = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf(y, x_, a_, b_);
  };
  stan::math::test::compare_cpu_opencl_prim_rev(f, x, alpha, beta);
}

TEST(ProbDistributionsMultinomialLogitGLM, opencl_matches_cpu_large) {
  int N = 153, M = 17, K = 11;
  vector<vector<int>> y(N, vector<int>(K));
  for (int n = 0; n < N; ++n)
    for (int k = 0; k < K; ++k)
      y[n][k] = (n * K + k) % 5;

  Matrix<double, Dynamic, Dynamic> x
      = Matrix<double, Dynamic, Dynamic>::Random(N, M);
  Matrix<double, 1, Dynamic> alpha = Matrix<double, 1, Dynamic>::Random(1, K);
  Matrix<double, Dynamic, Dynamic> beta
      = Matrix<double, Dynamic, Dynamic>::Random(M, K);

  auto f = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf(y, x_, a_, b_);
  };
  auto f_propto = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf<true>(y, x_, a_, b_);
  };
  stan::math::test::compare_cpu_opencl_prim_rev(f, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(f_propto, x, alpha, beta);
}

TEST(ProbDistributionsMultinomialLogitGLM,
     opencl_matches_cpu_large_matrix_alpha) {
  int N = 153, M = 17, K = 11;
  vector<vector<int>> y(N, vector<int>(K));
  for (int n = 0; n < N; ++n)
    for (int k = 0; k < K; ++k)
      y[n][k] = (n * K + k) % 5;

  Matrix<double, Dynamic, Dynamic> x
      = Matrix<double, Dynamic, Dynamic>::Random(N, M);
  Matrix<double, Dynamic, Dynamic> alpha
      = Matrix<double, Dynamic, Dynamic>::Random(N, K);
  Matrix<double, Dynamic, Dynamic> beta
      = Matrix<double, Dynamic, Dynamic>::Random(M, K);

  auto f = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf(y, x_, a_, b_);
  };
  auto f_propto = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf<true>(y, x_, a_, b_);
  };
  stan::math::test::compare_cpu_opencl_prim_rev(f, x, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(f_propto, x, alpha, beta);
}

TEST(ProbDistributionsMultinomialLogitGLM, opencl_neg_inf_alpha) {
  // alpha[n,k]=-inf forces softmax probability to 0; y[n,k]=0 for those
  // classes. Result must be finite and match the CPU prim result.
  int N = 2, M = 2, K = 3;
  vector<vector<int>> y{{2, 1, 0}, {0, 3, 2}};

  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << 1.0, 0.5, 0.3, -0.7;

  Matrix<double, Dynamic, Dynamic> beta(M, K);
  beta << 0.3, -0.2, 0.1, -0.1, 0.4, -0.3;

  Matrix<double, Dynamic, Dynamic> alpha(N, K);
  alpha << 0.2, -0.1, -stan::math::INFTY, -stan::math::INFTY, 0.4, 0.1;

  const double logp_cpu
      = stan::math::multinomial_logit_glm_lpmf(y, x, alpha, beta);
  ASSERT_TRUE(std::isfinite(logp_cpu));

  matrix_cl<double> x_cl(x), alpha_cl(alpha), beta_cl(beta);
  const double logp_cl
      = stan::math::multinomial_logit_glm_lpmf(y, x_cl, alpha_cl, beta_cl);

  EXPECT_TRUE(std::isfinite(logp_cl));
  EXPECT_FLOAT_EQ(logp_cpu, logp_cl);
}

// ---------------------------------------------------------------------------
// Signature 2: outcome counts supplied as an N x K matrix_cl<int> already on
// the device. The CPU reference uses the std::vector<std::vector<int>>
// signature (signature 1); the OpenCL side uses the matrix_cl<int> signature.
// ---------------------------------------------------------------------------
namespace {
matrix_cl<int> y_to_matrix_cl(const vector<vector<int>>& y, int N, int K) {
  Eigen::MatrixXi y_mat(N, K);
  for (int n = 0; n < N; ++n)
    for (int k = 0; k < K; ++k)
      y_mat(n, k) = y[n][k];
  return matrix_cl<int>(y_mat);
}
}  // namespace

TEST(ProbDistributionsMultinomialLogitGLM, opencl_matrix_cl_y_broadcast_alpha) {
  int N = 3, M = 2, K = 3;
  vector<vector<int>> y{{1, 2, 0}, {0, 3, 1}, {2, 0, 2}};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << 1.0, -0.5, 0.3, 0.7, -0.2, 1.1;
  Matrix<double, 1, Dynamic> alpha(1, K);
  alpha << 0.1, -0.3, 0.2;
  Matrix<double, Dynamic, Dynamic> beta(M, K);
  beta << 0.3, -0.2, 0.1, -0.1, 0.4, -0.3;

  matrix_cl<int> y_cl = y_to_matrix_cl(y, N, K);
  auto f_cpu = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf(y, x_, a_, b_);
  };
  auto f_cl = [&y_cl](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf(y_cl, x_, a_, b_);
  };
  auto f_cpu_propto = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf<true>(y, x_, a_, b_);
  };
  auto f_cl_propto = [&y_cl](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf<true>(y_cl, x_, a_, b_);
  };
  stan::math::test::compare_cpu_opencl_prim_rev_separate(f_cpu, f_cl, x, alpha,
                                                         beta);
  stan::math::test::compare_cpu_opencl_prim_rev_separate(
      f_cpu_propto, f_cl_propto, x, alpha, beta);
}

TEST(ProbDistributionsMultinomialLogitGLM, opencl_matrix_cl_y_matrix_alpha) {
  int N = 3, M = 2, K = 3;
  vector<vector<int>> y{{1, 2, 0}, {0, 3, 1}, {2, 0, 2}};
  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << 1.0, -0.5, 0.3, 0.7, -0.2, 1.1;
  Matrix<double, Dynamic, Dynamic> alpha(N, K);
  alpha << 0.1, -0.3, 0.2, -0.2, 0.1, 0.3, 0.0, 0.2, -0.1;
  Matrix<double, Dynamic, Dynamic> beta(M, K);
  beta << 0.3, -0.2, 0.1, -0.1, 0.4, -0.3;

  matrix_cl<int> y_cl = y_to_matrix_cl(y, N, K);
  auto f_cpu = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf(y, x_, a_, b_);
  };
  auto f_cl = [&y_cl](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf(y_cl, x_, a_, b_);
  };
  auto f_cpu_propto = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf<true>(y, x_, a_, b_);
  };
  auto f_cl_propto = [&y_cl](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf<true>(y_cl, x_, a_, b_);
  };
  stan::math::test::compare_cpu_opencl_prim_rev_separate(f_cpu, f_cl, x, alpha,
                                                         beta);
  stan::math::test::compare_cpu_opencl_prim_rev_separate(
      f_cpu_propto, f_cl_propto, x, alpha, beta);
}

TEST(ProbDistributionsMultinomialLogitGLM, opencl_matrix_cl_y_large) {
  int N = 153, M = 17, K = 11;
  vector<vector<int>> y(N, vector<int>(K));
  for (int n = 0; n < N; ++n)
    for (int k = 0; k < K; ++k)
      y[n][k] = (n * K + k) % 5;

  Matrix<double, Dynamic, Dynamic> x
      = Matrix<double, Dynamic, Dynamic>::Random(N, M);
  Matrix<double, 1, Dynamic> alpha = Matrix<double, 1, Dynamic>::Random(1, K);
  Matrix<double, Dynamic, Dynamic> beta
      = Matrix<double, Dynamic, Dynamic>::Random(M, K);

  matrix_cl<int> y_cl = y_to_matrix_cl(y, N, K);
  auto f_cpu = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf(y, x_, a_, b_);
  };
  auto f_cl = [&y_cl](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf(y_cl, x_, a_, b_);
  };
  stan::math::test::compare_cpu_opencl_prim_rev_separate(f_cpu, f_cl, x, alpha,
                                                         beta);
}

TEST(ProbDistributionsMultinomialLogitGLM, error_checking_matrix_cl_y) {
  int N = 3, M = 2, K = 3;
  vector<vector<int>> y{{1, 2, 0}, {0, 3, 1}, {2, 0, 2}};

  Matrix<double, Dynamic, Dynamic> x(N, M);
  x << 1.0, -0.5, 0.3, 0.7, -0.2, 1.1;
  Matrix<double, 1, Dynamic> alpha(1, K);
  alpha << 0.1, -0.3, 0.2;
  Matrix<double, Dynamic, Dynamic> beta(M, K);
  beta << 0.3, -0.2, 0.1, -0.1, 0.4, -0.3;

  matrix_cl<double> x_cl(x), alpha_cl(alpha), beta_cl(beta);
  matrix_cl<int> y_cl = y_to_matrix_cl(y, N, K);
  matrix_cl<int> y_bad_n_cl(Eigen::MatrixXi::Zero(N + 1, K));
  matrix_cl<int> y_bad_k_cl(Eigen::MatrixXi::Zero(N, K + 1));

  EXPECT_NO_THROW(
      stan::math::multinomial_logit_glm_lpmf(y_cl, x_cl, alpha_cl, beta_cl));
  // y rows must match the number of instances (rows of x).
  EXPECT_THROW(stan::math::multinomial_logit_glm_lpmf(y_bad_n_cl, x_cl,
                                                      alpha_cl, beta_cl),
               std::invalid_argument);
  // y cols must match the number of classes (cols of beta).
  EXPECT_THROW(stan::math::multinomial_logit_glm_lpmf(y_bad_k_cl, x_cl,
                                                      alpha_cl, beta_cl),
               std::invalid_argument);
}

TEST(ProbDistributionsMultinomialLogitGLM,
     empty_host_outcome_still_checks_dimensions) {
  vector<vector<int>> y;
  Matrix<double, Dynamic, Dynamic> x(1, 1);
  x << 0;
  Matrix<double, 1, Dynamic> alpha(1, 2);
  alpha << 0, 0;
  Matrix<double, Dynamic, Dynamic> beta(1, 2);
  beta << 0, 0;
  matrix_cl<double> x_cl(x), alpha_cl(alpha), beta_cl(beta);

  EXPECT_THROW(
      stan::math::multinomial_logit_glm_lpmf(y, x_cl, alpha_cl, beta_cl),
      std::invalid_argument);
}

TEST(ProbDistributionsMultinomialLogitGLM,
     extreme_logits_match_cpu_for_host_and_device_outcomes) {
  const int N = 2;
  const int M = 1;
  const int K = 3;
  vector<vector<int>> y{{0, 1, 0}, {2, 0, 3}};
  Matrix<double, Dynamic, Dynamic> x
      = Matrix<double, Dynamic, Dynamic>::Zero(N, M);
  Matrix<double, Dynamic, Dynamic> beta
      = Matrix<double, Dynamic, Dynamic>::Zero(M, K);
  Matrix<double, 1, Dynamic> alpha_row(1, K);
  alpha_row << 1000, 0, -1000;
  Matrix<double, Dynamic, Dynamic> alpha_mat(N, K);
  alpha_mat << 1000, 0, -1000, 0, -1000, -2000;

  matrix_cl<double> x_cl(x), beta_cl(beta), alpha_row_cl(alpha_row),
      alpha_mat_cl(alpha_mat);
  matrix_cl<int> y_cl = y_to_matrix_cl(y, N, K);

  const double expected_row
      = stan::math::multinomial_logit_glm_lpmf(y, x, alpha_row, beta);
  const double expected_mat
      = stan::math::multinomial_logit_glm_lpmf(y, x, alpha_mat, beta);
  EXPECT_NEAR(
      expected_row,
      stan::math::multinomial_logit_glm_lpmf(y, x_cl, alpha_row_cl, beta_cl),
      1e-10);
  EXPECT_NEAR(
      expected_row,
      stan::math::multinomial_logit_glm_lpmf(y_cl, x_cl, alpha_row_cl, beta_cl),
      1e-10);
  EXPECT_NEAR(
      expected_mat,
      stan::math::multinomial_logit_glm_lpmf(y, x_cl, alpha_mat_cl, beta_cl),
      1e-10);
  EXPECT_NEAR(
      expected_mat,
      stan::math::multinomial_logit_glm_lpmf(y_cl, x_cl, alpha_mat_cl, beta_cl),
      1e-10);

  Matrix<double, 1, Dynamic> alpha_shifted = alpha_row.array() + 1000000.0;
  matrix_cl<double> alpha_shifted_cl(alpha_shifted);
  EXPECT_NEAR(expected_row,
              stan::math::multinomial_logit_glm_lpmf(y, x_cl, alpha_shifted_cl,
                                                     beta_cl),
              1e-10);

  auto f_propto = [&y](const auto& x_, const auto& a_, const auto& b_) {
    return stan::math::multinomial_logit_glm_lpmf<true>(y, x_, a_, b_);
  };
  stan::math::test::compare_cpu_opencl_prim_rev(f_propto, x, alpha_row, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(f_propto, x, alpha_mat, beta);
}

TEST(ProbDistributionsMultinomialLogitGLM, neg_inf_and_zero_class_boundaries) {
  Matrix<double, Dynamic, Dynamic> x(1, 1);
  x << 0;
  Matrix<double, Dynamic, Dynamic> beta(1, 2);
  beta << 0, 0;
  Matrix<double, 1, Dynamic> alpha(1, 2);
  alpha << 0, -stan::math::INFTY;
  matrix_cl<double> x_cl(x), beta_cl(beta), alpha_cl(alpha);

  vector<vector<int>> y_positive{{0, 1}};
  matrix_cl<int> y_positive_cl = y_to_matrix_cl(y_positive, 1, 2);
  EXPECT_EQ(-stan::math::INFTY, stan::math::multinomial_logit_glm_lpmf(
                                    y_positive, x_cl, alpha_cl, beta_cl));
  EXPECT_EQ(-stan::math::INFTY, stan::math::multinomial_logit_glm_lpmf(
                                    y_positive_cl, x_cl, alpha_cl, beta_cl));

  alpha.setConstant(-stan::math::INFTY);
  alpha_cl = matrix_cl<double>(alpha);
  vector<vector<int>> y_zero{{0, 0}};
  matrix_cl<int> y_zero_cl = y_to_matrix_cl(y_zero, 1, 2);
  EXPECT_EQ(0, stan::math::multinomial_logit_glm_lpmf(y_zero, x_cl, alpha_cl,
                                                      beta_cl));
  EXPECT_EQ(0, stan::math::multinomial_logit_glm_lpmf(y_zero_cl, x_cl, alpha_cl,
                                                      beta_cl));

  stan::math::var_value<matrix_cl<double>> x_var(x_cl);
  stan::math::var_value<matrix_cl<double>> alpha_var(alpha_cl);
  stan::math::var_value<matrix_cl<double>> beta_var(beta_cl);
  var res = stan::math::multinomial_logit_glm_lpmf(y_zero, x_var, alpha_var,
                                                   beta_var);
  res.grad();
  EXPECT_EQ(0, res.val());
  EXPECT_TRUE(stan::math::from_matrix_cl(x_var.adj()).isZero());
  EXPECT_TRUE(stan::math::from_matrix_cl(alpha_var.adj()).isZero());
  EXPECT_TRUE(stan::math::from_matrix_cl(beta_var.adj()).isZero());

  Matrix<double, Dynamic, Dynamic> x0(2, 1);
  x0 << 0.5, -0.2;
  Matrix<double, 1, Dynamic> alpha0(1, 0);
  Matrix<double, Dynamic, Dynamic> beta0(1, 0);
  matrix_cl<double> x0_cl(x0), alpha0_cl(alpha0), beta0_cl(beta0);
  matrix_cl<int> y0_cl(2, 0);
  vector<vector<int>> y0(2);
  EXPECT_EQ(0, stan::math::multinomial_logit_glm_lpmf(y0, x0_cl, alpha0_cl,
                                                      beta0_cl));
  EXPECT_EQ(0, stan::math::multinomial_logit_glm_lpmf(y0_cl, x0_cl, alpha0_cl,
                                                      beta0_cl));
  Matrix<double, Dynamic, Dynamic> alpha0_mat(2, 0);
  matrix_cl<double> alpha0_mat_cl(alpha0_mat);
  EXPECT_EQ(0, stan::math::multinomial_logit_glm_lpmf(y0, x0_cl, alpha0_mat_cl,
                                                      beta0_cl));
  x0(0, 0) = stan::math::INFTY;
  matrix_cl<double> x0_inf_cl(x0);
  EXPECT_THROW(stan::math::multinomial_logit_glm_lpmf(y0, x0_inf_cl, alpha0_cl,
                                                      beta0_cl),
               std::domain_error);
}

TEST(ProbDistributionsMultinomialLogitGLM, count_total_above_int_max) {
  vector<vector<int>> y{{INT_MAX, INT_MAX - 1}};
  Matrix<double, Dynamic, Dynamic> x(1, 1);
  x << 0;
  Matrix<double, 1, Dynamic> alpha(1, 2);
  alpha << 0, 0;
  Matrix<double, Dynamic, Dynamic> beta(1, 2);
  beta << 0, 0;
  matrix_cl<double> x_cl(x), beta_cl(beta);
  matrix_cl<int> y_cl = y_to_matrix_cl(y, 1, 2);

  Eigen::Matrix<var, 1, Dynamic> alpha_cpu = alpha;
  stan::math::var_value<matrix_cl<double>> alpha_host{matrix_cl<double>(alpha)};
  stan::math::var_value<matrix_cl<double>> alpha_device{
      matrix_cl<double>(alpha)};
  var expected
      = stan::math::multinomial_logit_glm_lpmf<true>(y, x, alpha_cpu, beta);
  var actual_host = stan::math::multinomial_logit_glm_lpmf<true>(
      y, x_cl, alpha_host, beta_cl);
  var actual_device = stan::math::multinomial_logit_glm_lpmf<true>(
      y_cl, x_cl, alpha_device, beta_cl);
  (expected + actual_host + actual_device).grad();

  EXPECT_TRUE(std::isfinite(expected.val()));
  EXPECT_NEAR(expected.val(), actual_host.val(), 1e-4);
  EXPECT_NEAR(expected.val(), actual_device.val(), 1e-4);
  const Eigen::RowVectorXd host_adj
      = stan::math::from_matrix_cl(alpha_host.adj());
  const Eigen::RowVectorXd device_adj
      = stan::math::from_matrix_cl(alpha_device.adj());
  EXPECT_TRUE(alpha_cpu.adj().isApprox(host_adj, 1e-8));
  EXPECT_TRUE(alpha_cpu.adj().isApprox(device_adj, 1e-8));
}

TEST(ProbDistributionsMultinomialLogitGLM,
     integer_alpha_uses_double_custom_kernel_buffer) {
  vector<vector<int>> y{{1, 2}};
  Matrix<double, Dynamic, Dynamic> x(1, 1);
  x << 0.5;
  Matrix<int, 1, Dynamic> alpha(1, 2);
  alpha << 1, -2;
  Matrix<double, Dynamic, Dynamic> beta(1, 2);
  beta << 0.3, -0.1;
  matrix_cl<double> x_cl(x), beta_cl(beta);
  matrix_cl<int> alpha_cl(alpha);

  const double expected
      = stan::math::multinomial_logit_glm_lpmf(y, x, alpha, beta);
  EXPECT_NEAR(
      expected,
      stan::math::multinomial_logit_glm_lpmf(y, x_cl, alpha_cl, beta_cl),
      1e-10);
}

TEST(ProbDistributionsMultinomialLogitGLM,
     all_constant_propto_validates_device_inputs) {
  vector<vector<int>> y{{1, 0}};
  Matrix<double, Dynamic, Dynamic> x(1, 1);
  x << 0;
  Matrix<double, 1, Dynamic> alpha(1, 2);
  alpha << 0, 0;
  Matrix<double, Dynamic, Dynamic> beta(1, 2);
  beta << 0, 0;

  Matrix<double, 1, Dynamic> alpha_neg_inf(1, 2);
  alpha_neg_inf << 0, -stan::math::INFTY;
  matrix_cl<double> alpha_neg_inf_cl(alpha_neg_inf);
  matrix_cl<double> x_cl(x), alpha_cl(alpha), beta_cl(beta);
  EXPECT_EQ(0, stan::math::multinomial_logit_glm_lpmf<true>(
                   y, x_cl, alpha_neg_inf_cl, beta_cl));

  Matrix<double, Dynamic, Dynamic> x_nan = x;
  x_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();
  matrix_cl<double> x_nan_cl(x_nan);
  EXPECT_THROW(stan::math::multinomial_logit_glm_lpmf<true>(y, x_nan_cl,
                                                            alpha_cl, beta_cl),
               std::domain_error);

  beta(0, 0) = std::numeric_limits<double>::quiet_NaN();
  matrix_cl<double> beta_nan_cl(beta);
  EXPECT_THROW(stan::math::multinomial_logit_glm_lpmf<true>(y, x_cl, alpha_cl,
                                                            beta_nan_cl),
               std::domain_error);

  alpha(0, 0) = stan::math::INFTY;
  matrix_cl<double> alpha_inf_cl(alpha);
  EXPECT_THROW(stan::math::multinomial_logit_glm_lpmf<true>(
                   y, x_cl, alpha_inf_cl, beta_cl),
               std::domain_error);

  alpha(0, 0) = std::numeric_limits<double>::quiet_NaN();
  matrix_cl<double> alpha_nan_cl(alpha);
  EXPECT_THROW(stan::math::multinomial_logit_glm_lpmf<true>(
                   y, x_cl, alpha_nan_cl, beta_cl),
               std::domain_error);

  matrix_cl<int> y_negative_cl(Eigen::MatrixXi::Constant(1, 2, -1));
  EXPECT_THROW(stan::math::multinomial_logit_glm_lpmf<true>(y_negative_cl, x_cl,
                                                            alpha_cl, beta_cl),
               std::domain_error);
}

#endif
