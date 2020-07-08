#include <gtest/gtest.h>
#include <stan/math/prim.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <limits>
#include <vector>

TEST(ProbDistributionsBernoulliLogitGlm, vectorized) {
  using stan::math::bernoulli_logit_glm_rng;
  //  Test scalar/vector combinations.
  boost::random::mt19937 rng;

  Eigen::MatrixXd x(2, 3);
  x << 3.5, -1.5, 0.0, 2.0, 1.0, 3.0;

  double alpha_scalar = 1.0;
  std::vector<double> alpha{1.0, 3.0};
  Eigen::VectorXd alpha_vector(2);
  alpha_vector << 1.0, 3.0;
  Eigen::RowVectorXd alpha_vector_t(2);
  alpha_vector_t = alpha_vector;

  double beta_scalar = 2.0;
  std::vector<double> beta{2.0, 4.5, -1.0};
  Eigen::VectorXd beta_vector(3);
  beta_vector << 2.0, 4.5, -1.0;
  Eigen::RowVectorXd beta_vector_t(3);
  beta_vector_t = beta_vector;

  // Can't use VectorRNGTestRig since stan::math::size(alpha) !=
  // stan::math::size(beta) in general.

  EXPECT_NO_THROW(stan::math::bernoulli_logit_glm_rng(x, alpha, beta, rng));
  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_rng(x, alpha_scalar, beta, rng));
  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_rng(x, alpha_vector, beta, rng));
  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_rng(x, alpha_vector_t, beta, rng));

  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_rng(x, alpha, beta_scalar, rng));
  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_rng(x, alpha_scalar, beta_scalar, rng));
  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_rng(x, alpha_vector, beta_scalar, rng));
  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_rng(x, alpha_vector_t, beta_scalar, rng));

  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_rng(x, alpha, beta_vector, rng));
  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_rng(x, alpha_scalar, beta_vector, rng));
  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_rng(x, alpha_vector, beta_vector, rng));
  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_rng(x, alpha_vector_t, beta_vector, rng));

  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_rng(x, alpha, beta_vector_t, rng));
  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_rng(x, alpha_scalar, beta_vector_t, rng));
  EXPECT_NO_THROW(
      stan::math::bernoulli_logit_glm_rng(x, alpha_vector, beta_vector_t, rng));
  EXPECT_NO_THROW(stan::math::bernoulli_logit_glm_rng(x, alpha_vector_t,
                                                      beta_vector_t, rng));
}

TEST(ProbDistributionsBernoulliLogitGlm, errorCheck) {
  using stan::math::bernoulli_logit_glm_rng;
  // Check errors for nonfinite and wrong sizes.
  boost::random::mt19937 rng;

  int N = 3;
  int M = 2;
  int W = 4;

  Eigen::MatrixXd x = Eigen::MatrixXd::Random(N, M);
  Eigen::VectorXd alpha = Eigen::VectorXd::Random(N, 1);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(M, 1);

  EXPECT_NO_THROW(stan::math::bernoulli_logit_glm_rng(x, alpha, beta, rng));
  Eigen::MatrixXd xw1 = Eigen::MatrixXd::Random(W, M);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_rng(xw1, alpha, beta, rng),
               std::invalid_argument);
  Eigen::MatrixXd xw2 = Eigen::MatrixXd::Random(N, W);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_rng(xw2, alpha, beta, rng),
               std::invalid_argument);
  Eigen::MatrixXd xw3 = Eigen::MatrixXd::Random(N, M) * NAN;
  EXPECT_THROW(stan::math::bernoulli_logit_glm_rng(xw3, alpha, beta, rng),
               std::domain_error);
  Eigen::VectorXd alphaw1 = Eigen::VectorXd::Random(W, 1);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_rng(x, alphaw1, beta, rng),
               std::invalid_argument);
  Eigen::VectorXd alphaw2 = Eigen::VectorXd::Random(N, 1) * NAN;
  EXPECT_THROW(stan::math::bernoulli_logit_glm_rng(x, alphaw2, beta, rng),
               std::domain_error);
  Eigen::VectorXd betaw1 = Eigen::VectorXd::Random(W, 1);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_rng(x, alpha, betaw1, rng),
               std::invalid_argument);
  Eigen::VectorXd betaw2 = Eigen::VectorXd::Random(M, 1) * NAN;
  EXPECT_THROW(stan::math::bernoulli_logit_glm_rng(x, alpha, betaw2, rng),
               std::domain_error);
}

TEST(ProbDistributionsBernoulliLogitGlm, marginalChiSquareGoodnessFitTest) {
  using stan::math::bernoulli_logit_glm_rng;
  // Check distribution of result.
  boost::random::mt19937 rng;
  Eigen::MatrixXd x(2, 2);
  x << 3.5, -1.5, 2.0, -1.2;
  std::vector<double> alpha{2.0, 1.0};
  std::vector<double> beta{2.0, 4.5};

  //  sage: x = matrix([[3.5, -1.5], [2.0, -1.2]])
  //  sage: alpha = matrix([[2.0], [1.0]])
  //  sage: beta = matrix([[2.0], [4.5]])
  //  sage: z = alpha + x * beta
  //  sage: z
  //
  //  [  2.25000000000000]
  //  [-0.399999999999999]
  //  sage: p1 = 1 / (1 + e**(-z[0][0]))
  //  sage: p2 = 1 / (1 + e**(-z[1][0]))
  //  sage: p1
  //  0.904650535100891
  //  sage: p2
  //  0.401312339887548

  double p1 = 0.904650535100891;
  double p2 = 0.401312339887548;

  int N = 10000;

  // First bin is failures, second is successes. Take N samples, take
  // their first component. Should be (1 - p1) * N failures in the
  // first bin, or thereabouts. Now take their second
  // component. Should be (1 - p2) * N failures in the first bin.
  std::vector<double> bin_boundaries{0.1, 1.1};
  std::vector<double> proportions1{(1 - p1), p1};
  std::vector<double> proportions2{(1 - p2), p2};

  std::vector<double> samples1;
  std::vector<double> samples2;
  for (int i = 0; i < N; ++i) {
    std::vector<int> sample
        = stan::math::bernoulli_logit_glm_rng(x, alpha, beta, rng);
    samples1.push_back(sample[0]);
    samples2.push_back(sample[1]);
  }

  assert_matches_bins(samples1, bin_boundaries, proportions1, 1e-6);
  assert_matches_bins(samples2, bin_boundaries, proportions2, 1e-6);
}
