#include <gtest/gtest.h>
#include <stan/math/prim.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <exception>
#include <limits>
#include <vector>

class ProbDistributionsGaussianDLMInputsRng : public ::testing::Test {
 protected:
  virtual void SetUp() {
    FF = Eigen::MatrixXd::Random(2, 3);
    GG = Eigen::MatrixXd::Random(2, 2);
    V = Eigen::MatrixXd::Identity(3, 3);
    V_vec = Eigen::VectorXd::Constant(3, 1.0);
    W = Eigen::MatrixXd::Identity(2, 2);
    y = Eigen::MatrixXd::Random(3, 5);
    m0 = Eigen::VectorXd::Random(2);
    C0 = Eigen::MatrixXd::Identity(2, 2);
    T = 5;
  }

  Eigen::MatrixXd FF;
  Eigen::MatrixXd GG;
  Eigen::MatrixXd V;
  Eigen::VectorXd V_vec;
  Eigen::MatrixXd W;
  Eigen::MatrixXd y;
  Eigen::VectorXd m0;
  Eigen::MatrixXd C0;
  int T;
};

TEST_F(ProbDistributionsGaussianDLMInputsRng, PoliciesF) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::gaussian_dlm_obs_rng;
  boost::random::mt19937 rng;
  MatrixXd FF_sz1 = MatrixXd::Random(4, 3);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_sz1, GG, V, W, m0, C0, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_sz1, GG, V_vec, W, m0, C0, T, rng),
               std::invalid_argument);
  MatrixXd FF_sz2 = MatrixXd::Random(2, 4);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_sz2, GG, V, W, m0, C0, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_sz2, GG, V_vec, W, m0, C0, T, rng),
               std::invalid_argument);
  // finite and NaN
  MatrixXd FF_inf = FF;
  FF_inf(0, 0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_inf, GG, V, W, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_inf, GG, V_vec, W, m0, C0, T, rng),
               std::domain_error);
  MatrixXd FF_nan = FF;
  FF_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_nan, GG, V, W, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_nan, GG, V_vec, W, m0, C0, T, rng),
               std::domain_error);
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, PoliciesG) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::gaussian_dlm_obs_rng;
  boost::random::mt19937 rng;
  // size
  MatrixXd GG_sz1 = MatrixXd::Random(3, 3);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_sz1, V, W, m0, C0, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_sz1, V_vec, W, m0, C0, T, rng),
               std::invalid_argument);
  MatrixXd GG_sz2 = MatrixXd::Random(2, 3);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_sz2, V, W, m0, C0, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_sz2, V_vec, W, m0, C0, T, rng),
               std::invalid_argument);
  // finite and NaN
  MatrixXd GG_inf = GG;
  GG_inf(0, 0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_inf, V, W, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_inf, V_vec, W, m0, C0, T, rng),
               std::domain_error);
  MatrixXd GG_nan = GG;
  GG_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_nan, V, W, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_nan, V_vec, W, m0, C0, T, rng),
               std::domain_error);
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, PoliciesW) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::gaussian_dlm_obs_rng;
  boost::random::mt19937 rng;
  // Not symmetric
  MatrixXd W_asym = W;
  W_asym(0, 1) = 1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_asym, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_asym, m0, C0, T, rng),
               std::domain_error);
  // negative
  MatrixXd W_neg = W;
  W_neg(0, 0) = -1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_neg, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_neg, m0, C0, T, rng),
               std::domain_error);
  // finite and NaN
  MatrixXd W_infinite = W;
  W_infinite(0, 0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_infinite, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_infinite, m0, C0, T, rng),
               std::domain_error);
  MatrixXd W_nan = W;
  W_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_nan, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_nan, m0, C0, T, rng),
               std::domain_error);
  // wrong size
  MatrixXd W_sz = MatrixXd::Identity(4, 4);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_sz, m0, C0, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_sz, m0, C0, T, rng),
               std::invalid_argument);
  // not square
  MatrixXd W_notsq = MatrixXd::Identity(2, 3);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_notsq, m0, C0, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_notsq, m0, C0, T, rng),
               std::invalid_argument);
  // positive semi-definite is okay
  MatrixXd W_psd = MatrixXd::Zero(2, 2);
  W_psd(0, 0) = 1.0;
  EXPECT_NO_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_psd, m0, C0, T, rng));
  EXPECT_NO_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_psd, m0, C0, T, rng));
  // Non positive semi-definite is not.
  MatrixXd W_bad = MatrixXd::Zero(2, 2);
  W_bad << 1, 2, 2, 1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_bad, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_bad, m0, C0, T, rng),
               std::domain_error);
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, PoliciesVMatrix) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::gaussian_dlm_obs_rng;
  boost::random::mt19937 rng;
  // Not symmetric
  MatrixXd V_asym = V;
  V_asym(0, 2) = 1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_asym, W, m0, C0, T, rng),
               std::domain_error);
  // negative
  MatrixXd V_neg = V;
  V_neg(0, 2) = -1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_neg, W, m0, C0, T, rng),
               std::domain_error);
  // finite and NaN
  MatrixXd V_infinite = V;
  V_infinite(0, 0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_infinite, W, m0, C0, T, rng),
               std::domain_error);
  MatrixXd V_nan = V;
  V_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_nan, W, m0, C0, T, rng),
               std::domain_error);
  // wrong size
  MatrixXd V2 = MatrixXd::Identity(2, 2);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V2, W, m0, C0, T, rng),
               std::invalid_argument);
  // not square
  MatrixXd V3 = MatrixXd::Identity(2, 3);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V3, W, m0, C0, T, rng),
               std::invalid_argument);
  // positive semi-definite is okay
  MatrixXd V_psd = MatrixXd::Zero(3, 3);
  V_psd(0, 0) = 1.0;
  EXPECT_NO_THROW(gaussian_dlm_obs_rng(FF, GG, V_psd, W, m0, C0, T, rng));
  // Non positive semi-definite is not.
  MatrixXd V_bad = MatrixXd::Zero(3, 3);
  V_bad << 1, 2, 2, 2, 1, 2, 2, 2, 1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_bad, W, m0, C0, T, rng),
               std::domain_error);
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, PoliciesVVector) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::gaussian_dlm_obs_rng;
  boost::random::mt19937 rng;
  // negative
  MatrixXd V_neg = V_vec;
  V_neg(0) = -1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_neg, W, m0, C0, T, rng),
               std::invalid_argument);
  // finite and NaN
  MatrixXd V_infinite = V_vec;
  V_infinite(0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_infinite, W, m0, C0, T, rng),
               std::domain_error);
  MatrixXd V_nan = V_vec;
  V_nan(0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_nan, W, m0, C0, T, rng),
               std::domain_error);
  // wrong size
  VectorXd V_badsz = VectorXd::Constant(2, 1.0);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_badsz, W, m0, C0, T, rng),
               std::invalid_argument);
  // positive semi-definite is okay
  VectorXd V_psd = VectorXd::Zero(3);
  V_psd(0) = 1.0;
  EXPECT_NO_THROW(gaussian_dlm_obs_rng(FF, GG, V_psd, W, m0, C0, T, rng));
  // Non positive semi-definite is not, but with vector V, V is
  // diagonal so non psd means negative entries, and we already tested
  // that those give an error.
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, Policiesm0) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::gaussian_dlm_obs_rng;
  boost::random::mt19937 rng;
  // size
  VectorXd m0_sz = VectorXd::Zero(4, 1);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0_sz, C0, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W, m0_sz, C0, T, rng),
               std::invalid_argument);
  // finite and NaN
  VectorXd m0_inf = m0;
  m0_inf(0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0_inf, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W, m0_inf, C0, T, rng),
               std::domain_error);
  VectorXd m0_nan = m0;
  m0_nan(0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0_nan, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W, m0_nan, C0, T, rng),
               std::domain_error);
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, PoliciesC0) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::gaussian_dlm_obs_rng;
  boost::random::mt19937 rng;
  // size
  MatrixXd C0_sz = MatrixXd::Identity(3, 3);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0, C0_sz, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W, m0, C0_sz, T, rng),
               std::invalid_argument);
  // negative
  MatrixXd C0_neg = C0;
  C0_neg(0, 0) = -1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0, C0_neg, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W, m0, C0_neg, T, rng),
               std::domain_error);
  // asymmetric
  MatrixXd C0_asym = C0;
  C0_asym(0, 1) = 1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0, C0_asym, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W, m0, C0_neg, T, rng),
               std::domain_error);
  // not square
  MatrixXd C0_notsq = MatrixXd::Identity(3, 2);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0, C0_notsq, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W, m0, C0_notsq, T, rng),
               std::invalid_argument);
  // positive semi-definite okay.
  MatrixXd C0_psd = MatrixXd::Zero(2, 2);
  C0_psd(0, 0) = 1.0;
  EXPECT_NO_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0, C0_psd, T, rng));
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, PoliciesT) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::gaussian_dlm_obs_rng;
  boost::random::mt19937 rng;
  // Must be positive.
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0, C0, 0, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0, C0, -1, rng),
               std::domain_error);
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, chiSquaredGoodnessOfFit) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::gaussian_dlm_obs_rng;
  boost::random::mt19937 rng;

  // With identity state transition and initial mean 0 and identity
  // covariance matrices:
  // state at 0 = initial noise
  // state at t+1 = state at t + noise[t+1]
  // observed at T =
  //     initial noise
  //     + noise[1] + ... + noise[T]
  //     + observation noise.
  //   = sum of (T+2) variance 1 gaussians
  //   = variance T+2 gaussian noise
  // (Mean 0.)

  int n = 3;
  int r = 1;

  VectorXd m0_(n, 1);
  m0_ << 0, 0, 0;

  MatrixXd C0_ = MatrixXd::Identity(n, n);

  MatrixXd FF_(n, r);
  FF_ << 1, 0, 0;

  MatrixXd V_ = MatrixXd::Identity(r, r);

  MatrixXd GG_ = MatrixXd::Identity(n, n);

  MatrixXd W_ = MatrixXd::Identity(n, n);

  // Small enough that variance T_ + 2 and variance T_ + 1 can be
  // correctly distinguished. Just to keep the off-by-one errors under
  // control.
  int T_ = 3;

  int N = 10000;
  int K = boost::math::round(2 * std::pow(N, 0.4));
  boost::math::normal_distribution<> dist(0, sqrt(T_ + 2));

  std::vector<double> quantiles;
  for (int i = 1; i < K; i++)
    quantiles.push_back(quantile(dist, i * std::pow(K, -1.0)));
  quantiles.push_back(std::numeric_limits<double>::max());

  std::vector<double> samples;
  Eigen::MatrixXd y_;
  for (int n = 0; n < N; ++n) {
    y_ = gaussian_dlm_obs_rng(FF_, GG_, V_, W_, m0_, C0_, T_, rng);
    ASSERT_TRUE(y_.rows() == r && y_.cols() == T_);
    samples.push_back(y_(0, T_ - 1));
  }

  assert_matches_quantiles(samples, quantiles, 1e-6);
}
