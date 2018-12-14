#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <test/unit/math/prim/scal/prob/util.hpp>
#include <exception>
#include <limits>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;
using Eigen::MatrixXd;
using stan::math::gaussian_dlm_obs_rng;

class ProbDistributionsGaussianDLMInputsRng : public ::testing::Test {
 protected:
  virtual void SetUp() {
    FF = MatrixXd::Random(2, 3);
    GG = MatrixXd::Random(2, 2);
    V = MatrixXd::Identity(3, 3);
    V_vec = Matrix<double, Dynamic, 1>::Constant(3, 1.0);
    W = MatrixXd::Identity(2, 2);
    y = MatrixXd::Random(3, 5);
    m0 = Matrix<double, Dynamic, 1>::Random(2);
    C0 = MatrixXd::Identity(2, 2);
    T = 5;
  }

  Matrix<double, Dynamic, Dynamic> FF;
  Matrix<double, Dynamic, Dynamic> GG;
  Matrix<double, Dynamic, Dynamic> V;
  Matrix<double, Dynamic, 1> V_vec;
  Matrix<double, Dynamic, Dynamic> W;
  Matrix<double, Dynamic, Dynamic> y;
  Matrix<double, Dynamic, 1> m0;
  Matrix<double, Dynamic, Dynamic> C0;
  unsigned int T;
};

TEST_F(ProbDistributionsGaussianDLMInputsRng, PoliciesF) {
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, Dynamic> FF_sz1 = MatrixXd::Random(4, 3);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_sz1, GG, V, W, m0, C0, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_sz1, GG, V_vec, W, m0, C0, T, rng),
               std::invalid_argument);
  Matrix<double, Dynamic, Dynamic> FF_sz2 = MatrixXd::Random(2, 4);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_sz2, GG, V, W, m0, C0, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_sz2, GG, V_vec, W, m0, C0, T, rng),
               std::invalid_argument);
  // finite and NaN
  Matrix<double, Dynamic, Dynamic> FF_inf = FF;
  FF_inf(0, 0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_inf, GG, V, W, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_inf, GG, V_vec, W, m0, C0, T, rng),
               std::domain_error);
  Matrix<double, Dynamic, Dynamic> FF_nan = FF;
  FF_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_nan, GG, V, W, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF_nan, GG, V_vec, W, m0, C0, T, rng),
               std::domain_error);
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, PoliciesG) {
  boost::random::mt19937 rng;
  // size
  Matrix<double, Dynamic, Dynamic> GG_sz1 = MatrixXd::Random(3, 3);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_sz1, V, W, m0, C0, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_sz1, V_vec, W, m0, C0, T, rng),
               std::invalid_argument);
  Matrix<double, Dynamic, Dynamic> GG_sz2 = MatrixXd::Random(2, 3);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_sz2, V, W, m0, C0, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_sz2, V_vec, W, m0, C0, T, rng),
               std::invalid_argument);
  // finite and NaN
  Matrix<double, Dynamic, Dynamic> GG_inf = GG;
  GG_inf(0, 0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_inf, V, W, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_inf, V_vec, W, m0, C0, T, rng),
               std::domain_error);
  Matrix<double, Dynamic, Dynamic> GG_nan = GG;
  GG_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_nan, V, W, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG_nan, V_vec, W, m0, C0, T, rng),
               std::domain_error);
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, PoliciesW) {
  boost::random::mt19937 rng;
  // Not symmetric
  Matrix<double, Dynamic, Dynamic> W_asym = W;
  W_asym(0, 1) = 1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_asym, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_asym, m0, C0, T, rng),
               std::domain_error);
  // negative
  Matrix<double, Dynamic, Dynamic> W_neg = W;
  W_neg(0, 0) = -1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_neg, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_neg, m0, C0, T, rng),
               std::domain_error);
  // finite and NaN
  Matrix<double, Dynamic, Dynamic> W_infinite = W;
  W_infinite(0, 0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_infinite, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_infinite, m0, C0, T, rng),
               std::domain_error);
  Matrix<double, Dynamic, Dynamic> W_nan = W;
  W_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_nan, m0, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_nan, m0, C0, T, rng),
               std::domain_error);
  // wrong size
  Matrix<double, Dynamic, Dynamic> W_sz = MatrixXd::Identity(4, 4);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_sz, m0, C0, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_sz, m0, C0, T, rng),
               std::invalid_argument);
  // not square
  Matrix<double, Dynamic, Dynamic> W_notsq = MatrixXd::Identity(2, 3);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_notsq, m0, C0, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_notsq, m0, C0, T, rng),
               std::invalid_argument);
  // positive semi-definite is okay
  Matrix<double, Dynamic, Dynamic> W_psd = MatrixXd::Zero(2, 2);
  W_psd(0, 0) = 1.0;
  EXPECT_NO_THROW(gaussian_dlm_obs_rng(FF, GG, V, W_psd, m0, C0, T, rng));
  EXPECT_NO_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W_psd, m0, C0, T, rng));
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, PoliciesVMatrix) {
  boost::random::mt19937 rng;
  // Not symmetric
  Matrix<double, Dynamic, Dynamic> V_asym = V;
  V_asym(0, 2) = 1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_asym, W, m0, C0, T, rng),
               std::domain_error);
  // negative
  Matrix<double, Dynamic, Dynamic> V_neg = V;
  V_neg(0, 2) = -1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_neg, W, m0, C0, T, rng),
               std::domain_error);
  // finite and NaN
  Matrix<double, Dynamic, Dynamic> V_infinite = V;
  V_infinite(0, 0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_infinite, W, m0, C0, T, rng),
               std::domain_error);
  Matrix<double, Dynamic, Dynamic> V_nan = V;
  V_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_nan, W, m0, C0, T, rng),
               std::domain_error);
  // wrong size
  Matrix<double, Dynamic, Dynamic> V2 = MatrixXd::Identity(2, 2);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V2, W, m0, C0, T, rng),
               std::invalid_argument);
  // not square
  Matrix<double, Dynamic, Dynamic> V3 = MatrixXd::Identity(2, 3);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V3, W, m0, C0, T, rng),
               std::invalid_argument);
  // positive semi-definite is okay
  Matrix<double, Dynamic, Dynamic> V_psd = MatrixXd::Zero(3, 3);
  V_psd(0, 0) = 1.0;
  EXPECT_NO_THROW(gaussian_dlm_obs_rng(FF, GG, V_psd, W, m0, C0, T, rng));
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, PoliciesVVector) {
  boost::random::mt19937 rng;
  // negative
  Matrix<double, Dynamic, Dynamic> V_neg = V_vec;
  V_neg(0) = -1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_neg, W, m0, C0, T, rng),
               std::invalid_argument);
  // finite and NaN
  Matrix<double, Dynamic, Dynamic> V_infinite = V_vec;
  V_infinite(0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_infinite, W, m0, C0, T, rng),
               std::domain_error);
  Matrix<double, Dynamic, Dynamic> V_nan = V_vec;
  V_nan(0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_nan, W, m0, C0, T, rng),
               std::domain_error);
  // wrong size
  Matrix<double, Dynamic, 1> V_badsz
      = Matrix<double, Dynamic, 1>::Constant(2, 1.0);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_badsz, W, m0, C0, T, rng),
               std::invalid_argument);
  // positive semi-definite is okay
  Matrix<double, Dynamic, 1> V_psd = Matrix<double, Dynamic, 1>::Zero(3);
  V_psd(0) = 1.0;
  EXPECT_NO_THROW(gaussian_dlm_obs_rng(FF, GG, V_psd, W, m0, C0, T, rng));
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, Policiesm0) {
  boost::random::mt19937 rng;
  // m0
  // size
  Matrix<double, Dynamic, 1> m0_sz = Matrix<double, Dynamic, 1>::Zero(4, 1);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0_sz, C0, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W, m0_sz, C0, T, rng),
               std::invalid_argument);
  // finite and NaN
  Matrix<double, Dynamic, 1> m0_inf = m0;
  m0_inf(0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0_inf, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W, m0_inf, C0, T, rng),
               std::domain_error);
  Matrix<double, Dynamic, 1> m0_nan = m0;
  m0_nan(0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0_nan, C0, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W, m0_nan, C0, T, rng),
               std::domain_error);
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, PoliciesC0) {
  boost::random::mt19937 rng;
  // size
  Matrix<double, Dynamic, Dynamic> C0_sz = MatrixXd::Identity(3, 3);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0, C0_sz, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W, m0, C0_sz, T, rng),
               std::invalid_argument);
  // negative
  Matrix<double, Dynamic, Dynamic> C0_neg = C0;
  C0_neg(0, 0) = -1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0, C0_neg, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W, m0, C0_neg, T, rng),
               std::domain_error);
  // asymmetric
  Matrix<double, Dynamic, Dynamic> C0_asym = C0;
  C0_asym(0, 1) = 1;
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0, C0_asym, T, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W, m0, C0_neg, T, rng),
               std::domain_error);
  // not square
  Matrix<double, Dynamic, Dynamic> C0_notsq = MatrixXd::Identity(3, 2);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0, C0_notsq, T, rng),
               std::invalid_argument);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V_vec, W, m0, C0_notsq, T, rng),
               std::invalid_argument);
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, PoliciesT) {
  boost::random::mt19937 rng;
  // Must be positive.
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0, C0, 0, rng),
               std::domain_error);
  EXPECT_THROW(gaussian_dlm_obs_rng(FF, GG, V, W, m0, C0, -1, rng),
               std::domain_error);
}

TEST_F(ProbDistributionsGaussianDLMInputsRng, OutputCorrect) {
  boost::random::mt19937 rng;

  // Construct a Gaussian DLM process who's output is predictable, so
  // we can assert that the output is correct.

  // A robot starts at (0, 0) and wants to travel up and right, so
  // that at time t it is at position (t, t). The robot can control
  // its velocity, but with some noise. We observe the robot's
  // position with noise.

  // If the RNG is working correctly then we expect to get output like this:
  // [[0.1 , 0.9, 2.12, 3.01, ...],
  //  [-0.1, 1.1, 2.01, 2.98, ...]]

  // State space is:
  // (position x, position y,
  //  time t,
  //  velocity x, velocity y,
  //  constant component 1)
  const unsigned int n = 6;

  // Observation space is:
  // (observed position x, observed position y).
  const unsigned int r = 2;

  // Initially at (0, 0), at time 0, velocity (1, 1).
  // Constant component one is equal to 1.
  Matrix<double, Dynamic, 1> m0_(n, 1);
  m0_ << 0, 0, 0, 1, 1, 1;

  // Noise in initial state. Each component of the velocity is
  // affected by some Gaussian noise.
  double sigma2_vxx = 0.2;
  double sigma2_vxy = 0.1;
  double sigma2_vyy = 0.3;
  Matrix<double, Dynamic, Dynamic> C0_(n, n);
  C0_ <<
      // x, y, t, vx,         vy,         1
      0,
      0, 0, 0, 0, 0,                       // x
      0, 0, 0, 0, 0, 0,                    // y
      0, 0, 0, 0, 0, 0,                    // t
      0, 0, 0, sigma2_vxx, sigma2_vxy, 0,  // vx
      0, 0, 0, sigma2_vxy, sigma2_vyy, 0,  // vy
      0, 0, 0, 0, 0, 0;                    // 1

  // We observe just the position.
  Matrix<double, Dynamic, Dynamic> FF_(n, r);
  FF_ << 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;

  // We observe the position with noise.
  double sigma2_obs_xx = 0.3;
  double sigma2_obs_xy = 0.1;
  double sigma2_obs_yy = 0.2;
  Matrix<double, Dynamic, Dynamic> V_(r, r);
  V_ << sigma2_obs_xx, sigma2_obs_xy, sigma2_obs_xy, sigma2_obs_yy;

  // State transition
  // From physics:
  // x' = x + vx
  // y' = y + vy
  // t' = t + 1
  // From the controller:
  // vx' = (t + 2) - (x + vx) + noise
  // vy' = (t + 2) - (y + vy) + noise
  // one' = one
  // On the next step (t+1) the robot moves to x + vx, and due to the
  // time lag involved in controlling velocity while targeting
  // position, there is no changing that.
  // On the step after that (t+2) the robot will be at (x' + vx') = (x
  // + vx + (t + 2) - (x + vx) + noise) = t + 2 + noise.
  // So in general, x ~ t ~ y.
  Matrix<double, Dynamic, Dynamic> GG_(n, n);
  GG_ <<
      // x,  y, t, vx, vy, 1
      1,
      0, 0, 1, 0, 0,       // x
      0, 1, 0, 0, 1, 0,    // y
      0, 0, 1, 0, 0, 1,    // t
      -1, 0, 1, -1, 0, 2,  // vx
      0, -1, 1, 0, -1, 2,  // vx
      0, 0, 0, 0, 0, 1;    // 1

  // State transition noise.
  Matrix<double, Dynamic, Dynamic> W_(n, n);
  W_ <<
      // x, y, t,  vx,         vy,         1
      0,
      0, 0, 0, 0, 0,                       // x
      0, 0, 0, 0, 0, 0,                    // y
      0, 0, 0, 0, 0, 0,                    // t
      0, 0, 0, sigma2_vxx, sigma2_vxy, 0,  // vx
      0, 0, 0, sigma2_vxy, sigma2_vyy, 0,  // vy
      0, 0, 0, 0, 0, 0;                    // 1

  // The Gaussian DLM process is stable by construction. So we can run
  // it for a while.
  unsigned int T_ = 100;

  Matrix<double, Dynamic, Dynamic> y_
      = gaussian_dlm_obs_rng(FF_, GG_, V_, W_, m0_, C0_, T_, rng);

  assert(y_.rows() == r && r == 2 && y_.cols() == T_);

  double sum_squared_errors = 0;
  for (int t = 0; t < T_; ++t) {
    double observed_x_at_t = y_(0, t);
    double observed_y_at_t = y_(1, t);
    sum_squared_errors
        += pow(observed_x_at_t - t, 2) + pow(observed_y_at_t - t, 2);
  }

  // x[0] = 0
  // x[0] - 0 = 0
  // x[1] = x[0] + vx[0] = 0 + 1 + initial_noise_vx
  // x[1] - 1 = initial_noise_vx
  // For t>=1 ...
  // x[t+2] = x[t+1] + vx[t+1]
  //        = x[t] + vx[t] + ((t + 2) - (x[t] + vx[t])) + noise_vx[t]
  //        = (t + 2) + noise_vx[t]
  // x[t] - t = noise_vx[t - 2]
  // observed x[0] - 0 = observation_noise_x[0]
  // observed x[1] - 1 = initial_noise_vx + observation_noise_x[1]
  // observed x[t] - t = noise_vx[t - 2] + observation_noise_x[t]
  // (initial_noise_vx and noise_vx[0] are different noise, with, in
  // general, different variance. But here they are the same. (Note
  // that W_ = C0_.))
  // sum_squared_errors = sum (observed x[t] - t)^2 + (observed y[t] - t)^2

  // E[sum_squared_errors] ~=
  //    T * (sigma2_vxx + sigma2_vyy + sigma2_obs_xx + sigma2_obs_yy)

  const double expected_sse
      = T_ * (sigma2_vxx + sigma2_vyy + sigma2_obs_xx + sigma2_obs_yy);
  // expected_sse ~ 100.

  const double tolerance = 5;
  EXPECT_TRUE(fabs(sum_squared_errors - expected_sse) < tolerance);
}

TEST_F(ProbDistributionsGaussianDLMInputsRng,
       marginalChiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;

  // As in the previous test, but with simpler noise structure so we
  // can handle the distribution analytically. There is gaussian
  // noise of variance one in the state transition for the x-component
  // of velocity, and in the observation of the y component of
  // position. SSE will turn out to be chi-squared distributed.

  const unsigned int n = 6;
  const unsigned int r = 2;

  Matrix<double, Dynamic, 1> m0_(n, 1);
  m0_ << 0, 0, 0, 1, 1, 1;

  double sigma2_v = 1;
  Matrix<double, Dynamic, Dynamic> C0_(n, n);
  C0_ <<
      // x, y, t, vx,       vy, 1
      0,
      0, 0, 0, 0, 0,            // x
      0, 0, 0, 0, 0, 0,         // y
      0, 0, 0, 0, 0, 0,         // t
      0, 0, 0, sigma2_v, 0, 0,  // vx
      0, 0, 0, 0, 0, 0,         // vy
      0, 0, 0, 0, 0, 0;         // 1

  Matrix<double, Dynamic, Dynamic> FF_(n, r);
  FF_ << 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;

  double sigma2_obs = sigma2_v;
  Matrix<double, Dynamic, Dynamic> V_(r, r);
  V_ << 0, 0, 0, sigma2_obs;

  Matrix<double, Dynamic, Dynamic> GG_(n, n);
  GG_ <<
      // x,  y, t, vx, vy, 1
      1,
      0, 0, 1, 0, 0,       // x
      0, 1, 0, 0, 1, 0,    // y
      0, 0, 1, 0, 0, 1,    // t
      -1, 0, 1, -1, 0, 2,  // vx
      0, -1, 1, 0, -1, 2,  // vx
      0, 0, 0, 0, 0, 1;    // 1

  Matrix<double, Dynamic, Dynamic> W_(n, n);
  W_ <<
      // x, y, t,  vx,       vy, 1
      0,
      0, 0, 0, 0, 0,            // x
      0, 0, 0, 0, 0, 0,         // y
      0, 0, 0, 0, 0, 0,         // t
      0, 0, 0, sigma2_v, 0, 0,  // vx
      0, 0, 0, 0, 0, 0,         // vy
      0, 0, 0, 0, 0, 0;         // 1

  // We are sampling 10000 times, so T needs to be a bit small.
  unsigned int T_ = 10;

  // observed x[t] - t = noise_vx[t - 2] + observation_noise_x[t]
  //     = noise_vx[t - 2] + 0
  // observed y[t] - t = noise_vy[t - 2] + observation_noise_y[t]
  //     = 0 + observation_noise_y[t]
  // sum_squared_errors = sum (observed x[t] - t)^2 + (observed y[t] - t)^2
  //     = sum noise_vx[t - 2]^2 + observation_noise_y[t]^2
  // noise_vx[t], observation_noise_y[t] ~ normal(0, 1)
  // sum_squared_errors ~ chi_squared(2 * T - 1)
  // (2 * T - 1 because x[0] is noise free.)

  int N = 10000;
  int K = boost::math::round(2 * std::pow(N, 0.4));
  int df = 2 * T_ - 1;
  boost::math::chi_squared dist(df);

  std::vector<double> quantiles;
  for (int i = 1; i < K; i++)
    quantiles.push_back(quantile(dist, i * std::pow(K, -1.0)));
  quantiles.push_back(std::numeric_limits<double>::max());

  std::vector<double> samples;
  for (int count = 0; count < N; ++count) {
    Matrix<double, Dynamic, Dynamic> y_
        = gaussian_dlm_obs_rng(FF_, GG_, V_, W_, m0_, C0_, T_, rng);

    assert(y_.rows() == r && r == 2 && y_.cols() == T_);

    double sum_squared_errors = 0;
    for (int t = 0; t < T_; ++t) {
      double observed_x_at_t = y_(0, t);
      double observed_y_at_t = y_(1, t);
      sum_squared_errors
          += pow(observed_x_at_t - t, 2) + pow(observed_y_at_t - t, 2);
    }

    samples.push_back(sum_squared_errors);
  }

  assert_matches_quantiles(samples, quantiles, 1e-6);
}
