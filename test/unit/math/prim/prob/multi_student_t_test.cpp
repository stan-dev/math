#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <vector>
#include <limits>

TEST(ProbDistributionsMultiStudentT, NotVectorized) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_log;
  using stan::math::multi_student_t_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  double nu = 4.0;
  double lp = multi_student_t_log(y, nu, mu, Sigma);
  EXPECT_NO_THROW(multi_student_t_rng(nu, mu, Sigma, rng));
  // calc using R's mnormt package's dmt function
  EXPECT_NEAR(-10.1245958182, lp, 1e-9);
}
TEST(ProbDistributionsMultiStudentT, Vectorized) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_log;
  using stan::math::multi_student_t_rng;
  using std::vector;
  boost::random::mt19937 rng;
  vector<Matrix<double, Dynamic, 1> > vec_y(2);
  vector<Matrix<double, 1, Dynamic> > vec_y_t(2);
  Matrix<double, Dynamic, 1> y(3);
  Matrix<double, 1, Dynamic> y_t(3);
  y << 3.0, -2.0, 10.0;
  vec_y[0] = y;
  vec_y_t[0] = y;
  y << 3.0, -1.0, 5.0;
  vec_y[1] = y;
  vec_y_t[1] = y;
  y_t = y;

  vector<Matrix<double, Dynamic, 1> > vec_mu(2);
  vector<Matrix<double, 1, Dynamic> > vec_mu_t(2);
  Matrix<double, Dynamic, 1> mu(3);
  Matrix<double, 1, Dynamic> mu_t(3);
  mu << 2.0, -1.0, 4.0;
  vec_mu[0] = mu;
  vec_mu_t[0] = mu;
  mu << 1.0, -3.0, 4.0;
  vec_mu[1] = mu;
  vec_mu_t[1] = mu;
  mu_t = mu;

  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 10.0, -3.0, 0.0, -3.0, 5.0, 0.0, 0.0, 0.0, 5.0;

  double nu = 4.0;

  // y and mu vectorized
  EXPECT_FLOAT_EQ(-8.9286697030 - 6.8183896234,
                  stan::math::multi_student_t_log(vec_y, nu, vec_mu, Sigma));
  EXPECT_FLOAT_EQ(-8.9286697030 - 6.8183896234,
                  stan::math::multi_student_t_log(vec_y_t, nu, vec_mu, Sigma));
  EXPECT_FLOAT_EQ(-8.9286697030 - 6.8183896234,
                  stan::math::multi_student_t_log(vec_y, nu, vec_mu_t, Sigma));
  EXPECT_FLOAT_EQ(
      -8.9286697030 - 6.8183896234,
      stan::math::multi_student_t_log(vec_y_t, nu, vec_mu_t, Sigma));

  // y vectorized
  EXPECT_FLOAT_EQ(-9.1670535409 - 6.8183896234,
                  stan::math::multi_student_t_log(vec_y, nu, mu, Sigma));
  EXPECT_FLOAT_EQ(-9.1670535409 - 6.8183896234,
                  stan::math::multi_student_t_log(vec_y_t, nu, mu, Sigma));
  EXPECT_FLOAT_EQ(-9.1670535409 - 6.8183896234,
                  stan::math::multi_student_t_log(vec_y, nu, mu_t, Sigma));
  EXPECT_FLOAT_EQ(-9.1670535409 - 6.8183896234,
                  stan::math::multi_student_t_log(vec_y_t, nu, mu_t, Sigma));

  // mu vectorized
  EXPECT_FLOAT_EQ(-5.5280118939 - 6.8183896234,
                  stan::math::multi_student_t_log(y, nu, vec_mu, Sigma));
  EXPECT_FLOAT_EQ(-5.5280118939 - 6.8183896234,
                  stan::math::multi_student_t_log(y_t, nu, vec_mu, Sigma));
  EXPECT_FLOAT_EQ(-5.5280118939 - 6.8183896234,
                  stan::math::multi_student_t_log(y, nu, vec_mu_t, Sigma));
  EXPECT_FLOAT_EQ(-5.5280118939 - 6.8183896234,
                  stan::math::multi_student_t_log(y_t, nu, vec_mu_t, Sigma));
  EXPECT_NO_THROW(stan::math::multi_student_t_rng(nu, vec_mu, Sigma, rng));
  EXPECT_NO_THROW(stan::math::multi_student_t_rng(nu, vec_mu_t, Sigma, rng));
}

TEST(ProbDistributionsMultiStudentT, Sigma) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_log;
  using stan::math::multi_student_t_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  double nu = 4.0;
  EXPECT_NO_THROW(multi_student_t_log(y, nu, mu, Sigma));
  EXPECT_NO_THROW(multi_student_t_rng(nu, mu, Sigma, rng));

  // non-symmetric
  Sigma(0, 1) = 10;
  EXPECT_THROW(multi_student_t_log(y, nu, mu, Sigma), std::domain_error);
  EXPECT_THROW(multi_student_t_rng(nu, mu, Sigma, rng), std::domain_error);
}

TEST(ProbDistributionsMultiStudentT, Mu) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_log;
  using stan::math::multi_student_t_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  double nu = 4.0;
  EXPECT_NO_THROW(multi_student_t_log(y, nu, mu, Sigma));
  EXPECT_NO_THROW(multi_student_t_rng(nu, mu, Sigma, rng));

  mu(0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(multi_student_t_log(y, nu, mu, Sigma), std::domain_error);
  EXPECT_THROW(multi_student_t_rng(nu, mu, Sigma, rng), std::domain_error);

  mu(0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(multi_student_t_log(y, nu, mu, Sigma), std::domain_error);
  EXPECT_THROW(multi_student_t_rng(nu, mu, Sigma, rng), std::domain_error);

  mu(0) = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(multi_student_t_log(y, nu, mu, Sigma), std::domain_error);
  EXPECT_THROW(multi_student_t_rng(nu, mu, Sigma, rng), std::domain_error);
}

TEST(ProbDistributionsMultiStudentT, Y) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_log;
  using stan::math::multi_student_t_rng;
  using std::vector;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  double nu = 4.0;
  EXPECT_NO_THROW(multi_student_t_log(y, nu, mu, Sigma));

  y(0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(multi_student_t_log(y, nu, mu, Sigma), std::domain_error);

  y(0) = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(multi_student_t_log(y, nu, mu, Sigma));

  y(0) = -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(multi_student_t_log(y, nu, mu, Sigma));
}

TEST(ProbDistributionsMultiStudentT, Nu) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_log;
  using stan::math::multi_student_t_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  double nu = 4.0;
  EXPECT_NO_THROW(multi_student_t_log(y, nu, mu, Sigma));
  EXPECT_NO_THROW(multi_student_t_rng(nu, mu, Sigma, rng));

  nu = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(multi_student_t_log(y, nu, mu, Sigma), std::domain_error);
  EXPECT_THROW(multi_student_t_rng(nu, mu, Sigma, rng), std::domain_error);

  nu = 0.0;
  EXPECT_THROW(multi_student_t_log(y, nu, mu, Sigma), std::domain_error);
  EXPECT_THROW(multi_student_t_rng(nu, mu, Sigma, rng), std::domain_error);

  nu = -1.0;
  EXPECT_THROW(multi_student_t_log(y, nu, mu, Sigma), std::domain_error);
  EXPECT_THROW(multi_student_t_rng(nu, mu, Sigma, rng), std::domain_error);

  nu = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(multi_student_t_log(y, nu, mu, Sigma), std::domain_error);
  EXPECT_THROW(multi_student_t_rng(nu, mu, Sigma, rng), std::domain_error);

  // nu = infinity NOT OK
  nu = std::numeric_limits<double>::infinity();
  EXPECT_THROW(multi_student_t_log(y, nu, mu, Sigma), std::domain_error);
  EXPECT_THROW(multi_student_t_rng(nu, mu, Sigma, rng), std::domain_error);
}

TEST(ProbDistributionsMultiStudentT, ErrorSize1) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_log;
  using stan::math::multi_student_t_rng;
  using std::vector;
  Matrix<double, Dynamic, 1> y(2, 1);
  y << 2.0, -2.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  double nu = 4.0;

  EXPECT_THROW(multi_student_t_log(y, nu, mu, Sigma), std::invalid_argument);
}

TEST(ProbDistributionsMultiStudentT, ErrorSize2) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_log;
  using stan::math::multi_student_t_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 2);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0;
  double nu = 4.0;

  EXPECT_THROW(multi_student_t_log(y, nu, mu, Sigma), std::invalid_argument);
  EXPECT_THROW(multi_student_t_rng(nu, mu, Sigma, rng), std::invalid_argument);
}

TEST(ProbDistributionsMultiStudentT, ErrorSize3) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_log;
  using stan::math::multi_student_t_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(2, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0;
  double nu = 4.0;

  EXPECT_THROW(multi_student_t_log(y, nu, mu, Sigma), std::invalid_argument);
  EXPECT_THROW(multi_student_t_rng(nu, mu, Sigma, rng), std::invalid_argument);
}

TEST(ProbDistributionsMultiStudentT, ErrorSizeSigma) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_lpdf;
  using stan::math::multi_student_t_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(2, 2);
  Sigma << 1.0, 0.0, 0.0, 1.0;
  double nu = 4.0;

  EXPECT_THROW(multi_student_t_lpdf(y, nu, mu, Sigma), std::invalid_argument);
  EXPECT_THROW(multi_student_t_rng(nu, mu, Sigma, rng), std::invalid_argument);
}

TEST(ProbDistributionsMultiStudentT, ProptoAllDoublesZero) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_log;
  using stan::math::multi_student_t_rng;
  using std::vector;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  double nu = 4.0;

  EXPECT_FLOAT_EQ(0.0, multi_student_t_log<true>(y, nu, mu, Sigma));
}

TEST(ProbDistributionsMultiStudentT, marginalOneChiSquareGoodnessFitTest) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_log;
  using stan::math::multi_student_t_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 2.0, 3.0, 11.0;

  Matrix<double, Dynamic, Dynamic> s(3, 3);
  s << 10.0, 3.0, 11.0, 3.0, 9.0, 1.2, 11.0, 1.2, 16.0;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::students_t_distribution<> dist(3.0);
  boost::math::chi_squared mydist(K - 1);

  double loc[K - 1];
  for (int i = 1; i < K; i++)
    loc[i - 1] = quantile(dist, i * std::pow(K, -1.0));

  int count = 0;
  std::vector<int> bin(K);
  std::vector<double> expect(K, N / K);

  Eigen::VectorXd a(mu.rows());
  while (count < N) {
    a = stan::math::multi_student_t_rng(3.0, mu, s, rng);
    a(0) = (a(0) - mu(0, 0)) / std::sqrt(s(0, 0));
    int i = 0;
    while (i < K - 1 && a(0) > loc[i])
      ++i;
    ++bin[i];
    ++count;
  }

  double X = 0;
  for (int j = 0; j < K; j++)
    X += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);

  EXPECT_TRUE(X < quantile(complement(mydist, 1e-6)));
}

TEST(ProbDistributionsMultiStudentT, marginalTwoChiSquareGoodnessFitTest) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_log;
  using stan::math::multi_student_t_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 2.0, 3.0, 11.0;

  Matrix<double, Dynamic, Dynamic> s(3, 3);
  s << 10.0, 3.0, 11.0, 3.0, 9.0, 1.2, 11.0, 1.2, 16.0;
  int N = 10000;
  int K = stan::math::round(2 * std::pow(N, 0.4));
  boost::math::students_t_distribution<> dist(3.0);
  boost::math::chi_squared mydist(K - 1);

  double loc[K - 1];
  for (int i = 1; i < K; i++)
    loc[i - 1] = quantile(dist, i * std::pow(K, -1.0));

  int count = 0;
  std::vector<int> bin(K);
  std::vector<double> expect(K, N / K);

  Eigen::VectorXd a(mu.rows());
  while (count < N) {
    a = stan::math::multi_student_t_rng(3.0, mu, s, rng);
    a(1) = (a(1) - mu(1, 0)) / std::sqrt(s(1, 1));
    int i = 0;
    while (i < K - 1 && a(1) > loc[i])
      ++i;
    ++bin[i];
    ++count;
  }

  double X = 0;
  for (int j = 0; j < K; j++)
    X += ((bin[j] - expect[j]) * (bin[j] - expect[j]) / expect[j]);

  EXPECT_TRUE(X < quantile(complement(mydist, 1e-6)));
}

TEST(ProbDistributionsMultiStudentT, WrongSize) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_log;
  using stan::math::multi_student_t_rng;
  using std::vector;
  vector<Matrix<double, Dynamic, 1> > y_3_3(3);
  vector<Matrix<double, Dynamic, 1> > y_3_1(3);
  vector<Matrix<double, Dynamic, 1> > y_3_2(3);
  vector<Matrix<double, Dynamic, 1> > y_1_3(1);
  vector<Matrix<double, Dynamic, 1> > y_2_3(2);
  Matrix<double, Dynamic, 1> y_3(3);
  Matrix<double, Dynamic, 1> y_2(2);
  Matrix<double, Dynamic, 1> y_1(1);
  y_3 << 2.0, -2.0, 11.0;
  y_2 << 2.0, -2.0;
  y_1 << 2.0;
  y_3_3[0] = y_3;
  y_3_3[1] = y_3;
  y_3_3[2] = y_3;
  y_3_1[0] = y_1;
  y_3_1[1] = y_1;
  y_3_1[2] = y_1;
  y_3_2[0] = y_2;
  y_3_2[1] = y_2;
  y_3_2[2] = y_2;
  y_1_3[0] = y_3;
  y_2_3[0] = y_3;
  y_2_3[1] = y_3;

  vector<Matrix<double, Dynamic, 1> > mu_3_3(3);
  Matrix<double, Dynamic, 1> mu_3(3);
  mu_3 << 2.0, -2.0, 11.0;
  mu_3_3[0] = mu_3;
  mu_3_3[1] = mu_3;
  mu_3_3[2] = mu_3;

  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 10.0, -3.0, 0.0, -3.0, 5.0, 0.0, 0.0, 0.0, 5.0;

  double nu = 4.0;

  EXPECT_NO_THROW(stan::math::multi_student_t_lpdf(y_3_3, nu, mu_3_3, Sigma));
  EXPECT_NO_THROW(stan::math::multi_student_t_lpdf(y_3, nu, mu_3_3, Sigma));

  EXPECT_THROW(stan::math::multi_student_t_lpdf(y_1_3, nu, mu_3_3, Sigma),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_student_t_lpdf(y_2_3, nu, mu_3_3, Sigma),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_student_t_lpdf(y_3_1, nu, mu_3_3, Sigma),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_student_t_lpdf(y_3_2, nu, mu_3_3, Sigma),
               std::invalid_argument);
}
