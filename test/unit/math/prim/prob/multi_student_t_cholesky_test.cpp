#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <vector>
#include <limits>

TEST(ProbDistributionsMultiStudentTCholesky, NotVectorized) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_cholesky_lpdf;
  using stan::math::multi_student_t_cholesky_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;

  Matrix<double, Dynamic, Dynamic> L = Sigma.llt().matrixL();

  double nu = 4.0;
  double lp = multi_student_t_cholesky_lpdf(y, nu, mu, L);
  EXPECT_NO_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng));
  // calc using R's mnormt package's dmt function
  EXPECT_NEAR(-10.1245958182, lp, 1e-9);
}

TEST(ProbDistributionsMultiStudentTCholesky, Vectorized) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_cholesky_lpdf;
  using stan::math::multi_student_t_cholesky_rng;
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

  Matrix<double, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  double nu = 4.0;

  // y and mu vectorized
  EXPECT_FLOAT_EQ(
      -8.9286697030 - 6.8183896234,
      stan::math::multi_student_t_cholesky_lpdf(vec_y, nu, vec_mu, L));
  EXPECT_FLOAT_EQ(
      -8.9286697030 - 6.8183896234,
      stan::math::multi_student_t_cholesky_lpdf(vec_y_t, nu, vec_mu, L));
  EXPECT_FLOAT_EQ(
      -8.9286697030 - 6.8183896234,
      stan::math::multi_student_t_cholesky_lpdf(vec_y, nu, vec_mu_t, L));
  EXPECT_FLOAT_EQ(
      -8.9286697030 - 6.8183896234,
      stan::math::multi_student_t_cholesky_lpdf(vec_y_t, nu, vec_mu_t, L));

  // y vectorized
  EXPECT_FLOAT_EQ(-9.1670535409 - 6.8183896234,
                  stan::math::multi_student_t_cholesky_lpdf(vec_y, nu, mu, L));
  EXPECT_FLOAT_EQ(
      -9.1670535409 - 6.8183896234,
      stan::math::multi_student_t_cholesky_lpdf(vec_y_t, nu, mu, L));
  EXPECT_FLOAT_EQ(
      -9.1670535409 - 6.8183896234,
      stan::math::multi_student_t_cholesky_lpdf(vec_y, nu, mu_t, L));
  EXPECT_FLOAT_EQ(
      -9.1670535409 - 6.8183896234,
      stan::math::multi_student_t_cholesky_lpdf(vec_y_t, nu, mu_t, L));

  // mu vectorized
  EXPECT_FLOAT_EQ(-5.5280118939 - 6.8183896234,
                  stan::math::multi_student_t_cholesky_lpdf(y, nu, vec_mu, L));
  EXPECT_FLOAT_EQ(
      -5.5280118939 - 6.8183896234,
      stan::math::multi_student_t_cholesky_lpdf(y_t, nu, vec_mu, L));
  EXPECT_FLOAT_EQ(
      -5.5280118939 - 6.8183896234,
      stan::math::multi_student_t_cholesky_lpdf(y, nu, vec_mu_t, L));
  EXPECT_FLOAT_EQ(
      -5.5280118939 - 6.8183896234,
      stan::math::multi_student_t_cholesky_lpdf(y_t, nu, vec_mu_t, L));
  EXPECT_NO_THROW(stan::math::multi_student_t_cholesky_rng(nu, vec_mu, L, rng));
  EXPECT_NO_THROW(
      stan::math::multi_student_t_cholesky_rng(nu, vec_mu_t, L, rng));
}

TEST(ProbDistributionsMultiStudentTCholesky, Sigma) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_cholesky_lpdf;
  using stan::math::multi_student_t_cholesky_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  Matrix<double, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  double nu = 4.0;
  EXPECT_NO_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L));
  EXPECT_NO_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng));

  // non-symmetric
  L(0, 1) = 10;
  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L), std::domain_error);
  EXPECT_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng), std::domain_error);
}

TEST(ProbDistributionsMultiStudentTCholesky, Mu) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_cholesky_lpdf;
  using stan::math::multi_student_t_cholesky_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  Matrix<double, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  double nu = 4.0;
  EXPECT_NO_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L));
  EXPECT_NO_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng));

  mu(0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L), std::domain_error);
  EXPECT_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng), std::domain_error);

  mu(0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L), std::domain_error);
  EXPECT_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng), std::domain_error);

  mu(0) = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L), std::domain_error);
  EXPECT_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng), std::domain_error);
}

TEST(ProbDistributionsMultiStudentTCholesky, Y) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_cholesky_lpdf;
  using std::vector;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  Matrix<double, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  double nu = 4.0;
  EXPECT_NO_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L));

  y(0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L), std::domain_error);

  y(0) = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L));

  y(0) = -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L));
}

TEST(ProbDistributionsMultiStudentTCholesky, Nu) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_cholesky_lpdf;
  using stan::math::multi_student_t_cholesky_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  Matrix<double, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  double nu = 4.0;
  EXPECT_NO_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L));
  EXPECT_NO_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng));

  nu = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L), std::domain_error);
  EXPECT_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng), std::domain_error);

  nu = 0.0;
  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L), std::domain_error);
  EXPECT_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng), std::domain_error);

  nu = -1.0;
  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L), std::domain_error);
  EXPECT_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng), std::domain_error);

  nu = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L), std::domain_error);
  EXPECT_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng), std::domain_error);

  // nu = infinity is NOT OK
  nu = std::numeric_limits<double>::infinity();
  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L), std::domain_error);
  EXPECT_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng), std::domain_error);
}

TEST(ProbDistributionsMultiStudentTCholesky, ErrorSize1) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_cholesky_lpdf;
  using std::vector;
  Matrix<double, Dynamic, 1> y(2, 1);
  y << 2.0, -2.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  Matrix<double, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  double nu = 4.0;

  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L),
               std::invalid_argument);
}

TEST(ProbDistributionsMultiStudentTCholesky, ErrorSize0) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_cholesky_lpdf;
  using stan::math::multi_student_t_cholesky_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y_empty(0, 1);
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;

  Matrix<double, Dynamic, 1> mu_empty(0, 1);
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;

  Matrix<double, Dynamic, Dynamic> Sigma_empty(0, 0);
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  Matrix<double, Dynamic, Dynamic> L = Sigma.llt().matrixL();

  double nu_empty = std::numeric_limits<double>::quiet_NaN();
  double nu = 4;

  EXPECT_THROW(multi_student_t_cholesky_lpdf(y_empty, nu, mu, L),
               std::invalid_argument);
  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu_empty, mu, L),
               std::domain_error);
  EXPECT_THROW(multi_student_t_cholesky_lpdf(y_empty, nu, mu_empty, L),
               std::invalid_argument);
  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, Sigma_empty),
               std::domain_error);

  EXPECT_THROW(multi_student_t_cholesky_rng(nu_empty, mu, L, rng),
               std::domain_error);
  EXPECT_THROW(multi_student_t_cholesky_rng(nu, mu_empty, L, rng),
               std::invalid_argument);
  EXPECT_THROW(multi_student_t_cholesky_rng(nu, mu, Sigma_empty, rng),
               std::domain_error);
}

TEST(ProbDistributionsMultiStudentTCholesky, ErrorSize2) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_cholesky_lpdf;
  using stan::math::multi_student_t_cholesky_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(2, 2);
  Sigma << 1.0, 0.0, 0.0, 1.0;
  Matrix<double, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  double nu = 4.0;

  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L),
               std::invalid_argument);
  EXPECT_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng),
               std::invalid_argument);
}

TEST(ProbDistributionsMultiStudentTCholesky, ErrorSize3) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_cholesky_lpdf;
  using stan::math::multi_student_t_cholesky_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(4, 4);
  Sigma << 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
      0.0, 1.0;
  Matrix<double, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  double nu = 4.0;

  EXPECT_THROW(multi_student_t_cholesky_lpdf(y, nu, mu, L),
               std::invalid_argument);
  EXPECT_THROW(multi_student_t_cholesky_rng(nu, mu, L, rng),
               std::invalid_argument);
}

TEST(ProbDistributionsMultiStudentTCholesky, ProptoAllDoublesZero) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_cholesky_lpdf;
  using stan::math::multi_student_t_rng;
  using std::vector;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  Matrix<double, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  double nu = 4.0;

  EXPECT_FLOAT_EQ(0.0, multi_student_t_cholesky_lpdf<true>(y, nu, mu, L));
}

TEST(ProbDistributionsMultiStudentTCholesky,
     marginalOneChiSquareGoodnessFitTest) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_cholesky_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 2.0, 3.0, 11.0;

  Matrix<double, Dynamic, Dynamic> s(3, 3);
  s << 10.0, 3.0, 11.0, 3.0, 9.0, 1.2, 11.0, 1.2, 16.0;
  Matrix<double, Dynamic, Dynamic> L = s.llt().matrixL();
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
    a = multi_student_t_cholesky_rng(3.0, mu, L, rng);
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

TEST(ProbDistributionsMultiStudentTCholesky,
     marginalTwoChiSquareGoodnessFitTest) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_cholesky_lpdf;
  using stan::math::multi_student_t_cholesky_rng;
  using std::vector;
  boost::random::mt19937 rng;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 2.0, 3.0, 11.0;

  Matrix<double, Dynamic, Dynamic> s(3, 3);
  s << 10.0, 3.0, 11.0, 3.0, 9.0, 1.2, 11.0, 1.2, 16.0;
  Matrix<double, Dynamic, Dynamic> L = s.llt().matrixL();
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
    a = multi_student_t_cholesky_rng(3.0, mu, L, rng);
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

TEST(ProbDistributionsMultiStudentTCholesky, WrongSize) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::multi_student_t_cholesky_lpdf;
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
  Matrix<double, Dynamic, Dynamic> L = Sigma.llt().matrixL();

  double nu = 4.0;

  EXPECT_NO_THROW(
      stan::math::multi_student_t_cholesky_lpdf(y_3_3, nu, mu_3_3, L));
  EXPECT_NO_THROW(
      stan::math::multi_student_t_cholesky_lpdf(y_3, nu, mu_3_3, L));

  EXPECT_THROW(stan::math::multi_student_t_cholesky_lpdf(y_1_3, nu, mu_3_3, L),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_student_t_cholesky_lpdf(y_2_3, nu, mu_3_3, L),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_student_t_cholesky_lpdf(y_3_1, nu, mu_3_3, L),
               std::invalid_argument);
  EXPECT_THROW(stan::math::multi_student_t_cholesky_lpdf(y_3_2, nu, mu_3_3, L),
               std::invalid_argument);
}
