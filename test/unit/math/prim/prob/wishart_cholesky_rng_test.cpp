#include <stan/math/prim.hpp>
#include <test/unit/math/prim/util.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

TEST(ProbDistributionsWishartCholesky, rng) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;

  using stan::math::wishart_cholesky_rng;

  boost::random::mt19937 rng;

  MatrixXd omega(3, 4);
  EXPECT_THROW(wishart_cholesky_rng(3.0, omega, rng), std::domain_error);

  MatrixXd sigma(3, 3);
  sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 1.0, 0.0, 1.0, 3.0;

  Matrix<double, Dynamic, Dynamic> LS = sigma.llt().matrixL();

  EXPECT_NO_THROW(wishart_cholesky_rng(3.0, LS, rng));
  EXPECT_THROW(wishart_cholesky_rng(2, LS, rng), std::domain_error);
  EXPECT_THROW(wishart_cholesky_rng(-1, LS, rng), std::domain_error);
  LS(2, 2) = -1.0;
  EXPECT_THROW(wishart_cholesky_rng(3.0, LS, rng), std::domain_error);
}

TEST(ProbDistributionsWishartCholesky, rng_pos_def) {
  using Eigen::MatrixXd;
  using stan::math::wishart_cholesky_rng;

  boost::random::mt19937 rng;

  MatrixXd Sigma(2, 2);
  MatrixXd Sigma_non_pos_def(2, 2);

  Sigma << 1, 0, 0, 1;
  Sigma_non_pos_def << -1, 0, 0, 1;

  unsigned int dof = 5;

  EXPECT_NO_THROW(wishart_cholesky_rng(dof, Sigma, rng));
  EXPECT_THROW(wishart_cholesky_rng(dof, Sigma_non_pos_def, rng),
               std::domain_error);
}

TEST(ProbDistributionsWishartCholesky, marginalTwoChiSquareGoodnessFitTest) {
  using boost::math::chi_squared;
  using boost::math::digamma;
  using Eigen::MatrixXd;
  using stan::math::determinant;
  using stan::math::wishart_cholesky_rng;
  using std::log;

  boost::random::mt19937 rng;
  MatrixXd sigma(3, 3);
  sigma << 9.0, -3.0, 2.0, -3.0, 4.0, 0.0, 2.0, 0.0, 3.0;
  int N = 10000;

  double avg = 0;
  double expect = sigma.rows() * log(2.0) + log(determinant(sigma))
                  + digamma(5.0 / 2.0) + digamma(4.0 / 2.0)
                  + digamma(3.0 / 2.0);

  MatrixXd a(sigma.rows(), sigma.rows());
  for (int count = 0; count < N; ++count) {
    a = wishart_cholesky_rng(5.0, stan::math::cholesky_decompose(sigma), rng);
    avg += stan::math::sum(stan::math::log(a.diagonal()));
  }
  avg /= N;
  double chi = (expect - avg) * (expect - avg) / expect;
  chi_squared mydist(1);
  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}

TEST(ProbDistributionsWishartCholesky, SpecialRNGTest) {
  // For any vector C != 0
  // (C' * W * C) / (C' * S * C)
  // must be chi-square distributed with df = k
  // which has mean = k and variance = 2k

  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::multiply_lower_tri_self_transpose;
  using stan::math::wishart_cholesky_rng;

  boost::random::mt19937 rng(1234);

  MatrixXd sigma(3, 3);

  sigma << 9.0, 2.0, 2.0, 2.0, 4.0, 1.0, 2.0, 1.0, 3.0;

  Matrix<double, Dynamic, Dynamic> LS = sigma.llt().matrixL();

  VectorXd C(3);
  C << 2, 1, 3;

  size_t N = 1e4;
  int k = 20;
  // tolerance for variance
  double tol = 0.2;
  std::vector<double> acum;
  acum.reserve(N);
  for (size_t i = 0; i < N; i++)
    acum.push_back(
        (C.transpose()
         * multiply_lower_tri_self_transpose(wishart_cholesky_rng(k, LS, rng))
         * C)(0)
        / (C.transpose() * sigma * C)(0));

  EXPECT_NEAR(1, stan::math::mean(acum) / k, tol * tol);
  EXPECT_NEAR(1, stan::math::variance(acum) / (2 * k), tol);
}
