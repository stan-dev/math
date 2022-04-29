#include <stan/math/prim.hpp>
#include <test/unit/math/prim/util.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

TEST(ProbDistributionsInvWishartCholesky, rng) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;

  using stan::math::inv_wishart_cholesky_rng;

  boost::random::mt19937 rng;

  MatrixXd omega(3, 4);
  EXPECT_THROW(inv_wishart_cholesky_rng(3.0, omega, rng), std::domain_error);

  MatrixXd sigma(3, 3);
  sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 1.0, 0.0, 1.0, 3.0;

  Matrix<double, Dynamic, Dynamic> LS = sigma.llt().matrixL();

  EXPECT_NO_THROW(inv_wishart_cholesky_rng(3.0, LS, rng));
  EXPECT_THROW(inv_wishart_cholesky_rng(2, LS, rng), std::domain_error);
  EXPECT_THROW(inv_wishart_cholesky_rng(-1, LS, rng), std::domain_error);
  LS(2, 2) = -1;
  EXPECT_THROW(inv_wishart_cholesky_rng(3.0, LS, rng), std::domain_error);
}

TEST(ProbDistributionsInvWishartCholesky, rng_pos_def) {
  using Eigen::MatrixXd;
  using stan::math::inv_wishart_cholesky_rng;

  boost::random::mt19937 rng;

  MatrixXd Sigma(2, 2);
  MatrixXd Sigma_non_pos_def(2, 2);

  Sigma << 1, 0, 0, 1;
  Sigma_non_pos_def << -1, 0, 0, 1;

  unsigned int dof = 5;

  EXPECT_NO_THROW(inv_wishart_cholesky_rng(dof, Sigma, rng));
  EXPECT_THROW(inv_wishart_cholesky_rng(dof, Sigma_non_pos_def, rng),
               std::domain_error);
}

TEST(ProbDistributionsInvWishartCholesky, marginalTwoChiSquareGoodnessFitTest) {
  using boost::math::chi_squared;
  using boost::math::digamma;
  using Eigen::MatrixXd;
  using stan::math::determinant;
  using stan::math::inv_wishart_cholesky_rng;
  using std::log;

  boost::random::mt19937 rng;
  MatrixXd sigma(3, 3);
  sigma << 9.0, -3.0, 2.0, -3.0, 4.0, 0.0, 2.0, 0.0, 3.0;

  MatrixXd siginv(3, 3);
  siginv = sigma.inverse();
  int N = 10000;

  double avg = 0;
  double expect = sigma.rows() * log(2.0) + log(determinant(siginv))
                  + digamma(5.0 / 2.0) + digamma(4.0 / 2.0)
                  + digamma(3.0 / 2.0);

  MatrixXd a(sigma.rows(), sigma.rows());
  for (int count = 0; count < N; ++count) {
    a = inv_wishart_cholesky_rng(5.0, stan::math::cholesky_decompose(sigma),
                                 rng);
    avg += stan::math::sum(stan::math::log(a.diagonal()));
  }
  avg /= N;
  double chi = (expect - avg) * (expect - avg) / expect;
  chi_squared mydist(1);
  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}

TEST(ProbDistributionsInvWishartCholesky, SpecialRNGTest) {
  // When the scale matrix is an identity matrix and df = k + 2
  // The avg of the samples should also be an identity matrix

  using Eigen::MatrixXd;
  using stan::math::inv_wishart_cholesky_rng;
  using stan::math::multiply_lower_tri_self_transpose;

  boost::random::mt19937 rng(1234U);
  int N = 1e5;
  double tol = 0.1;
  for (int k = 1; k < 5; k++) {
    MatrixXd sigma = MatrixXd::Identity(k, k);
    MatrixXd Z = MatrixXd::Zero(k, k);
    for (int i = 0; i < N; i++) {
      Z += stan::math::crossprod(inv_wishart_cholesky_rng(k + 2, sigma, rng));
    }
    Z /= N;
    for (int j = 0; j < k; j++) {
      for (int i = 0; i < k; i++) {
        if (j == i)
          EXPECT_NEAR(Z(i, j), 1.0, tol);
        else
          EXPECT_NEAR(Z(i, j), 0.0, tol);
      }
    }
  }
}
