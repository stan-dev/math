#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/util.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <stdexcept>

TEST(ProbDistributionsInvWishartRng, rng) {
  using Eigen::MatrixXd;
  using stan::math::inv_wishart_rng;
  boost::random::mt19937 rng;

  MatrixXd omega(3, 4);
  EXPECT_THROW(inv_wishart_rng(3.0, omega, rng), std::invalid_argument);

  MatrixXd sigma(3, 3);
  sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 2.0, 1.0, 3.0;
  EXPECT_NO_THROW(inv_wishart_rng(3.0, sigma, rng));
  EXPECT_THROW(inv_wishart_rng(2, sigma, rng), std::domain_error);
  EXPECT_THROW(inv_wishart_rng(-1, sigma, rng), std::domain_error);
}
TEST(probdistributionsInvWishartRng, symmetry) {
  using Eigen::MatrixXd;
  using stan::math::inv_wishart_rng;
  using stan::test::unit::expect_symmetric;
  using stan::test::unit::spd_rng;

  boost::random::mt19937 rng;
  for (int k = 1; k < 20; ++k)
    for (double nu = k - 0.5; nu < k + 20; ++nu)
      for (int n = 0; n < 10; ++n)
        expect_symmetric(inv_wishart_rng(nu, spd_rng(k, rng), rng));
}

TEST(ProbDistributionsInvWishart, chiSquareGoodnessFitTest) {
  using Eigen::MatrixXd;
  using boost::math::chi_squared;
  using boost::math::digamma;
  using stan::math::determinant;
  using stan::math::inv_wishart_rng;
  using std::log;

  boost::random::mt19937 rng;
  MatrixXd sigma(3, 3);
  sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 1.0, 0.0, 1.0, 3.0;
  int N = 10000;

  MatrixXd siginv(3, 3);
  siginv = sigma.inverse();
  int count = 0;
  double avg = 0;
  double expect = sigma.rows() * log(2.0) + log(determinant(siginv))
                  + digamma(5.0 / 2.0) + digamma(4.0 / 2.0)
                  + digamma(3.0 / 2.0);

  MatrixXd a(sigma.rows(), sigma.rows());
  while (count < N) {
    a = inv_wishart_rng(5.0, sigma, rng);
    avg += log(determinant(a)) / N;
    count++;
  }
  double chi = (expect - avg) * (expect - avg) / expect;
  chi_squared mydist(1);
  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}

TEST(ProbDistributionsInvWishart, SpecialRNGTest) {
  // When the scale matrix is an identity matrix and df = k + 2
  // The avg of the samples should also be an identity matrix

  using Eigen::MatrixXd;
  using stan::math::inv_wishart_rng;

  boost::random::mt19937 rng(1234U);
  int N = 1e5;
  double tol = 0.1;
  for (int k = 1; k < 5; k++) {
    MatrixXd sigma = MatrixXd::Identity(k, k);
    MatrixXd Z = MatrixXd::Zero(k, k);
    for (int i = 0; i < N; i++)
      Z += inv_wishart_rng(k + 2, sigma, rng);
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
