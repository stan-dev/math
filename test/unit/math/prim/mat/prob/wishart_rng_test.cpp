#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/util.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

TEST(ProbDistributionsWishartRng, rng) {
  using Eigen::MatrixXd;
  using stan::math::wishart_rng;

  boost::random::mt19937 rng;

  MatrixXd omega(3, 4);
  EXPECT_THROW(wishart_rng(3.0, omega, rng), std::invalid_argument);

  MatrixXd sigma(3, 3);
  sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 2.0, 1.0, 3.0;
  EXPECT_NO_THROW(wishart_rng(3.0, sigma, rng));
  EXPECT_THROW(wishart_rng(2, sigma, rng), std::domain_error);
  EXPECT_THROW(wishart_rng(-1, sigma, rng), std::domain_error);
}

TEST(probdistributionsWishartRng, symmetry) {
  using Eigen::MatrixXd;
  using stan::math::wishart_rng;
  using stan::test::unit::expect_symmetric;
  using stan::test::unit::spd_rng;

  boost::random::mt19937 rng;
  for (int k = 1; k < 20; ++k)
    for (double nu = k - 0.9; nu < k + 10; ++nu)
      for (int n = 0; n < 10; ++n)
        expect_symmetric(wishart_rng(nu, spd_rng(k, rng), rng));
}

TEST(ProbDistributionsWishart, marginalTwoChiSquareGoodnessFitTest) {
  using Eigen::MatrixXd;
  using boost::math::chi_squared;
  using boost::math::digamma;
  using stan::math::determinant;
  using stan::math::wishart_rng;
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
    a = wishart_rng(5.0, sigma, rng);
    avg += log(determinant(a));
  }
  avg /= N;
  double chi = (expect - avg) * (expect - avg) / expect;
  chi_squared mydist(1);
  EXPECT_TRUE(chi < quantile(complement(mydist, 1e-6)));
}

TEST(ProbDistributionsWishart, SpecialRNGTest) {
  // For any vector C != 0
  // (C' * W * C) / (C' * S * C)
  // must be chi-square distributed with df = k
  // which has mean = k and variance = 2k

  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::wishart_rng;

  boost::random::mt19937 rng(1234);

  MatrixXd sigma(3, 3);
  MatrixXd sigma_sym(3, 3);

  // wishart_rng should take only the lower part
  sigma << 9.0, -3.0, 1.0, 2.0, 4.0, -1.0, 2.0, 1.0, 3.0;

  sigma_sym << 9.0, 2.0, 2.0, 2.0, 4.0, 1.0, 2.0, 1.0, 3.0;

  VectorXd C(3);
  C << 2, 1, 3;

  size_t N = 1e4;
  int k = 20;
  // tolerance for variance
  double tol = 0.2;
  std::vector<double> acum;
  acum.reserve(N);
  for (size_t i = 0; i < N; i++)
    acum.push_back((C.transpose() * wishart_rng(k, sigma, rng) * C)(0)
                   / (C.transpose() * sigma_sym * C)(0));

  EXPECT_NEAR(1, stan::math::mean(acum) / k, tol * tol);
  EXPECT_NEAR(1, stan::math::variance(acum) / (2 * k), tol);
}
