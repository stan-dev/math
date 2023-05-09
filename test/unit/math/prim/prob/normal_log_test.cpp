#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbNormal, log_matches_lpdf) {
  double y = 0.8;
  double mu = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::normal_lpdf(y, mu, sigma)),
                  (stan::math::normal_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::normal_lpdf<true>(y, mu, sigma)),
                  (stan::math::normal_log<true>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::normal_lpdf<false>(y, mu, sigma)),
                  (stan::math::normal_log<false>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::normal_lpdf<true>(y, mu, sigma)),
                  (stan::math::normal_log<true>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::normal_lpdf<false>(y, mu, sigma)),
                  (stan::math::normal_log<false>(y, mu, sigma)));
  EXPECT_FLOAT_EQ((stan::math::normal_lpdf(y, mu, sigma)),
                  (stan::math::normal_log(y, mu, sigma)));
}

TEST(ProbNormal, test_vlpdf) {
  Eigen::Matrix<double, -1, 1> Y = Eigen::Matrix<double, -1, 1>::Random(5);
  Eigen::Matrix<double, -1, 1> Mu = Eigen::Matrix<double, -1, 1>::Random(5);
  Eigen::Matrix<double, -1, 1> Sigma
      = stan::math::abs(Eigen::Matrix<double, -1, 1>::Random(5));
  Eigen::Matrix<double, -1, 1> A
      = stan::math::normal_lpdf<false, stan::math::ProbReturnType::Vector>(
          Y, Mu, Sigma);
}
