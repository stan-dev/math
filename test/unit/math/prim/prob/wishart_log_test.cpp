#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbWishart, log_matches_lpmf) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Sigma(2, 2);
  Sigma << 1.848220, 1.899623, 1.899623, 12.751941;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Y(2, 2);
  Y << 2.011108, -11.206611, -11.206611, 112.94139;

  unsigned int dof = 3;

  EXPECT_FLOAT_EQ((stan::math::wishart_lpdf(Y, dof, Sigma)),
                  (stan::math::wishart_log(Y, dof, Sigma)));
  EXPECT_FLOAT_EQ((stan::math::wishart_lpdf<true>(Y, dof, Sigma)),
                  (stan::math::wishart_log<true>(Y, dof, Sigma)));
  EXPECT_FLOAT_EQ((stan::math::wishart_lpdf<false>(Y, dof, Sigma)),
                  (stan::math::wishart_log<false>(Y, dof, Sigma)));
}
