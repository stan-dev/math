#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbInvWishart, log_matches_lpmf) {
  Eigen::MatrixXd Y(4, 4);
  Eigen::MatrixXd Sigma(4, 4);
  Y << 7.988168, -9.555605, -14.47483, 4.395895, -9.555605, 44.750570, 49.21577,
      -18.454186, -14.474830, 49.21577, 60.08987, -21.48108, 4.395895,
      -18.454186, -21.48108, 7.885833;
  Sigma << 2.9983662, 0.2898776, -2.650523, 0.1055911, 0.2898776, 11.4803610,
      7.1579931, -3.1129955, -2.650523, 7.1579931, 11.676181, -3.5866852,
      0.1055911, -3.1129955, -3.5866852, 1.4482736;
  unsigned int dof = 5;

  EXPECT_FLOAT_EQ((stan::math::inv_wishart_lpdf(Y, dof, Sigma)),
                  (stan::math::inv_wishart_log(Y, dof, Sigma)));
  EXPECT_FLOAT_EQ((stan::math::inv_wishart_lpdf<true>(Y, dof, Sigma)),
                  (stan::math::inv_wishart_log<true>(Y, dof, Sigma)));
  EXPECT_FLOAT_EQ((stan::math::inv_wishart_lpdf<false>(Y, dof, Sigma)),
                  (stan::math::inv_wishart_log<false>(Y, dof, Sigma)));
  EXPECT_FLOAT_EQ((stan::math::inv_wishart_lpdf<true>(Y, dof, Sigma)),
                  (stan::math::inv_wishart_log<true>(Y, dof, Sigma)));
  EXPECT_FLOAT_EQ((stan::math::inv_wishart_lpdf<false>(Y, dof, Sigma)),
                  (stan::math::inv_wishart_log<false>(Y, dof, Sigma)));
}
