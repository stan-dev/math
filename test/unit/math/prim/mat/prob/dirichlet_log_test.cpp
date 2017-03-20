#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(ProbDirichlet, log_matches_lpmf) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> theta(3,1);
  theta << 0.2, 0.3, 0.5;
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha(3,1);
  alpha << 1.0, 1.0, 1.0;
  EXPECT_FLOAT_EQ((stan::math::dirichlet_lpmf<true, double, double>(theta, alpha)),
                  (stan::math::dirichlet_log<true, double, double>(theta,alpha)));
  EXPECT_FLOAT_EQ((stan::math::dirichlet_lpmf<false, double, double>(theta, alpha)),
                  (stan::math::dirichlet_log<false, double, double>(theta,alpha)));
  EXPECT_FLOAT_EQ((stan::math::dirichlet_lpmf<double, double>(theta, alpha)),
                  (stan::math::dirichlet_log<double, double>(theta,alpha)));
  EXPECT_FLOAT_EQ((stan::math::dirichlet_lpmf(theta, alpha)),
                  (stan::math::dirichlet_log(theta,alpha)));
  EXPECT_FLOAT_EQ((stan::math::dirichlet_lpmf<true>(theta, alpha)),
                  (stan::math::dirichlet_log<true>(theta,alpha)));
  EXPECT_FLOAT_EQ((stan::math::dirichlet_lpmf<false>(theta, alpha)),
                  (stan::math::dirichlet_log<false>(theta,alpha)));
}
