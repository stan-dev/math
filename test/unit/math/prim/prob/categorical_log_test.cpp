#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsCategorical, log_matches_lpmf) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> theta(3, 1);
  theta << 0.3, 0.5, 0.2;

  EXPECT_FLOAT_EQ((stan::math::categorical_lpmf<true>(1, theta)),
                  (stan::math::categorical_log<true>(1, theta)));
  EXPECT_FLOAT_EQ((stan::math::categorical_lpmf<false>(1, theta)),
                  (stan::math::categorical_log<false>(1, theta)));
  EXPECT_FLOAT_EQ((stan::math::categorical_lpmf(1, theta)),
                  (stan::math::categorical_log(1, theta)));

  std::vector<int> ns(5);
  ns[0] = 1;
  ns[1] = 2;
  ns[2] = 2;
  ns[3] = 1;
  ns[4] = 3;
  EXPECT_FLOAT_EQ((stan::math::categorical_lpmf<true>(ns, theta)),
                  (stan::math::categorical_log<true>(ns, theta)));
  EXPECT_FLOAT_EQ((stan::math::categorical_lpmf<false>(ns, theta)),
                  (stan::math::categorical_log<false>(ns, theta)));
  EXPECT_FLOAT_EQ((stan::math::categorical_lpmf(ns, theta)),
                  (stan::math::categorical_log(ns, theta)));
}
