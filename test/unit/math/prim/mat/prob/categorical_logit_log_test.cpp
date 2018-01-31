#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbCategoricalLogit, log_matches_lpmf) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> theta(3, 1);
  theta << -1, 2, -10;

  EXPECT_FLOAT_EQ((stan::math::categorical_logit_lpmf<true, double>(1, theta)),
                  (stan::math::categorical_logit_log<true, double>(1, theta)));
  EXPECT_FLOAT_EQ((stan::math::categorical_logit_lpmf<false, double>(1, theta)),
                  (stan::math::categorical_logit_log<false, double>(1, theta)));
  EXPECT_FLOAT_EQ((stan::math::categorical_logit_lpmf<double>(1, theta)),
                  (stan::math::categorical_logit_log<double>(1, theta)));
  EXPECT_FLOAT_EQ((stan::math::categorical_logit_lpmf(1, theta)),
                  (stan::math::categorical_logit_log(1, theta)));
  EXPECT_FLOAT_EQ((stan::math::categorical_logit_lpmf<true>(1, theta)),
                  (stan::math::categorical_logit_log<true>(1, theta)));
  EXPECT_FLOAT_EQ((stan::math::categorical_logit_lpmf<false>(1, theta)),
                  (stan::math::categorical_logit_log<false>(1, theta)));

  std::vector<int> ns(5);
  ns[0] = 1;
  ns[1] = 1;
  ns[2] = 3;
  ns[3] = 2;
  ns[4] = 3;

  EXPECT_FLOAT_EQ((stan::math::categorical_logit_lpmf<true, double>(ns, theta)),
                  (stan::math::categorical_logit_log<true, double>(ns, theta)));
  EXPECT_FLOAT_EQ(
      (stan::math::categorical_logit_lpmf<false, double>(ns, theta)),
      (stan::math::categorical_logit_log<false, double>(ns, theta)));
  EXPECT_FLOAT_EQ((stan::math::categorical_logit_lpmf<double>(ns, theta)),
                  (stan::math::categorical_logit_log<double>(ns, theta)));
  EXPECT_FLOAT_EQ((stan::math::categorical_logit_lpmf(ns, theta)),
                  (stan::math::categorical_logit_log(ns, theta)));
  EXPECT_FLOAT_EQ((stan::math::categorical_logit_lpmf<true>(ns, theta)),
                  (stan::math::categorical_logit_log<true>(ns, theta)));
  EXPECT_FLOAT_EQ((stan::math::categorical_logit_lpmf<false>(ns, theta)),
                  (stan::math::categorical_logit_log<false>(ns, theta)));
}
