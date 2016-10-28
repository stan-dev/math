#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(ProbMultinomial, log_matches_lpmf) {
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Eigen::Matrix<double,Eigen::Dynamic,1> theta(3,1);
  theta << 0.2, 0.3, 0.5;
  EXPECT_FLOAT_EQ((stan::math::multinomial_lpmf(ns,theta)),
                  (stan::math::multinomial_log(ns,theta)));
  EXPECT_FLOAT_EQ((stan::math::multinomial_lpmf<true>(ns,theta)),
                  (stan::math::multinomial_log<true>(ns,theta)));
  EXPECT_FLOAT_EQ((stan::math::multinomial_lpmf<false>(ns,theta)),
                  (stan::math::multinomial_log<false>(ns,theta)));
  EXPECT_FLOAT_EQ((stan::math::multinomial_lpmf<true, double>(ns,theta)),
                  (stan::math::multinomial_log<true, double>(ns,theta)));
  EXPECT_FLOAT_EQ((stan::math::multinomial_lpmf<false, double>(ns,theta)),
                  (stan::math::multinomial_log<false, double>(ns,theta)));
  EXPECT_FLOAT_EQ((stan::math::multinomial_lpmf<double>(ns,theta)),
                  (stan::math::multinomial_log<double>(ns,theta)));
}
