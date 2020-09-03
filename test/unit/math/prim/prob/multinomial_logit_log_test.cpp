#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbMultinomialLogit, log_matches_lpmf) {
  using stan::math::multinomial_logit_log;
  using stan::math::multinomial_logit_lpmf;

  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Eigen::Matrix<double, Eigen::Dynamic, 1> theta(3, 1);
  theta << log(0.2), log(0.3), log(0.5);
  EXPECT_FLOAT_EQ((multinomial_logit_lpmf(ns, theta)),
                  (multinomial_logit_log(ns, theta)));
  EXPECT_FLOAT_EQ((multinomial_logit_lpmf<true>(ns, theta)),
                  (multinomial_logit_log<true>(ns, theta)));
  EXPECT_FLOAT_EQ((multinomial_logit_lpmf<false>(ns, theta)),
                  (multinomial_logit_log<false>(ns, theta)));
}
