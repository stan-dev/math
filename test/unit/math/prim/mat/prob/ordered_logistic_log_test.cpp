#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(ProbOrderedLogistic, log_matches_lpmf) {
  int K = 5;
  Eigen::Matrix<double,Eigen::Dynamic,1> c(K-1);
  c << -1.7, -0.3, 1.2, 2.6;
  double lambda = 1.1;
  
  EXPECT_FLOAT_EQ((stan::math::ordered_logistic_lpmf(3, lambda, c)),
                  (stan::math::ordered_logistic_log(3, lambda, c)));
  EXPECT_FLOAT_EQ((stan::math::ordered_logistic_lpmf<true>(3, lambda, c)),
                  (stan::math::ordered_logistic_log<true>(3, lambda, c)));
  EXPECT_FLOAT_EQ((stan::math::ordered_logistic_lpmf<false>(3, lambda, c)),
                  (stan::math::ordered_logistic_log<false>(3, lambda, c)));
  EXPECT_FLOAT_EQ((stan::math::ordered_logistic_lpmf<true, double, double>(3, lambda, c)),
                  (stan::math::ordered_logistic_log<true, double, double>(3, lambda, c)));
  EXPECT_FLOAT_EQ((stan::math::ordered_logistic_lpmf<false, double, double>(3, lambda, c)),
                  (stan::math::ordered_logistic_log<false, double, double>(3, lambda, c)));
  EXPECT_FLOAT_EQ((stan::math::ordered_logistic_lpmf<double, double>(3, lambda, c)),
                  (stan::math::ordered_logistic_log<double, double>(3, lambda, c)));

}
