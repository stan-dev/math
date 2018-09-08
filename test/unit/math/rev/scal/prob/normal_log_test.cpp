#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbDistributionsNormal, intVsDouble) {
  using stan::math::var;
  for (double thetaval = -5.0; thetaval < 6.0; thetaval += 0.5) {
    var theta(thetaval);
    var lp1(0.0);
    lp1 += stan::math::normal_log<true>(0, theta, 1);
    double lp1val = lp1.val();
    stan::math::grad(lp1.vi_);
    double lp1adj = lp1.adj();

    var theta2(thetaval);
    var lp2(0.0);
    lp2 += stan::math::normal_log<true>(theta2, 0, 1);
    double lp2val = lp2.val();
    stan::math::grad(lp2.vi_);
    double lp2adj = lp2.adj();
    EXPECT_FLOAT_EQ(lp1val, lp2val);
    EXPECT_FLOAT_EQ(lp1adj, lp2adj);
  }
}

TEST(ProbDistributionsNormal, propto) {
  using stan::math::var;

  double mu1 = 1.0;
  double mu2 = 2.0;

  double d1 = stan::math::normal_log<true>(0, mu1, 1)
              - stan::math::normal_log<false>(0, mu1, 1);
  double d2 = stan::math::normal_log<true>(0, mu2, 1)
              - stan::math::normal_log<false>(0, mu2, 1);

  double v1
      = stan::math::value_of(stan::math::normal_log<true>(0, var(mu1), 1)
                             - stan::math::normal_log<false>(0, var(mu1), 1));
  double v2
      = stan::math::value_of(stan::math::normal_log<true>(0, var(mu2), 1)
                             - stan::math::normal_log<false>(0, var(mu2), 1));

  EXPECT_FLOAT_EQ(v1, v2);
  EXPECT_FLOAT_EQ(d1, d2);
}
