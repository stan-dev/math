#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbNormal, ccdf_log_matches_lccdf) {
  double y = 0.8;
  double mu = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::normal_lccdf(y, mu, sigma)),
                  (stan::math::normal_ccdf_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::normal_lccdf<double, double, double>(y, mu, sigma)),
      (stan::math::normal_ccdf_log<double, double, double>(y, mu, sigma)));
}

TEST(ProbNormal, lccdf_tail) {
  using stan::math::normal_lccdf;

  EXPECT_FLOAT_EQ(-6.661338147750941214694e-16, normal_lccdf(-8.0, 0, 1));
  EXPECT_FLOAT_EQ(-3.186340080674249758114e-14, normal_lccdf(-7.5, 0, 1));
  EXPECT_FLOAT_EQ(-1.279865102788699562477e-12, normal_lccdf(-7.0, 0, 1));
  EXPECT_FLOAT_EQ(-4.015998644826973564545e-11, normal_lccdf(-6.5, 0, 1));
  EXPECT_FLOAT_EQ(-9.865877009111571184118e-10, normal_lccdf(-6.0, 0, 1));
  EXPECT_FLOAT_EQ(-1.898956265833514866414e-08, normal_lccdf(-5.5, 0, 1));
  EXPECT_FLOAT_EQ(-2.866516130081049047962e-07, normal_lccdf(-5.0, 0, 1));
  EXPECT_FLOAT_EQ(-3.397678896843115195074e-06, normal_lccdf(-4.5, 0, 1));
  EXPECT_FLOAT_EQ(-3.167174337748932124543e-05, normal_lccdf(-4.0, 0, 1));
  EXPECT_FLOAT_EQ(-0.0002326561413768195969113, normal_lccdf(-3.5, 0, 1));
  EXPECT_FLOAT_EQ(-0.001350809964748202673598, normal_lccdf(-3.0, 0, 1));
  EXPECT_FLOAT_EQ(-0.0062290254858600267035, normal_lccdf(-2.5, 0, 1));
  EXPECT_FLOAT_EQ(-0.02301290932896348992442, normal_lccdf(-2.0, 0, 1));
  EXPECT_FLOAT_EQ(-0.06914345561223400604689, normal_lccdf(-1.5, 0, 1));
  EXPECT_FLOAT_EQ(-0.1727537790234499048836, normal_lccdf(-1.0, 0, 1));
  EXPECT_FLOAT_EQ(-0.3689464152886565151412, normal_lccdf(-0.5, 0, 1));
  EXPECT_FLOAT_EQ(-0.6931471805599452862268, normal_lccdf(0, 0, 1));
  EXPECT_FLOAT_EQ(-1.175911761593618320987, normal_lccdf(0.5, 0, 1));
  EXPECT_FLOAT_EQ(-1.841021645009263352222, normal_lccdf(1.0, 0, 1));
  EXPECT_FLOAT_EQ(-2.705944400823889317564, normal_lccdf(1.5, 0, 1));
  EXPECT_FLOAT_EQ(-3.78318433368203210776, normal_lccdf(2.0, 0, 1));
  EXPECT_FLOAT_EQ(-5.081648277278686620662, normal_lccdf(2.5, 0, 1));
  EXPECT_FLOAT_EQ(-6.607726221510342945464, normal_lccdf(3.0, 0, 1));
  EXPECT_FLOAT_EQ(-8.366065308344028395027, normal_lccdf(3.5, 0, 1));
  EXPECT_FLOAT_EQ(-10.36010148652728979357, normal_lccdf(4.0, 0, 1));
  EXPECT_FLOAT_EQ(-12.59241973571053385683, normal_lccdf(4.5, 0, 1));
  EXPECT_FLOAT_EQ(-15.06499839383403838156, normal_lccdf(5.0, 0, 1));
  EXPECT_FLOAT_EQ(-17.77937635198566113104, normal_lccdf(5.5, 0, 1));
  EXPECT_FLOAT_EQ(-20.73676889383495947072, normal_lccdf(6.0, 0, 1));
  EXPECT_FLOAT_EQ(-23.93814997800869548428, normal_lccdf(6.5, 0, 1));
}
