#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbStdNormal, ccdf_log_matches_lccdf) {
  double y = 0.8;

  EXPECT_FLOAT_EQ((stan::math::std_normal_lccdf(y)),
                  (stan::math::std_normal_ccdf_log(y)));
  EXPECT_FLOAT_EQ((stan::math::std_normal_lccdf<double>(y)),
                  (stan::math::std_normal_ccdf_log<double>(y)));
}

TEST(ProbStdNormal, lccdf_tail) {
  using stan::math::std_normal_lccdf;

  EXPECT_FLOAT_EQ(1, -6.661338147750941214694e-16 / std_normal_lccdf(-8.0));
  EXPECT_FLOAT_EQ(1, -3.186340080674249758114e-14 / std_normal_lccdf(-7.5));
  EXPECT_FLOAT_EQ(1, -1.279865102788699562477e-12 / std_normal_lccdf(-7.0));
  EXPECT_FLOAT_EQ(1, -4.015998644826973564545e-11 / std_normal_lccdf(-6.5));
  EXPECT_FLOAT_EQ(1, -9.865877009111571184118e-10 / std_normal_lccdf(-6.0));
  EXPECT_FLOAT_EQ(1, -1.898956265833514866414e-08 / std_normal_lccdf(-5.5));
  EXPECT_FLOAT_EQ(1, -2.866516130081049047962e-07 / std_normal_lccdf(-5.0));
  EXPECT_FLOAT_EQ(1, -3.397678896843115195074e-06 / std_normal_lccdf(-4.5));
  EXPECT_FLOAT_EQ(1, -3.167174337748932124543e-05 / std_normal_lccdf(-4.0));
  EXPECT_FLOAT_EQ(1, -0.0002326561413768195969113 / std_normal_lccdf(-3.5));
  EXPECT_FLOAT_EQ(1, -0.001350809964748202673598 / std_normal_lccdf(-3.0));
  EXPECT_FLOAT_EQ(1, -0.0062290254858600267035 / std_normal_lccdf(-2.5));
  EXPECT_FLOAT_EQ(1, -0.02301290932896348992442 / std_normal_lccdf(-2.0));
  EXPECT_FLOAT_EQ(1, -0.06914345561223400604689 / std_normal_lccdf(-1.5));
  EXPECT_FLOAT_EQ(1, -0.1727537790234499048836 / std_normal_lccdf(-1.0));
  EXPECT_FLOAT_EQ(1, -0.3689464152886565151412 / std_normal_lccdf(-0.5));
  EXPECT_FLOAT_EQ(1, -0.6931471805599452862268 / std_normal_lccdf(0));
  EXPECT_FLOAT_EQ(1, -1.175911761593618320987 / std_normal_lccdf(0.5));
  EXPECT_FLOAT_EQ(1, -1.841021645009263352222 / std_normal_lccdf(1.0));
  EXPECT_FLOAT_EQ(1, -2.705944400823889317564 / std_normal_lccdf(1.5));
  EXPECT_FLOAT_EQ(1, -3.78318433368203210776 / std_normal_lccdf(2.0));
  EXPECT_FLOAT_EQ(1, -5.081648277278686620662 / std_normal_lccdf(2.5));
  EXPECT_FLOAT_EQ(1, -6.607726221510342945464 / std_normal_lccdf(3.0));
  EXPECT_FLOAT_EQ(1, -8.366065308344028395027 / std_normal_lccdf(3.5));
  EXPECT_FLOAT_EQ(1, -10.36010148652728979357 / std_normal_lccdf(4.0));
  EXPECT_FLOAT_EQ(1, -12.59241973571053385683 / std_normal_lccdf(4.5));
  EXPECT_FLOAT_EQ(1, -15.06499839383403838156 / std_normal_lccdf(5.0));
  EXPECT_FLOAT_EQ(1, -17.77937635198566113104 / std_normal_lccdf(5.5));
  EXPECT_FLOAT_EQ(1, -20.73676889383495947072 / std_normal_lccdf(6.0));
  EXPECT_FLOAT_EQ(1, -23.93814997800869548428 / std_normal_lccdf(6.5));
}
