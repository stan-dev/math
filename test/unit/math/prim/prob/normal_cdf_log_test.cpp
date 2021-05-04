#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST(ProbNormal, cdf_log_matches_lcdf) {
  double y = 0.8;
  double mu = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::normal_lcdf(y, mu, sigma)),
                  (stan::math::normal_cdf_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::normal_lcdf<double, double, double>(y, mu, sigma)),
      (stan::math::normal_cdf_log<double, double, double>(y, mu, sigma)));
}

TEST(ProbNormal, lcdf_tails) {
  using stan::math::normal_lcdf;
  using std::exp;

  EXPECT_FLOAT_EQ(4.60535300958196e-308, exp(normal_lcdf(-37.5, 0, 1)));
  EXPECT_FLOAT_EQ(5.72557122252458e-300, exp(normal_lcdf(-37, 0, 1)));
  EXPECT_FLOAT_EQ(5.54472571307484e-292, exp(normal_lcdf(-36.5, 0, 1)));
  EXPECT_FLOAT_EQ(4.18262406579728e-284, exp(normal_lcdf(-36, 0, 1)));
  EXPECT_FLOAT_EQ(2.45769154066194e-276, exp(normal_lcdf(-35.5, 0, 1)));
  EXPECT_FLOAT_EQ(1.12491070647241e-268, exp(normal_lcdf(-35, 0, 1)));
  EXPECT_FLOAT_EQ(4.01072896657726e-261, exp(normal_lcdf(-34.5, 0, 1)));
  EXPECT_FLOAT_EQ(1.11389878557438e-253, exp(normal_lcdf(-34, 0, 1)));
  EXPECT_FLOAT_EQ(2.40983869512039e-246, exp(normal_lcdf(-33.5, 0, 1)));
  EXPECT_FLOAT_EQ(4.06118562091586e-239, exp(normal_lcdf(-33, 0, 1)));
  EXPECT_FLOAT_EQ(5.33142435967881e-232, exp(normal_lcdf(-32.5, 0, 1)));
  EXPECT_FLOAT_EQ(5.4520806035124e-225, exp(normal_lcdf(-32, 0, 1)));
  EXPECT_FLOAT_EQ(4.34323260103177e-218, exp(normal_lcdf(-31.5, 0, 1)));
  EXPECT_FLOAT_EQ(2.6952500812005e-211, exp(normal_lcdf(-31, 0, 1)));
  EXPECT_FLOAT_EQ(1.30293791317808e-204, exp(normal_lcdf(-30.5, 0, 1)));
  EXPECT_FLOAT_EQ(4.90671392714819e-198, exp(normal_lcdf(-30, 0, 1)));
  EXPECT_FLOAT_EQ(1.43947455222918e-191, exp(normal_lcdf(-29.5, 0, 1)));
  EXPECT_FLOAT_EQ(3.28978526670438e-185, exp(normal_lcdf(-29, 0, 1)));
  EXPECT_FLOAT_EQ(5.85714125380634e-179, exp(normal_lcdf(-28.5, 0, 1)));
  EXPECT_FLOAT_EQ(8.12386946965943e-173, exp(normal_lcdf(-28, 0, 1)));
  EXPECT_FLOAT_EQ(8.77817055687808e-167, exp(normal_lcdf(-27.5, 0, 1)));
  EXPECT_FLOAT_EQ(7.38948100688502e-161, exp(normal_lcdf(-27, 0, 1)));
  EXPECT_FLOAT_EQ(4.84616266030332e-155, exp(normal_lcdf(-26.5, 0, 1)));
  EXPECT_FLOAT_EQ(2.47606331550339e-149, exp(normal_lcdf(-26, 0, 1)));
  EXPECT_FLOAT_EQ(9.85623651896393e-144, exp(normal_lcdf(-25.5, 0, 1)));
  EXPECT_FLOAT_EQ(3.05669670638256e-138, exp(normal_lcdf(-25, 0, 1)));
  EXPECT_FLOAT_EQ(7.38570686148941e-133, exp(normal_lcdf(-24.5, 0, 1)));
  EXPECT_FLOAT_EQ(1.3903921185497e-127, exp(normal_lcdf(-24, 0, 1)));
  EXPECT_FLOAT_EQ(2.03936756324998e-122, exp(normal_lcdf(-23.5, 0, 1)));
  EXPECT_FLOAT_EQ(2.33063700622065e-117, exp(normal_lcdf(-23, 0, 1)));
  EXPECT_FLOAT_EQ(2.07531079906636e-112, exp(normal_lcdf(-22.5, 0, 1)));
  EXPECT_FLOAT_EQ(1.43989243514508e-107, exp(normal_lcdf(-22, 0, 1)));
  EXPECT_FLOAT_EQ(7.78439707718263e-103, exp(normal_lcdf(-21.5, 0, 1)));
  EXPECT_FLOAT_EQ(3.27927801897904e-98, exp(normal_lcdf(-21, 0, 1)));
  EXPECT_FLOAT_EQ(1.0764673258791e-93, exp(normal_lcdf(-20.5, 0, 1)));
  EXPECT_FLOAT_EQ(2.75362411860623e-89, exp(normal_lcdf(-20, 0, 1)));
  EXPECT_FLOAT_EQ(5.48911547566041e-85, exp(normal_lcdf(-19.5, 0, 1)));
  EXPECT_FLOAT_EQ(8.52722395263098e-81, exp(normal_lcdf(-19, 0, 1)));
  EXPECT_FLOAT_EQ(1.03236986895633e-76, exp(normal_lcdf(-18.5, 0, 1)));
  EXPECT_FLOAT_EQ(9.74094891893715e-73, exp(normal_lcdf(-18, 0, 1)));
  EXPECT_FLOAT_EQ(7.16345876623504e-69, exp(normal_lcdf(-17.5, 0, 1)));
  EXPECT_FLOAT_EQ(4.10599620209891e-65, exp(normal_lcdf(-17, 0, 1)));
  EXPECT_FLOAT_EQ(1.83446300316473e-61, exp(normal_lcdf(-16.5, 0, 1)));
  EXPECT_FLOAT_EQ(6.38875440053809e-58, exp(normal_lcdf(-16, 0, 1)));
  EXPECT_FLOAT_EQ(1.73446079179387e-54, exp(normal_lcdf(-15.5, 0, 1)));
  EXPECT_FLOAT_EQ(3.67096619931275e-51, exp(normal_lcdf(-15, 0, 1)));
  EXPECT_FLOAT_EQ(6.05749476441522e-48, exp(normal_lcdf(-14.5, 0, 1)));
  EXPECT_FLOAT_EQ(7.7935368191928e-45, exp(normal_lcdf(-14, 0, 1)));
  EXPECT_FLOAT_EQ(7.81880730565789e-42, exp(normal_lcdf(-13.5, 0, 1)));
  EXPECT_FLOAT_EQ(6.11716439954988e-39, exp(normal_lcdf(-13, 0, 1)));
  EXPECT_FLOAT_EQ(3.73256429887771e-36, exp(normal_lcdf(-12.5, 0, 1)));
  EXPECT_FLOAT_EQ(1.77648211207768e-33, exp(normal_lcdf(-12, 0, 1)));
  EXPECT_FLOAT_EQ(6.59577144611367e-31, exp(normal_lcdf(-11.5, 0, 1)));
  EXPECT_FLOAT_EQ(1.91065957449868e-28, exp(normal_lcdf(-11, 0, 1)));
  EXPECT_FLOAT_EQ(4.31900631780923e-26, exp(normal_lcdf(-10.5, 0, 1)));
  EXPECT_FLOAT_EQ(7.61985302416053e-24, exp(normal_lcdf(-10, 0, 1)));
  EXPECT_FLOAT_EQ(1.04945150753626e-21, exp(normal_lcdf(-9.5, 0, 1)));
  EXPECT_FLOAT_EQ(1.12858840595384e-19, exp(normal_lcdf(-9, 0, 1)));
  EXPECT_FLOAT_EQ(9.47953482220332e-18, exp(normal_lcdf(-8.5, 0, 1)));
  EXPECT_FLOAT_EQ(6.22096057427178e-16, exp(normal_lcdf(-8, 0, 1)));
  EXPECT_FLOAT_EQ(3.1908916729109e-14, exp(normal_lcdf(-7.5, 0, 1)));
  EXPECT_FLOAT_EQ(1.27981254388584e-12, exp(normal_lcdf(-7, 0, 1)));
  EXPECT_FLOAT_EQ(4.01600058385912e-11, exp(normal_lcdf(-6.5, 0, 1)));
  EXPECT_FLOAT_EQ(9.86587645037698e-10, exp(normal_lcdf(-6, 0, 1)));
  EXPECT_FLOAT_EQ(1.89895624658877e-08, exp(normal_lcdf(-5.5, 0, 1)));
  EXPECT_FLOAT_EQ(2.86651571879194e-07, exp(normal_lcdf(-5, 0, 1)));
  EXPECT_FLOAT_EQ(3.39767312473006e-06, exp(normal_lcdf(-4.5, 0, 1)));
  EXPECT_FLOAT_EQ(3.16712418331199e-05, exp(normal_lcdf(-4, 0, 1)));
  EXPECT_FLOAT_EQ(0.000232629079035525, exp(normal_lcdf(-3.5, 0, 1)));
  EXPECT_FLOAT_EQ(0.00134989803163009, exp(normal_lcdf(-3, 0, 1)));
  EXPECT_FLOAT_EQ(0.00620966532577613, exp(normal_lcdf(-2.5, 0, 1)));
  EXPECT_FLOAT_EQ(0.0227501319481792, exp(normal_lcdf(-2, 0, 1)));
  EXPECT_FLOAT_EQ(0.0668072012688581, exp(normal_lcdf(-1.5, 0, 1)));
  EXPECT_FLOAT_EQ(0.158655253931457, exp(normal_lcdf(-1, 0, 1)));
  EXPECT_FLOAT_EQ(0.308537538725987, exp(normal_lcdf(-0.5, 0, 1)));
  EXPECT_FLOAT_EQ(0.5, exp(normal_lcdf(0, 0, 1)));
  EXPECT_FLOAT_EQ(0.691462461274013, exp(normal_lcdf(0.5, 0, 1)));
  EXPECT_FLOAT_EQ(0.841344746068543, exp(normal_lcdf(1, 0, 1)));
  EXPECT_FLOAT_EQ(0.933192798731142, exp(normal_lcdf(1.5, 0, 1)));
  EXPECT_FLOAT_EQ(0.977249868051821, exp(normal_lcdf(2, 0, 1)));
  EXPECT_FLOAT_EQ(0.993790334674224, exp(normal_lcdf(2.5, 0, 1)));
  EXPECT_FLOAT_EQ(0.99865010196837, exp(normal_lcdf(3, 0, 1)));
  EXPECT_FLOAT_EQ(0.999767370920964, exp(normal_lcdf(3.5, 0, 1)));
  EXPECT_FLOAT_EQ(0.999968328758167, exp(normal_lcdf(4, 0, 1)));
  EXPECT_FLOAT_EQ(0.999996602326875, exp(normal_lcdf(4.5, 0, 1)));
  EXPECT_FLOAT_EQ(0.999999713348428, exp(normal_lcdf(5, 0, 1)));
  EXPECT_FLOAT_EQ(0.999999981010438, exp(normal_lcdf(5.5, 0, 1)));
  EXPECT_FLOAT_EQ(0.999999999013412, exp(normal_lcdf(6, 0, 1)));
  EXPECT_FLOAT_EQ(0.99999999995984, exp(normal_lcdf(6.5, 0, 1)));
  EXPECT_FLOAT_EQ(0.99999999999872, exp(normal_lcdf(7, 0, 1)));
  EXPECT_FLOAT_EQ(0.999999999999968, exp(normal_lcdf(7.5, 0, 1)));
  EXPECT_FLOAT_EQ(0.999999999999999, exp(normal_lcdf(8, 0, 1)));
  EXPECT_FLOAT_EQ(1, exp(normal_lcdf(8.5, 0, 1)));
  EXPECT_FLOAT_EQ(1, exp(normal_lcdf(9, 0, 1)));
  EXPECT_FLOAT_EQ(1, exp(normal_lcdf(9.5, 0, 1)));
  EXPECT_FLOAT_EQ(1, exp(normal_lcdf(10, 0, 1)));
}
