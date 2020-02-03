#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbStdNormal, cdf_log_matches_lcdf) {
  double y = 0.8;

  EXPECT_FLOAT_EQ((stan::math::std_normal_lcdf(y)),
                  (stan::math::std_normal_cdf_log(y)));
  EXPECT_FLOAT_EQ((stan::math::std_normal_lcdf<double>(y)),
                  (stan::math::std_normal_cdf_log<double>(y)));
}

TEST(ProbStdNormal, lcdf_tails) {
  using stan::math::std_normal_lcdf;
  using std::exp;

  EXPECT_FLOAT_EQ(4.60535300958196e-308, exp(std_normal_lcdf(-37.5)));
  EXPECT_FLOAT_EQ(5.72557122252458e-300, exp(std_normal_lcdf(-37)));
  EXPECT_FLOAT_EQ(5.54472571307484e-292, exp(std_normal_lcdf(-36.5)));
  EXPECT_FLOAT_EQ(4.18262406579728e-284, exp(std_normal_lcdf(-36)));
  EXPECT_FLOAT_EQ(2.45769154066194e-276, exp(std_normal_lcdf(-35.5)));
  EXPECT_FLOAT_EQ(1.12491070647241e-268, exp(std_normal_lcdf(-35)));
  EXPECT_FLOAT_EQ(4.01072896657726e-261, exp(std_normal_lcdf(-34.5)));
  EXPECT_FLOAT_EQ(1.11389878557438e-253, exp(std_normal_lcdf(-34)));
  EXPECT_FLOAT_EQ(2.40983869512039e-246, exp(std_normal_lcdf(-33.5)));
  EXPECT_FLOAT_EQ(4.06118562091586e-239, exp(std_normal_lcdf(-33)));
  EXPECT_FLOAT_EQ(5.33142435967881e-232, exp(std_normal_lcdf(-32.5)));
  EXPECT_FLOAT_EQ(5.4520806035124e-225, exp(std_normal_lcdf(-32)));
  EXPECT_FLOAT_EQ(4.34323260103177e-218, exp(std_normal_lcdf(-31.5)));
  EXPECT_FLOAT_EQ(2.6952500812005e-211, exp(std_normal_lcdf(-31)));
  EXPECT_FLOAT_EQ(1.30293791317808e-204, exp(std_normal_lcdf(-30.5)));
  EXPECT_FLOAT_EQ(4.90671392714819e-198, exp(std_normal_lcdf(-30)));
  EXPECT_FLOAT_EQ(1.43947455222918e-191, exp(std_normal_lcdf(-29.5)));
  EXPECT_FLOAT_EQ(3.28978526670438e-185, exp(std_normal_lcdf(-29)));
  EXPECT_FLOAT_EQ(5.85714125380634e-179, exp(std_normal_lcdf(-28.5)));
  EXPECT_FLOAT_EQ(8.12386946965943e-173, exp(std_normal_lcdf(-28)));
  EXPECT_FLOAT_EQ(8.77817055687808e-167, exp(std_normal_lcdf(-27.5)));
  EXPECT_FLOAT_EQ(7.38948100688502e-161, exp(std_normal_lcdf(-27)));
  EXPECT_FLOAT_EQ(4.84616266030332e-155, exp(std_normal_lcdf(-26.5)));
  EXPECT_FLOAT_EQ(2.47606331550339e-149, exp(std_normal_lcdf(-26)));
  EXPECT_FLOAT_EQ(9.85623651896393e-144, exp(std_normal_lcdf(-25.5)));
  EXPECT_FLOAT_EQ(3.05669670638256e-138, exp(std_normal_lcdf(-25)));
  EXPECT_FLOAT_EQ(7.38570686148941e-133, exp(std_normal_lcdf(-24.5)));
  EXPECT_FLOAT_EQ(1.3903921185497e-127, exp(std_normal_lcdf(-24)));
  EXPECT_FLOAT_EQ(2.03936756324998e-122, exp(std_normal_lcdf(-23.5)));
  EXPECT_FLOAT_EQ(2.33063700622065e-117, exp(std_normal_lcdf(-23)));
  EXPECT_FLOAT_EQ(2.07531079906636e-112, exp(std_normal_lcdf(-22.5)));
  EXPECT_FLOAT_EQ(1.43989243514508e-107, exp(std_normal_lcdf(-22)));
  EXPECT_FLOAT_EQ(7.78439707718263e-103, exp(std_normal_lcdf(-21.5)));
  EXPECT_FLOAT_EQ(3.27927801897904e-98, exp(std_normal_lcdf(-21)));
  EXPECT_FLOAT_EQ(1.0764673258791e-93, exp(std_normal_lcdf(-20.5)));
  EXPECT_FLOAT_EQ(2.75362411860623e-89, exp(std_normal_lcdf(-20)));
  EXPECT_FLOAT_EQ(5.48911547566041e-85, exp(std_normal_lcdf(-19.5)));
  EXPECT_FLOAT_EQ(8.52722395263098e-81, exp(std_normal_lcdf(-19)));
  EXPECT_FLOAT_EQ(1.03236986895633e-76, exp(std_normal_lcdf(-18.5)));
  EXPECT_FLOAT_EQ(9.74094891893715e-73, exp(std_normal_lcdf(-18)));
  EXPECT_FLOAT_EQ(7.16345876623504e-69, exp(std_normal_lcdf(-17.5)));
  EXPECT_FLOAT_EQ(4.10599620209891e-65, exp(std_normal_lcdf(-17)));
  EXPECT_FLOAT_EQ(1.83446300316473e-61, exp(std_normal_lcdf(-16.5)));
  EXPECT_FLOAT_EQ(6.38875440053809e-58, exp(std_normal_lcdf(-16)));
  EXPECT_FLOAT_EQ(1.73446079179387e-54, exp(std_normal_lcdf(-15.5)));
  EXPECT_FLOAT_EQ(3.67096619931275e-51, exp(std_normal_lcdf(-15)));
  EXPECT_FLOAT_EQ(6.05749476441522e-48, exp(std_normal_lcdf(-14.5)));
  EXPECT_FLOAT_EQ(7.7935368191928e-45, exp(std_normal_lcdf(-14)));
  EXPECT_FLOAT_EQ(7.81880730565789e-42, exp(std_normal_lcdf(-13.5)));
  EXPECT_FLOAT_EQ(6.11716439954988e-39, exp(std_normal_lcdf(-13)));
  EXPECT_FLOAT_EQ(3.73256429887771e-36, exp(std_normal_lcdf(-12.5)));
  EXPECT_FLOAT_EQ(1.77648211207768e-33, exp(std_normal_lcdf(-12)));
  EXPECT_FLOAT_EQ(6.59577144611367e-31, exp(std_normal_lcdf(-11.5)));
  EXPECT_FLOAT_EQ(1.91065957449868e-28, exp(std_normal_lcdf(-11)));
  EXPECT_FLOAT_EQ(4.31900631780923e-26, exp(std_normal_lcdf(-10.5)));
  EXPECT_FLOAT_EQ(7.61985302416053e-24, exp(std_normal_lcdf(-10)));
  EXPECT_FLOAT_EQ(1.04945150753626e-21, exp(std_normal_lcdf(-9.5)));
  EXPECT_FLOAT_EQ(1.12858840595384e-19, exp(std_normal_lcdf(-9)));
  EXPECT_FLOAT_EQ(9.47953482220332e-18, exp(std_normal_lcdf(-8.5)));
  EXPECT_FLOAT_EQ(6.22096057427178e-16, exp(std_normal_lcdf(-8)));
  EXPECT_FLOAT_EQ(3.1908916729109e-14, exp(std_normal_lcdf(-7.5)));
  EXPECT_FLOAT_EQ(1.27981254388584e-12, exp(std_normal_lcdf(-7)));
  EXPECT_FLOAT_EQ(4.01600058385912e-11, exp(std_normal_lcdf(-6.5)));
  EXPECT_FLOAT_EQ(9.86587645037698e-10, exp(std_normal_lcdf(-6)));
  EXPECT_FLOAT_EQ(1.89895624658877e-08, exp(std_normal_lcdf(-5.5)));
  EXPECT_FLOAT_EQ(2.86651571879194e-07, exp(std_normal_lcdf(-5)));
  EXPECT_FLOAT_EQ(3.39767312473006e-06, exp(std_normal_lcdf(-4.5)));
  EXPECT_FLOAT_EQ(3.16712418331199e-05, exp(std_normal_lcdf(-4)));
  EXPECT_FLOAT_EQ(0.000232629079035525, exp(std_normal_lcdf(-3.5)));
  EXPECT_FLOAT_EQ(0.00134989803163009, exp(std_normal_lcdf(-3)));
  EXPECT_FLOAT_EQ(0.00620966532577613, exp(std_normal_lcdf(-2.5)));
  EXPECT_FLOAT_EQ(0.0227501319481792, exp(std_normal_lcdf(-2)));
  EXPECT_FLOAT_EQ(0.0668072012688581, exp(std_normal_lcdf(-1.5)));
  EXPECT_FLOAT_EQ(0.158655253931457, exp(std_normal_lcdf(-1)));
  EXPECT_FLOAT_EQ(0.308537538725987, exp(std_normal_lcdf(-0.5)));
  EXPECT_FLOAT_EQ(0.5, exp(std_normal_lcdf(0)));
  EXPECT_FLOAT_EQ(0.691462461274013, exp(std_normal_lcdf(0.5)));
  EXPECT_FLOAT_EQ(0.841344746068543, exp(std_normal_lcdf(1)));
  EXPECT_FLOAT_EQ(0.933192798731142, exp(std_normal_lcdf(1.5)));
  EXPECT_FLOAT_EQ(0.977249868051821, exp(std_normal_lcdf(2)));
  EXPECT_FLOAT_EQ(0.993790334674224, exp(std_normal_lcdf(2.5)));
  EXPECT_FLOAT_EQ(0.99865010196837, exp(std_normal_lcdf(3)));
  EXPECT_FLOAT_EQ(0.999767370920964, exp(std_normal_lcdf(3.5)));
  EXPECT_FLOAT_EQ(0.999968328758167, exp(std_normal_lcdf(4)));
  EXPECT_FLOAT_EQ(0.999996602326875, exp(std_normal_lcdf(4.5)));
  EXPECT_FLOAT_EQ(0.999999713348428, exp(std_normal_lcdf(5)));
  EXPECT_FLOAT_EQ(0.999999981010438, exp(std_normal_lcdf(5.5)));
  EXPECT_FLOAT_EQ(0.999999999013412, exp(std_normal_lcdf(6)));
  EXPECT_FLOAT_EQ(0.99999999995984, exp(std_normal_lcdf(6.5)));
  EXPECT_FLOAT_EQ(0.99999999999872, exp(std_normal_lcdf(7)));
  EXPECT_FLOAT_EQ(0.999999999999968, exp(std_normal_lcdf(7.5)));
  EXPECT_FLOAT_EQ(0.999999999999999, exp(std_normal_lcdf(8)));
  EXPECT_FLOAT_EQ(1, exp(std_normal_lcdf(8.5)));
  EXPECT_FLOAT_EQ(1, exp(std_normal_lcdf(9)));
  EXPECT_FLOAT_EQ(1, exp(std_normal_lcdf(9.5)));
  EXPECT_FLOAT_EQ(1, exp(std_normal_lcdf(10)));
}
