#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <string>
#include <vector>

TEST(std_normal_cdf, tails) {
  using stan::math::std_normal_cdf;
  using stan::math::var;

  EXPECT_FLOAT_EQ(1, 4.60535300958196e-308 / std_normal_cdf(var(-37.5)).val());
  EXPECT_FLOAT_EQ(1, 5.72557122252458e-300 / std_normal_cdf(var(-37)).val());
  EXPECT_FLOAT_EQ(1, 5.54472571307484e-292 / std_normal_cdf(var(-36.5)).val());
  EXPECT_FLOAT_EQ(1, 4.18262406579728e-284 / std_normal_cdf(var(-36)).val());
  EXPECT_FLOAT_EQ(1, 2.45769154066194e-276 / std_normal_cdf(var(-35.5)).val());
  EXPECT_FLOAT_EQ(1, 1.12491070647241e-268 / std_normal_cdf(var(-35)).val());
  EXPECT_FLOAT_EQ(1, 4.01072896657726e-261 / std_normal_cdf(var(-34.5)).val());
  EXPECT_FLOAT_EQ(1, 1.11389878557438e-253 / std_normal_cdf(var(-34)).val());
  EXPECT_FLOAT_EQ(1, 2.40983869512039e-246 / std_normal_cdf(var(-33.5)).val());
  EXPECT_FLOAT_EQ(1, 4.06118562091586e-239 / std_normal_cdf(var(-33)).val());
  EXPECT_FLOAT_EQ(1, 5.33142435967881e-232 / std_normal_cdf(var(-32.5)).val());
  EXPECT_FLOAT_EQ(1, 5.4520806035124e-225 / std_normal_cdf(var(-32)).val());
  EXPECT_FLOAT_EQ(1, 4.34323260103177e-218 / std_normal_cdf(var(-31.5)).val());
  EXPECT_FLOAT_EQ(1, 2.6952500812005e-211 / std_normal_cdf(var(-31)).val());
  EXPECT_FLOAT_EQ(1, 1.30293791317808e-204 / std_normal_cdf(var(-30.5)).val());
  EXPECT_FLOAT_EQ(1, 4.90671392714819e-198 / std_normal_cdf(var(-30)).val());
  EXPECT_FLOAT_EQ(1, 1.43947455222918e-191 / std_normal_cdf(var(-29.5)).val());
  EXPECT_FLOAT_EQ(1, 3.28978526670438e-185 / std_normal_cdf(var(-29)).val());
  EXPECT_FLOAT_EQ(1, 5.85714125380634e-179 / std_normal_cdf(var(-28.5)).val());
  EXPECT_FLOAT_EQ(1, 8.12386946965943e-173 / std_normal_cdf(var(-28)).val());
  EXPECT_FLOAT_EQ(1, 8.77817055687808e-167 / std_normal_cdf(var(-27.5)).val());
  EXPECT_FLOAT_EQ(1, 7.38948100688502e-161 / std_normal_cdf(var(-27)).val());
  EXPECT_FLOAT_EQ(1, 4.84616266030332e-155 / std_normal_cdf(var(-26.5)).val());
  EXPECT_FLOAT_EQ(1, 2.47606331550339e-149 / std_normal_cdf(var(-26)).val());
  EXPECT_FLOAT_EQ(1, 9.85623651896393e-144 / std_normal_cdf(var(-25.5)).val());
  EXPECT_FLOAT_EQ(1, 3.05669670638256e-138 / std_normal_cdf(var(-25)).val());
  EXPECT_FLOAT_EQ(1, 7.38570686148941e-133 / std_normal_cdf(var(-24.5)).val());
  EXPECT_FLOAT_EQ(1, 1.3903921185497e-127 / std_normal_cdf(var(-24)).val());
  EXPECT_FLOAT_EQ(1, 2.03936756324998e-122 / std_normal_cdf(var(-23.5)).val());
  EXPECT_FLOAT_EQ(1, 2.33063700622065e-117 / std_normal_cdf(var(-23)).val());
  EXPECT_FLOAT_EQ(1, 2.07531079906636e-112 / std_normal_cdf(var(-22.5)).val());
  EXPECT_FLOAT_EQ(1, 1.43989243514508e-107 / std_normal_cdf(var(-22)).val());
  EXPECT_FLOAT_EQ(1, 7.78439707718263e-103 / std_normal_cdf(var(-21.5)).val());
  EXPECT_FLOAT_EQ(1, 3.27927801897904e-98 / std_normal_cdf(var(-21)).val());
  EXPECT_FLOAT_EQ(1, 1.0764673258791e-93 / std_normal_cdf(var(-20.5)).val());
  EXPECT_FLOAT_EQ(1, 2.75362411860623e-89 / std_normal_cdf(var(-20)).val());
  EXPECT_FLOAT_EQ(1, 5.48911547566041e-85 / std_normal_cdf(var(-19.5)).val());
  EXPECT_FLOAT_EQ(1, 8.52722395263098e-81 / std_normal_cdf(var(-19)).val());
  EXPECT_FLOAT_EQ(1, 1.03236986895633e-76 / std_normal_cdf(var(-18.5)).val());
  EXPECT_FLOAT_EQ(1, 9.74094891893715e-73 / std_normal_cdf(var(-18)).val());
  EXPECT_FLOAT_EQ(1, 7.16345876623504e-69 / std_normal_cdf(var(-17.5)).val());
  EXPECT_FLOAT_EQ(1, 4.10599620209891e-65 / std_normal_cdf(var(-17)).val());
  EXPECT_FLOAT_EQ(1, 1.83446300316473e-61 / std_normal_cdf(var(-16.5)).val());
  EXPECT_FLOAT_EQ(1, 6.38875440053809e-58 / std_normal_cdf(var(-16)).val());
  EXPECT_FLOAT_EQ(1, 1.73446079179387e-54 / std_normal_cdf(var(-15.5)).val());
  EXPECT_FLOAT_EQ(1, 3.67096619931275e-51 / std_normal_cdf(var(-15)).val());
  EXPECT_FLOAT_EQ(1, 6.05749476441522e-48 / std_normal_cdf(var(-14.5)).val());
  EXPECT_FLOAT_EQ(1, 7.7935368191928e-45 / std_normal_cdf(var(-14)).val());
  EXPECT_FLOAT_EQ(1, 7.81880730565789e-42 / std_normal_cdf(var(-13.5)).val());
  EXPECT_FLOAT_EQ(1, 6.11716439954988e-39 / std_normal_cdf(var(-13)).val());
  EXPECT_FLOAT_EQ(1, 3.73256429887771e-36 / std_normal_cdf(var(-12.5)).val());
  EXPECT_FLOAT_EQ(1, 1.77648211207768e-33 / std_normal_cdf(var(-12)).val());
  EXPECT_FLOAT_EQ(1, 6.59577144611367e-31 / std_normal_cdf(var(-11.5)).val());
  EXPECT_FLOAT_EQ(1, 1.91065957449868e-28 / std_normal_cdf(var(-11)).val());
  EXPECT_FLOAT_EQ(1, 4.31900631780923e-26 / std_normal_cdf(var(-10.5)).val());
  EXPECT_FLOAT_EQ(1, 7.61985302416053e-24 / std_normal_cdf(var(-10)).val());
  EXPECT_FLOAT_EQ(1, 1.04945150753626e-21 / std_normal_cdf(var(-9.5)).val());
  EXPECT_FLOAT_EQ(1, 1.12858840595384e-19 / std_normal_cdf(var(-9)).val());
  EXPECT_FLOAT_EQ(1, 9.47953482220332e-18 / std_normal_cdf(var(-8.5)).val());
  EXPECT_FLOAT_EQ(1, 6.22096057427178e-16 / std_normal_cdf(var(-8)).val());
  EXPECT_FLOAT_EQ(1, 3.1908916729109e-14 / std_normal_cdf(var(-7.5)).val());
  EXPECT_FLOAT_EQ(1, 1.27981254388584e-12 / std_normal_cdf(var(-7)).val());
  EXPECT_FLOAT_EQ(1, 4.01600058385912e-11 / std_normal_cdf(var(-6.5)).val());
  EXPECT_FLOAT_EQ(1, 9.86587645037698e-10 / std_normal_cdf(var(-6)).val());
  EXPECT_FLOAT_EQ(1, 1.89895624658877e-08 / std_normal_cdf(var(-5.5)).val());
  EXPECT_FLOAT_EQ(1, 2.86651571879194e-07 / std_normal_cdf(var(-5)).val());
  EXPECT_FLOAT_EQ(1, 3.39767312473006e-06 / std_normal_cdf(var(-4.5)).val());
  EXPECT_FLOAT_EQ(1, 3.16712418331199e-05 / std_normal_cdf(var(-4)).val());
  EXPECT_FLOAT_EQ(1, 0.000232629079035525 / std_normal_cdf(var(-3.5)).val());
  EXPECT_FLOAT_EQ(1, 0.00134989803163009 / std_normal_cdf(var(-3)).val());
  EXPECT_FLOAT_EQ(1, 0.00620966532577613 / std_normal_cdf(var(-2.5)).val());
  EXPECT_FLOAT_EQ(1, 0.0227501319481792 / std_normal_cdf(var(-2)).val());
  EXPECT_FLOAT_EQ(1, 0.0668072012688581 / std_normal_cdf(var(-1.5)).val());
  EXPECT_FLOAT_EQ(1, 0.158655253931457 / std_normal_cdf(var(-1)).val());
  EXPECT_FLOAT_EQ(1, 0.308537538725987 / std_normal_cdf(var(-0.5)).val());
  EXPECT_FLOAT_EQ(1, 0.5 / std_normal_cdf(var(0)).val());
  EXPECT_FLOAT_EQ(1, 0.691462461274013 / std_normal_cdf(var(0.5)).val());
  EXPECT_FLOAT_EQ(1, 0.841344746068543 / std_normal_cdf(var(1)).val());
  EXPECT_FLOAT_EQ(1, 0.933192798731142 / std_normal_cdf(var(1.5)).val());
  EXPECT_FLOAT_EQ(1, 0.977249868051821 / std_normal_cdf(var(2)).val());
  EXPECT_FLOAT_EQ(1, 0.993790334674224 / std_normal_cdf(var(2.5)).val());
  EXPECT_FLOAT_EQ(1, 0.99865010196837 / std_normal_cdf(var(3)).val());
  EXPECT_FLOAT_EQ(1, 0.999767370920964 / std_normal_cdf(var(3.5)).val());
  EXPECT_FLOAT_EQ(1, 0.999968328758167 / std_normal_cdf(var(4)).val());
  EXPECT_FLOAT_EQ(1, 0.999996602326875 / std_normal_cdf(var(4.5)).val());
  EXPECT_FLOAT_EQ(1, 0.999999713348428 / std_normal_cdf(var(5)).val());
  EXPECT_FLOAT_EQ(1, 0.999999981010438 / std_normal_cdf(var(5.5)).val());
  EXPECT_FLOAT_EQ(1, 0.999999999013412 / std_normal_cdf(var(6)).val());
  EXPECT_FLOAT_EQ(1, 0.99999999995984 / std_normal_cdf(var(6.5)).val());
  EXPECT_FLOAT_EQ(1, 0.99999999999872 / std_normal_cdf(var(7)).val());
  EXPECT_FLOAT_EQ(1, 0.999999999999968 / std_normal_cdf(var(7.5)).val());
  EXPECT_FLOAT_EQ(1, 0.999999999999999 / std_normal_cdf(var(8)).val());
  EXPECT_FLOAT_EQ(1, 1 / std_normal_cdf(var(8.5)).val());
  EXPECT_FLOAT_EQ(1, 1 / std_normal_cdf(var(9)).val());
  EXPECT_FLOAT_EQ(1, 1 / std_normal_cdf(var(9.5)).val());
  EXPECT_FLOAT_EQ(1, 1 / std_normal_cdf(var(10)).val());

  stan::math::recover_memory();
}

void test_value_and_derivatives(double expected_val, double y_dbl) {
  using stan::math::is_nan;
  using stan::math::std_normal_cdf;
  using stan::math::var;

  std::stringstream msg_ss;
  msg_ss << "parameters: (" << y_dbl << ")";
  std::string msg = msg_ss.str();
  SCOPED_TRACE(msg);

  var y(y_dbl);
  std::vector<double> gradients;
  var val = std_normal_cdf(y);
  std::vector<var> x;
  x.push_back(y);
  gradients.clear();
  val.grad(x, gradients);

  double e = 1e-10;
  double inv2e = 0.5 / e;
  double finite_diff
      = (std_normal_cdf(y_dbl + e) - std_normal_cdf(y_dbl - e)) * inv2e;

  EXPECT_FLOAT_EQ(expected_val, val.val());
  EXPECT_FALSE(is_nan(gradients[0]));
  if (!is_nan(gradients[0])) {
    if (!is_nan(finite_diff))
      EXPECT_NEAR(finite_diff, gradients[0], 1e-2);
    else
      EXPECT_FLOAT_EQ(0.0, gradients[0]);
  }
}

TEST(std_normal_cdf, derivatives) {
  test_value_and_derivatives(0.5, 0.0);
  test_value_and_derivatives(0.0, -20.0);
  test_value_and_derivatives(0.0, -30.0);
  test_value_and_derivatives(0.0, -50.0);
  test_value_and_derivatives(1.0, 30.0);
}
