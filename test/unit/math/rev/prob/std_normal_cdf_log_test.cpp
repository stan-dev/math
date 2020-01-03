#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

TEST(std_normal_lcdf, tails) {
  using stan::math::std_normal_lcdf;
  using stan::math::var;
  using std::exp;

  EXPECT_FLOAT_EQ(
      1, 4.60535300958196e-308 / exp(std_normal_lcdf(var(-37.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  5.72557122252458e-300 / exp(std_normal_lcdf(var(-37)).val()));
  EXPECT_FLOAT_EQ(
      1, 5.54472571307484e-292 / exp(std_normal_lcdf(var(-36.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  4.18262406579728e-284 / exp(std_normal_lcdf(var(-36)).val()));
  EXPECT_FLOAT_EQ(
      1, 2.45769154066194e-276 / exp(std_normal_lcdf(var(-35.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  1.12491070647241e-268 / exp(std_normal_lcdf(var(-35)).val()));
  EXPECT_FLOAT_EQ(
      1, 4.01072896657726e-261 / exp(std_normal_lcdf(var(-34.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  1.11389878557438e-253 / exp(std_normal_lcdf(var(-34)).val()));
  EXPECT_FLOAT_EQ(
      1, 2.40983869512039e-246 / exp(std_normal_lcdf(var(-33.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  4.06118562091586e-239 / exp(std_normal_lcdf(var(-33)).val()));
  EXPECT_FLOAT_EQ(
      1, 5.33142435967881e-232 / exp(std_normal_lcdf(var(-32.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  5.4520806035124e-225 / exp(std_normal_lcdf(var(-32)).val()));
  EXPECT_FLOAT_EQ(
      1, 4.34323260103177e-218 / exp(std_normal_lcdf(var(-31.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  2.6952500812005e-211 / exp(std_normal_lcdf(var(-31)).val()));
  EXPECT_FLOAT_EQ(
      1, 1.30293791317808e-204 / exp(std_normal_lcdf(var(-30.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  4.90671392714819e-198 / exp(std_normal_lcdf(var(-30)).val()));
  EXPECT_FLOAT_EQ(
      1, 1.43947455222918e-191 / exp(std_normal_lcdf(var(-29.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  3.28978526670438e-185 / exp(std_normal_lcdf(var(-29)).val()));
  EXPECT_FLOAT_EQ(
      1, 5.85714125380634e-179 / exp(std_normal_lcdf(var(-28.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  8.12386946965943e-173 / exp(std_normal_lcdf(var(-28)).val()));
  EXPECT_FLOAT_EQ(
      1, 8.77817055687808e-167 / exp(std_normal_lcdf(var(-27.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  7.38948100688502e-161 / exp(std_normal_lcdf(var(-27)).val()));
  EXPECT_FLOAT_EQ(
      1, 4.84616266030332e-155 / exp(std_normal_lcdf(var(-26.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  2.47606331550339e-149 / exp(std_normal_lcdf(var(-26)).val()));
  EXPECT_FLOAT_EQ(
      1, 9.85623651896393e-144 / exp(std_normal_lcdf(var(-25.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  3.05669670638256e-138 / exp(std_normal_lcdf(var(-25)).val()));
  EXPECT_FLOAT_EQ(
      1, 7.38570686148941e-133 / exp(std_normal_lcdf(var(-24.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  1.3903921185497e-127 / exp(std_normal_lcdf(var(-24)).val()));
  EXPECT_FLOAT_EQ(
      1, 2.03936756324998e-122 / exp(std_normal_lcdf(var(-23.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  2.33063700622065e-117 / exp(std_normal_lcdf(var(-23)).val()));
  EXPECT_FLOAT_EQ(
      1, 2.07531079906636e-112 / exp(std_normal_lcdf(var(-22.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  1.43989243514508e-107 / exp(std_normal_lcdf(var(-22)).val()));
  EXPECT_FLOAT_EQ(
      1, 7.78439707718263e-103 / exp(std_normal_lcdf(var(-21.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  3.27927801897904e-98 / exp(std_normal_lcdf(var(-21)).val()));
  EXPECT_FLOAT_EQ(1,
                  1.0764673258791e-93 / exp(std_normal_lcdf(var(-20.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  2.75362411860623e-89 / exp(std_normal_lcdf(var(-20)).val()));
  EXPECT_FLOAT_EQ(
      1, 5.48911547566041e-85 / exp(std_normal_lcdf(var(-19.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  8.52722395263098e-81 / exp(std_normal_lcdf(var(-19)).val()));
  EXPECT_FLOAT_EQ(
      1, 1.03236986895633e-76 / exp(std_normal_lcdf(var(-18.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  9.74094891893715e-73 / exp(std_normal_lcdf(var(-18)).val()));
  EXPECT_FLOAT_EQ(
      1, 7.16345876623504e-69 / exp(std_normal_lcdf(var(-17.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  4.10599620209891e-65 / exp(std_normal_lcdf(var(-17)).val()));
  EXPECT_FLOAT_EQ(
      1, 1.83446300316473e-61 / exp(std_normal_lcdf(var(-16.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  6.38875440053809e-58 / exp(std_normal_lcdf(var(-16)).val()));
  EXPECT_FLOAT_EQ(
      1, 1.73446079179387e-54 / exp(std_normal_lcdf(var(-15.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  3.67096619931275e-51 / exp(std_normal_lcdf(var(-15)).val()));
  EXPECT_FLOAT_EQ(
      1, 6.05749476441522e-48 / exp(std_normal_lcdf(var(-14.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  7.7935368191928e-45 / exp(std_normal_lcdf(var(-14)).val()));
  EXPECT_FLOAT_EQ(
      1, 7.81880730565789e-42 / exp(std_normal_lcdf(var(-13.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  6.11716439954988e-39 / exp(std_normal_lcdf(var(-13)).val()));
  EXPECT_FLOAT_EQ(
      1, 3.73256429887771e-36 / exp(std_normal_lcdf(var(-12.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  1.77648211207768e-33 / exp(std_normal_lcdf(var(-12)).val()));
  EXPECT_FLOAT_EQ(
      1, 6.59577144611367e-31 / exp(std_normal_lcdf(var(-11.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  1.91065957449868e-28 / exp(std_normal_lcdf(var(-11)).val()));
  EXPECT_FLOAT_EQ(
      1, 4.31900631780923e-26 / exp(std_normal_lcdf(var(-10.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  7.61985302416053e-24 / exp(std_normal_lcdf(var(-10)).val()));
  EXPECT_FLOAT_EQ(1,
                  1.04945150753626e-21 / exp(std_normal_lcdf(var(-9.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  1.12858840595384e-19 / exp(std_normal_lcdf(var(-9)).val()));
  EXPECT_FLOAT_EQ(1,
                  9.47953482220332e-18 / exp(std_normal_lcdf(var(-8.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  6.22096057427178e-16 / exp(std_normal_lcdf(var(-8)).val()));
  EXPECT_FLOAT_EQ(1,
                  3.1908916729109e-14 / exp(std_normal_lcdf(var(-7.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  1.27981254388584e-12 / exp(std_normal_lcdf(var(-7)).val()));
  EXPECT_FLOAT_EQ(1,
                  4.01600058385912e-11 / exp(std_normal_lcdf(var(-6.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  9.86587645037698e-10 / exp(std_normal_lcdf(var(-6)).val()));
  EXPECT_FLOAT_EQ(1,
                  1.89895624658877e-08 / exp(std_normal_lcdf(var(-5.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  2.86651571879194e-07 / exp(std_normal_lcdf(var(-5)).val()));
  EXPECT_FLOAT_EQ(1,
                  3.39767312473006e-06 / exp(std_normal_lcdf(var(-4.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  3.16712418331199e-05 / exp(std_normal_lcdf(var(-4)).val()));
  EXPECT_FLOAT_EQ(1,
                  0.000232629079035525 / exp(std_normal_lcdf(var(-3.5)).val()));
  EXPECT_FLOAT_EQ(1, 0.00134989803163009 / exp(std_normal_lcdf(var(-3)).val()));
  EXPECT_FLOAT_EQ(1,
                  0.00620966532577613 / exp(std_normal_lcdf(var(-2.5)).val()));
  EXPECT_FLOAT_EQ(1, 0.0227501319481792 / exp(std_normal_lcdf(var(-2)).val()));
  EXPECT_FLOAT_EQ(1,
                  0.0668072012688581 / exp(std_normal_lcdf(var(-1.5)).val()));
  EXPECT_FLOAT_EQ(1, 0.158655253931457 / exp(std_normal_lcdf(var(-1)).val()));
  EXPECT_FLOAT_EQ(1, 0.308537538725987 / exp(std_normal_lcdf(var(-0.5)).val()));
  EXPECT_FLOAT_EQ(1, 0.5 / exp(std_normal_lcdf(var(0)).val()));
  EXPECT_FLOAT_EQ(1, 0.691462461274013 / exp(std_normal_lcdf(var(0.5)).val()));
  EXPECT_FLOAT_EQ(1, 0.841344746068543 / exp(std_normal_lcdf(var(1)).val()));
  EXPECT_FLOAT_EQ(1, 0.933192798731142 / exp(std_normal_lcdf(var(1.5)).val()));
  EXPECT_FLOAT_EQ(1, 0.977249868051821 / exp(std_normal_lcdf(var(2)).val()));
  EXPECT_FLOAT_EQ(1, 0.993790334674224 / exp(std_normal_lcdf(var(2.5)).val()));
  EXPECT_FLOAT_EQ(1, 0.99865010196837 / exp(std_normal_lcdf(var(3)).val()));
  EXPECT_FLOAT_EQ(1, 0.999767370920964 / exp(std_normal_lcdf(var(3.5)).val()));
  EXPECT_FLOAT_EQ(1, 0.999968328758167 / exp(std_normal_lcdf(var(4)).val()));
  EXPECT_FLOAT_EQ(1, 0.999996602326875 / exp(std_normal_lcdf(var(4.5)).val()));
  EXPECT_FLOAT_EQ(1, 0.999999713348428 / exp(std_normal_lcdf(var(5)).val()));
  EXPECT_FLOAT_EQ(1, 0.999999981010438 / exp(std_normal_lcdf(var(5.5)).val()));
  EXPECT_FLOAT_EQ(1, 0.999999999013412 / exp(std_normal_lcdf(var(6)).val()));
  EXPECT_FLOAT_EQ(1, 0.99999999995984 / exp(std_normal_lcdf(var(6.5)).val()));
  EXPECT_FLOAT_EQ(1, 0.99999999999872 / exp(std_normal_lcdf(var(7)).val()));
  EXPECT_FLOAT_EQ(1, 0.999999999999968 / exp(std_normal_lcdf(var(7.5)).val()));
  EXPECT_FLOAT_EQ(1, 0.999999999999999 / exp(std_normal_lcdf(var(8)).val()));
  EXPECT_FLOAT_EQ(1, 1 / exp(std_normal_lcdf(var(8.5)).val()));
  EXPECT_FLOAT_EQ(1, 1 / exp(std_normal_lcdf(var(9)).val()));
  EXPECT_FLOAT_EQ(1, 1 / exp(std_normal_lcdf(var(9.5)).val()));
  EXPECT_FLOAT_EQ(1, 1 / exp(std_normal_lcdf(var(10)).val()));

  stan::math::recover_memory();
}

void test_value_and_derivatives(double expected_val, double y_dbl) {
  using stan::math::is_nan;
  using stan::math::std_normal_lcdf;
  using stan::math::var;

  std::stringstream msg_ss;
  msg_ss << "for y = " << y_dbl;
  std::string msg = msg_ss.str();
  SCOPED_TRACE(msg);

  var y(y_dbl);

  std::vector<double> gradients;
  var val = std_normal_lcdf(y);
  std::vector<var> x;
  x.push_back(y);
  gradients.clear();
  val.grad(x, gradients);

  double e = 1e-10;
  double inv2e = 0.5 / e;
  double finite_diff
      = (std_normal_lcdf(y_dbl + e) - std_normal_lcdf(y_dbl - e)) * inv2e;

  EXPECT_FLOAT_EQ(expected_val, val.val());
  EXPECT_FALSE(is_nan(gradients[0]));

  // use relative error check in first case as logs of functions do not
  // necessarily converge to asymptotic values
  if (!is_nan(gradients[0])) {
    if (!is_nan(finite_diff)) {
      if (abs((finite_diff - gradients[0]) / finite_diff) > 0.01) {
        EXPECT_NEAR(finite_diff, gradients[0], 1e-2);
      }
    } else {
      EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), gradients[0]);
    }
  }
}

TEST(std_normal_lcdf, derivatives) {
  test_value_and_derivatives(log(0.5), 0.0);
  test_value_and_derivatives(-std::numeric_limits<double>::infinity(),
                             -1.18292e31);
  test_value_and_derivatives(0.0, 30.0);
}
