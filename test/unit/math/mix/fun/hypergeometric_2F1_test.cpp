#include <stan/math/rev/fun/hypergeometric_2F1.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, hypergeometric2F1_1) {
  using stan::math::fvar;
  fvar<double> a1 = 3.70975;
  fvar<double> a2 = 1;
  fvar<double> b = 2.70975;
  fvar<double> z = -0.2;

  a1.d_ = 1;
  a2.d_ = 1;
  b.d_ = 1;
  z.d_ = 1;

  fvar<double> res = stan::math::hypergeometric_2F1(a1, a2, b, z);

  EXPECT_FLOAT_EQ(-0.0488658806159776 - 0.193844936204681 + 0.0677809985598383
                      + 0.865295247272367,
                  res.d_);
}

TEST(mathMixScalFun, hypergeometric2F1_2) {
  using stan::math::fvar;
  using stan::math::var;
  fvar<var> a1 = 2;
  fvar<var> a2 = 1;
  fvar<var> b = 2;
  fvar<var> z = 0.4;

  fvar<var> res = stan::math::hypergeometric_2F1(a1, a2, b, z);
  res.val_.grad();

  EXPECT_FLOAT_EQ(0.461773432358295, a1.val().adj());
  EXPECT_FLOAT_EQ(0.851376039609984, a2.val().adj());
  EXPECT_FLOAT_EQ(-0.461773432358295, b.val().adj());
  EXPECT_FLOAT_EQ(2.77777777777778, z.val().adj());
}

TEST(mathMixScalFun, hypergeometric2F1_3_euler) {
  using stan::math::fvar;
  fvar<double> a1 = 1;
  fvar<double> a2 = 1;
  fvar<double> b = 2;
  fvar<double> z = -5;

  a1.d_ = 1;
  a2.d_ = 1;
  b.d_ = 1;
  z.d_ = 1;

  fvar<double> res = stan::math::hypergeometric_2F1(a1, a2, b, z);

  EXPECT_FLOAT_EQ(-0.321040199556840 - 0.321040199556840 + 0.129536268190289
                      + 0.0383370454357889,
                  res.d_);
}
