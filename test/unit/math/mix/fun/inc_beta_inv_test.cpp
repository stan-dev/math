#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>

TEST(ProbInternalMath, inc_beta_inv_fv1) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::inc_beta_inv;
  double a_d = 1;
  double b_d = 2;
  double p_d = 0.5;
  fvar<var> a_v = a_d;
  fvar<var> b_v = b_d;
  fvar<var> p_v = p_d;
  a_v.d_ = 1.0;
  b_v.d_ = 1.0;
  p_v.d_ = 1.0;

  fvar<var> res = inc_beta_inv(a_v, b_v, p_v);
  res.val_.grad();

  EXPECT_FLOAT_EQ(a_v.val_.adj(), 0.287698278597);
  EXPECT_FLOAT_EQ(b_v.val_.adj(), -0.122532267934);
  EXPECT_FLOAT_EQ(p_v.val_.adj(), 0.707106781187);

  a_v = a_d;
  b_v = b_d;
  p_v = p_d;
  a_v.d_ = 1.0;
  b_v.d_ = 1.0;
  p_v.d_ = 1.0;

  res = inc_beta_inv(a_d, b_v, p_v);
  res.val_.grad();

  EXPECT_FLOAT_EQ(b_v.val_.adj(), -0.122532267934);
  EXPECT_FLOAT_EQ(p_v.val_.adj(), 0.707106781187);

  b_v = b_d;
  p_v = p_d;
  b_v.d_ = 1.0;
  p_v.d_ = 1.0;

  res = inc_beta_inv(a_v, b_d, p_v);
  res.val_.grad();

  EXPECT_FLOAT_EQ(a_v.val_.adj(), 0.287698278597);
  EXPECT_FLOAT_EQ(p_v.val_.adj(), 0.707106781187);

  a_v = a_d;
  p_v = p_d;
  a_v.d_ = 1.0;
  p_v.d_ = 1.0;

  res = inc_beta_inv(a_v, b_v, p_d);
  res.val_.grad();

  EXPECT_FLOAT_EQ(a_v.val_.adj(), 0.287698278597);
  EXPECT_FLOAT_EQ(b_v.val_.adj(), -0.122532267934);
}

TEST(ProbInternalMath, inc_beta_inv_fv2) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::inc_beta_inv;
  fvar<fvar<var>> a = 2;
  fvar<fvar<var>> b = 5;
  fvar<fvar<var>> p = 0.1;
  a.d_ = 1.0;
  b.d_ = 1.0;
  p.d_ = 1.0;

  fvar<fvar<var>> res = inc_beta_inv(a, b, p);
  res.val_.val_.grad();

  EXPECT_FLOAT_EQ(a.val_.val_.adj(), 0.0783025374798);
  EXPECT_FLOAT_EQ(b.val_.val_.adj(), -0.0161882044585);
  EXPECT_FLOAT_EQ(p.val_.val_.adj(), 0.530989359806);
}
