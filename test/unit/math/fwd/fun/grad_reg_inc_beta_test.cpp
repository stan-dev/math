#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>

TEST(ProbInternalMath, grad_reg_inc_beta_fd) {
  using stan::math::digamma;
  using stan::math::exp;
  using stan::math::fvar;
  using stan::math::lbeta;

  fvar<double> a = 1.0;
  fvar<double> b = 1.0;
  fvar<double> g = 0.4;
  a.d_ = 1.0;
  b.d_ = 1.0;
  g.d_ = 1.0;
  fvar<double> dig_a = digamma(a);
  fvar<double> dig_b = digamma(b);
  fvar<double> dig_sum = digamma(a + b);
  fvar<double> beta_ab = exp(lbeta(a, b));
  fvar<double> g_a;
  fvar<double> g_b;

  stan::math::grad_reg_inc_beta(g_a, g_b, a, b, g, dig_a, dig_b, dig_sum,
                                beta_ab);
  EXPECT_FLOAT_EQ(-0.36651629442883944183907601651838247842001142107486495485,
                  g_a.val_);
  EXPECT_NEAR(0.306495375042422864944011633197968575202046200428315551199,
              g_b.val_, 1e-6);
}

TEST(ProbInternalMath, grad_reg_inc_beta_ffd) {
  using stan::math::digamma;
  using stan::math::exp;
  using stan::math::fvar;
  using stan::math::lbeta;

  fvar<fvar<double> > a = 1.0;
  fvar<fvar<double> > b = 1.0;
  fvar<fvar<double> > g = 0.4;
  a.d_ = 1.0;
  b.d_ = 1.0;
  g.d_ = 1.0;
  fvar<fvar<double> > dig_a = digamma(a);
  fvar<fvar<double> > dig_b = digamma(b);
  fvar<fvar<double> > dig_sum = digamma(a + b);
  fvar<fvar<double> > beta_ab = exp(lbeta(a, b));
  fvar<fvar<double> > g_a;
  fvar<fvar<double> > g_b;

  stan::math::grad_reg_inc_beta(g_a, g_b, a, b, g, dig_a, dig_b, dig_sum,
                                beta_ab);
  EXPECT_FLOAT_EQ(-0.36651629442883944183907601651838247842001142107486495485,
                  g_a.val_.val_);
  EXPECT_NEAR(0.306495375042422864944011633197968575202046200428315551199,
              g_b.val_.val_, 1e-6);
}
