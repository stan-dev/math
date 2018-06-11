#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradFwdLogInvLogitDiff, FvarVar) {
  using stan::math::fvar;
  using stan::math::log_inv_logit_diff;
  using stan::math::var;

  fvar<var> x(0.5, 1.0);
  fvar<var> y(-1.0, 1.0);
  fvar<var> z = log_inv_logit_diff(x, y);
  z.d_.grad();

  EXPECT_FLOAT_EQ(0.664757585587, x.d_.adj());
  EXPECT_FLOAT_EQ(-0.556158338159, y.d_.adj());
  EXPECT_FLOAT_EQ(0.108599247428, z.d_.val());
}

TEST(AgradFwdLogInvLogitDiff, FvarVar_Dbl) {
  using stan::math::fvar;
  using stan::math::log_inv_logit_diff;
  using stan::math::var;

  fvar<var> x(0.5, 1.0);
  double y = -1.0;
  fvar<var> z = log_inv_logit_diff(x, y);
  z.d_.grad();

  EXPECT_FLOAT_EQ(0.664757585587, x.d_.adj());
  EXPECT_FLOAT_EQ(x.d_.adj(), z.d_.val());

  double a = 0.5;
  fvar<var> b(-1.0, 1.0);
  fvar<var> c = log_inv_logit_diff(a, b);
  c.d_.grad();

  EXPECT_FLOAT_EQ(-0.556158338159, b.d_.adj());
  EXPECT_FLOAT_EQ(b.d_.adj(), c.d_.val());
}

TEST(AgradFwdLogInvLogitDiff, FvarFvarVar) {
  using stan::math::fvar;
  using stan::math::log_inv_logit_diff;
  using stan::math::var;

  fvar<fvar<var>> x(0.5, 1.0);
  fvar<fvar<var>> y(-1.0, 1.0);
  fvar<fvar<var>> z = log_inv_logit_diff(x, y);
  z.d_.val_.grad();

  EXPECT_FLOAT_EQ(0.664757585587, x.d_.val_.adj());
  EXPECT_FLOAT_EQ(-0.556158338159, y.d_.val_.adj());
  EXPECT_FLOAT_EQ(0.108599247428, z.d_.val_.val());
}

TEST(AgradFwdLogInvLogitDiff, FvarFVarVar_Dbl) {
  using stan::math::fvar;
  using stan::math::log_inv_logit_diff;
  using stan::math::var;

  fvar<fvar<var>> x(0.5, 1.0);
  double y = -1.0;
  fvar<fvar<var>> z = log_inv_logit_diff(x, y);
  z.d_.val_.grad();

  EXPECT_FLOAT_EQ(0.664757585587, x.d_.val_.adj());
  EXPECT_FLOAT_EQ(x.d_.val_.adj(), z.d_.val_.val());

  double a = 0.5;
  fvar<fvar<var>> b(-1.0, 1.0);
  fvar<fvar<var>> c = log_inv_logit_diff(a, b);
  c.d_.val_.grad();

  EXPECT_FLOAT_EQ(-0.556158338159, b.d_.val_.adj());
  EXPECT_FLOAT_EQ(b.d_.val_.adj(), c.d_.val_.val());
}
