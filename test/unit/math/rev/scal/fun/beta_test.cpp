#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRevBeta, beta_vv) {
  using stan::math::beta;
  using stan::math::lbeta;
  using stan::math::var;
  using std::exp;

  AVAR x = 5.2;
  AVAR y = 6.7;

  AVAR a = beta(x, y);
  a.grad();

  EXPECT_FLOAT_EQ(a.val(), 0.00043056878890093439172760842);
  EXPECT_FLOAT_EQ(x.adj(), -0.0003808378386711703475058416);
  EXPECT_FLOAT_EQ(y.adj(), -0.0002619172217554305098971095);

  AVAR i = 5.2;
  AVAR j = 6.7;

  AVAR k = exp(lbeta(i, j));
  k.grad();

  EXPECT_FLOAT_EQ(k.val(), 0.00043056878890093439172760842);
  EXPECT_FLOAT_EQ(i.adj(), -0.0003808378386711703475058416);
  EXPECT_FLOAT_EQ(j.adj(), -0.0002619172217554305098971095);
}

TEST(AgradRevBeta, beta_vd) {
  using stan::math::beta;
  using stan::math::lbeta;
  using stan::math::var;
  using std::exp;

  AVAR x = 3.1;
  double y = 2.5;

  AVAR a = beta(x, y);
  a.grad();
  EXPECT_FLOAT_EQ(a.val(), 0.04746059440167382140512361959);
  EXPECT_FLOAT_EQ(x.adj(), -0.0317656356629060817425389983);

  AVAR i = 3.1;
  double j = 2.5;

  AVAR k = exp(lbeta(i, j));
  k.grad();
  EXPECT_FLOAT_EQ(k.val(), 0.04746059440167382140512361959);
  EXPECT_FLOAT_EQ(i.adj(), -0.0317656356629060817425389983);
}

TEST(AgradRevBeta, beta_dv) {
  using stan::math::beta;
  using stan::math::lbeta;
  using stan::math::var;
  using std::exp;

  double x = 7.5;
  AVAR y = 1.8;
  AVAR a = beta(x, y);
  a.grad();
  EXPECT_FLOAT_EQ(a.val(), 0.02262404450446503812849214638);
  EXPECT_FLOAT_EQ(y.adj(), -0.0427661664280179724695518476);

  double i = 7.5;
  AVAR j = 1.8;
  AVAR k = exp(lbeta(i, j));
  k.grad();
  EXPECT_FLOAT_EQ(k.val(), 0.02262404450446503812849214638);
  EXPECT_FLOAT_EQ(j.adj(), -0.0427661664280179724695518476);
}
