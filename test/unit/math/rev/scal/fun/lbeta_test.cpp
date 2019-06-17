#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRevLbeta, lbeta_vv) {
  using stan::math::lbeta;
  using stan::math::var;

  AVAR x = 0.5;
  AVAR y = 1.2;

  AVAR a = lbeta(x, y);
  a.grad();
  EXPECT_FLOAT_EQ(a.val(), 0.5827985503284501018789171306);
  EXPECT_FLOAT_EQ(x.adj(), -2.1720579008949174361209405170);
  EXPECT_FLOAT_EQ(y.adj(), -0.4975877714656822522271721464);
}

TEST(AgradRevLbeta, lbeta_vd) {
  using stan::math::lbeta;
  using stan::math::var;

  AVAR x = 6.7;
  double y = 3.1;

  AVAR a = lbeta(x, y);
  a.grad();
  EXPECT_FLOAT_EQ(a.val(), -5.5417862655832854195801518779);
  EXPECT_FLOAT_EQ(x.adj(), -0.4048668189669833227198141126);
}

TEST(AgradRevLbeta, lbeta_dv) {
  using stan::math::lbeta;
  using stan::math::var;

  double x = 12.3;
  AVAR y = 4.8;

  AVAR a = lbeta(x, y);
  a.grad();
  EXPECT_FLOAT_EQ(a.val(), -9.8322070838104728837715023046);
  EXPECT_FLOAT_EQ(y.adj(), -1.3487060660277426588120358549);
}
