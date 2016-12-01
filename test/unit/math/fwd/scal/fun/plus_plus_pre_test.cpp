#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>

typedef stan::math::fvar<double> fd_t;
typedef stan::math::fvar<fd_t> ffd_t;

TEST(MathFunctions, plusPlusPreFvarValue) {
  using stan::math::plus_plus_pre;
  fd_t x = 4;
  fd_t y = plus_plus_pre(x);
  EXPECT_FLOAT_EQ(5, y.val_);
  EXPECT_FLOAT_EQ(5, x.val_);
}

TEST(MathFunctions, plusPlusPreFvarDeriv) {
  using stan::math::plus_plus_pre;
  fd_t x(4, 2.3);
  fd_t y = plus_plus_pre(x);
  EXPECT_FLOAT_EQ(2.3, y.d_);
}

TEST(MathFunctions, plusPlusPreFvarRef) {
  using stan::math::plus_plus_pre;
  fd_t x(4, 2.3);
  EXPECT_EQ(&x, &plus_plus_pre(x));
}


TEST(MathFunctions, plusPlusPreF2varValue) {
  using stan::math::plus_plus_pre;
  ffd_t x(4, 2.3);
  ffd_t y = plus_plus_pre(x);
  EXPECT_FLOAT_EQ(5, y.val_.val_);
  EXPECT_FLOAT_EQ(5, x.val_.val_);
}

TEST(MathFunctions, plusPlusPreF2varDeriv) {
  using stan::math::plus_plus_pre;
  ffd_t x(4, 2.3);
  ffd_t y = plus_plus_pre(x);
  EXPECT_FLOAT_EQ(2.3, y.d_.val_);
  EXPECT_FLOAT_EQ(2.3, x.d_.val_);
}

TEST(MathFunctions, plusPlusPreF2varRef) {
  using stan::math::plus_plus_pre;
  ffd_t x(4, 2.3);
  EXPECT_EQ(&x, &plus_plus_pre(x));
}
