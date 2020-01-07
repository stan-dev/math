#include <gtest/gtest.h>
#include <stan/math/fwd/core.hpp>
#include <sstream>

// thorough tests for operators with fvar<double>, fvar<fvar<double>>,
// fvar<var>, and fvar<fvar<var>> are in test/unit/math/core/operator_*

TEST(mathFwdCoreFvar, copyCtor) {
  using stan::math::fvar;
  fvar<double> a(2, 1);
  fvar<double> b(a);
  EXPECT_FLOAT_EQ(a.val_, b.val_);
  EXPECT_FLOAT_EQ(a.d_, b.d_);
}

TEST(mathFwdCoreFvar, ctor) {
  using stan::math::fvar;

  fvar<double> a;
  EXPECT_FLOAT_EQ(0.0, a.val_);
  EXPECT_FLOAT_EQ(0.0, a.d_);

  fvar<double> b(1.9);
  EXPECT_FLOAT_EQ(1.9, b.val_);
  EXPECT_FLOAT_EQ(0.0, b.d_);

  fvar<double> c(1.93, -27.832);
  EXPECT_FLOAT_EQ(1.93, c.val_);
  EXPECT_FLOAT_EQ(-27.832, c.d_);
}

TEST(AgradFwdFvar, insertionOoperator) {
  using stan::math::fvar;
  fvar<double> a(5.0);
  std::stringstream ss;
  ss << a;
  EXPECT_EQ("5", ss.str());
}
