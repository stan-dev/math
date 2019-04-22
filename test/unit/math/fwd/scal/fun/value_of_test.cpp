#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/fpclassify.hpp>
#include <limits>

TEST(AgradFwd, value_of) {
  using stan::math::fvar;
  using stan::math::value_of;
  using stan::math::value_of;

  fvar<double> a = 5.0;
  EXPECT_FLOAT_EQ(5.0, value_of(a));
  // make sure all work together
  EXPECT_FLOAT_EQ(5.0, value_of(5.0));
  EXPECT_FLOAT_EQ(5.0, value_of(5));
}

TEST(AgradFwd, value_of_nan) {
  using stan::math::fvar;
  using stan::math::value_of;
  using stan::math::value_of;
  double nan = std::numeric_limits<double>::quiet_NaN();

  fvar<double> a = nan;
  EXPECT_TRUE(stan::math::is_nan(value_of(a)));
  EXPECT_TRUE(stan::math::is_nan(value_of(nan)));
}
