#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/fpclassify.hpp>
#include <test/unit/math/rev/scal/fun/util.hpp>
#include <limits>

TEST(AgradFwd, primitiveValueRevNested) {
  using stan::math::fvar;
  using stan::math::primitive_value;
  using stan::math::var;

  fvar<var> a = 5.2;
  EXPECT_FLOAT_EQ(5.2, primitive_value(a));

  // make sure all work together
  EXPECT_FLOAT_EQ(5.3, primitive_value(5.3));
  EXPECT_EQ(3, primitive_value(3));
}

TEST(AgradFwd, primitiveValueNanRevNested) {
  using stan::math::fvar;
  using stan::math::primitive_value;
  using stan::math::var;
  double nan = std::numeric_limits<double>::quiet_NaN();

  fvar<var> a = nan;
  EXPECT_TRUE(stan::math::is_nan(primitive_value(a)));
  EXPECT_TRUE(stan::math::is_nan(primitive_value(nan)));
}
