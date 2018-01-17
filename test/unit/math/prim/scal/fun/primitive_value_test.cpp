#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, primitiveValue) {
  using stan::math::primitive_value;
  // int
  EXPECT_EQ(5, primitive_value(5));
  // uint
  EXPECT_EQ(5U, primitive_value(5U));
  // long >> int
  EXPECT_EQ(10000000000L, primitive_value(10000000000L));
  // char
  EXPECT_EQ('a', primitive_value('a'));

  // double
  EXPECT_EQ(7.3, primitive_value(7.3));
  // float
  EXPECT_EQ(7.3f, primitive_value(7.3f));
}

TEST(MathFunctions, primiviteValueNaN) {
  using boost::math::isnan;
  using stan::math::primitive_value;
  using std::numeric_limits;

  EXPECT_TRUE(
      isnan<double>(primitive_value(numeric_limits<double>::quiet_NaN())));
}
