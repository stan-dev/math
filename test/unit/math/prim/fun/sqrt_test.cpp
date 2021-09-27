#include <stan/math.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, sqrtInt) {
  using stan::math::sqrt;
  EXPECT_FLOAT_EQ(std::sqrt(3.0), sqrt(3));
  EXPECT_FLOAT_EQ(std::sqrt(3.1), sqrt(3.1));
  EXPECT_TRUE(stan::math::is_nan(sqrt(-2)));

  uint32_t ulong = 1;
  uint64_t ulonglong = 1;
  long double ldouble = 1.5;
  EXPECT_FLOAT_EQ(std::sqrt(ulong), sqrt(ulong));
  EXPECT_FLOAT_EQ(std::sqrt(ulonglong), sqrt(ulonglong));
  EXPECT_FLOAT_EQ(std::sqrt(ldouble), sqrt(ldouble));
}
