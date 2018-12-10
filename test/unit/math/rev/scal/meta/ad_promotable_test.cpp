#include <stan/math/rev/scal/meta/ad_promotable.hpp>
#include <gtest/gtest.h>
#include <string>

TEST(MathMeta, primitive_to_var) {
  EXPECT_TRUE((stan::math::ad_promotable<bool, stan::math::var>::value));
  EXPECT_TRUE((stan::math::ad_promotable<char, stan::math::var>::value));
  EXPECT_TRUE(
      (stan::math::ad_promotable<unsigned char, stan::math::var>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_TRUE((stan::math::ad_promotable<short, stan::math::var>::value));
  EXPECT_TRUE(
      // NOLINTNEXTLINE(runtime/int)
      (stan::math::ad_promotable<unsigned short, stan::math::var>::value));
  EXPECT_TRUE((stan::math::ad_promotable<int, stan::math::var>::value));
  EXPECT_TRUE(
      (stan::math::ad_promotable<unsigned int, stan::math::var>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_TRUE((stan::math::ad_promotable<long, stan::math::var>::value));
  EXPECT_TRUE(
      // NOLINTNEXTLINE(runtime/int)
      (stan::math::ad_promotable<unsigned long, stan::math::var>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_TRUE((stan::math::ad_promotable<long long, stan::math::var>::value));
  EXPECT_TRUE(
      // NOLINTNEXTLINE(runtime/int)
      (stan::math::ad_promotable<unsigned long long, stan::math::var>::value));
  EXPECT_TRUE((stan::math::ad_promotable<float, stan::math::var>::value));
  EXPECT_TRUE((stan::math::ad_promotable<double, stan::math::var>::value));
  EXPECT_TRUE((stan::math::ad_promotable<long double, stan::math::var>::value));
}

TEST(MathMeta, nonprimitive_to_var) {
  EXPECT_FALSE(
      (stan::math::ad_promotable<std::string, stan::math::var>::value));
}
