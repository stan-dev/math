#include <stan/math/prim/scal/meta/ad_promotable.hpp>
#include <gtest/gtest.h>

TEST(MathMeta, primitive_to_double) {
  EXPECT_TRUE((stan::math::ad_promotable<bool, double>::value));
  EXPECT_TRUE((stan::math::ad_promotable<char, double>::value));
  EXPECT_TRUE((stan::math::ad_promotable<unsigned char, double>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_TRUE((stan::math::ad_promotable<short, double>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_TRUE((stan::math::ad_promotable<unsigned short, double>::value));
  EXPECT_TRUE((stan::math::ad_promotable<int, double>::value));
  EXPECT_TRUE((stan::math::ad_promotable<unsigned int, double>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_TRUE((stan::math::ad_promotable<long, double>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_TRUE((stan::math::ad_promotable<unsigned long, double>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_TRUE((stan::math::ad_promotable<long long, double>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_TRUE((stan::math::ad_promotable<unsigned long long, double>::value));
  EXPECT_TRUE((stan::math::ad_promotable<float, double>::value));
  EXPECT_TRUE((stan::math::ad_promotable<double, double>::value));
  EXPECT_TRUE((stan::math::ad_promotable<long double, double>::value));
}

TEST(MathMeta, primitive_to_float) {
  EXPECT_FALSE((stan::math::ad_promotable<bool, float>::value));
  EXPECT_FALSE((stan::math::ad_promotable<char, float>::value));
  EXPECT_FALSE((stan::math::ad_promotable<unsigned char, float>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_FALSE((stan::math::ad_promotable<short, float>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_FALSE((stan::math::ad_promotable<unsigned short, float>::value));
  EXPECT_FALSE((stan::math::ad_promotable<int, float>::value));
  EXPECT_FALSE((stan::math::ad_promotable<unsigned int, float>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_FALSE((stan::math::ad_promotable<long, float>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_FALSE((stan::math::ad_promotable<unsigned long, float>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_FALSE((stan::math::ad_promotable<long long, float>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_FALSE((stan::math::ad_promotable<unsigned long long, float>::value));
  EXPECT_FALSE((stan::math::ad_promotable<double, float>::value));
  EXPECT_FALSE((stan::math::ad_promotable<long double, float>::value));

  EXPECT_TRUE((stan::math::ad_promotable<float, float>::value))
      << "All primitive types should be promotable to the same type";
}
