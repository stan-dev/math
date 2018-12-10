#include <stan/math/fwd/scal/meta/ad_promotable.hpp>
#include <gtest/gtest.h>
#include <string>

TEST(MathMeta, primitive_to_fvar) {
  EXPECT_TRUE(
      (stan::math::ad_promotable<bool, stan::math::fvar<double>>::value));
  EXPECT_TRUE(
      (stan::math::ad_promotable<char, stan::math::fvar<double>>::value));
  EXPECT_TRUE((stan::math::ad_promotable<unsigned char,
                                         stan::math::fvar<double>>::value));
  EXPECT_TRUE(
      // NOLINTNEXTLINE(runtime/int)
      (stan::math::ad_promotable<short, stan::math::fvar<double>>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_TRUE((stan::math::ad_promotable<unsigned short,
                                         stan::math::fvar<double>>::value));
  EXPECT_TRUE(
      (stan::math::ad_promotable<int, stan::math::fvar<double>>::value));
  EXPECT_TRUE((stan::math::ad_promotable<unsigned int,
                                         stan::math::fvar<double>>::value));
  EXPECT_TRUE(
      // NOLINTNEXTLINE(runtime/int)
      (stan::math::ad_promotable<long, stan::math::fvar<double>>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_TRUE((stan::math::ad_promotable<unsigned long,
                                         stan::math::fvar<double>>::value));
  EXPECT_TRUE(
      // NOLINTNEXTLINE(runtime/int)
      (stan::math::ad_promotable<long long, stan::math::fvar<double>>::value));
  // NOLINTNEXTLINE(runtime/int)
  EXPECT_TRUE((stan::math::ad_promotable<unsigned long long,
                                         stan::math::fvar<double>>::value));
  EXPECT_TRUE(
      (stan::math::ad_promotable<float, stan::math::fvar<double>>::value));
  EXPECT_TRUE(
      (stan::math::ad_promotable<double, stan::math::fvar<double>>::value));
  EXPECT_TRUE((
      stan::math::ad_promotable<long double, stan::math::fvar<double>>::value));
}

TEST(MathMeta, nonprimitive_to_fvar) {
  EXPECT_FALSE((
      stan::math::ad_promotable<std::string, stan::math::fvar<double>>::value));
}
