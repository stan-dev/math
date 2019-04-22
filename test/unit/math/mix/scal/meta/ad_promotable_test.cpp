#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <string>

TEST(MathMeta, primitive_to_mix) {
  EXPECT_TRUE(
      (stan::math::ad_promotable<bool,
                                 stan::math::fvar<stan::math::var>>::value));
  EXPECT_TRUE(
      (stan::math::ad_promotable<char,
                                 stan::math::fvar<stan::math::var>>::value));
  EXPECT_TRUE(
      (stan::math::ad_promotable<unsigned char,
                                 stan::math::fvar<stan::math::var>>::value));
  EXPECT_TRUE(
      // NOLINTNEXTLINE(runtime/int)
      (stan::math::ad_promotable<short,
                                 stan::math::fvar<stan::math::var>>::value));
  EXPECT_TRUE(
      // NOLINTNEXTLINE(runtime/int)
      (stan::math::ad_promotable<unsigned short,
                                 stan::math::fvar<stan::math::var>>::value));
  EXPECT_TRUE(
      (stan::math::ad_promotable<int,
                                 stan::math::fvar<stan::math::var>>::value));
  EXPECT_TRUE(
      (stan::math::ad_promotable<unsigned int,
                                 stan::math::fvar<stan::math::var>>::value));
  EXPECT_TRUE(
      // NOLINTNEXTLINE(runtime/int)
      (stan::math::ad_promotable<long,
                                 stan::math::fvar<stan::math::var>>::value));
  EXPECT_TRUE(
      // NOLINTNEXTLINE(runtime/int)
      (stan::math::ad_promotable<unsigned long,
                                 stan::math::fvar<stan::math::var>>::value));
  EXPECT_TRUE(
      // NOLINTNEXTLINE(runtime/int)
      (stan::math::ad_promotable<long long,
                                 stan::math::fvar<stan::math::var>>::value));
  EXPECT_TRUE(
      // NOLINTNEXTLINE(runtime/int)
      (stan::math::ad_promotable<unsigned long long,
                                 stan::math::fvar<stan::math::var>>::value));
  EXPECT_TRUE(
      (stan::math::ad_promotable<float,
                                 stan::math::fvar<stan::math::var>>::value));
  EXPECT_TRUE(
      (stan::math::ad_promotable<double,
                                 stan::math::fvar<stan::math::var>>::value));
  EXPECT_TRUE(
      (stan::math::ad_promotable<long double,
                                 stan::math::fvar<stan::math::var>>::value));
}

TEST(MathMeta, rev_to_mix) {
  EXPECT_TRUE(
      (stan::math::ad_promotable<stan::math::var,
                                 stan::math::fvar<stan::math::var>>::value));
}

TEST(MathMeta, fwd_to_mix) {
  EXPECT_FALSE(
      (stan::math::ad_promotable<stan::math::fvar<double>,
                                 stan::math::fvar<stan::math::var>>::value));
}

TEST(MathMeta, nonprimitive_to_mix) {
  EXPECT_FALSE(
      (stan::math::ad_promotable<std::string,
                                 stan::math::fvar<stan::math::var>>::value));
}
