#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, modulus) {
  using stan::math::modulus;
  EXPECT_EQ(0, modulus(4, 2));
  EXPECT_EQ(0, modulus(6, 3));
  EXPECT_EQ(0, modulus(16, 4));
  EXPECT_EQ(0, modulus(34, 17));
  EXPECT_EQ(0, modulus(-4, 2));
  EXPECT_EQ(0, modulus(-6, 3));
  EXPECT_EQ(0, modulus(-16, 4));
  EXPECT_EQ(0, modulus(34, -17));
  EXPECT_EQ(0, modulus(4, -2));
  EXPECT_EQ(0, modulus(6, -3));
  EXPECT_EQ(0, modulus(16, -4));
  EXPECT_EQ(0, modulus(-34, 17));

  EXPECT_EQ(1, modulus(17, 4));
  EXPECT_EQ(1, modulus(22, 3));
  EXPECT_EQ(2, modulus(22, 4));
  EXPECT_EQ(6, modulus(34, 7));
  EXPECT_EQ(10, modulus(44, 17));

  EXPECT_EQ(1, modulus(17, -4));
  EXPECT_EQ(1, modulus(22, -3));
  EXPECT_EQ(2, modulus(22, -4));
  EXPECT_EQ(6, modulus(34, -7));
  EXPECT_EQ(10, modulus(44, -17));

  EXPECT_EQ(-1, modulus(-17, 4));
  EXPECT_EQ(-1, modulus(-22, 3));
  EXPECT_EQ(-2, modulus(-22, 4));
  EXPECT_EQ(-6, modulus(-34, 7));
  EXPECT_EQ(-10, modulus(-44, 17));
}

TEST(MathFunctions, int_modulus_0) {
  int x = 1;
  int y = 0;
  EXPECT_THROW(stan::math::modulus(x, y), std::domain_error);
}
