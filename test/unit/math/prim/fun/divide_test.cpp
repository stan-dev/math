
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>




TEST(MathFunctions, divide) {
  using stan::math::divide;
  EXPECT_EQ(2, divide(4, 2));
  EXPECT_EQ(2, divide(6, 3));
  EXPECT_EQ(4, divide(16, 4));
  EXPECT_EQ(2, divide(34, 17));
  EXPECT_EQ(-2, divide(-4, 2));
  EXPECT_EQ(-2, divide(-6, 3));
  EXPECT_EQ(-4, divide(-16, 4));
  EXPECT_EQ(-2, divide(34, -17));
  EXPECT_EQ(-2, divide(4, -2));
  EXPECT_EQ(-2, divide(6, -3));
  EXPECT_EQ(-4, divide(16, -4));
  EXPECT_EQ(-2, divide(-34, 17));

  EXPECT_EQ(4, divide(17, 4));
  EXPECT_EQ(7, divide(22, 3));
  EXPECT_EQ(5, divide(22, 4));
  EXPECT_EQ(4, divide(34, 7));
  EXPECT_EQ(2, divide(44, 17));

  EXPECT_EQ(-4, divide(17, -4));
  EXPECT_EQ(-7, divide(22, -3));
  EXPECT_EQ(-5, divide(22, -4));
  EXPECT_EQ(-4, divide(34, -7));
  EXPECT_EQ(-2, divide(44, -17));

  EXPECT_EQ(-4, divide(-17, 4));
  EXPECT_EQ(-7, divide(-22, 3));
  EXPECT_EQ(-5, divide(-22, 4));
  EXPECT_EQ(-4, divide(-34, 7));
  EXPECT_EQ(-2, divide(-44, 17));
}

void test_divide_modulus(int a, int b) {
  using stan::math::divide;
  using stan::math::modulus;
  EXPECT_EQ(a, divide(a, b) * b + modulus(a, b));
}

TEST(MathFunctions, divide_modulus) {
  for (int i = 1; i < 50; i++)
    for (int j = 1; j < 50; j++)
      test_divide_modulus(i, j);
}

TEST(MathFunctions, int_divide_by_0) {
  int x = 1;
  int y = 0;
  EXPECT_THROW(stan::math::divide(x, y), std::domain_error);
}



TEST(MathMatrix_mat, divide) {
  stan::math::vector_d v0;
  stan::math::row_vector_d rv0;
  stan::math::matrix_d m0;

  using stan::math::divide;
  EXPECT_NO_THROW(divide(v0, 2.0));
  EXPECT_NO_THROW(divide(rv0, 2.0));
  EXPECT_NO_THROW(divide(m0, 2.0));
}
