#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <stdexcept>

TEST(MathMatrixPrimMat, singular_values) {
  using stan::math::matrix_d;
  using stan::math::singular_values;
  using stan::math::vector_d;

  matrix_d m0(0, 0);
  EXPECT_THROW(singular_values(m0), std::invalid_argument);

  matrix_d m1(1, 1);
  m1 << 1.0;
  EXPECT_NO_THROW(singular_values(m1));

  matrix_d m22(2, 2);
  m22 << 1, 9, -4, 2;
  vector_d m22_D(2);
  m22_D << 9.220341788859560, 4.121322275266775;
  EXPECT_MATRIX_FLOAT_EQ(m22_D, singular_values(m22));

  matrix_d m23(2, 3);
  m23 << 1, 3, -5, 7, 9, -11;
  vector_d m23_D(2);
  m23_D << 16.821011215675149, 1.747450051398016;
  EXPECT_MATRIX_FLOAT_EQ(m23_D, singular_values(m23));

  matrix_d m32(3, 2);
  m32 << 1, 3, -5, 7, 9, -11;
  vector_d m32_D(2);
  m32_D << 16.698998232964481, 2.672724829728879;
  EXPECT_MATRIX_FLOAT_EQ(m32_D, singular_values(m32));
}
