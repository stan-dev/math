#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, log_determinant_ldlt) {
  using stan::math::determinant;
  using std::fabs;
  using std::log;

  stan::math::matrix_d x(2, 2);
  {
    x << 2, 1, 1, 3;
    auto ldlt_x = stan::math::make_ldlt_factor(x);
    EXPECT_FLOAT_EQ(log(fabs(determinant(x))),
                    stan::math::log_determinant_ldlt(ldlt_x));
  }

  {
    x << 1, 0, 0, 3;
    auto ldlt_x = stan::math::make_ldlt_factor(x);
    EXPECT_FLOAT_EQ(log(3.0), stan::math::log_determinant_ldlt(ldlt_x));
  }
}

TEST(MathMatrixPrimMat, log_determinant_ldlt_0x0) {
  using stan::math::determinant;
  using std::fabs;
  using std::log;

  stan::math::matrix_d x(0, 0);
  auto ldlt_x = stan::math::make_ldlt_factor(x);

  EXPECT_FLOAT_EQ(log(fabs(determinant(x))),
                  stan::math::log_determinant_ldlt(ldlt_x));
}
