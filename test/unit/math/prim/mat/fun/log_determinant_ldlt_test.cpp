#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, log_deterimant_ldlt) {
  using stan::math::determinant;
  using std::log;
  using std::fabs;
  
  stan::math::matrix_d x(2,2);  
  stan::math::LDLT_factor<double,-1,-1> ldlt_x;

  x << 2, 1, 1, 3;
  ldlt_x.compute(x);
  ASSERT_TRUE(ldlt_x.success());
  
  EXPECT_FLOAT_EQ(log(fabs(determinant(x))),
                  stan::math::log_determinant_ldlt(ldlt_x));

  x << 1, 0, 0, 3;
  ldlt_x.compute(x);
  ASSERT_TRUE(ldlt_x.success());
  EXPECT_FLOAT_EQ(log(3.0),
                  stan::math::log_determinant_ldlt(ldlt_x));
}
