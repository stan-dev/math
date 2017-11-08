#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix,identity_m_exception_pass) {
  stan::math::matrix_gpu m;
  
  EXPECT_NO_THROW(m = stan::math::identity(1));
  EXPECT_NO_THROW(m = stan::math::identity(1));
  EXPECT_NO_THROW(m = stan::math::identity(1));

  stan::math::matrix_gpu m0;
  EXPECT_NO_THROW(m0 = stan::math::identity(0));
  
}

