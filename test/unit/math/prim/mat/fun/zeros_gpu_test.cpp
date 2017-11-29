#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixGPU,zero_m_exception_pass) {
  stan::math::matrix_gpu m(1,1);
  
  EXPECT_NO_THROW(m = stan::math::identity(1));
  EXPECT_NO_THROW(stan::math::zeros(m, stan::math::UPPER));
  EXPECT_NO_THROW(stan::math::zeros(m, stan::math::LOWER));

  stan::math::matrix_gpu m0;
  EXPECT_NO_THROW(m0 = stan::math::identity(0));
  
}
