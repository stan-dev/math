#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>


TEST(MathMatrixGPU, copy_triangular) {
  stan::math::matrix_d d1;

  d1.resize(3,3);
  stan::math::matrix_gpu d11(d1);
  EXPECT_NO_THROW(stan::math::matrix_gpu d22 = 
   stan::math::copy_triangular(d11, stan::math::LOWER));
     EXPECT_NO_THROW(stan::math::matrix_gpu d33 = 
   stan::math::copy_triangular(d11, stan::math::UPPER));
}


TEST(MathMatrixGPU, copy_triangular_transpose) {
  stan::math::matrix_d d1;

  d1.resize(3,3);
  stan::math::matrix_gpu d11(d1);
  EXPECT_NO_THROW(stan::math::copy_triangular_transposed(d11,
   stan::math::LOWER_TO_UPPER_TRIANGULAR));
  EXPECT_NO_THROW(stan::math::copy_triangular_transposed(d11,
   stan::math::UPPER_TO_LOWER_TRIANGULAR));
}

TEST(MathMatrixGPU, copy_submatrix) {
  stan::math::matrix_d d1;

  d1.resize(3,3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(2,2);
  EXPECT_NO_THROW(stan::math::copy_submatrix(d11, d22, 1, 1, 0, 0, 2, 2));
}
