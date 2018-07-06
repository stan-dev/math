#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/gpu/copy_submatrix_opencl.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixGPU, copy_submatrix_exception) {
  stan::math::matrix_d d1;

  d1.resize(3, 3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(2, 2);
  EXPECT_NO_THROW(stan::math::copy_submatrix(d11, d22, 1, 1, 0, 0, 2, 2));

  EXPECT_THROW(stan::math::copy_submatrix(d11, d22, 1, 1, 0, 0, 4, 4),
               std::domain_error);
  EXPECT_THROW(stan::math::copy_submatrix(d11, d22, 4, 4, 0, 0, 2, 2),
               std::domain_error);
  EXPECT_THROW(stan::math::copy_submatrix(d11, d22, 1, 1, 3, 3, 4, 4),
               std::domain_error);
  EXPECT_THROW(stan::math::copy_submatrix(d11, d22, 1, 1, 3, 3, 3, 3),
               std::domain_error);
}

#endif
