#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/gpu/opencl_context.hpp>
#include <stan/math/gpu/matrix_gpu.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

TEST(MathMatrixGPU, matrix_gpu_creation) {
  stan::math::vector_d d1;
  stan::math::matrix_d d2;

  d1.resize(3);
  d2.resize(2, 3);
  EXPECT_NO_THROW(stan::math::matrix_gpu A(1, 1));
  EXPECT_NO_THROW(stan::math::matrix_gpu d11(d1));
  EXPECT_NO_THROW(stan::math::matrix_gpu d22(d2));
}

#endif
