#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

TEST(MathMatrixGPU, matrix_cl_creation) {
  stan::math::vector_d d1;
  stan::math::matrix_d d2;
  stan::math::matrix_d d3;

  d1.resize(3);
  d2.resize(2, 3);
  EXPECT_NO_THROW(stan::math::matrix_cl A(1, 1));
  EXPECT_NO_THROW(stan::math::matrix_cl d11(d1));
  EXPECT_NO_THROW(stan::math::matrix_cl d22(d2));
  EXPECT_NO_THROW(stan::math::matrix_cl d33(d3));
}

#endif
