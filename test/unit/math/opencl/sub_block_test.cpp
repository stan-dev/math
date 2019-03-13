#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixGPU, sub_block_pass) {
  stan::math::matrix_d d1;
  stan::math::matrix_d d2;

  d1.resize(3, 3);
  d2.resize(4, 4);

  d1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  d2 << 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1;

  stan::math::matrix_cl d11(d1);
  stan::math::matrix_cl d22(d2);
  d22.sub_block(d11, 0, 0, 0, 0, 2, 2);
  stan::math::copy(d2, d22);
  EXPECT_EQ(1, d2(0, 0));
  EXPECT_EQ(2, d2(0, 1));
  EXPECT_EQ(4, d2(1, 0));
  EXPECT_EQ(5, d2(1, 1));
}

TEST(MathMatrixGPU, sub_block_exception) {
  stan::math::matrix_d d1;
  stan::math::matrix_d d2;

  d1.resize(3, 3);
  d2.resize(4, 4);
  stan::math::matrix_cl d11(d1);
  stan::math::matrix_cl d22(d2);
  EXPECT_THROW(d22.sub_block(d11, 1, 1, 0, 0, 4, 4), std::domain_error);
  EXPECT_THROW(d22.sub_block(d11, 4, 4, 0, 0, 2, 2), std::domain_error);
}

TEST(MathMatrixGPU, sub_block_lower_pass) {
  stan::math::matrix_d d1;
  stan::math::matrix_d d2;

  d1.resize(3, 3);
  d2.resize(4, 4);

  d1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  d2 << 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1;

  stan::math::matrix_cl d11(d1);
  stan::math::matrix_cl d22(d2);
  d22.sub_block<stan::math::TriangularViewCL::Lower>(d11, 0, 0, 0, 0, 3, 3);
  stan::math::copy(d2, d22);
  EXPECT_EQ(1, d2(0, 0));
  EXPECT_EQ(15, d2(0, 1));
  EXPECT_EQ(14, d2(0, 2));
  EXPECT_EQ(10, d2(1, 2));
  EXPECT_EQ(7, d2(2, 0));
}

#endif
