#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/transpose.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathGpu, transpose_size_zero) {
  stan::math::vector_d v0;
  stan::math::row_vector_d rv0;
  stan::math::matrix_d m0;

  stan::math::matrix_cl v00(v0);
  stan::math::matrix_cl rv00(rv0);
  stan::math::matrix_cl m00(m0);

  stan::math::matrix_cl v00_dst(v0.cols(), v0.rows());
  stan::math::matrix_cl rv00_dst(rv0.cols(), rv0.rows());
  stan::math::matrix_cl m00_dst(m0.cols(), m0.rows());

  using stan::math::transpose;
  EXPECT_NO_THROW(transpose(v00));
  EXPECT_NO_THROW(transpose(rv00));
  EXPECT_NO_THROW(transpose(m00));
  EXPECT_NO_THROW(v00_dst = transpose(v00));
  EXPECT_NO_THROW(rv00_dst = transpose(rv00));
  EXPECT_NO_THROW(m00_dst = transpose(m00));
}

TEST(MathGpu, transpose) {
  stan::math::vector_d v0(3);
  stan::math::row_vector_d v0_dst(3);
  stan::math::row_vector_d rv0(3);
  stan::math::vector_d rv0_dst(3);
  stan::math::matrix_d m0(3, 3);
  stan::math::matrix_d m0_dst(m0.cols(), m0.rows());

  v0 << 1, 2, 3;
  rv0 << 1, 2, 3;
  m0 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::math::matrix_cl v00(v0);
  stan::math::matrix_cl rv00(rv0);
  stan::math::matrix_cl m00(m0);

  stan::math::matrix_cl v00_dst(v0.cols(), v0.rows());
  stan::math::matrix_cl rv00_dst(rv0.cols(), rv0.rows());
  stan::math::matrix_cl m00_dst(m0.cols(), m0.rows());

  using stan::math::transpose;
  v00_dst = transpose(v00);
  EXPECT_NO_THROW(rv00_dst = transpose(rv00));
  EXPECT_NO_THROW(m00_dst = transpose(m00));

  EXPECT_NO_THROW(stan::math::copy(v0_dst, v00_dst));
  EXPECT_EQ(v0(0), v0_dst(0));
  EXPECT_EQ(v0(1), v0_dst(1));
  EXPECT_EQ(v0(2), v0_dst(2));
  EXPECT_NO_THROW(stan::math::copy(rv0_dst, rv00_dst));
  EXPECT_EQ(rv0(0), rv0_dst(0));
  EXPECT_EQ(rv0(1), rv0_dst(1));
  EXPECT_EQ(rv0(2), rv0_dst(2));
  EXPECT_NO_THROW(stan::math::copy(m0_dst, m00_dst));
  EXPECT_EQ(m0(0, 0), m0_dst(0, 0));
  EXPECT_EQ(m0(0, 1), m0_dst(1, 0));
  EXPECT_EQ(m0(0, 2), m0_dst(2, 0));
  EXPECT_EQ(m0(1, 0), m0_dst(0, 1));
  EXPECT_EQ(m0(1, 1), m0_dst(1, 1));
  EXPECT_EQ(m0(1, 2), m0_dst(2, 1));
  EXPECT_EQ(m0(2, 0), m0_dst(0, 2));
  EXPECT_EQ(m0(2, 1), m0_dst(1, 2));
  EXPECT_EQ(m0(2, 2), m0_dst(2, 2));
}
#endif
