#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, transpose) {
  stan::math::vector_d v0;
  stan::math::row_vector_d rv0;
  stan::math::matrix_d m0;

  stan::math::matrix_gpu v00(v0);
  stan::math::matrix_gpu rv00(rv0);
  stan::math::matrix_gpu m00(m0);

  stan::math::matrix_gpu v00_dst(v0.cols(), v0.rows());
  stan::math::matrix_gpu rv00_dst(rv0.cols(), rv0.rows());
  stan::math::matrix_gpu m00_dst(m0.cols(), m0. rows());

  using stan::math::transpose;
  EXPECT_NO_THROW(transpose(v00));
  EXPECT_NO_THROW(transpose(rv00));
  EXPECT_NO_THROW(transpose(m00));
  EXPECT_NO_THROW(v00_dst = transpose(v00));
  EXPECT_NO_THROW(rv00_dst = transpose(rv00));
  EXPECT_NO_THROW(m00_dst = transpose(m00));
}
