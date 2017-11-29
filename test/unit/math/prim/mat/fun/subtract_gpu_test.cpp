#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixGPU,subtract_v_exception_pass) {
  stan::math::vector_d d1, d2;

  d1.resize(3);
  d2.resize(3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(3,1);
  EXPECT_NO_THROW(d33 = stan::math::subtract(d11, d22));
}

// TODO(Steve): This should probably throw expection? 
//The dimensions do not match
TEST(MathMatrixGPU,subtract_v_exception_pass_zero) {
  stan::math::vector_d d1, d2;
  d1.resize(0);
  d2.resize(0);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(0,1);
  EXPECT_NO_THROW(d33 = stan::math::subtract(d11, d22));
}

TEST(MathMatrixGPU,subtract_v_exception_fail_zero) {
  stan::math::vector_d d1, d2;
  d1.resize(2);
  d2.resize(3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(3,3);
  EXPECT_THROW(d33 = stan::math::subtract(d11, d22), std::invalid_argument);
}

TEST(MathMatrixGPU,subtract_rv_exception_pass) {
  stan::math::row_vector_d d1, d2;

  d1.resize(3);
  d2.resize(3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(1,3);
  EXPECT_NO_THROW(d33 = stan::math::subtract(d11, d22));
}

TEST(MathMatrixGPU,subtract_rv_exception_pass_zero) {
  stan::math::row_vector_d d1, d2;
  d1.resize(0);
  d2.resize(0);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(1,0);
  EXPECT_NO_THROW(d33 = stan::math::subtract(d11, d22));
}

TEST(MathMatrixGPU,subtract_rv_exception_fail_zero) {
  stan::math::row_vector_d d1, d2;
  d1.resize(2);
  d2.resize(3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(3,1);
  EXPECT_THROW(d33 = stan::math::subtract(d11, d22), std::invalid_argument);
}

TEST(MathMatrixGPU,subtract_m_exception) {
  stan::math::matrix_d d1, d2;
  d1.resize(2,3);
  d2.resize(2,3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(2,3);
  EXPECT_NO_THROW(d33 = stan::math::subtract(d11, d22));
}

TEST(MathMatrixGPU,subtract_m_exception_pass_zero) {
  stan::math::matrix_d d1, d2;
  d1.resize(0,0);
  d2.resize(0,0);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(0,0);
  EXPECT_NO_THROW(d33 = stan::math::subtract(d11, d22));
}

TEST(MathMatrixGPU,subtract_m_exception_fail) {
  stan::math::matrix_d d1, d2;
  d1.resize(2,3);
  d2.resize(3,3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(3,2);
  EXPECT_THROW(d33 = stan::math::subtract(d11, d22), std::invalid_argument);
}

TEST(MathMatrixGPU, subtract_exception) {
  stan::math::vector_d v1(2);
  v1 << 1, 2;
  stan::math::vector_d v2(3);
  v2 << 10, 100, 1000;

  stan::math::row_vector_d rv1(2);
  v1 << 1, 2;
  stan::math::row_vector_d rv2(3);
  v2 << 10, 100, 1000;

  stan::math::matrix_d m1(2,3);
  m1 << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_d m2(3,2);
  m2 << 10, 100, 1000, 0, -10, -12;

  using stan::math::subtract;
  using stan::math::matrix_gpu;
  matrix_gpu v11(v1);
  matrix_gpu v22(v2);
  matrix_gpu v33(v1);
  matrix_gpu rv11(rv1);
  matrix_gpu rv22(rv2);
  matrix_gpu rv33(rv1);
  matrix_gpu m11(m1);
  matrix_gpu m22(m2);
  matrix_gpu m33(m1);
  EXPECT_THROW(v33 = subtract(v11, v22),std::invalid_argument);
  EXPECT_THROW(rv33 = subtract(rv11, rv22),std::invalid_argument);
  EXPECT_THROW(m33 = subtract(m11, m22),std::invalid_argument);
}
