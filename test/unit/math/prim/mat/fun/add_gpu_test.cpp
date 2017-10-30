#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix,add_v_exception_pass) {
  stan::math::vector_d d1, d2;

  d1.resize(3);
  d2.resize(3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(3,1);
  EXPECT_NO_THROW(stan::math::add(d11, d22, d33));
}

TEST(MathMatrix,add_v_exception_pass_zero) {
  stan::math::vector_d d1, d2;
  d1.resize(0);
  d2.resize(0);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(1,1);
  EXPECT_NO_THROW(stan::math::add(d11, d22, d33));
}

TEST(MathMatrix,add_v_exception_pass_invalid_arg) {
  stan::math::row_vector_d d1, d2;

  d1.resize(2);
  d2.resize(3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(3,0);
  EXPECT_THROW(stan::math::add(d11, d22, d33), std::invalid_argument);
}

TEST(MathMatrix,add_rv_exception_pass) {
  stan::math::row_vector_d d1, d2;

  d1.resize(3);
  d2.resize(3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(1,1);
  EXPECT_NO_THROW(stan::math::add(d11, d22, d33));
}

TEST(MathMatrix,add_rv_exception_pass_zero) {
  stan::math::row_vector_d d1, d2;

  d1.resize(0);
  d2.resize(0);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(1,1);
  EXPECT_NO_THROW(stan::math::add(d11, d22, d33));
}

TEST(MathMatrix,add_rv_exception_fail_invalid_arg) {
  stan::math::row_vector_d d1, d2;

  d1.resize(2);
  d2.resize(3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(3,1);
  EXPECT_THROW(stan::math::add(d11, d22, d33), std::invalid_argument);
}

TEST(MathMatrix,add_m_exception_pass_simple) {
  stan::math::matrix_d d1, d2;

  d1.resize(2,3);
  d2.resize(2,3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(2,3);
  EXPECT_NO_THROW(stan::math::add(d1, d2));
}

TEST(MathMatrix,add_m_exception_pass_zero) {
  stan::math::matrix_d d1, d2;
  d1.resize(0,0);
  d2.resize(0,0);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(0,0);
  EXPECT_NO_THROW(stan::math::add(d1, d2));
}

TEST(MathMatrix,add_m_exception_fail_invalid_arg) {
  stan::math::matrix_d d1, d2;
  d1.resize(2,3);
  d2.resize(3,3);
  stan::math::matrix_gpu d11(d1);
  stan::math::matrix_gpu d22(d2);
  stan::math::matrix_gpu d33(2,3);
  EXPECT_THROW(stan::math::add(d1, d2), std::invalid_argument);
}

//TODO(Steve): Matrix GPU function for adding constant
/*
TEST(MathMatrix,add_c_m) {
  stan::math::matrix_d v(2,2);
  v << 1, 2, 3, 4;
  stan::math::matrix_d result;

  result = stan::math::add(2.0,v);
  EXPECT_FLOAT_EQ(3.0,result(0,0));
  EXPECT_FLOAT_EQ(4.0,result(0,1));
  EXPECT_FLOAT_EQ(5.0,result(1,0));
  EXPECT_FLOAT_EQ(6.0,result(1,1));

  result = stan::math::add(v,2.0);
  EXPECT_FLOAT_EQ(3.0,result(0,0));
  EXPECT_FLOAT_EQ(4.0,result(0,1));
  EXPECT_FLOAT_EQ(5.0,result(1,0));
  EXPECT_FLOAT_EQ(6.0,result(1,1));
}

TEST(MathMatrix,add_c_rv) {
  stan::math::row_vector_d v(3);
  v << 1, 2, 3;
  stan::math::row_vector_d result;

  result = stan::math::add(2.0,v);
  EXPECT_FLOAT_EQ(3.0,result(0));
  EXPECT_FLOAT_EQ(4.0,result(1));
  EXPECT_FLOAT_EQ(5.0,result(2));

  result = stan::math::add(v,2.0);
  EXPECT_FLOAT_EQ(3.0,result(0));
  EXPECT_FLOAT_EQ(4.0,result(1));
  EXPECT_FLOAT_EQ(5.0,result(2));
}


TEST(MathMatrix,add_c_v) {
  stan::math::vector_d v(3);
  v << 1, 2, 3;
  stan::math::vector_d result;

  result = stan::math::add(2.0,v);
  EXPECT_FLOAT_EQ(3.0,result(0));
  EXPECT_FLOAT_EQ(4.0,result(1));
  EXPECT_FLOAT_EQ(5.0,result(2));

  result = stan::math::add(v,2.0);
  EXPECT_FLOAT_EQ(3.0,result(0));
  EXPECT_FLOAT_EQ(4.0,result(1));
  EXPECT_FLOAT_EQ(5.0,result(2));
}
*/

TEST(MathMatrix, add) {
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

  using stan::math::add;
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

  EXPECT_THROW(add(v11, v22, v33),std::invalid_argument);
  EXPECT_THROW(add(rv11, rv22, rv33),std::invalid_argument);
  EXPECT_THROW(add(m11, m22, m33),std::invalid_argument);
}
