#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix,multiply_c_v) {
  stan::math::vector_d v(3);
  v << 1, 2, 3;
  stan::math::matrix_gpu vv(v);
  vv = stan::math::multiply(vv, 2.0);
  stan::math::copy(vv, v);  
  EXPECT_FLOAT_EQ(2.0, v(0));
  EXPECT_FLOAT_EQ(4.0, v(1));
  EXPECT_FLOAT_EQ(6.0, v(2));
}

TEST(MathMatrix,multiply_c_rv) {
  stan::math::row_vector_d rv(3);
  rv << 1, 2, 3;
  stan::math::matrix_gpu vv(rv);
  vv = stan::math::multiply(vv, 2.0);
  stan::math::copy(vv, rv);  
  EXPECT_FLOAT_EQ(2.0, rv(0));
  EXPECT_FLOAT_EQ(4.0, rv(1));
  EXPECT_FLOAT_EQ(6.0, rv(2));
}

TEST(MathMatrix, multiply_c_m) {
  stan::math::matrix_d m(2,3);
  m << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_gpu mm(m);
  mm = stan::math::multiply(mm, 2.0);
  stan::math::copy(mm, m);  
  EXPECT_FLOAT_EQ(2.0, m(0,0));
  EXPECT_FLOAT_EQ(4.0, m(0,1));
  EXPECT_FLOAT_EQ(6.0, m(0,2));
  EXPECT_FLOAT_EQ(8.0, m(1,0));
  EXPECT_FLOAT_EQ(10.0, m(1,1));
  EXPECT_FLOAT_EQ(12.0, m(1,2));
}

TEST(MathMatrix,multiply_rv_v_exception_pass_vec) {
  stan::math::row_vector_d rv;
  stan::math::vector_d v;
  
  rv.resize(3);
  v.resize(3);
  stan::math::matrix_gpu rvv(rv);
  stan::math::matrix_gpu vv(v);
  stan::math::matrix_gpu ans_v(1,1);
  EXPECT_NO_THROW(ans_v = stan::math::multiply(rvv, vv));
}

TEST(MathMatrix,multiply_rv_v_exception_pass_empty) {
  stan::math::row_vector_d rv;
  stan::math::vector_d v;
  rv.resize(0);
  v.resize(0);
  stan::math::matrix_gpu rvv(rv);
  stan::math::matrix_gpu vv(v);
  stan::math::matrix_gpu ans_vv(1,1);
  EXPECT_NO_THROW(ans_vv = stan::math::multiply(rvv, vv));
}

TEST(MathMatrix,multiply_rv_v_exception_fail) {
  stan::math::row_vector_d rv;
  stan::math::vector_d v;
  rv.resize(2);
  v.resize(3);
  stan::math::matrix_gpu rvv(rv);
  stan::math::matrix_gpu vv(v);
  stan::math::matrix_gpu ans_vv(2,3);
  EXPECT_THROW(ans_vv = stan::math::multiply(rvv, vv), std::invalid_argument);
}

TEST(MathMatrix,multiply_m_v_exception_pass) {
  stan::math::matrix_d m;
  stan::math::vector_d v;
  
  m.resize(3, 5);
  v.resize(5);
  stan::math::matrix_gpu mm(m);
  stan::math::matrix_gpu vv(v);
  stan::math::matrix_gpu ans_mm(3,1);
  EXPECT_NO_THROW(ans_mm = stan::math::multiply(mm, vv));
}
//TODO(Steve): Should this really throw an exception?
//The dimensions match, so I think it should not
//I made it no_throw
TEST(MathMatrix,multiply_m_v_exception_fail_zero) {
  stan::math::matrix_d m;
  stan::math::vector_d v;
  m.resize(3, 0);
  v.resize(0);
  stan::math::matrix_gpu mm(m);
  stan::math::matrix_gpu vv(v);
  stan::math::matrix_gpu ans_mm(1, 0);
  EXPECT_NO_THROW(ans_mm = stan::math::multiply(mm, vv));
}

TEST(MathMatrix,multiply_m_v_exception_pass_pass) {
  stan::math::matrix_d m;
  stan::math::vector_d v;
  m.resize(2, 3);
  v.resize(2);
  stan::math::matrix_gpu mm(m);
  stan::math::matrix_gpu vv(v);
  stan::math::matrix_gpu ans_mm(2, 1);
  EXPECT_THROW(ans_mm = stan::math::multiply(mm, vv), std::invalid_argument);  
}

TEST(MathMatrix,multiply_rv_m_exception_pass_vec1) {
  stan::math::row_vector_d rv;
  stan::math::matrix_d m;
    
  rv.resize(3);
  m.resize(3, 5);
  stan::math::matrix_gpu mm(m);
  stan::math::matrix_gpu rvv(rv);
  stan::math::matrix_gpu ans_mm(1, 5);
  EXPECT_NO_THROW(ans_mm = stan::math::multiply(rvv, mm));
}

//TODO(Steve): Should this really throw an exception?
//The dimensions match, so I think it should not
//I made it no_throw
TEST(MathMatrix,multiply_rv_m_exception_fail_zero1) {
  stan::math::row_vector_d rv;
  stan::math::matrix_d m;
  rv.resize(0);
  m.resize(0, 3);
  stan::math::matrix_gpu mm(m);
  stan::math::matrix_gpu rvv(rv);
  stan::math::matrix_gpu ans_mm(1, 0);
  EXPECT_NO_THROW(ans_mm = stan::math::multiply(rvv, mm));
}

TEST(MathMatrix,multiply_rv_m_exception_fail_dims) {
  stan::math::row_vector_d rv;
  stan::math::matrix_d m;
  rv.resize(3);
  m.resize(2, 3);
  stan::math::matrix_gpu mm(m);
  stan::math::matrix_gpu rvv(rv);
  stan::math::matrix_gpu ans_mm(1,3);
  EXPECT_THROW(ans_mm = stan::math::multiply(rvv, mm), std::invalid_argument);
}

TEST(MathMatrix,multiply_m_exception_pass_diagonal_mul) {
  stan::math::matrix_d m0;
  m0.resize(3,2);
  m0 << 1, 1, 1,
  1, 1, 1;
  stan::math::matrix_gpu mm(m0);
  EXPECT_NO_THROW(stan::math::diagonal_multiply(mm, 1.0));
}

TEST(MathMatrix,multiply_m_m_exception_pass_dim) {
  stan::math::matrix_d m1, m2;
  
  m1.resize(1, 3);
  m2.resize(3, 5);
  stan::math::matrix_gpu mm1(m1);
  stan::math::matrix_gpu mm2(m2);
  stan::math::matrix_gpu mm3(1, 5);
  EXPECT_NO_THROW(mm3 = stan::math::multiply(mm1, mm2));
}

//TODO(Steve): Should this really throw an exception?
//The dimensions match, so I think it should not
//I made it no_throw
TEST(MathMatrix,multiply_m_m_exception_fail_dim_zero) {
  stan::math::matrix_d m1, m2; 
  m1.resize(2, 0);
  m2.resize(0, 3);
  stan::math::matrix_gpu mm1(m1);
  stan::math::matrix_gpu mm2(m2);
  stan::math::matrix_gpu mm3(3, 3);
  EXPECT_NO_THROW(mm3 = stan::math::multiply(mm1, mm2));
}

TEST(MathMatrix,multiply_m_m_exception_fail_dim) {
  stan::math::matrix_d m1, m2;
  m1.resize(4, 3);
  m2.resize(2, 3);
  stan::math::matrix_gpu mm1(m1);
  stan::math::matrix_gpu mm2(m2);
  stan::math::matrix_gpu mm3(3, 3);
  EXPECT_THROW(mm3 = stan::math::multiply(mm1, mm2), std::invalid_argument);
}

TEST(MathMatrix, multiply) {
  stan::math::vector_d v0;
  stan::math::row_vector_d rv0;
  stan::math::matrix_d m0;
  
  stan::math::matrix_gpu v00(v0);
  
  
  stan::math::matrix_gpu rv00(rv0);
  stan::math::matrix_gpu m00(m0);
  using stan::math::multiply;
  EXPECT_NO_THROW(multiply(v00,2.0));
  EXPECT_NO_THROW(multiply(rv00,2.0));
  EXPECT_NO_THROW(multiply(m00,2.0));
  EXPECT_NO_THROW(multiply(2.0,v00));
  EXPECT_NO_THROW(multiply(2.0,rv00));
  EXPECT_NO_THROW(multiply(2.0,m00));
}


TEST(AgradRevMatrix,multiply_int) {
  using stan::math::multiply;
  using stan::math::assign;

  typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
  
  int d_int = 2;
  vector_d vec(4);
  vec << 1, 2, 3, 4;
  vector_d t_vec(4);
  assign(t_vec, multiply(vec,d_int));

}

TEST(AgradRevMatrix,multiply_vector_int) {
  using stan::math::multiply;
  using stan::math::vector_d;

  vector_d dvec(3);
  dvec << 1, 2, 3;
  int a = 2;
  vector_d prod_vec = multiply(dvec,a);
  EXPECT_EQ(3,prod_vec.size());
  EXPECT_EQ(2.0, prod_vec[0]);
  EXPECT_EQ(4.0, prod_vec[1]);
  EXPECT_EQ(6.0, prod_vec[2]);
}
