#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <gtest/gtest.h>
#include <algorithm>

#define EXPECT_MATRIX_EQ(A, B)       \
  for (int i = 0; i < A.size(); i++) \
    EXPECT_EQ(A(i), B(i));

TEST(MathMatrixCL, rep_matrix_size_one_exception_pass) {
  stan::math::matrix_cl<double> a(2, 2);
  EXPECT_THROW(stan::math::rep_matrix(a, 5, 5), std::invalid_argument);
  EXPECT_THROW(stan::math::rep_matrix(a, 5, -1), std::domain_error);
  stan::math::matrix_cl<double> b(1, 1);
  EXPECT_THROW(stan::math::rep_matrix(b, -1, 2), std::domain_error);
  stan::math::matrix_cl<double> c(1, 0);
  EXPECT_THROW(stan::math::rep_matrix(c, 5, 2), std::invalid_argument);

  EXPECT_NO_THROW(stan::math::rep_matrix(b, 0, 2));
  EXPECT_NO_THROW(stan::math::rep_matrix(b, 5, 0));

  EXPECT_THROW(stan::math::rep_matrix(a, 5, 0), std::invalid_argument);
  EXPECT_THROW(stan::math::rep_matrix(c, 5, 0), std::invalid_argument);
}

TEST(MathMatrixCL, rep_matrix_size_one_value_check) {
  stan::math::matrix_d m0(1, 1);
  m0 << 2;
  stan::math::matrix_cl<double> m0_cl(m0);

  stan::math::matrix_cl<double> m1_cl = stan::math::rep_matrix(m0_cl, 1, 1);
  stan::math::matrix_d m1 = stan::math::from_matrix_cl(m1_cl);
  EXPECT_EQ(1, m1.size());
  EXPECT_EQ(2, m1(0, 0));

  stan::math::matrix_cl<double> m2_cl = stan::math::rep_matrix(m0_cl, 2, 2);
  stan::math::matrix_d m2 = stan::math::from_matrix_cl(m2_cl);
  EXPECT_EQ(2, m2.rows());
  EXPECT_EQ(2, m2.cols());
  EXPECT_EQ(2, m2(0, 0));
  EXPECT_EQ(2, m2(0, 1));
  EXPECT_EQ(2, m2(1, 0));
  EXPECT_EQ(2, m2(1, 1));

  stan::math::matrix_d m00(1, 1);
  m00 << -5;
  stan::math::matrix_cl<double> m00_cl(m00);
  stan::math::matrix_cl<double> m3_cl = stan::math::rep_matrix(m00_cl, 1, 4);
  stan::math::matrix_d m3 = stan::math::from_matrix_cl(m3_cl);
  EXPECT_EQ(1, m3.rows());
  EXPECT_EQ(4, m3.cols());
  EXPECT_EQ(-5, m3(0, 0));
  EXPECT_EQ(-5, m3(0, 1));
  EXPECT_EQ(-5, m3(0, 2));
  EXPECT_EQ(-5, m3(0, 3));

  stan::math::matrix_cl<double> m4_cl = stan::math::rep_matrix(m00_cl, 5, 1);
  stan::math::matrix_d m4 = stan::math::from_matrix_cl(m4_cl);
  EXPECT_EQ(5, m4.rows());
  EXPECT_EQ(1, m4.cols());
  EXPECT_EQ(-5, m4(0, 0));
  EXPECT_EQ(-5, m4(1, 0));
  EXPECT_EQ(-5, m4(2, 0));
  EXPECT_EQ(-5, m4(3, 0));
  EXPECT_EQ(-5, m4(4, 0));
}

TEST(MathMatrixCL, rep_matrix_v_rv_exception_pass) {
  stan::math::matrix_cl<double> a(2, 2);
  EXPECT_THROW(stan::math::rep_matrix(a, 5), std::invalid_argument);
  EXPECT_THROW(stan::math::rep_matrix(a, -1), std::domain_error);

  stan::math::matrix_cl<double> b(1, 1);
  EXPECT_THROW(stan::math::rep_matrix(b, 5), std::invalid_argument);
  stan::math::matrix_cl<double> c(1, 0);
  EXPECT_NO_THROW(stan::math::rep_matrix(c, 5));
  stan::math::matrix_cl<double> d(0, 5);
  EXPECT_THROW(stan::math::rep_matrix(d, 6), std::invalid_argument);
}

TEST(MathMatrixCL, rep_matrix_v_value_check) {
  stan::math::vector_d m0(6);
  m0 << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_cl<double> m0_cl(m0);

  stan::math::matrix_d m1 = stan::math::rep_matrix(m0, 12);
  stan::math::matrix_cl<double> m1_cl = stan::math::rep_matrix(m0_cl, 12);
  stan::math::matrix_d m1_cl_res = stan::math::from_matrix_cl(m1_cl);
  EXPECT_EQ(m1.rows(), m1_cl_res.rows());
  EXPECT_EQ(m1.cols(), m1_cl_res.cols());
  EXPECT_MATRIX_EQ(m1, m1_cl_res);
}

TEST(MathMatrixCL, rep_matrix_rv_value_check) {
  stan::math::row_vector_d m0(5);
  m0 << 1, 2, 3, 4, 5;
  stan::math::matrix_cl<double> m0_cl(m0);

  stan::math::matrix_d m1 = stan::math::rep_matrix(m0, 7);
  stan::math::matrix_cl<double> m1_cl = stan::math::rep_matrix(m0_cl, 7);
  stan::math::matrix_d m1_cl_res = stan::math::from_matrix_cl(m1_cl);
  EXPECT_EQ(m1.rows(), m1_cl_res.rows());
  EXPECT_EQ(m1.cols(), m1_cl_res.cols());
  EXPECT_MATRIX_EQ(m1, m1_cl_res);
}
#endif
