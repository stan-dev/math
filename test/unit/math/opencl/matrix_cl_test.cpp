#ifdef STAN_OPENCL
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/sub_block.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

template <typename T>
inline void test_matrix_creation() {
  std::vector<T> vec_1({0, 1, 2, 3});
  Eigen::Matrix<T, 2, 2> mat_1;
  Eigen::Matrix<T, 2, 2> mat_2;
  mat_1 << 1, 2, 3, 4;
  EXPECT_NO_THROW(stan::math::matrix_cl<T> A(1, 1));
  EXPECT_NO_THROW(stan::math::matrix_cl<T> vec_1cl(vec_1, 2, 2));
  EXPECT_NO_THROW(stan::math::matrix_cl<T> mat_1cl(mat_1));
  EXPECT_NO_THROW(stan::math::matrix_cl<T> mat_2cl(mat_2));
}

TEST(MathMatrixCL, matrix_cl_types_creation) {
  test_matrix_creation<double>();
  test_matrix_creation<float>();
  test_matrix_creation<int>();
  test_matrix_creation<long double>();
}

#ifndef STAN_OPENCL_NOCACHE
TEST(MathMatrixCL, matrix_cl_cache) {
  using stan::math::matrix_cl;
  Eigen::MatrixXd m(2, 2);
  m << 1, 2, 3, 4;

  const matrix_cl<double> m1_cl = matrix_cl<double>::constant(m);

  cl_mem mem_handle = m.opencl_buffer_();
  EXPECT_EQ(mem_handle, m1_cl.buffer()());

  const matrix_cl<double> m2_cl = matrix_cl<double>::constant(m);

  EXPECT_EQ(mem_handle, m.opencl_buffer_());
  EXPECT_EQ(mem_handle, m1_cl.buffer()());
  EXPECT_EQ(mem_handle, m2_cl.buffer()());
}
#else
TEST(MathMatrixCL, matrix_cl_nocache) {
  using stan::math::matrix_cl;
  Eigen::MatrixXd m(2, 2);
  m << 1, 2, 3, 4;

  const matrix_cl<double> m1_cl = matrix_cl<double>::constant(m);
  const matrix_cl<double> m2_cl = matrix_cl<double>::constant(m);

  EXPECT_NE(m1_cl.buffer()(), m2_cl.buffer()());
}
#endif

TEST(MathMatrixCL, matrix_cl_constructor_nocache) {
  using stan::math::matrix_cl;
  Eigen::MatrixXd m(2, 2);
  m << 1, 2, 3, 4;

  const matrix_cl<double> m1_cl(m);
  const matrix_cl<double> m2_cl(m);

  EXPECT_NE(m1_cl.buffer()(), m2_cl.buffer()());
}
#endif
