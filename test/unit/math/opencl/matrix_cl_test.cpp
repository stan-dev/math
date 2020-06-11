#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/zeros.hpp>
#include <stan/math/opencl/sub_block.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
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

TEST(MathMatrixCL, matrix_cl_value) {
  Eigen::Matrix<double, -1, -1> col_major(3, 3);
  col_major << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::math::matrix_cl<double> cl_from_col_major(col_major);
  Eigen::MatrixXd res = stan::math::from_matrix_cl(cl_from_col_major);
  expect_matrix_eq(col_major, res);

  Eigen::Matrix<double, -1, -1, Eigen::RowMajor> row_major(3, 3);
  row_major << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::math::matrix_cl<double> cl_from_row_major(row_major);
  res = stan::math::from_matrix_cl(cl_from_row_major);
  expect_matrix_eq(row_major, res);
}

TEST(MathMatrixCL, assignment) {
  using stan::math::matrix_cl;
  Eigen::Matrix<double, 2, 2> mat_1;
  Eigen::Matrix<double, 2, 2> mat_2;
  mat_1 << 1, 2, 3, 4;
  matrix_cl<double> mat1_cl(mat_1);
  matrix_cl<double> mat2_cl(2, 2);
  {
    Eigen::Matrix<double, 2, 2> mat_local_1;
    mat_local_1 << 1, 2, 3, 4;
    matrix_cl<double> mat1_local_cl(mat_local_1);
    matrix_cl<double> mat2_local_cl(mat_local_1);
    mat2_cl = mat1_local_cl;
    mat1_local_cl.template zeros<stan::math::matrix_cl_view::Entire>();
  }
  Eigen::Matrix<double, 2, 2> mat_2_fromcl
      = stan::math::from_matrix_cl(mat2_cl);
  // Make sure mat_2_from_cl matches the local scoped values
  EXPECT_EQ(mat_2_fromcl(0), 1);
  EXPECT_EQ(mat_2_fromcl(1), 3);
  EXPECT_EQ(mat_2_fromcl(2), 2);
  EXPECT_EQ(mat_2_fromcl(3), 4);

  EXPECT_NO_THROW({ mat2_cl = std::move(mat1_cl); });
  // mat1_cl should be null
  EXPECT_NE(mat2_cl.buffer()(), mat1_cl.buffer()());
  EXPECT_EQ(nullptr, mat1_cl.buffer()());
}

#endif
