#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

// Test that values from a dense matrix in sparse format are extracted.
TEST(SparseStuff, csr_extract_w_dense) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();
  stan::math::vector_d result = stan::math::csr_extract_w(a);
  EXPECT_FLOAT_EQ( 2.0, result(0));
  EXPECT_FLOAT_EQ( 4.0, result(1));
  EXPECT_FLOAT_EQ( 6.0, result(2));
  EXPECT_FLOAT_EQ( 8.0, result(3));
  EXPECT_FLOAT_EQ(10.0, result(4));
  EXPECT_FLOAT_EQ(12.0, result(5));
}

// Test that values from a dense matrix are extracted.
TEST(SparseStuff, csr_extract_w_dense_dense) {
  stan::math::matrix_d m(2, 3);
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  stan::math::vector_d result = stan::math::csr_extract_w(m);
  EXPECT_FLOAT_EQ( 2.0, result(0));
  EXPECT_FLOAT_EQ( 4.0, result(1));
  EXPECT_FLOAT_EQ( 6.0, result(2));
  EXPECT_FLOAT_EQ( 8.0, result(3));
  EXPECT_FLOAT_EQ(10.0, result(4));
  EXPECT_FLOAT_EQ(12.0, result(5));
}

// Test that values from a dense matrix in sparse format are extracted
// after A.makeCompressed();
TEST(SparseStuff, csr_extract_w_dense_compressed) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();
  a.makeCompressed();
  stan::math::vector_d result = stan::math::csr_extract_w(a);
  EXPECT_FLOAT_EQ( 2.0, result(0));
  EXPECT_FLOAT_EQ( 4.0, result(1));
  EXPECT_FLOAT_EQ( 6.0, result(2));
  EXPECT_FLOAT_EQ( 8.0, result(3));
  EXPECT_FLOAT_EQ(10.0, result(4));
  EXPECT_FLOAT_EQ(12.0, result(5));
}

// Test that values from a sparse matrix in sparse format are extracted
TEST(SparseStuff, csr_extract_w_sparse) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 0.0, 0.0, 0.0;
  a = m.sparseView();
  stan::math::vector_d result = stan::math::csr_extract_w(a);
  EXPECT_FLOAT_EQ(2.0, result(0));
  EXPECT_FLOAT_EQ(4.0, result(1));
  EXPECT_FLOAT_EQ(6.0, result(2));
}


// Test that values from a sparse matrix in sparse format are extracted
// after A.makeCompressed()
TEST(SparseStuff, csr_extract_w_sparse_compressed) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 0.0, 0.0, 0.0;
  a = m.sparseView();
  stan::math::vector_d result = stan::math::csr_extract_w(a);
  EXPECT_FLOAT_EQ( 2.0, result(0));
  EXPECT_FLOAT_EQ( 4.0, result(1));
  EXPECT_FLOAT_EQ( 6.0, result(2));
}

