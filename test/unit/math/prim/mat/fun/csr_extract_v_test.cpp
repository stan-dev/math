#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <vector>

// Test that column indexes of values from dense matrix in sparse
// format are extracted.
TEST(SparseStuff, csr_extract_v_dense) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();
  std::vector<int> result = stan::math::csr_extract_v(a);
  EXPECT_EQ(1, result[0]);
  EXPECT_EQ(2, result[1]);
  EXPECT_EQ(3, result[2]);
  EXPECT_EQ(1, result[3]);
  EXPECT_EQ(2, result[4]);
  EXPECT_EQ(3, result[5]);
}

// Test that column indexes of values from dense matrix in dense
// format are extracted.
TEST(SparseStuff, csr_extract_v_dense_dense) {
  stan::math::matrix_d m(2, 3);
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  std::vector<int> result = stan::math::csr_extract_v(m);
  EXPECT_EQ(1, result[0]);
  EXPECT_EQ(2, result[1]);
  EXPECT_EQ(3, result[2]);
  EXPECT_EQ(1, result[3]);
  EXPECT_EQ(2, result[4]);
  EXPECT_EQ(3, result[5]);
}

// Test that column indexes of values from dense matrix in sparse
// format are extracted after makeCompressed().
TEST(SparseStuff, csr_extract_v_dense_compressed) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();
  a.makeCompressed();
  std::vector<int> result = stan::math::csr_extract_v(a);
  EXPECT_EQ(1, result[0]);
  EXPECT_EQ(2, result[1]);
  EXPECT_EQ(3, result[2]);
  EXPECT_EQ(1, result[3]);
  EXPECT_EQ(2, result[4]);
  EXPECT_EQ(3, result[5]);
}

// Test that column indexes of values from sparse matrix in sparse
// format are extracted.
TEST(SparseStuff, csr_extract_v_sparse) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 0.0, 0.0, 0.0;
  a = m.sparseView();
  std::vector<int> result = stan::math::csr_extract_v(a);
  EXPECT_EQ(1, result[0]);
  EXPECT_EQ(2, result[1]);
  EXPECT_EQ(3, result[2]);
  EXPECT_EQ(3, result.size());
}

// Test that column indexes of values from sparse matrix in sparse
// format are extracted after makeCompressed().
TEST(SparseStuff, csr_extract_v_sparse_compressed) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 0.0, 0.0, 0.0;
  a = m.sparseView();
  a.makeCompressed();
  std::vector<int> result = stan::math::csr_extract_v(a);
  EXPECT_EQ(1, result[0]);
  EXPECT_EQ(2, result[1]);
  EXPECT_EQ(3, result[2]);
  EXPECT_EQ(3, result.size());
}



