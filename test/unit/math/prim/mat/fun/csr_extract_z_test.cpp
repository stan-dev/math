#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <vector>

// Test that non-zero entry counts of values from dense matrix in
// sparse format are extracted.
TEST(SparseStuff, csr_extract_z_dense) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();
  std::vector<int> result = stan::math::csr_extract_z(a);
  EXPECT_EQ(3, result[0]);
  EXPECT_EQ(3, result[1]);
}

// Test that non-zero entry counts of values from dense matrix in
// sparse format are extracted.
TEST(SparseStuff, csr_extract_z_dense_dense) {
  stan::math::matrix_d m(2, 3);
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  std::vector<int> result = stan::math::csr_extract_z(m);
  EXPECT_EQ(3, result[0]);
  EXPECT_EQ(3, result[1]);
}

// Test that non-zero entry counts of values from dense matrix in
// sparse format are extracted after makeCompressed();
TEST(SparseStuff, csr_extract_z_dense_compressed) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();
  a.makeCompressed();
  std::vector<int> result = stan::math::csr_extract_z(a);
  EXPECT_EQ(3, result[0]);
  EXPECT_EQ(3, result[1]);
}


// Test that non-zero entry counts of values from sparse matrix in
// sparse format are extracted.
TEST(SparseStuff, csr_extract_z_sparse) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 0.0, 0.0, 0.0;
  a = m.sparseView();
  std::vector<int> result = stan::math::csr_extract_z(a);
  EXPECT_EQ(3, result[0]);
  EXPECT_EQ(0, result[1]);
}

// Test that non-zero entry counts of values from sparse matrix in
// sparse format are extracted after makeCompressed();
TEST(SparseStuff, csr_extract_z_sparse_compressed) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 0.0, 0.0, 0.0;
  a = m.sparseView();
  a.makeCompressed();
  std::vector<int> result = stan::math::csr_extract_z(a);
  EXPECT_EQ(3, result[0]);
  EXPECT_EQ(0, result[1]);
}

