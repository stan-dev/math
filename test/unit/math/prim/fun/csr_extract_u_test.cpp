#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <vector>

// Test that start of each row's values in NZE vector (w) is correctly
// extracted from a dense matrix in sparse format.
TEST(SparseStuff, csr_extract_u_dense) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();
  std::vector<int> result = stan::math::csr_extract_u(a);
  EXPECT_EQ(1, result[0]);
  EXPECT_EQ(4, result[1]);
  EXPECT_EQ(7, result[2]);
}

// Test that start of each row's values in NZE vector (w) is correctly
// extracted from a dense matrix in sparse format.
TEST(SparseStuff, csr_extract_u_dense_dense) {
  stan::math::matrix_d m(2, 3);
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  std::vector<int> result = stan::math::csr_extract_u(m);
  EXPECT_EQ(1, result[0]);
  EXPECT_EQ(4, result[1]);
  EXPECT_EQ(7, result[2]);
}

// Test that start of each row's values in NZE vector (w) is correctly
// extracted from a dense matrix in sparse format after
// makeCompressed().
TEST(SparseStuff, csr_extract_u_dense_compressed) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();
  a.makeCompressed();
  std::vector<int> result = stan::math::csr_extract_u(a);
  EXPECT_EQ(1, result[0]);
  EXPECT_EQ(4, result[1]);
  EXPECT_EQ(7, result[2]);
}

// Test that start of each row's values in NZE vector (w) is correctly
// extracted from a sparse matrix in sparse format.
TEST(SparseStuff, csr_extract_u_sparse) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 0.0, 0.0, 0.0;
  a = m.sparseView();
  std::vector<int> result = stan::math::csr_extract_u(a);
  EXPECT_EQ(1, result[0]);
  EXPECT_EQ(4, result[1]);
  EXPECT_EQ(4, result[2]);
  EXPECT_EQ(3, result.size());
}

// Test that start of each row's values in NZE vector (w) is correctly
// extracted from a dense matrix in sparse format after
// makeCompressed().
TEST(SparseStuff, csr_extract_u_sparse_compressed) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 0.0, 0.0, 0.0;
  a = m.sparseView();
  a.makeCompressed();
  std::vector<int> result = stan::math::csr_extract_u(a);
  EXPECT_EQ(1, result[0]);
  EXPECT_EQ(4, result[1]);
  EXPECT_EQ(4, result[2]);
  EXPECT_EQ(3, result.size());
}
