#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <vector>

// Test that dense multiplication results is correct (CSR).
TEST(SparseStuff, csr_matrix_times_vector_dense_multiply) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();

  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);

  stan::math::vector_d b(3);
  b << 22, 33, 44;

  stan::math::vector_d result
      = stan::math::csr_matrix_times_vector(2, 3, X_w, X_v, X_u, b);
  EXPECT_FLOAT_EQ(440.0, result(0));
  EXPECT_FLOAT_EQ(1034.0, result(1));
}

// Test that sparse multiplication is correct when columns are all zero.
TEST(SparseStuff, csr_matrix_times_vector_empty_row_multiply) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 0.0, 8.0, 10.0, 0.0;
  a = m.sparseView();

  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);

  stan::math::vector_d b(3);
  b << 22, 33, 44;

  stan::math::vector_d result
      = stan::math::csr_matrix_times_vector(2, 3, X_w, X_v, X_u, b);
  EXPECT_FLOAT_EQ(176.0, result(0));
  EXPECT_FLOAT_EQ(506.0, result(1));
}

// Test that sparse multiplication is correct when columns are all zero.
TEST(SparseStuff, csr_matrix_times_vector_empty_column_multiply) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 0.0, 6.0, 8.0, 0.0, 12.0;
  a = m.sparseView();

  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);

  stan::math::vector_d b(3);
  b << 22, 33, 44;

  stan::math::vector_d result
      = stan::math::csr_matrix_times_vector(2, 3, X_w, X_v, X_u, b);

  EXPECT_FLOAT_EQ(308.0, result(0));
  EXPECT_FLOAT_EQ(704.0, result(1));
}

// Test that m=0 throws (CSR).
TEST(SparseStuff, csr_matrix_times_vector_m0) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();

  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);

  stan::math::vector_d b(3);
  b << 22, 33, 44;

  EXPECT_THROW(stan::math::csr_matrix_times_vector(0, 3, X_w, X_v, X_u, b),
               std::domain_error);
}

// Test that n=0 throws (CSR).
TEST(SparseStuff, csr_matrix_times_vector_n0) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();

  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);

  stan::math::vector_d b(3);
  b << 22, 33, 44;
  EXPECT_THROW(stan::math::csr_matrix_times_vector(2, 0, X_w, X_v, X_u, b),
               std::domain_error);
}

// Test that short b throws (CSR).
TEST(SparseStuff, csr_matrix_times_vector_b_short) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();

  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);

  stan::math::vector_d b(2);  // short b
  b << 22, 33;

  EXPECT_THROW(stan::math::csr_matrix_times_vector(2, 3, X_w, X_v, X_u, b),
               std::invalid_argument);
}

// Test that short u throws (CSR).
TEST(SparseStuff, csr_matrix_times_vector_u_short) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();

  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);

  X_u.erase(X_u.begin());  // make a short u:

  stan::math::vector_d b(3);
  b << 22, 33, 44;

  EXPECT_THROW(stan::math::csr_matrix_times_vector(2, 3, X_w, X_v, X_u, b),
               std::invalid_argument);
}

// Test that short v throws (CSR).
TEST(SparseStuff, csr_matrix_times_vector_v_short) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();

  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);
  X_v.erase(X_v.begin() + 4);  // make a short v:

  stan::math::vector_d b(3);
  b << 22, 33, 44;
  EXPECT_THROW(stan::math::csr_matrix_times_vector(2, 3, X_w, X_v, X_u, b),
               std::invalid_argument);
}
