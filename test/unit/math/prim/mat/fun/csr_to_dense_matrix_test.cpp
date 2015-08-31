#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <vector>


// Test that sparse and dense multiplication results is the same after
// plumbing through csr_extract_*.
TEST(SparseStuff, csr_to_dense_matrix_two_route) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();
  
  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);
  std::vector<int> X_z = stan::math::csr_extract_z(a);
  Eigen::Matrix<double, 2, 3> A = stan::math::csr_to_dense_matrix(2, 3, X_w, X_v, X_u, X_z);

  stan::math::vector_d b(3);
  b << 22, 33, 44;

  stan::math::vector_d result_sparse = stan::math::csr_matrix_times_vector(
      2, 3, X_w, X_v, X_u, X_z, b);
  stan::math::vector_d result_dense = A * b;
  stan::math::vector_d result_straight_dense = m * b;

  EXPECT_FLOAT_EQ(result_dense(0), result_sparse(0));
  EXPECT_FLOAT_EQ(result_dense(1), result_sparse(1));
  EXPECT_FLOAT_EQ(result_dense(0), result_straight_dense(0));
  EXPECT_FLOAT_EQ(result_dense(1), result_straight_dense(1));
}

// Test that sparse and dense multiplication results is the same after
// plumbing through csr_extract_* with sparse test matrix.
TEST(SparseStuff, csr_to_dense_matrix_two_route_sparse) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 0.0, 0.0, 0.0, 0.0;
  a = m.sparseView();
  
  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);
  std::vector<int> X_z = stan::math::csr_extract_z(a);
  Eigen::Matrix<double, 2, 3> A = stan::math::csr_to_dense_matrix(2, 3, X_w, X_v, X_u, X_z);

  stan::math::vector_d b(3);
  b << 22, 33, 44;

  stan::math::vector_d result_sparse = stan::math::csr_matrix_times_vector(
      2, 3, X_w, X_v, X_u, X_z, b);
  stan::math::vector_d result_dense = A * b;
  stan::math::vector_d result_straight_dense = m * b;

  EXPECT_FLOAT_EQ(result_dense(0), result_sparse(0));
  EXPECT_FLOAT_EQ(result_dense(1), result_sparse(1));
  EXPECT_FLOAT_EQ(result_dense(0), result_straight_dense(0));
  EXPECT_FLOAT_EQ(result_dense(1), result_straight_dense(1));
}

// Test that m=0 throws (CSR).
TEST(SparseStuff, csr_to_dense_matrix_m0) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();

  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);
  std::vector<int> X_z = stan::math::csr_extract_z(a);

  Eigen::Matrix<double, 2, 3> result;

  EXPECT_THROW({
  result = stan::math::csr_to_dense_matrix(0, 3, X_w, X_v, X_u, X_z);},
  std::domain_error);
}

// Test that n=0 throws (CSR).
TEST(SparseStuff, csr_to_dense_matrix_n0) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();

  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);
  std::vector<int> X_z = stan::math::csr_extract_z(a);

  Eigen::Matrix<double, 2, 3> result;

  EXPECT_THROW({
  result = stan::math::csr_to_dense_matrix(2, 0, X_w, X_v, X_u, X_z);},
  std::domain_error);
}

// Test that short u throws (CSR).
TEST(SparseStuff, csr_to_dense_matrix_u_short) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();

  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);
  std::vector<int> X_z = stan::math::csr_extract_z(a);

  X_u.erase(X_u.begin());  // make a short u:

  Eigen::Matrix<double, 2, 3> result;

  EXPECT_THROW({
  result = stan::math::csr_to_dense_matrix(2, 3, X_w, X_v, X_u, X_z);},
  std::invalid_argument);
}

// Test that short z throws (CSR).
TEST(SparseStuff, csr_to_dense_matrix_z_short) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();

  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);
  std::vector<int> X_z = stan::math::csr_extract_z(a);

  X_z.erase(X_z.begin());  // make a short z:

  Eigen::Matrix<double, 2, 3> result;

  EXPECT_THROW({
  result = stan::math::csr_to_dense_matrix(2, 3, X_w, X_v, X_u, X_z);},
  std::invalid_argument);
}

// Test that short v throws (CSR).
TEST(SparseStuff, csr_to_dense_matrix_v_short) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();

  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);
  std::vector<int> X_z = stan::math::csr_extract_z(a);

  X_v.erase(X_v.begin()+4);  // make a short v:

  Eigen::Matrix<double, 2, 3> result;

  EXPECT_THROW({
  result = stan::math::csr_to_dense_matrix(2, 3, X_w, X_v, X_u, X_z);},
  std::invalid_argument);
}

// Test that wrong z throws (CSR).
TEST(SparseStuff, csr_to_dense_matrix_z_wrong) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();

  stan::math::vector_d X_w = stan::math::csr_extract_w(a);
  std::vector<int> X_v = stan::math::csr_extract_v(a);
  std::vector<int> X_u = stan::math::csr_extract_u(a);
  std::vector<int> X_z = stan::math::csr_extract_z(a);

  X_z[X_z.size()-1] += 1;  // make a wrong z:

  Eigen::Matrix<double, 2, 3> result;

  EXPECT_THROW({
  result = stan::math::csr_to_dense_matrix(2, 3, X_w, X_v, X_u, X_z);},
  std::invalid_argument);
}


