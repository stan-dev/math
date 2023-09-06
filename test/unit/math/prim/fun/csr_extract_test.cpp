#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

// Test that values from a dense matrix in sparse format are extracted.
TEST(SparseStuff, csr_extract_dense) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();

  stan::math::vector_d w;
  std::vector<int> v, u;
  std::tie(w, v, u) = stan::math::csr_extract(a);

  EXPECT_FLOAT_EQ(2.0, w(0));
  EXPECT_FLOAT_EQ(4.0, w(1));
  EXPECT_FLOAT_EQ(6.0, w(2));
  EXPECT_FLOAT_EQ(8.0, w(3));
  EXPECT_FLOAT_EQ(10.0, w(4));
  EXPECT_FLOAT_EQ(12.0, w(5));

  EXPECT_EQ(1, v[0]);
  EXPECT_EQ(2, v[1]);
  EXPECT_EQ(3, v[2]);
  EXPECT_EQ(1, v[3]);
  EXPECT_EQ(2, v[4]);
  EXPECT_EQ(3, v[5]);

  EXPECT_EQ(1, u[0]);
  EXPECT_EQ(4, u[1]);
  EXPECT_EQ(7, u[2]);
}

// Test that values from a dense matrix are extracted.
TEST(SparseStuff, csr_extract_dense_dense) {
  stan::math::matrix_d m(2, 3);
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  stan::math::vector_d w;
  std::vector<int> v, u;
  std::tie(w, v, u) = stan::math::csr_extract(m);

  EXPECT_FLOAT_EQ(2.0, w(0));
  EXPECT_FLOAT_EQ(4.0, w(1));
  EXPECT_FLOAT_EQ(6.0, w(2));
  EXPECT_FLOAT_EQ(8.0, w(3));
  EXPECT_FLOAT_EQ(10.0, w(4));
  EXPECT_FLOAT_EQ(12.0, w(5));

  EXPECT_EQ(1, v[0]);
  EXPECT_EQ(2, v[1]);
  EXPECT_EQ(3, v[2]);
  EXPECT_EQ(1, v[3]);
  EXPECT_EQ(2, v[4]);
  EXPECT_EQ(3, v[5]);

  EXPECT_EQ(1, u[0]);
  EXPECT_EQ(4, u[1]);
  EXPECT_EQ(7, u[2]);
}

// Test that values from a dense matrix in sparse format are extracted
// after A.makeCompressed();
TEST(SparseStuff, csr_extract_dense_compressed) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 8.0, 10.0, 12.0;
  a = m.sparseView();
  a.makeCompressed();
  stan::math::vector_d w;
  std::vector<int> v, u;
  std::tie(w, v, u) = stan::math::csr_extract(a);

  EXPECT_FLOAT_EQ(2.0, w(0));
  EXPECT_FLOAT_EQ(4.0, w(1));
  EXPECT_FLOAT_EQ(6.0, w(2));
  EXPECT_FLOAT_EQ(8.0, w(3));
  EXPECT_FLOAT_EQ(10.0, w(4));
  EXPECT_FLOAT_EQ(12.0, w(5));

  EXPECT_EQ(1, v[0]);
  EXPECT_EQ(2, v[1]);
  EXPECT_EQ(3, v[2]);
  EXPECT_EQ(1, v[3]);
  EXPECT_EQ(2, v[4]);
  EXPECT_EQ(3, v[5]);

  EXPECT_EQ(1, u[0]);
  EXPECT_EQ(4, u[1]);
  EXPECT_EQ(7, u[2]);
}

// Test that values from a sparse matrix in sparse format are extracted
TEST(SparseStuff, csr_extract_sparse) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 0.0, 0.0, 0.0;
  a = m.sparseView();
  stan::math::vector_d w;
  std::vector<int> v, u;
  std::tie(w, v, u) = stan::math::csr_extract(a);

  EXPECT_FLOAT_EQ(2.0, w(0));
  EXPECT_FLOAT_EQ(4.0, w(1));
  EXPECT_FLOAT_EQ(6.0, w(2));

  EXPECT_EQ(1, v[0]);
  EXPECT_EQ(2, v[1]);
  EXPECT_EQ(3, v[2]);
  EXPECT_EQ(3, v.size());

  EXPECT_EQ(1, u[0]);
  EXPECT_EQ(4, u[1]);
  EXPECT_EQ(4, u[2]);
  EXPECT_EQ(3, u.size());
}

// Test that values from a sparse matrix in sparse format are extracted
// after A.makeCompressed()
TEST(SparseStuff, csr_extract_sparse_compressed) {
  stan::math::matrix_d m(2, 3);
  Eigen::SparseMatrix<double, Eigen::RowMajor> a;
  m << 2.0, 4.0, 6.0, 0.0, 0.0, 0.0;
  a = m.sparseView();
  stan::math::vector_d w;
  std::vector<int> v, u;
  std::tie(w, v, u) = stan::math::csr_extract(a);
  EXPECT_FLOAT_EQ(2.0, w(0));
  EXPECT_FLOAT_EQ(4.0, w(1));
  EXPECT_FLOAT_EQ(6.0, w(2));

  EXPECT_EQ(1, v[0]);
  EXPECT_EQ(2, v[1]);
  EXPECT_EQ(3, v[2]);
  EXPECT_EQ(3, v.size());

  EXPECT_EQ(1, u[0]);
  EXPECT_EQ(4, u[1]);
  EXPECT_EQ(4, u[2]);
  EXPECT_EQ(3, u.size());
}
