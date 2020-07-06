#include <stan/math/rev/core.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <sstream>

TEST(AgradRev, insertion_operator) {
  stan::math::vari v(5);
  std::stringstream ss;
  ss << &v;
  EXPECT_EQ("5:0", ss.str());
}

TEST(AgradRev, long_double_test) {
  stan::math::vari_value<long double> v(5);
  std::stringstream ss;
  ss << &v;
  EXPECT_EQ("5:0", ss.str());
}

TEST(AgradRev, dense_matrix_vari) {
  using stan::math::vari_value;
  using eig_mat = Eigen::MatrixXd;
  vari_value<eig_mat> A_vari(eig_mat::Random(3, 3));
  eig_mat B(eig_mat::Random(3, 3));
  vari_value<eig_mat> B_vari(B);
  EXPECT_MATRIX_FLOAT_EQ(B, B_vari.val_);
}
TEST(AgradRev, dense_vector_vari) {
  using stan::math::vari_value;
  using eig_vec = Eigen::Matrix<double, -1, 1>;
  vari_value<eig_vec> A_vari(eig_vec::Random(3));
  eig_vec B(eig_vec::Random(3));
  vari_value<eig_vec> B_vari(B);
  EXPECT_MATRIX_FLOAT_EQ(B, B_vari.val_);
}

TEST(AgradRev, dense_row_vector_vari) {
  using stan::math::vari_value;
  using eig_row_vec = Eigen::Matrix<double, 1, -1>;
  vari_value<eig_row_vec> A_vari(eig_row_vec::Random(3));
  eig_row_vec B(eig_row_vec::Random(3));
  vari_value<eig_row_vec> B_vari(B);
  EXPECT_MATRIX_FLOAT_EQ(B, B_vari.val_);
}

TEST(AgradRev, sparse_matrix_vari) {
  using stan::math::vari_value;
  using eig_mat = Eigen::SparseMatrix<double>;
  using inner_iterator = typename eig_mat::InnerIterator;
  using stan::test::make_sparse_matrix_random;
  vari_value<eig_mat> A_vari(make_sparse_matrix_random(10, 10));
  eig_mat B = make_sparse_matrix_random(10, 10);
  vari_value<eig_mat> B_vari(B);
  for (int k = 0; k < B.outerSize(); ++k) {
    for (inner_iterator it(B, k), iz(B_vari.val_, k); it; ++it, ++iz) {
      EXPECT_FLOAT_EQ(iz.value(), it.value());
    }
  }
}
