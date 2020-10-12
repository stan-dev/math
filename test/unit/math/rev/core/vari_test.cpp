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

TEST(AgradRev, arena_matrix_matrix_vari) {
  using stan::math::arena_matrix;
  using stan::math::vari_value;
  arena_matrix<Eigen::MatrixXd> x(Eigen::MatrixXd::Random(5, 5));
  const auto& x_ref = x;
  auto* A = new vari_value<Eigen::MatrixXd>(x);
  EXPECT_MATRIX_FLOAT_EQ((*A).val_, x);
  auto* B = new vari_value<Eigen::MatrixXd>(x_ref);
  EXPECT_MATRIX_FLOAT_EQ((*B).val_, x);
  auto* C = new vari_value<Eigen::MatrixXd>(x, true);
  EXPECT_MATRIX_FLOAT_EQ((*C).val_, x);
  auto* D = new vari_value<Eigen::MatrixXd>(x_ref, true);
  EXPECT_MATRIX_FLOAT_EQ((*D).val_, x);
}

TEST(AgradRev, dense_vari_matrix_views) {
  using stan::math::vari_value;
  using eig_mat = Eigen::MatrixXd;
  eig_mat A(5, 5);
  for (int i = 0; i < A.size(); ++i) {
    A(i) = i;
  }
  stan::math::vari_value<eig_mat> A_v(A);
  auto A_head = A_v.block(1, 1, 3, 3);
  EXPECT_MATRIX_FLOAT_EQ(A_head.val_, A_v.val_.block(1, 1, 3, 3));
  EXPECT_MATRIX_FLOAT_EQ(A, A_v.val_);

  auto A_row = A_v.row(3);
  EXPECT_MATRIX_FLOAT_EQ(A_row.val_, A_v.val_.row(3));
  EXPECT_MATRIX_FLOAT_EQ(A, A_v.val_);

  auto A_col = A_v.col(3);
  EXPECT_MATRIX_FLOAT_EQ(A_col.val_, A_v.val_.col(3));
  EXPECT_MATRIX_FLOAT_EQ(A, A_v.val_);

  auto A_op_par = A_v(3);
  EXPECT_FLOAT_EQ(A_op_par.val_, A_v.val_(3));
  EXPECT_MATRIX_FLOAT_EQ(A, A_v.val_);

  auto A_op_par2 = A_v(3, 3);
  EXPECT_FLOAT_EQ(A_op_par2.val_, A_v.val_(3, 3));
  EXPECT_MATRIX_FLOAT_EQ(A, A_v.val_);

  auto A_op_coeff = A_v.coeff(3, 3);
  EXPECT_FLOAT_EQ(A_op_par2.val_, A_v.val_.coeff(3, 3));
  EXPECT_MATRIX_FLOAT_EQ(A, A_v.val_);
}

TEST(AgradRev, dense_vari_vector_views) {
  using stan::math::vari_value;
  using eig_vec = Eigen::VectorXd;
  eig_vec A(10);
  for (int i = 0; i < A.size(); ++i) {
    A(i) = i;
  }
  stan::math::vari_value<eig_vec> A_v(A);
  auto A_sub = A_v.head(3);
  EXPECT_MATRIX_FLOAT_EQ(A_sub.val_, A_v.val_.head(3));
  EXPECT_MATRIX_FLOAT_EQ(A, A_v.val_);

  auto A_sub_tail = A_v.tail(3);
  EXPECT_MATRIX_FLOAT_EQ(A_sub_tail.val_, A_v.val_.tail(3));
  EXPECT_MATRIX_FLOAT_EQ(A, A_v.val_);

  auto A_segment = A_v.segment(3, 5);
  EXPECT_MATRIX_FLOAT_EQ(A_segment.val_, A_v.val_.segment(3, 5));
  EXPECT_MATRIX_FLOAT_EQ(A, A_v.val_);
}
