#include <stan/math.hpp>
#include <stan/math/prim.hpp>
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
  stan::math::vari_value<Eigen::MatrixXd> A(Eigen::MatrixXd::Random(3, 3).eval());
}

TEST(AgradRev, dense_matrix_views) {
  using stan::math::vari_value;
  using eig_mat = Eigen::MatrixXd;
  eig_mat A(10, 10);
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
}

TEST(AgradRev, dense_vector_views) {
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


TEST(AgradRev, sparse_matrix_views) {
  using stan::math::vari_value;
  using eig_mat = Eigen::SparseMatrix<double>;
  using inner_iterator = typename eig_mat::InnerIterator;
  eig_mat A = test::make_sparse_matrix_random(10, 10);
  stan::math::vari_value<eig_mat> A_v(A);
  eig_mat A_block = A.block(1, 1, 3, 3);
  auto A_v_block = A_v.block(1, 1, 3, 3);
  for (int k = 0; k < A_block.outerSize(); ++k) {
    for (inner_iterator it(A_block, k), iz(A_v_block.val_, k); it; ++it, ++iz) {
      EXPECT_FLOAT_EQ(iz.value(), it.value());
    }
  }
}
