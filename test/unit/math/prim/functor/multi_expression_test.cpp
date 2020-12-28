#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <tuple>

TEST(MathFunctions, multi_expression_mismatched_sizes) {
  using RowMat
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  Eigen::MatrixXd a(2, 2);
  Eigen::MatrixXd b(2, 3);
  RowMat c(2, 3);
  RowMat d(3, 2);
  EXPECT_THROW(stan::math::eigen_expressions(a, b), std::invalid_argument);
  EXPECT_THROW(stan::math::eigen_expressions(b, c), std::invalid_argument);
  EXPECT_NO_THROW(stan::math::eigen_expressions(b, d));
}

TEST(MathFunctions, multi_expression_all_col_major) {
  Eigen::MatrixXd a(2, 3);
  a << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd b(2, 3);
  b << 1, 2, 3, 5, 6, 2;
  Eigen::MatrixXd c, d;
  stan::math::eigen_results(c, d) = stan::math::eigen_expressions(a * 2, b + a);

  Eigen::MatrixXd c_res = a * 2;
  Eigen::MatrixXd d_res = b + a;
  EXPECT_MATRIX_EQ(c, c_res);
  EXPECT_MATRIX_EQ(d, d_res);
}

TEST(MathFunctions, multi_expression_all_row_major) {
  using RowMat
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  RowMat a(2, 3);
  a << 1, 2, 3, 4, 4, 5;
  RowMat b(2, 3);
  b << 1, 2, 3, 5, 7, 8;
  RowMat c, d;
  stan::math::eigen_results(c, d) = stan::math::eigen_expressions(a * 2, b + a);

  RowMat c_res = a * 2;
  RowMat d_res = b + a;
  EXPECT_MATRIX_EQ(c, c_res);
  EXPECT_MATRIX_EQ(d, d_res);
}

TEST(MathFunctions, multi_expression_more_row_major) {
  using RowMat
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  Eigen::MatrixXd a(2, 2);
  a << 1, 2, 3, 4;
  RowMat b(2, 2);
  b << 1, 2, 3, 5;
  RowMat c, d;
  stan::math::eigen_results(c, d) = stan::math::eigen_expressions(a * 2, b + a);

  RowMat c_res = a * 2;
  RowMat d_res = b + a;
  EXPECT_MATRIX_EQ(c, c_res);
  EXPECT_MATRIX_EQ(d, d_res);
}

TEST(MathFunctions, multi_expression_more_col_major) {
  using RowMat
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  Eigen::MatrixXd a(2, 2);
  a << 1, 2, 3, 4;
  RowMat b(2, 2);
  b << 1, 2, 3, 5;
  Eigen::MatrixXd c, d;
  stan::math::eigen_results(c, d) = stan::math::eigen_expressions(a * 2, b + a);

  Eigen::MatrixXd c_res = a * 2;
  Eigen::MatrixXd d_res = b + a;
  EXPECT_MATRIX_EQ(c, c_res);
  EXPECT_MATRIX_EQ(d, d_res);
}

TEST(MathFunctions, multi_expression_with_block_more_row_major) {
  using RowMat
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  Eigen::MatrixXd a(2, 2);
  a << 1, 2, 3, 4;
  RowMat b(3, 3);
  b << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  RowMat c, d;
  stan::math::eigen_results(c, d)
      = stan::math::eigen_expressions(a * 2, b.block(1, 0, 2, 2) + a);

  RowMat c_res = a * 2;
  RowMat d_res = b.block(1, 0, 2, 2) + a;
  EXPECT_MATRIX_EQ(c, c_res);
  EXPECT_MATRIX_EQ(d, d_res);
}

TEST(MathFunctions, multi_expression_with_block_more_col_major) {
  using RowMat
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  Eigen::MatrixXd a(2, 2);
  a << 1, 2, 3, 4;
  RowMat b(3, 3);
  b << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  Eigen::MatrixXd c, d;
  stan::math::eigen_results(c, d)
      = stan::math::eigen_expressions(a * 2, b.block(1, 0, 2, 2) + a);

  Eigen::MatrixXd c_res = a * 2;
  Eigen::MatrixXd d_res = b.block(1, 0, 2, 2) + a;
  EXPECT_MATRIX_EQ(c, c_res);
  EXPECT_MATRIX_EQ(d, d_res);
}

TEST(MathFunctions, multi_expression_transpose_more_row_major) {
  using RowMat
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  RowMat a(2, 2);
  a << 1, 2, 3, 4;
  RowMat c, d;
  stan::math::eigen_results(c, d)
      = stan::math::eigen_expressions(a * 2, a.transpose());

  RowMat c_res = a * 2;
  RowMat d_res = a.transpose();
  EXPECT_MATRIX_EQ(c, c_res);
  EXPECT_MATRIX_EQ(d, d_res);
}

TEST(MathFunctions, multi_expression_transpose_more_col_major) {
  using RowMat
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  Eigen::MatrixXd a(2, 2);
  a << 1, 2, 3, 4;
  Eigen::MatrixXd c, d;
  stan::math::eigen_results(c, d)
      = stan::math::eigen_expressions(a * 2, a.transpose());

  Eigen::MatrixXd c_res = a * 2;
  Eigen::MatrixXd d_res = a.transpose();
  EXPECT_MATRIX_EQ(c, c_res);
  EXPECT_MATRIX_EQ(d, d_res);
}

TEST(MathFunctions, multi_expression_many_expressions) {
  using RowMat
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  Eigen::MatrixXd a(2, 3);
  a << 1, 2, 3, 4, 5, 6;
  RowMat b(3, 3);
  b << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  Eigen::MatrixXd c(4, 4);
  c << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
  Eigen::MatrixXd d, e(c);
  RowMat f, g;
  stan::math::eigen_results(d, e.block(2, 1, 2, 3), f, g)
      = stan::math::eigen_expressions(
          a * 2, b.block(0, 1, 3, 2).transpose() + a, c.block(2, 1, 2, 3) + a,
          b.block(0, 1, 3, 2).transpose().array() - 2);

  Eigen::MatrixXd d_res = a * 2;
  Eigen::MatrixXd e_res = c;
  e_res.block(2, 1, 2, 3) = b.block(0, 1, 3, 2).transpose() + a;
  RowMat f_res = c.block(2, 1, 2, 3) + a;
  RowMat g_res = b.block(0, 1, 3, 2).transpose().array() - 2;
  EXPECT_MATRIX_EQ(d, d_res);
  EXPECT_MATRIX_EQ(e, e_res);
  EXPECT_MATRIX_EQ(f, f_res);
  EXPECT_MATRIX_EQ(g, g_res);
}

TEST(MathFunctions, multi_expression_many_expressions2) {
  using RowMat
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  Eigen::MatrixXd a(2, 2);
  a << 1, 2, 3, 4;
  RowMat b(3, 3);
  b << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  Eigen::MatrixXd c(4, 4);
  c << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
  Eigen::MatrixXd d, e(c);
  RowMat f, g;
  stan::math::eigen_results(d, e.block(2, 1, 2, 2), f, g)
      = stan::math::eigen_expressions(
          a * 2, b.block(0, 1, 2, 2).transpose() + a,
          c.block(2, 1, 2, 2) + a.transpose(), b.block(1, 0, 2, 2).array() - 2);

  Eigen::MatrixXd d_res = a * 2;
  Eigen::MatrixXd e_res = c;
  e_res.block(2, 1, 2, 2) = b.block(0, 1, 2, 2).transpose() + a;
  RowMat f_res = c.block(2, 1, 2, 2) + a.transpose();
  RowMat g_res = b.block(1, 0, 2, 2).array() - 2;
  EXPECT_MATRIX_EQ(d, d_res);
  EXPECT_MATRIX_EQ(e, e_res);
  EXPECT_MATRIX_EQ(f, f_res);
  EXPECT_MATRIX_EQ(g, g_res);
}

TEST(MathFunctions, multi_expression_compound_addition_assignment) {
  Eigen::MatrixXd a(2, 3);
  a << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd b(2, 3);
  b << 1, 2, 3, 5, 6, 2;
  Eigen::MatrixXd c(2, 3);
  a << 1, 2, -3, 4, -5, 6;
  Eigen::MatrixXd d(2, 3);
  b << 1, -2, 3, -5, 6, -2;
  Eigen::MatrixXd c_res = c;
  Eigen::MatrixXd d_res = d;
  stan::math::eigen_results(c, d)
      += stan::math::eigen_expressions(a * 2, b + a);

  c_res += a * 2;
  d_res += b + a;
  EXPECT_MATRIX_EQ(c, c_res);
  EXPECT_MATRIX_EQ(d, d_res);
}

TEST(MathFunctions, multi_expression_matrix_product) {
  Eigen::MatrixXd a(2, 3);
  a << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd b(2, 3);
  b << 1, 2, 3, 5, 6, 2;
  Eigen::MatrixXd b2(3, 3);
  b2 << 1, 2, 3, 5, 1, 6, 9, 8, 7;
  Eigen::MatrixXd c, d;
  stan::math::eigen_results(c, d)
      = stan::math::eigen_expressions(a * 2, b + a * b2);

  Eigen::MatrixXd c_res = a * 2;
  Eigen::MatrixXd d_res = b + a * b2;
  EXPECT_MATRIX_EQ(c, c_res);
  EXPECT_MATRIX_EQ(d, d_res);
}
