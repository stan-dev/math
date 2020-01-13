#include <gtest/gtest.h>
#include <stan/math/prim.hpp>
#include <limits>
#include <string>
#include <vector>

TEST(MathPrimMat, zero_sizes_add_diag) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat(0, 0);
  Eigen::Matrix<double, Eigen::Dynamic, 1> to_add1(0);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> output;

  EXPECT_NO_THROW(output = stan::math::add_diag(mat, to_add1));
  EXPECT_EQ(0, output.rows());
  EXPECT_EQ(0, output.cols());
}

TEST(MathPrimMat, vector_correct_size_add_diag) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat(2, 3);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat2(0, 0);
  Eigen::Matrix<double, Eigen::Dynamic, 1> to_add1(1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> to_add2(2);
  Eigen::Matrix<double, Eigen::Dynamic, 1> to_add3(3);

  EXPECT_THROW(stan::math::add_diag(mat, to_add1), std::invalid_argument);
  EXPECT_NO_THROW(stan::math::add_diag(mat, to_add2));
  EXPECT_THROW(stan::math::add_diag(mat, to_add3), std::invalid_argument);
  EXPECT_THROW(stan::math::add_diag(mat2, to_add2), std::invalid_argument);
}

TEST(MathPrimMat, row_vector_correct_size_add_diag) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat(2, 3);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat2(0, 0);
  Eigen::Matrix<double, 1, Eigen::Dynamic> to_add1(1);
  Eigen::Matrix<double, 1, Eigen::Dynamic> to_add2(2);
  Eigen::Matrix<double, 1, Eigen::Dynamic> to_add3(3);

  EXPECT_THROW(stan::math::add_diag(mat, to_add1), std::invalid_argument);
  EXPECT_NO_THROW(stan::math::add_diag(mat, to_add2));
  EXPECT_THROW(stan::math::add_diag(mat, to_add3), std::invalid_argument);
  EXPECT_THROW(stan::math::add_diag(mat2, to_add2), std::invalid_argument);
}

TEST(MathPrimMat, double_mat_double_add_diag) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat(2, 3);
  mat << 1, 1, 1, 1, 1, 1;

  double jitter = 1e-10;

  Eigen::MatrixXd out_mat;
  EXPECT_NO_THROW(out_mat = stan::math::add_diag(mat, jitter));
  for (int i = 0; i < 2; ++i)
    EXPECT_FLOAT_EQ(1.0 + jitter, out_mat(i, i))
        << "index: ( " << i << ", " << i << ")";
}

TEST(MathPrimMat, double_mat_double_vec_add_diag) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat(2, 3);
  mat << 1, 1, 1, 1, 1, 1;

  Eigen::Matrix<double, Eigen::Dynamic, 1> to_add(2);
  to_add << 0, 1;

  Eigen::MatrixXd out_mat;
  EXPECT_NO_THROW(out_mat = stan::math::add_diag(mat, to_add));
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (i == j)
        EXPECT_FLOAT_EQ(1 + to_add(i), out_mat(i, i))
            << "index: ( " << i << ", " << i << ")";
      else
        EXPECT_FLOAT_EQ(1, out_mat(i, j))
            << "index: ( " << i << ", " << i << ")";
    }
  }
}

TEST(MathPrimMat, double_mat_double_rvec_add_diag) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat(2, 3);
  mat << 1, 1, 1, 1, 1, 1;

  Eigen::Matrix<double, 1, Eigen::Dynamic> to_add(2);
  to_add << 0, 1;

  Eigen::MatrixXd out_mat;
  EXPECT_NO_THROW(out_mat = stan::math::add_diag(mat, to_add));
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (i == j)
        EXPECT_FLOAT_EQ(1 + to_add(i), out_mat(i, j))
            << "index: ( " << i << ", " << i << ")";
      else
        EXPECT_FLOAT_EQ(1, out_mat(i, j))
            << "index: ( " << i << ", " << i << ")";
    }
  }
}
