#include <gtest/gtest.h>
#include <stan/math/rev/mat.hpp>
#include <limits>
#include <string>
#include <vector>

TEST(MathPrimMat, var_mat_double_add_diag) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> mat(2, 3);
  mat << 1, 1, 1, 1, 1, 1;

  double jitter = 1e-10;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> out_mat;
  EXPECT_NO_THROW(out_mat = stan::math::add_diag(mat, jitter));
  for (int i = 0; i < 2; ++i)
    EXPECT_FLOAT_EQ(1.0 + jitter, value_of(out_mat(i, i)))
        << "index: ( " << i << ", " << i << ")";
}

TEST(MathPrimMat, double_mat_var_add_diag) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat(2, 3);
  mat << 1, 1, 1, 1, 1, 1;

  stan::math::var jitter = 1e-10;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> out_mat;
  EXPECT_NO_THROW(out_mat = stan::math::add_diag(mat, jitter));
  for (int i = 0; i < 2; ++i)
    EXPECT_FLOAT_EQ(1.0 + value_of(jitter), value_of(out_mat(i, i)))
        << "index: ( " << i << ", " << i << ")";
}

TEST(MathPrimMat, var_mat_var_add_diag) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> mat(2, 3);
  mat << 1, 1, 1, 1, 1, 1;

  stan::math::var jitter = 1e-10;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> out_mat;
  EXPECT_NO_THROW(out_mat = stan::math::add_diag(mat, jitter));
  for (int i = 0; i < 2; ++i)
    EXPECT_FLOAT_EQ(1.0 + value_of(jitter), value_of(out_mat(i, i)))
        << "index: ( " << i << ", " << i << ")";
}

TEST(MathPrimMat, double_mat_var_vec_add_diag) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat(2, 3);
  mat << 1, 1, 1, 1, 1, 1;

  Eigen::Matrix<stan::math::var, 1, -1> to_add(2);
  to_add << 0, 1;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> out_mat;
  EXPECT_NO_THROW(out_mat = stan::math::add_diag(mat, to_add));
  for (int i = 0; i < 2; ++i)
    EXPECT_FLOAT_EQ(1 + value_of(to_add[i]), value_of(out_mat(i, i)))
        << "index: ( " << i << ", " << i << ")";
}

TEST(MathPrimMat, double_mat_var_rvec_add_diag) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat(2, 3);
  mat << 1, 1, 1, 1, 1, 1;

  Eigen::Matrix<stan::math::var, -1, 1> to_add(2);
  to_add << 0, 1;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> out_mat;
  EXPECT_NO_THROW(out_mat = stan::math::add_diag(mat, to_add));
  for (int i = 0; i < 2; ++i)
    EXPECT_FLOAT_EQ(1 + value_of(to_add[i]), value_of(out_mat(i, i)))
        << "index: ( " << i << ", " << i << ")";
}
