#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <test/unit/util.hpp>

// this is just testing the nan behavior of the built-in fma
// there is no longer a stan::math::fma, just the agrad versions
// instead, the top-level ::fma should be used by including <cmath>

TEST(MathFunctions, fma) {
  using stan::math::fma;
  EXPECT_FLOAT_EQ(5.0, fma(1.0, 2.0, 3.0));
  EXPECT_FLOAT_EQ(10.0, fma(2.0, 3.0, 4.0));
  EXPECT_FLOAT_EQ(
      11.0, fma(static_cast<int>(3), static_cast<int>(2), static_cast<int>(5)));
}

TEST(MathFunctions, fma_nan) {
  using stan::math::fma;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(fma(1.0, 2.0, nan)));

  EXPECT_TRUE(std::isnan(fma(1.0, nan, 3.0)));

  EXPECT_TRUE(std::isnan(fma(1.0, nan, nan)));

  EXPECT_TRUE(std::isnan(fma(nan, 2.0, 3.0)));

  EXPECT_TRUE(std::isnan(fma(nan, 2.0, nan)));

  EXPECT_TRUE(std::isnan(fma(nan, nan, 3.0)));

  EXPECT_TRUE(std::isnan(fma(nan, nan, nan)));
}

TEST(MathFunctions, fma_matrix) {
  using stan::math::add;
  using stan::math::elt_multiply;
  using stan::math::fma;

  double x_scalar = 1.0;
  Eigen::VectorXd x_vec(2);
  x_vec << 1.0, 2.0;
  Eigen::RowVectorXd x_rowvec(2);
  x_rowvec << 1.0, 2.0;
  Eigen::MatrixXd x_mat(2, 2);
  x_mat << 1.0, 2.0, -1.0, 1.1;

  double y_scalar = 2.0;
  Eigen::VectorXd y_vec(2);
  y_vec << 2.0, -3.0;
  Eigen::RowVectorXd y_rowvec(2);
  y_rowvec << 2.0, -3.0;
  Eigen::MatrixXd y_mat(2, 2);
  x_mat << 1.0, 2.0, -1.0, 1.1;

  double z_scalar = 3.0;
  Eigen::VectorXd z_vec(2);
  z_vec << -3.0, 4.0;
  Eigen::RowVectorXd z_rowvec(2);
  z_rowvec << -3.0, 4.0;
  Eigen::MatrixXd zm(2, 2);
  x_mat << 3.0, 4.0, -1.0, 1.1;

  EXPECT_MATRIX_EQ(add(elt_multiply(x_scalar, y_scalar), z_vec),
                   fma(x_scalar, y_scalar, z_vec));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_scalar, y_vec), z_scalar),
                   fma(x_scalar, y_vec, z_scalar));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_scalar, y_vec), z_vec),
                   fma(x_scalar, y_vec, z_vec));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_vec, y_scalar), z_scalar),
                   fma(x_vec, y_scalar, z_scalar));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_vec, y_scalar), z_vec),
                   fma(x_vec, y_scalar, z_vec));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_vec, y_vec), z_scalar),
                   fma(x_vec, y_vec, z_scalar));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_vec, y_vec), z_vec),
                   fma(x_vec, y_vec, z_vec));

  EXPECT_MATRIX_EQ(add(elt_multiply(x_scalar, y_scalar), z_rowvec),
                   fma(x_scalar, y_scalar, z_rowvec));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_scalar, y_rowvec), z_scalar),
                   fma(x_scalar, y_rowvec, z_scalar));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_scalar, y_rowvec), z_rowvec),
                   fma(x_scalar, y_rowvec, z_rowvec));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_rowvec, y_scalar), z_scalar),
                   fma(x_rowvec, y_scalar, z_scalar));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_rowvec, y_scalar), z_rowvec),
                   fma(x_rowvec, y_scalar, z_rowvec));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_rowvec, y_rowvec), z_scalar),
                   fma(x_rowvec, y_rowvec, z_scalar));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_rowvec, y_rowvec), z_rowvec),
                   fma(x_rowvec, y_rowvec, z_rowvec));

  EXPECT_MATRIX_EQ(add(elt_multiply(x_scalar, y_scalar), zm),
                   fma(x_scalar, y_scalar, zm));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_scalar, y_mat), z_scalar),
                   fma(x_scalar, y_mat, z_scalar));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_scalar, y_mat), zm),
                   fma(x_scalar, y_mat, zm));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_mat, y_scalar), z_scalar),
                   fma(x_mat, y_scalar, z_scalar));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_mat, y_scalar), zm),
                   fma(x_mat, y_scalar, zm));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_mat, y_mat), z_scalar),
                   fma(x_mat, y_mat, z_scalar));
  EXPECT_MATRIX_EQ(add(elt_multiply(x_mat, y_mat), zm), fma(x_mat, y_mat, zm));
}

TEST(MathFunctions, fma_matrix_error) {
  using stan::math::fma;
  double x_scalar = 1.0;
  Eigen::RowVectorXd x_rowvec(2);
  x_rowvec << 1.0, 2.0;

  double y_scalar = 2.0;
  Eigen::VectorXd y_vec(2);
  y_vec << 2.0, -3.0;
  Eigen::RowVectorXd y_rowvec(2);
  y_rowvec << 2.0, -3.0;

  double z_scalar = 3.0;
  Eigen::VectorXd z_vec(2);
  z_vec << -3.0, 4.0;

  EXPECT_THROW(fma(x_scalar, y_rowvec, z_vec), std::invalid_argument);
  EXPECT_THROW(fma(x_rowvec, y_scalar, z_vec), std::invalid_argument);
  EXPECT_THROW(fma(x_rowvec, y_vec, z_scalar), std::invalid_argument);
}
