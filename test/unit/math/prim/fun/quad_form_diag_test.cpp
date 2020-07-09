#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrim, quadFormDiag) {
  using stan::math::quad_form_diag;
  Eigen::MatrixXd m(1, 1);
  m << 3;

  Eigen::VectorXd v(1);
  v << 9;

  Eigen::MatrixXd v_m = v.asDiagonal();

  EXPECT_MATRIX_FLOAT_EQ(v_m * m * v_m, quad_form_diag(m, v));
}

TEST(MathMatrixPrim, quadFormDiag2) {
  using stan::math::quad_form_diag;
  Eigen::MatrixXd m(3, 3);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9;

  Eigen::VectorXd v(3);
  v << 1, 2, 3;

  Eigen::MatrixXd v_m = v.asDiagonal();

  EXPECT_MATRIX_FLOAT_EQ(v_m * m * v_m, quad_form_diag(m, v));

  Eigen::RowVectorXd rv(3);
  rv << 1, 2, 3;
  EXPECT_MATRIX_FLOAT_EQ(v_m * m * v_m, quad_form_diag(m, rv));
}

TEST(MathMatrixPrim, quadFormDiagException) {
  using stan::math::quad_form_diag;
  Eigen::MatrixXd m(2, 2);
  m << 2, 3, 4, 5;
  Eigen::VectorXd v(3);
  v << 1, 2, 3;
  EXPECT_THROW(quad_form_diag(m, v), std::invalid_argument);

  Eigen::MatrixXd m2(3, 2);
  m2 << 2, 3, 4, 5, 6, 7;

  Eigen::VectorXd v2(2);
  v2 << 1, 2;

  EXPECT_THROW(quad_form_diag(m2, v), std::invalid_argument);
  EXPECT_THROW(quad_form_diag(m2, v2), std::invalid_argument);
}
