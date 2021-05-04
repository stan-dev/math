#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrim, mdivide_left_tri_val) {
  using stan::math::mdivide_left_tri;
  stan::math::matrix_d I = Eigen::MatrixXd::Identity(2, 2);

  stan::math::matrix_d Ad(2, 2);
  Ad << 2.0, 0.0, 5.0, 7.0;
  EXPECT_MATRIX_FLOAT_EQ(I, mdivide_left_tri<Eigen::Lower>(Ad, Ad));

  stan::math::matrix_d A_Ainv = Ad * mdivide_left_tri<Eigen::Lower>(Ad);
  EXPECT_MATRIX_NEAR(I, A_Ainv, 1e-15);

  Ad << 2.0, 3.0, 0.0, 7.0;
  EXPECT_MATRIX_FLOAT_EQ(I, mdivide_left_tri<Eigen::Upper>(Ad, Ad));
}

TEST(MathMatrixPrim, mdivide_left_tri_size_zero) {
  using stan::math::mdivide_left_tri;
  stan::math::matrix_d Ad(0, 0);
  stan::math::matrix_d b0(0, 2);
  stan::math::matrix_d I;

  I = mdivide_left_tri<Eigen::Lower>(Ad, Ad);
  EXPECT_EQ(0, I.rows());
  EXPECT_EQ(0, I.cols());

  I = mdivide_left_tri<Eigen::Upper>(Ad, Ad);
  EXPECT_EQ(0, I.rows());
  EXPECT_EQ(0, I.cols());

  I = mdivide_left_tri<Eigen::Lower>(Ad);
  EXPECT_EQ(0, I.rows());
  EXPECT_EQ(0, I.cols());

  I = mdivide_left_tri<Eigen::Upper>(Ad);
  EXPECT_EQ(0, I.rows());
  EXPECT_EQ(0, I.cols());

  I = mdivide_left_tri<Eigen::Lower>(Ad, b0);
  EXPECT_EQ(0, I.rows());
  EXPECT_EQ(b0.cols(), I.cols());

  I = mdivide_left_tri<Eigen::Upper>(Ad, b0);
  EXPECT_EQ(0, I.rows());
  EXPECT_EQ(b0.cols(), I.cols());
}
