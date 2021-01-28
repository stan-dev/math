#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <stdexcept>

TEST(MathMatrixPrimMat, svd_U) {
  using stan::math::matrix_d;
  using stan::math::svd_U;

  // Values generated using R base::svd

  matrix_d m00(0, 0);
  EXPECT_THROW(svd_U(m00), std::invalid_argument);

  matrix_d m11(1, 1);
  m11 << 5;
  matrix_d m11_U(1, 1);
  m11_U << 1;
  EXPECT_MATRIX_FLOAT_EQ(m11_U, svd_U(m11));

  matrix_d m22(2, 2);
  m22 << 1, 9, -4, 2;
  matrix_d m22_U(2, 2);
  m22_U << 0.97759158130035961, 0.21051057021124275, 0.21051057021124264,
      -0.97759158130035972;
  EXPECT_MATRIX_FLOAT_EQ(m22_U, svd_U(m22));

  matrix_d m23(2, 3);
  m23 << 1, 3, -5, 7, 9, -11;
  matrix_d m23_U(2, 2);
  m23_U << -0.33784321032557019, -0.94120240396894062, -0.94120240396894039,
      0.33784321032557019;
  EXPECT_MATRIX_FLOAT_EQ(m23_U, svd_U(m23));

  matrix_d m32(3, 2);
  m32 << 1, 3, -5, 7, 9, -11;
  matrix_d m32_U(3, 2);
  m32_U << 0.10657276921942949, 0.978015199679778457, 0.51489182605641137,
      0.099934023569151917, -0.85060487438128174, 0.183028577355020677;
  EXPECT_MATRIX_FLOAT_EQ(
      m32_U, svd_U(m32));  // R's SVD returns different signs than Eigen.
}
