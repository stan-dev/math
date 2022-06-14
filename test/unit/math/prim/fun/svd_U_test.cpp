#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <stdexcept>
#include <stan/math/prim/fun/Eigen.hpp>
#include <complex>

TEST(MathMatrixPrimMat, svd_U) {
  using stan::math::matrix_d;
  using stan::math::svd_U;
  using compl_t = std::complex<double>;
  using matrix_c = Eigen::Matrix<compl_t, Eigen::Dynamic, Eigen::Dynamic>;

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

  matrix_c c32(3, 2);
  c32 << compl_t(0.86636546, 0.34306449), compl_t(0.28267243, 0.52462912),
      compl_t(0.12104914, 0.2533793), compl_t(0.66889264, 0.39276455),
      compl_t(0.02184348, 0.0614428), compl_t(0.96599692, 0.16180684);
  matrix_c c32_U(3, 2);
  c32_U << compl_t(0.50789057, 0.35782384), compl_t(0.74507868, 0.14261495),
      compl_t(0.4823205, 0.19300139), compl_t(-0.29600299, 0.17466116),
      compl_t(0.58489862, -0.04494745), compl_t(-0.53765935, 0.13159357);

  EXPECT_MATRIX_COMPLEX_FLOAT_EQ(c32_U, svd_U(c32));
}
