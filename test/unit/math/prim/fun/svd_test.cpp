#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <stdexcept>
#include <stan/math/prim/fun/Eigen.hpp>
#include <complex>

TEST(MathMatrixPrimMat, svd) {
  using stan::math::matrix_d;
  using stan::math::svd;
  using stan::math::vector_d;
  using compl_t = std::complex<double>;
  using matrix_c = Eigen::Matrix<compl_t, Eigen::Dynamic, Eigen::Dynamic>;

  // Values generated using R base::svd

  matrix_d m00(0, 0);

  EXPECT_NO_THROW(svd(m00));

  matrix_d m11(1, 1);
  m11 << 5;
  matrix_d m11_U(1, 1);
  m11_U << 1;
  matrix_d m11_V(1, 1);
  m11_V << 1;

  matrix_d m11_U_out, m11_V_out;
  std::tie(m11_U_out, std::ignore, m11_V_out) = svd(m11);
  EXPECT_MATRIX_FLOAT_EQ(m11_U, m11_U_out);
  EXPECT_MATRIX_FLOAT_EQ(m11_V, m11_V_out);

  matrix_d m22(2, 2);
  m22 << 1, 9, -4, 2;
  matrix_d m22_U(2, 2);
  m22_U << 0.97759158130035961, 0.21051057021124275, 0.21051057021124264,
      -0.97759158130035972;
  matrix_d m22_V(2, 2);
  m22_V << 0.014701114509569043, 0.999891932776825976, 0.999891932776825976,
      -0.014701114509569043;
  vector_d m22_D(2);
  m22_D << 9.220341788859560, 4.121322275266775;

  matrix_d m22_U_out, m22_V_out;
  vector_d m22_D_out;
  std::tie(m22_U_out, m22_D_out, m22_V_out) = svd(m22);

  EXPECT_MATRIX_FLOAT_EQ(m22_U, m22_U_out);
  EXPECT_MATRIX_FLOAT_EQ(m22_V, m22_V_out);
  EXPECT_MATRIX_FLOAT_EQ(m22_D, m22_D_out);

  matrix_d m23(2, 3);
  m23 << 1, 3, -5, 7, 9, -11;
  matrix_d m23_U(2, 2);
  m23_U << -0.33784321032557019, -0.94120240396894062, -0.94120240396894039,
      0.33784321032557019;
  matrix_d m23_V(3, 2);
  m23_V << -0.41176240532160857, 0.81473005032163681, -0.56383954240865775,
      0.12417046246885260, 0.71591667949570703, 0.56638912538393127;
  vector_d m23_D(2);
  m23_D << 16.821011215675149, 1.747450051398016;

  matrix_d m23_U_out, m23_V_out;
  vector_d m23_D_out;
  std::tie(m23_U_out, m23_D_out, m23_V_out) = svd(m23);

  EXPECT_MATRIX_FLOAT_EQ(m23_U, m23_U_out);
  EXPECT_MATRIX_FLOAT_EQ(m23_V, m23_V_out);
  EXPECT_MATRIX_FLOAT_EQ(m23_D, m23_D_out);

  matrix_d m32(3, 2);
  m32 << 1, 3, -5, 7, 9, -11;
  matrix_d m32_U(3, 2);
  m32_U << 0.10657276921942949, 0.978015199679778457, 0.51489182605641137,
      0.099934023569151917, -0.85060487438128174, 0.183028577355020677;
  matrix_d m32_V(2, 2);
  m32_V << -0.60622380392317887, 0.79529409626685355, 0.79529409626685355,
      0.60622380392317887;
  vector_d m32_D(2);
  m32_D << 16.698998232964481, 2.672724829728879;

  matrix_d m32_U_out, m32_V_out;
  vector_d m32_D_out;
  std::tie(m32_U_out, m32_D_out, m32_V_out) = svd(m32);
  // R's SVD returns different signs than Eigen.
  EXPECT_MATRIX_FLOAT_EQ(m32_U, m32_U_out);
  EXPECT_MATRIX_FLOAT_EQ(m32_V, m32_V_out);
  EXPECT_MATRIX_FLOAT_EQ(m32_D, m32_D_out);

  matrix_c c32(3, 2);
  c32 << compl_t(0.86636546, 0.34306449), compl_t(0.28267243, 0.52462912),
      compl_t(0.12104914, 0.2533793), compl_t(0.66889264, 0.39276455),
      compl_t(0.02184348, 0.0614428), compl_t(0.96599692, 0.16180684);
  matrix_c c32_U(3, 2);
  c32_U << compl_t(0.50789057, 0.35782384), compl_t(0.74507868, 0.14261495),
      compl_t(0.4823205, 0.19300139), compl_t(-0.29600299, 0.17466116),
      compl_t(0.58489862, -0.04494745), compl_t(-0.53765935, 0.13159357);
  matrix_c c32_V(2, 2);
  c32_V << compl_t(0.45315061, 0.), compl_t(0.89143395, 0.),
      compl_t(0.85785925, -0.2423469), compl_t(-0.43608328, 0.12319437);
  vector_d c32_D(2);
  c32_D << 1.50077492, 0.78435681;

  matrix_c c32_U_out, c32_V_out;
  vector_d c32_D_out;
  std::tie(c32_U_out, c32_D_out, c32_V_out) = svd(c32);
  EXPECT_MATRIX_COMPLEX_FLOAT_EQ(c32_U, c32_U_out);
}
