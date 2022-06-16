#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <stdexcept>
#include <stan/math/prim/fun/Eigen.hpp>
#include <complex>

TEST(MathMatrixPrimMat, svd_V) {
  using stan::math::matrix_d;
  using stan::math::svd_V;
  using compl_t = std::complex<double>;
  using matrix_c = Eigen::Matrix<compl_t, Eigen::Dynamic, Eigen::Dynamic>;

  // Values generated using R base::svd

  matrix_d m00(0, 0);
  EXPECT_THROW(svd_V(m00), std::invalid_argument);

  matrix_d m11(1, 1);
  m11 << 5;
  matrix_d m11_V(1, 1);
  m11_V << 1;
  EXPECT_MATRIX_FLOAT_EQ(m11_V, svd_V(m11));

  matrix_d m22(2, 2);
  m22 << 1, 9, -4, 2;
  matrix_d m22_V(2, 2);
  m22_V << 0.014701114509569043, 0.999891932776825976, 0.999891932776825976,
      -0.014701114509569043;
  EXPECT_MATRIX_FLOAT_EQ(m22_V, svd_V(m22));

  matrix_d m23(2, 3);
  m23 << 1, 3, -5, 7, 9, -11;
  matrix_d m23_V(3, 2);
  m23_V << -0.41176240532160857, 0.81473005032163681, -0.56383954240865775,
      0.12417046246885260, 0.71591667949570703, 0.56638912538393127;
  EXPECT_MATRIX_FLOAT_EQ(m23_V, svd_V(m23));

  matrix_d m32(3, 2);
  m32 << 1, 3, -5, 7, 9, -11;
  matrix_d m32_V(2, 2);
  m32_V << -0.60622380392317887, 0.79529409626685355, 0.79529409626685355,
      0.60622380392317887;
  EXPECT_MATRIX_FLOAT_EQ(
      m32_V, svd_V(m32));  // R's SVD returns different signs than Eigen.

  matrix_c c32(3, 2);
  c32 << compl_t(0.86636546, 0.34306449), compl_t(0.28267243, 0.52462912),
      compl_t(0.12104914, 0.2533793), compl_t(0.66889264, 0.39276455),
      compl_t(0.02184348, 0.0614428), compl_t(0.96599692, 0.16180684);
  matrix_c c32_V(2, 2);
  c32_V << compl_t(0.45315061, 0.), compl_t(0.89143395, 0.),
      compl_t(0.85785925, -0.2423469), compl_t(-0.43608328, 0.12319437);

  EXPECT_MATRIX_COMPLEX_FLOAT_EQ(c32_V, svd_V(c32));
}
