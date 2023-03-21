#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(primFun, rfft) {
  using v_t = Eigen::VectorXd;
  using cv_t = Eigen::Matrix<std::complex<double>, -1, 1>;
  using stan::math::rfft;

  // reference answers calculated using Scipy.fft with double precision

  v_t x0(0);
  cv_t y0 = rfft(x0);
  EXPECT_EQ(0, y0.size());

  v_t x1(1);
  x1 << -3.247;
  cv_t y1 = rfft(x1);
  EXPECT_EQ(1, y1.size());
  EXPECT_EQ(real(y1[0]), -3.247);
  EXPECT_EQ(imag(y1[0]), 0);

  cv_t y1b = rfft(x1 + x1);
  EXPECT_NEAR(real(y1b[0]), -3.247 * 2, 1e-6);
  EXPECT_NEAR(imag(y1b[0]), 0, 1e-6);

  v_t x(3);
  x << 1, -3, 7;
  cv_t y = rfft(x);
  EXPECT_EQ(3, y.size());
  EXPECT_NEAR(real(y[0]), 5, 1e-6);
  EXPECT_NEAR(imag(y[0]), 0, 1e-6);
  EXPECT_NEAR(real(y[1]), -1, 1e-6);
  EXPECT_NEAR(imag(y[1]), 8.6602540, 1e-6);
  EXPECT_NEAR(real(y[2]), -1, 1e-6);
  EXPECT_NEAR(imag(y[2]), -8.6602540, 1e-6);

  cv_t yb = rfft(x + x);
  EXPECT_EQ(3, yb.size());
  EXPECT_NEAR(real(yb[0]), 2 * 5, 1e-6);
  EXPECT_NEAR(imag(yb[0]), 0, 1e-6);
  EXPECT_NEAR(real(yb[1]), 2 * -1, 1e-6);
  EXPECT_NEAR(imag(yb[1]), 2 * 8.6602540, 1e-6);
  EXPECT_NEAR(real(yb[2]), 2 * -1, 1e-6);
  EXPECT_NEAR(imag(yb[2]), 2 * -8.6602540, 1e-6);
}

TEST(primFun, inv_rfft) {
  using v_t = Eigen::VectorXd;
  using c_t = std::complex<double>;
  using cv_t = Eigen::Matrix<std::complex<double>, -1, 1>;
  using stan::math::inv_rfft;

  // reference answers calculated using Scipy.fft with double precision

  cv_t y0(0);
  v_t x0 = inv_rfft(y0);
  EXPECT_EQ(0, x0.size());

  cv_t y1(1);
  y1 << c_t(-3.247, 0);
  v_t x1 = inv_rfft(y1);

  EXPECT_EQ(1, x1.size());
  EXPECT_EQ(x1[0], -3.247);

  cv_t y(3);
  y << c_t(5, 0), c_t(-1, 8.6602540), c_t(-1, -8.6602540);
  v_t x = inv_rfft(y);
  EXPECT_EQ(3, x.size());
  v_t x_expected(3);
  x_expected << 1, -3, 7;
  EXPECT_MATRIX_NEAR(x_expected, x, 1e-7);
}

TEST(primFun, rfft2) {
  using c_t = std::complex<double>;
  using cm_t = Eigen::Matrix<std::complex<double>, -1, -1>;
  using m_t = Eigen::MatrixXd;
  using stan::math::rfft2;

  // reference answers calculated using Scipy.fft with double
  // precision, and verified with R

  m_t x(0, 0);
  cm_t y = rfft2(x);
  EXPECT_EQ(0, y.rows());
  EXPECT_EQ(0, y.cols());

  m_t x11(1, 1);
  x11 << -3.9;
  cm_t y11 = rfft2(x11);
  EXPECT_NEAR(real(y11(0, 0)), -3.9, 1e-6);
  EXPECT_NEAR(imag(y11(0, 0)), 0, 1e-6);

  m_t x12(1, 2);
  x12 << -3.9, 0.2;
  cm_t y12 = rfft2(x12);
  cm_t y12_expected(1, 2);
  y12_expected << c_t(-3.7, 0), c_t(-4.1, 0);
  EXPECT_MATRIX_COMPLEX_NEAR(y12_expected, y12, 1e-8);

  m_t x33(3, 3);
  x33 << 2, -1.4, 1, -9, 2, 3.9, 13, 1.3, -2.2;
  cm_t y33 = rfft2(x33);
  cm_t y33_expected(3, 3);
  y33_expected << c_t(10.6, -0.), c_t(3.7, +0.69282032), c_t(3.7, -0.69282032),
      c_t(-2.9, +13.16358614), c_t(5.5, +24.76832655), c_t(-2.6, +19.22576396),
      c_t(-2.9, -13.16358614), c_t(-2.6, -19.22576396), c_t(5.5, -24.76832655);
  EXPECT_MATRIX_COMPLEX_NEAR(y33_expected, y33, 1e-8);
}

TEST(primFunFFT, invfft2) {
  using c_t = std::complex<double>;
  using cm_t = Eigen::Matrix<std::complex<double>, -1, -1>;
  using m_t = Eigen::MatrixXd;
  using stan::math::inv_rfft2;
  using stan::math::rfft2;

  // reference answers calculated using Scipy.fft with double
  // precision, and verified with R

  cm_t x(0, 0);
  m_t y = inv_rfft2(x);
  EXPECT_EQ(0, y.rows());
  EXPECT_EQ(0, y.cols());

  cm_t x11(1, 1);
  x11 << c_t(-3.9,0);
  m_t y11 = inv_rfft2(x11);
  EXPECT_NEAR(y11(0, 0), -3.9, 1e-6);

  cm_t x12(1, 2);
  x12 << c_t(-3.7, 0), c_t(-4.1, 0);
  m_t y12 = inv_rfft2(x12);
  m_t y12_expected(1, 2);
  y12_expected << -3.9, 0.2;
  EXPECT_MATRIX_NEAR(y12_expected, y12, 1e-8);

  cm_t x33(3, 3);
  x33 << c_t(10.6, -0.), c_t(3.7, +0.69282032), c_t(3.7, -0.69282032),
      c_t(-2.9, +13.16358614), c_t(5.5, +24.76832655), c_t(-2.6, +19.22576396),
      c_t(-2.9, -13.16358614), c_t(-2.6, -19.22576396), c_t(5.5, -24.76832655);
  m_t y33 = inv_rfft2(x33);
  m_t y33_expected(3, 3);
  y33_expected << 2, -1.4, 1, -9, 2, 3.9, 13, 1.3, -2.2;
  EXPECT_MATRIX_NEAR(y33_expected, y33, 1e-8);
}
