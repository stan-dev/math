#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(primFun, fft) {
  using c_t = std::complex<double>;
  using cv_t = Eigen::Matrix<std::complex<double>, -1, 1>;
  using stan::math::fft;

  // reference answers calculated using Scipy.fft with double precision

  cv_t x0(0);
  cv_t y0 = fft(x0);
  EXPECT_EQ(0, y0.size());

  cv_t x1(1);
  x1 << c_t(-3.247, 1.98555);
  cv_t y1 = fft(x1);
  EXPECT_EQ(1, y1.size());
  EXPECT_EQ(real(y1[0]), -3.247);
  EXPECT_EQ(imag(y1[0]), 1.98555);

  cv_t y1b = fft(x1 + x1);
  EXPECT_NEAR(real(y1b[0]), -3.247 * 2, 1e-6);
  EXPECT_NEAR(imag(y1b[0]), 1.98555 * 2, 1e-6);

  cv_t x(3);
  x << c_t(1, -2), c_t(-3, 5), c_t(-7, 11);
  cv_t y = fft(x);
  EXPECT_EQ(3, y.size());
  EXPECT_NEAR(real(y[0]), -9, 1e-6);
  EXPECT_NEAR(imag(y[0]), 14, 1e-6);
  EXPECT_NEAR(real(y[1]), 0.80384758, 1e-6);
  EXPECT_NEAR(imag(y[1]), -13.46410162, 1e-6);
  EXPECT_NEAR(real(y[2]), 11.19615242, 1e-6);
  EXPECT_NEAR(imag(y[2]), -6.53589838, 1e-6);

  cv_t yb = fft(x + x);
  EXPECT_EQ(3, yb.size());
  EXPECT_NEAR(real(yb[0]), 2 * -9, 1e-6);
  EXPECT_NEAR(imag(yb[0]), 2 * 14, 1e-6);
  EXPECT_NEAR(real(yb[1]), 2 * 0.80384758, 1e-6);
  EXPECT_NEAR(imag(yb[1]), 2 * -13.46410162, 1e-6);
  EXPECT_NEAR(real(yb[2]), 2 * 11.19615242, 1e-6);
  EXPECT_NEAR(imag(yb[2]), 2 * -6.53589838, 1e-6);
}

TEST(primFun, inv_fft) {
  using c_t = std::complex<double>;
  using cv_t = Eigen::Matrix<std::complex<double>, -1, 1>;
  using stan::math::inv_fft;

  // reference answers calculated using Scipy.fft with double precision

  cv_t y0(0);
  cv_t x0 = inv_fft(y0);
  cv_t x0_expected(0);
  EXPECT_EQ(x0_expected, x0);

  cv_t y1(1);
  y1 << c_t(-3.247, 1.98555);
  cv_t x1 = inv_fft(y1);
  cv_t x1_expected(1);
  x1_expected << c_t(-3.247, 1.98555);
  EXPECT_MATRIX_COMPLEX_NEAR(x1_expected, x1, 1e-8);

  EXPECT_EQ(1, x1.size());
  EXPECT_EQ(real(x1[0]), -3.247);
  EXPECT_EQ(imag(x1[0]), 1.98555);

  cv_t y(3);
  y << c_t(-9, 14), c_t(0.80384758, -13.46410162),
      c_t(11.19615242, -6.53589838);
  cv_t x = inv_fft(y);
  EXPECT_EQ(3, y.size());
  Eigen::VectorXcd x_expected(3);
  x_expected << c_t(1, -2), c_t(-3, 5), c_t(-7, 11);
  EXPECT_MATRIX_COMPLEX_NEAR(x_expected, x, 1e-8);
}

TEST(primFun, fft2) {
  using c_t = std::complex<double>;
  using cm_t = Eigen::Matrix<std::complex<double>, -1, -1>;
  using stan::math::fft2;
  using stan::math::inv_fft2;

  // reference answers calculated using Scipy.fft with double
  // precision, and verified with R

  cm_t x(0, 0);
  cm_t y = fft2(x);
  EXPECT_EQ(0, y.rows());
  EXPECT_EQ(0, y.cols());

  cm_t x11(1, 1);
  x11 << c_t(1.0, -3.9);
  cm_t y11 = fft2(x11);
  EXPECT_NEAR(1.0, std::real(y11(0, 0)), 1e-6);
  EXPECT_NEAR(-3.9, std::imag(y11(0, 0)), 1e-6);

  cm_t x12(1, 2);
  x12 << c_t(1.0, -3.9), c_t(-8.6, 0.2);
  cm_t y12 = fft2(x12);
  cm_t y12_expected(1, 2);
  y12_expected << c_t(-7.6, -3.7), c_t(9.6, -4.1);
  EXPECT_MATRIX_COMPLEX_NEAR(y12_expected, y12, 1e-8);

  cm_t x33(3, 3);
  x33 << c_t(1, 2), c_t(3, -1.4), c_t(2, 1), c_t(3, -9), c_t(2, -1.3),
      c_t(3.9, -1.8), c_t(13, -1.8), c_t(1.3, 1.9), c_t(-2.2, -2.2);
  cm_t y33 = fft2(x33);
  cm_t y33_expected(3, 3);
  y33_expected << c_t(27, -12.6), c_t(13.90525589, -9.15166605),
      c_t(10.09474411, -4.64833395), c_t(-13.160254038, 11.471281292),
      c_t(-13.29326674, 20.88153533), c_t(-13.25262794, 15.82794549),
      c_t(4.160254038, 5.928718708), c_t(-11.34737206, -7.72794549),
      c_t(4.89326674, -1.98153533);
  EXPECT_MATRIX_COMPLEX_NEAR(y33_expected, y33, 1e-8);
}

TEST(primFunFFT, invfft2) {
  using c_t = std::complex<double>;
  using cm_t = Eigen::Matrix<std::complex<double>, -1, -1>;
  using stan::math::fft2;
  using stan::math::inv_fft;
  using stan::math::inv_fft2;

  // reference answers calculated using R
  cm_t x(0, 0);
  cm_t y = inv_fft2(x);
  EXPECT_EQ(0, y.rows());
  EXPECT_EQ(0, y.cols());

  cm_t x11(1, 1);
  x11 << c_t(1.0, -3.9);
  cm_t y11 = inv_fft2(x11);
  EXPECT_NEAR(1.0, std::real(y11(0, 0)), 1e-6);
  EXPECT_NEAR(-3.9, std::imag(y11(0, 0)), 1e-6);

  cm_t x13(1, 3);
  x13 << c_t(-2.3, 1.82), c_t(1.18, 9.32), c_t(1.15, -14.1);
  cm_t y13 = inv_fft2(x13);
  cm_t y13copy = inv_fft(x13.row(0));
  EXPECT_MATRIX_COMPLEX_NEAR(y13, y13copy.transpose(), 1e-8);

  cm_t x33(3, 3);
  x33 << c_t(1, 2), c_t(3, -1.4), c_t(2, 1), c_t(3, -9), c_t(2, -1.3),
      c_t(3.9, -1.8), c_t(13, -1.8), c_t(1.3, 1.9), c_t(-2.2, -2.2);
  cm_t y33 = fft2(x33);

  // check versus results from R
  EXPECT_NEAR(27.000000, real(y33(0, 0)), 1e-5);
  EXPECT_NEAR(-12.600000, imag(y33(0, 0)), 1e-5);
  EXPECT_NEAR(13.90526, real(y33(0, 1)), 1e-5);
  EXPECT_NEAR(-9.15167, imag(y33(0, 1)), 1e-5);
  EXPECT_NEAR(10.094744, real(y33(0, 2)), 1e-5);
  EXPECT_NEAR(-4.648334, imag(y33(0, 2)), 1e-5);
  EXPECT_NEAR(-13.160254, real(y33(1, 0)), 1e-5);
  EXPECT_NEAR(11.471281, imag(y33(1, 0)), 1e-5);
  EXPECT_NEAR(-13.29327, real(y33(1, 1)), 1e-5);
  EXPECT_NEAR(20.88154, imag(y33(1, 1)), 1e-5);
  EXPECT_NEAR(-13.252628, real(y33(1, 2)), 1e-5);
  EXPECT_NEAR(15.827945, imag(y33(1, 2)), 1e-5);
  EXPECT_NEAR(4.160254, real(y33(2, 0)), 1e-5);
  EXPECT_NEAR(5.928719, imag(y33(2, 0)), 1e-5);
  EXPECT_NEAR(-11.34737, real(y33(2, 1)), 1e-5);
  EXPECT_NEAR(-7.72795, imag(y33(2, 1)), 1e-5);
  EXPECT_NEAR(4.893267, real(y33(2, 2)), 1e-5);
  EXPECT_NEAR(-1.981535, imag(y33(2, 2)), 1e-5);

  // check round trips inv_fft(fft(x))
  cm_t x33copy = inv_fft2(y33);
  EXPECT_MATRIX_COMPLEX_NEAR(x33, x33copy, 1e-8);

  // check round trip fft(inv_fft(x))
  cm_t z33 = inv_fft2(x33);
  cm_t x33copy2 = fft2(z33);
  EXPECT_MATRIX_COMPLEX_NEAR(x33, x33copy2, 1e-8);
}
