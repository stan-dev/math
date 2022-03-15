#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(primFun, fft) {
  typedef std::complex<double> c_t;
  typedef Eigen::Matrix<std::complex<double>, -1, 1> cv_t;
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
  x <<  c_t(1, -2), c_t(-3, 5), c_t(-7, 11);
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
  typedef std::complex<double> c_t;
  typedef Eigen::Matrix<std::complex<double>, -1, 1> cv_t;
  using stan::math::inv_fft;

  // reference answers calculated using Scipy.fft with double precision

  cv_t y0(0);
  cv_t x0 = inv_fft(y0);
  EXPECT_EQ(0, x0.size());

  cv_t y1(1);
  y1 << c_t(-3.247, 1.98555);
  cv_t x1 = inv_fft(y1);
  EXPECT_EQ(1, x1.size());
  EXPECT_EQ(real(x1[0]), -3.247);
  EXPECT_EQ(imag(x1[0]), 1.98555);
  
  cv_t y(3);
  y <<  c_t(-9, 14), c_t(0.80384758, -13.46410162), c_t(11.19615242, -6.53589838);
  cv_t x = inv_fft(y);
  EXPECT_EQ(3, y.size());
  EXPECT_NEAR(real(x[0]), 1, 1e-6);
  EXPECT_NEAR(imag(x[0]), -2, 1e-6);
  EXPECT_NEAR(real(x[1]), -3, 1e-6);
  EXPECT_NEAR(imag(x[1]), 5, 1e-6);
  EXPECT_NEAR(real(x[2]), -7, 1e-6);
  EXPECT_NEAR(imag(x[2]), 11, 1e-6);
}

TEST(primFun, fft2) {
  typedef std::complex<double> c_t;
  typedef Eigen::Matrix<std::complex<double>, -1, -1> cm_t;
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
  EXPECT_EQ(1, y12.rows());
  EXPECT_EQ(2, y12.cols());
  EXPECT_NEAR(-7.6, real(y12(0, 0)), 1e-6);
  EXPECT_NEAR(-3.7, imag(y12(0, 0)), 1e-6);
  EXPECT_NEAR(9.6, real(y12(0, 1)), 1e-6);
  EXPECT_NEAR(-4.1, imag(y12(0, 1)), 1e-6);

  cm_t x33(3, 3);
  x33 <<
    c_t(1, 2), c_t(3, -1.4), c_t(2, 1),
    c_t(3, -9), c_t(2, -1.3), c_t(3.9, -1.8),
    c_t(13, -1.8), c_t(1.3, 1.9), c_t(-2.2, -2.2);
  cm_t y33 = fft2(x33);
  EXPECT_EQ(3, y33.rows());
  EXPECT_EQ(3, y33.cols());

  // check vs. results from R (matches SciPy)
  EXPECT_NEAR(27, real(y33(0, 0)), 1e-6);
  EXPECT_NEAR(-12.6, imag(y33(0, 0)), 1e-6);
  EXPECT_NEAR(13.90525589, real(y33(0, 1)), 1e-6);
  EXPECT_NEAR(-9.15166605, imag(y33(0, 1)), 1e-6);
  EXPECT_NEAR(10.09474411, real(y33(0, 2)), 1e-1);
  EXPECT_NEAR(-4.64833395, imag(y33(0, 2)), 2e-1);
  EXPECT_NEAR(-13.160254038, real(y33(1, 0)), 1e-1);
  EXPECT_NEAR(11.471281292, imag(y33(1, 0)), 2e-1);
  EXPECT_NEAR(-13.29326674, real(y33(1, 1)), 1e-1);
  EXPECT_NEAR(20.88153533, imag(y33(1, 1)), 1e-1);
  EXPECT_NEAR(-13.25262794, real(y33(1, 2)), 1e-1);
  EXPECT_NEAR(15.82794549, imag(y33(1, 2)), 1e-1);
  EXPECT_NEAR(4.160254038, real(y33(2, 0)), 1e-1);
  EXPECT_NEAR(5.928718708, imag(y33(2, 0)), 1e-1);
  EXPECT_NEAR(-11.34737206, real(y33(2, 1)), 1e-1);
  EXPECT_NEAR(-7.72794549, imag(y33(2, 1)), 1e-1);
  EXPECT_NEAR(4.89326674, real(y33(2, 2)), 1e-1);
  EXPECT_NEAR(-1.98153533, imag(y33(2, 2)), 1e-1);
}

void expect_complex_matrix_float_eq(
        const Eigen::Matrix<std::complex<double>, -1, -1>& x,
	const Eigen::Matrix<std::complex<double>, -1, -1>& y) {
  EXPECT_EQ(x.rows(), y.rows());
  EXPECT_EQ(x.cols(), y.cols());
  for (int j = 0; j < y.cols(); ++j) {
    for (int i = 0; i < y.rows(); ++i) {
      EXPECT_FLOAT_EQ(real(x(i, j)), real(y(i, j)));
      EXPECT_FLOAT_EQ(imag(x(i, j)), imag(y(i, j)));
    }
  }
}

TEST(primFunFFT, invfft2) {
  typedef std::complex<double> c_t;
  typedef Eigen::Matrix<std::complex<double>, -1, -1> cm_t;
  using stan::math::fft2;
  using stan::math::inv_fft2;
  using stan::math::inv_fft;

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
  expect_complex_matrix_float_eq(y13, y13copy.transpose());
  
  cm_t x33(3, 3);
  x33 <<
    c_t(1, 2), c_t(3, -1.4), c_t(2, 1),
    c_t(3, -9), c_t(2, -1.3), c_t(3.9, -1.8),
    c_t(13, -1.8), c_t(1.3, 1.9), c_t(-2.2, -2.2);
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
  expect_complex_matrix_float_eq(x33, x33copy);

  // check round trip fft(inv_fft(x))
  cm_t z33 = inv_fft2(x33);
  cm_t x33copy2 = fft2(z33);
  expect_complex_matrix_float_eq(x33, x33copy2);
}
