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
