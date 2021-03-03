#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <test/unit/util.hpp>

// this is just testing the nan behavior of the built-in fma
// there is no longer a stan::math::fma, just the agrad versions
// instead, the top-level ::fma should be used by including <cmath>

TEST(MathFunctions, fma) {
  using stan::math::fma;
  EXPECT_FLOAT_EQ(5.0, fma(1.0, 2.0, 3.0));
  EXPECT_FLOAT_EQ(10.0, fma(2.0, 3.0, 4.0));
  EXPECT_FLOAT_EQ(
      11.0, fma(static_cast<int>(3), static_cast<int>(2), static_cast<int>(5)));
}

TEST(MathFunctions, fma_nan) {
  using stan::math::fma;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(fma(1.0, 2.0, nan)));

  EXPECT_TRUE(std::isnan(fma(1.0, nan, 3.0)));

  EXPECT_TRUE(std::isnan(fma(1.0, nan, nan)));

  EXPECT_TRUE(std::isnan(fma(nan, 2.0, 3.0)));

  EXPECT_TRUE(std::isnan(fma(nan, 2.0, nan)));

  EXPECT_TRUE(std::isnan(fma(nan, nan, 3.0)));

  EXPECT_TRUE(std::isnan(fma(nan, nan, nan)));
}

TEST(MathFunctions, fma_matrix) {
  using stan::math::add;
  using stan::math::elt_multiply;
  using stan::math::fma;

  double xd = 1.0;
  Eigen::VectorXd xv(2);
  xv << 1.0, 2.0;
  Eigen::RowVectorXd xr(2);
  xr << 1.0, 2.0;
  Eigen::MatrixXd xm(2, 2);
  xm << 1.0, 2.0, -1.0, 1.1;

  double yd = 2.0;
  Eigen::VectorXd yv(2);
  yv << 2.0, -3.0;
  Eigen::RowVectorXd yr(2);
  yr << 2.0, -3.0;
  Eigen::MatrixXd ym(2, 2);
  xm << 1.0, 2.0, -1.0, 1.1;

  double zd = 3.0;
  Eigen::VectorXd zv(2);
  zv << -3.0, 4.0;
  Eigen::RowVectorXd zr(2);
  zr << -3.0, 4.0;
  Eigen::MatrixXd zm(2, 2);
  xm << 3.0, 4.0, -1.0, 1.1;

  EXPECT_MATRIX_EQ(add(elt_multiply(xd, yd), zv), fma(xd, yd, zv));
  EXPECT_MATRIX_EQ(add(elt_multiply(xd, yv), zd), fma(xd, yv, zd));
  EXPECT_MATRIX_EQ(add(elt_multiply(xd, yv), zv), fma(xd, yv, zv));
  EXPECT_MATRIX_EQ(add(elt_multiply(xv, yd), zd), fma(xv, yd, zd));
  EXPECT_MATRIX_EQ(add(elt_multiply(xv, yd), zv), fma(xv, yd, zv));
  EXPECT_MATRIX_EQ(add(elt_multiply(xv, yv), zd), fma(xv, yv, zd));
  EXPECT_MATRIX_EQ(add(elt_multiply(xv, yv), zv), fma(xv, yv, zv));

  EXPECT_MATRIX_EQ(add(elt_multiply(xd, yd), zr), fma(xd, yd, zr));
  EXPECT_MATRIX_EQ(add(elt_multiply(xd, yr), zd), fma(xd, yr, zd));
  EXPECT_MATRIX_EQ(add(elt_multiply(xd, yr), zr), fma(xd, yr, zr));
  EXPECT_MATRIX_EQ(add(elt_multiply(xr, yd), zd), fma(xr, yd, zd));
  EXPECT_MATRIX_EQ(add(elt_multiply(xr, yd), zr), fma(xr, yd, zr));
  EXPECT_MATRIX_EQ(add(elt_multiply(xr, yr), zd), fma(xr, yr, zd));
  EXPECT_MATRIX_EQ(add(elt_multiply(xr, yr), zr), fma(xr, yr, zr));

  EXPECT_MATRIX_EQ(add(elt_multiply(xd, yd), zm), fma(xd, yd, zm));
  EXPECT_MATRIX_EQ(add(elt_multiply(xd, ym), zd), fma(xd, ym, zd));
  EXPECT_MATRIX_EQ(add(elt_multiply(xd, ym), zm), fma(xd, ym, zm));
  EXPECT_MATRIX_EQ(add(elt_multiply(xm, yd), zd), fma(xm, yd, zd));
  EXPECT_MATRIX_EQ(add(elt_multiply(xm, yd), zm), fma(xm, yd, zm));
  EXPECT_MATRIX_EQ(add(elt_multiply(xm, ym), zd), fma(xm, ym, zd));
  EXPECT_MATRIX_EQ(add(elt_multiply(xm, ym), zm), fma(xm, ym, zm));
}
