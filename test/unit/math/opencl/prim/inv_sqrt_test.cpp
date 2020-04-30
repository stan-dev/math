#include <stan/math/opencl/prim/inv_sqrt.hpp>
#include <stan/math/prim/fun/inv_sqrt.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

TEST(MathMatrixCL, inv_sqrt) {
  double y = 4.0;
  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> z_cl = stan::math::inv_sqrt(y_cl);
  stan::math::matrix_d z = stan::math::from_matrix_cl(z_cl);
  EXPECT_FLOAT_EQ(1 / 2.0, z(0,0));

  y_cl = stan::math::to_matrix_cl(25.0);
  z_cl = stan::math::inv_sqrt(y_cl);
  z = stan::math::from_matrix_cl(z_cl);
  EXPECT_FLOAT_EQ(1 / 5.0, z(0,0));

  y_cl = stan::math::to_matrix_cl(0.0);
  z_cl = stan::math::inv_sqrt(y_cl);
  z = stan::math::from_matrix_cl(z_cl);
  EXPECT_FLOAT_EQ(stan::math::positive_infinity(), z(0,0));


  y_cl = stan::math::to_matrix_cl(-50.0);
  z_cl = stan::math::inv_sqrt(y_cl);
  z = stan::math::from_matrix_cl(z_cl);
  EXPECT_TRUE(std::isnan(z(0,0)));

  y_cl = stan::math::to_matrix_cl(std::numeric_limits<double>::quiet_NaN());
  z_cl = stan::math::inv_sqrt(y_cl);
  z = stan::math::from_matrix_cl(z_cl);
  EXPECT_TRUE(std::isnan(z(0,0)));

  stan::math::matrix_d a(3,3);
  a << 1, 2, 3, 4, 22, 0.5, 5, 18, 99;

  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> b_cl = stan::math::inv_sqrt(a_cl);
  stan::math::matrix_d b = stan::math::from_matrix_cl(b_cl);
  EXPECT_MATRIX_NEAR(b, stan::math::inv_sqrt(a), 1E-08);
}
