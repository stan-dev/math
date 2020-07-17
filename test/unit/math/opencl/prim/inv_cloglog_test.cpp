#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/inv_cloglog.hpp>
#include <stan/math/prim/fun/inv_cloglog.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathMatrixCL, inv_cloglog) {
  stan::math::matrix_d y(1, 4);
  y << 3.7, 0.0, -2.93, std::numeric_limits<double>::quiet_NaN();

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> z_cl = stan::math::inv_cloglog(y_cl);
  stan::math::matrix_d z = stan::math::from_matrix_cl(z_cl);
  EXPECT_FLOAT_EQ(1 - std::exp(-std::exp(3.7)), z(0, 0));
  EXPECT_FLOAT_EQ(1 - std::exp(-std::exp(0.0)), z(0, 1));
  EXPECT_FLOAT_EQ(1 - std::exp(-std::exp(-2.93)), z(0, 2));
  EXPECT_TRUE(std::isnan(stan::math::inv_cloglog(z(0, 3))));

  stan::math::vector_d av(9);
  av << 1, 2, 3, 4, 22, 0.5, 5, 18, 99;
  stan::math::matrix_cl<double> av_cl(av);
  stan::math::matrix_cl<double> bv_cl = stan::math::inv_cloglog(av_cl);
  stan::math::vector_d bv = stan::math::from_matrix_cl(bv_cl);
  EXPECT_MATRIX_NEAR(bv, stan::math::inv_cloglog(av), 1E-08);

  stan::math::row_vector_d arv(9);
  arv << 1, 2, 3, 4, 22, 0.5, 5, 18, 99;
  stan::math::matrix_cl<double> arv_cl(arv);
  stan::math::matrix_cl<double> brv_cl = stan::math::inv_cloglog(arv_cl);
  stan::math::row_vector_d brv = stan::math::from_matrix_cl(brv_cl);
  EXPECT_MATRIX_NEAR(brv, stan::math::inv_cloglog(arv), 1E-08);
}
#endif
