#ifdef STAN_OPENCL
#include <stan/math/opencl/prim.hpp>
#include <stan/math/prim/fun/logit.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathMatrixCL, logit) {
  stan::math::matrix_d a(5, 1);
  a << 0.5, 0.4, 0.3, 0.2, 0.1;
  using stan::math::logit;
  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> res_cl = logit(a_cl);
  stan::math::matrix_d res = stan::math::from_matrix_cl(res_cl);
  EXPECT_MATRIX_NEAR(logit(a), res, 1E-08);
}

TEST(MathMatrixCL, logit_nan) {
  stan::math::matrix_d a(1, 1);
  a << std::numeric_limits<double>::quiet_NaN();
  using stan::math::logit;
  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> res_cl = logit(a_cl);
  stan::math::matrix_d res = stan::math::from_matrix_cl(res_cl);

  EXPECT_TRUE(std::isnan(logit(res(0))));
}

#endif
