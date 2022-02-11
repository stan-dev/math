#ifdef STAN_OPENCL
#include <stan/math/opencl/prim.hpp>
#include <stan/math/prim/fun/log1m_inv_logit.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathMatrixCL, log1m_inv_logit) {
  stan::math::matrix_d a(5, 1);
  a << -7.2, 0.0, 1.9, -5.0, 0.1;
  using stan::math::log1m_inv_logit;
  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> res_cl = log1m_inv_logit(a_cl);
  stan::math::matrix_d res = stan::math::from_matrix_cl(res_cl);
  EXPECT_MATRIX_NEAR(log1m_inv_logit(a), res, 1E-08);
}

TEST(MathMatrixCL, log1m_inv_logit_nan) {
  stan::math::matrix_d a(1, 1);
  a << std::numeric_limits<double>::quiet_NaN();
  using stan::math::log1m_inv_logit;
  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> res_cl = log1m_inv_logit(a_cl);
  stan::math::matrix_d res = stan::math::from_matrix_cl(res_cl);

  EXPECT_TRUE(std::isnan(log1m_inv_logit(res(0))));
}
#endif
