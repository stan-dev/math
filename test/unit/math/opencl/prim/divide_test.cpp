#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixCL, divide_exception_pass) {
  stan::math::matrix_cl<double> m;
  EXPECT_NO_THROW(stan::math::divide(m, 2));
}

TEST(MathMatrixCL, divide_val) {
  using stan::math::divide;
  stan::math::matrix_d a(10, 1);
  a << 4, 6, 8, 10, -21, -55, -5, 0.5, 0, 1;

  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> res_cl = stan::math::divide(a_cl, 1.5);
  stan::math::matrix_d res = stan::math::from_matrix_cl(res_cl);
  stan::math::matrix_d res_cpu = stan::math::divide(a, 1.5);
  EXPECT_MATRIX_NEAR(res_cpu, res, 1E-8);
}
#endif
