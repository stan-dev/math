#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixCL, col_exception) {
  stan::math::matrix_cl<double> m1_cl(3, 3);
  EXPECT_THROW(col(m1_cl, 0), std::invalid_argument);
  EXPECT_THROW(col(m1_cl, 4), std::invalid_argument);

  stan::math::var_value<stan::math::matrix_cl<double>> m2_cl = m1_cl;
  EXPECT_THROW(col(m2_cl, 0), std::invalid_argument);
  EXPECT_THROW(col(m2_cl, 4), std::invalid_argument);
}

auto col_functor = [](const auto& a) { return stan::math::col(a, 2); };

TEST(MathMatrixCL, col_value_check) {
  stan::math::matrix_d m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;

  stan::math::test::compare_cpu_opencl_prim_rev(col_functor, m1);
}

#endif
