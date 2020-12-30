#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixCL, segment_exception) {
  stan::math::matrix_cl<double> m1_cl(3, 3);
  EXPECT_THROW(segment(m1_cl, 1, 5), std::invalid_argument);

  stan::math::matrix_cl<double> m2_cl(3, 1);
  EXPECT_THROW(segment(m2_cl, 0, 2), std::invalid_argument);
  EXPECT_THROW(segment(m2_cl, 1, 5), std::invalid_argument);

  stan::math::matrix_cl<double> m3_cl(1, 3);
  EXPECT_THROW(segment(m3_cl, 0, 2), std::invalid_argument);
  EXPECT_THROW(segment(m3_cl, 1, 5), std::invalid_argument);

  stan::math::var_value<stan::math::matrix_cl<double>> m4_cl = m1_cl;
  EXPECT_THROW(segment(m1_cl, 1, 5), std::invalid_argument);

  stan::math::var_value<stan::math::matrix_cl<double>> m5_cl = m2_cl;
  EXPECT_THROW(segment(m2_cl, 0, 2), std::invalid_argument);
  EXPECT_THROW(segment(m2_cl, 1, 5), std::invalid_argument);

  stan::math::var_value<stan::math::matrix_cl<double>> m6_cl = m3_cl;
  EXPECT_THROW(segment(m3_cl, 0, 2), std::invalid_argument);
  EXPECT_THROW(segment(m3_cl, 1, 5), std::invalid_argument);
}

auto segment_functor = [](const auto& a, size_t i, size_t n) {
  return stan::math::segment(a, i, n);
};

TEST(MathMatrixCL, segment_value_check) {
  stan::math::vector_d m1(16);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

  stan::math::test::compare_cpu_opencl_prim_rev(segment_functor, m1, 5, 10);
}

#endif
