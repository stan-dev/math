#ifdef STAN_OPENCL

#include <stan/math/prim.hpp>
#include <stan/math/opencl/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

// Test we can handle a semi-complicated but normal queue
TEST(async_opencl, thrash_opencl) {
  auto size = 100;
  stan::math::matrix_d m1 = stan::math::matrix_d::Random(size, size);
  stan::math::matrix_d m1_result = m1 + m1 * m1 - m1;
  stan::math::matrix_cl<double> m1_cl(size, size);
  m1_cl = stan::math::to_matrix_cl(m1);
  stan::math::matrix_cl<double> m1_result_cl = m1_cl + m1_cl * m1_cl - m1_cl;
  stan::math::matrix_d m1_result_test(size, size);
  m1_result_test = stan::math::from_matrix_cl(m1_result_cl);
  EXPECT_MATRIX_NEAR(m1_result, m1_result_test, 1e-12)
}

// Test we handle writes to self
TEST(async_opencl, assign_miss) {
  auto size = 100;
  stan::math::matrix_d m1 = stan::math::matrix_d::Random(size, size);
  stan::math::matrix_d m1_result = m1 + m1 * m1 - m1;
  m1_result = m1_result + m1_result * m1_result - m1_result;
  stan::math::matrix_cl<double> m1_cl(size, size);
  m1_cl = stan::math::to_matrix_cl(m1);
  stan::math::matrix_cl<double> m1_result_cl = m1_cl + m1_cl * m1_cl - m1_cl;
  m1_result_cl = m1_result_cl + m1_result_cl * m1_result_cl - m1_result_cl;
  stan::math::matrix_d m1_result_test(size, size);
  m1_result_test = stan::math::from_matrix_cl(m1_result_cl);
  EXPECT_MATRIX_NEAR(m1_result, m1_result_test, 1e-12)
}

// test we handle reads correctly
TEST(async_opencl, read_miss) {
  auto size = 100;
  stan::math::matrix_d m1 = stan::math::matrix_d::Random(size, size);
  stan::math::matrix_d m1_result = m1 + m1 * m1 - m1;
  stan::math::matrix_cl<double> m1_cl(size, size);
  m1_cl = stan::math::to_matrix_cl(m1);
  stan::math::matrix_cl<double> m1_result_cl = m1_cl + m1_cl * m1_cl - m1_cl;
  m1_cl = m1_cl * 2;
  stan::math::matrix_d m1_result_test(size, size);
  m1_result_test = stan::math::from_matrix_cl(m1_result_cl);
  EXPECT_MATRIX_NEAR(m1_result, m1_result_test, 1e-12)
  m1 = m1 * 2;
  stan::math::matrix_d m1_multiply_test(size, size);
  m1_multiply_test = stan::math::from_matrix_cl(m1_cl);
  EXPECT_MATRIX_NEAR(m1, m1_multiply_test, 1e-12)
}

#endif
