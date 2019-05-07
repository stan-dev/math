#ifdef STAN_OPENCL

#include <stan/math/prim/mat.hpp>
#include <stan/math/opencl/opencl.hpp>
#include <gtest/gtest.h>

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

// Test we can handle a semi-complicated but normal queue
TEST(async_opencl, thrash_opencl) {
  auto size = 100;
  stan::math::matrix_d m1 = stan::math::matrix_d::Random(size, size);
  stan::math::matrix_d m1_result = m1 + m1 * m1 - m1;
  stan::math::matrix_cl m1_cl(size, size);
  stan::math::copy(m1_cl, m1);  // NOLINT
  stan::math::matrix_cl m1_result_cl = m1_cl + m1_cl * m1_cl - m1_cl;
  stan::math::matrix_d m1_result_test(size, size);
  stan::math::copy(m1_result_test, m1_result_cl);  // NOLINT
  EXPECT_MATRIX_NEAR(m1_result, m1_result_test, 1e-12)
}

// Test we handle writes to self
TEST(async_opencl, assign_miss) {
  auto size = 100;
  stan::math::matrix_d m1 = stan::math::matrix_d::Random(size, size);
  stan::math::matrix_d m1_result = m1 + m1 * m1 - m1;
  m1_result = m1_result + m1_result * m1_result - m1_result;
  stan::math::matrix_cl m1_cl(size, size);
  stan::math::copy(m1_cl, m1);  // NOLINT
  stan::math::matrix_cl m1_result_cl = m1_cl + m1_cl * m1_cl - m1_cl;
  m1_result_cl = m1_result_cl + m1_result_cl * m1_result_cl - m1_result_cl;
  stan::math::matrix_d m1_result_test(size, size);
  stan::math::copy(m1_result_test, m1_result_cl);  // NOLINT
  EXPECT_MATRIX_NEAR(m1_result, m1_result_test, 1e-12)
}

// test we handle reads correctly
TEST(async_opencl, read_miss) {
  auto size = 100;
  stan::math::matrix_d m1 = stan::math::matrix_d::Random(size, size);
  stan::math::matrix_d m1_result = m1 + m1 * m1 - m1;
  stan::math::matrix_cl m1_cl(size, size);
  stan::math::copy(m1_cl, m1);  // NOLINT
  stan::math::matrix_cl m1_result_cl = m1_cl + m1_cl * m1_cl - m1_cl;
  m1_cl = m1_cl * 2;
  stan::math::matrix_d m1_result_test(size, size);
  stan::math::copy(m1_result_test, m1_result_cl);  // NOLINT
  EXPECT_MATRIX_NEAR(m1_result, m1_result_test, 1e-12)
  m1 = m1 * 2;
  stan::math::matrix_d m1_multiply_test(size, size);
  stan::math::copy(m1_multiply_test, m1_cl);  // NOLINT
  EXPECT_MATRIX_NEAR(m1, m1_multiply_test, 1e-12)
}

#endif
