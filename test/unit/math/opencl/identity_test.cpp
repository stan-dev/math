#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/identity.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixGPU, identity_m_exception_pass) {
  stan::math::matrix_cl m(1, 1);
  stan::math::matrix_cl m0;

  EXPECT_NO_THROW(m = stan::math::identity(1));
  EXPECT_NO_THROW(m0 = stan::math::identity(0));
}

TEST(MathMatrixGPU, identity_m_value_check) {
  stan::math::matrix_d m0(2, 2);
  m0 << 2, 2, 2, 2;
  stan::math::matrix_cl m(m0);

  EXPECT_NO_THROW(m = stan::math::identity(2));

  stan::math::copy(m0, m);
  EXPECT_EQ(1, m0(0, 0));
  EXPECT_EQ(0, m0(0, 1));
  EXPECT_EQ(0, m0(1, 0));
  EXPECT_EQ(1, m0(1, 1));
}

#endif
