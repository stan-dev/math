#ifdef STAN_OPENCL
#include <stan/math/gpu/opencl_context.hpp>
#include <gtest/gtest.h>

TEST(MathGpu, setup) {
  cl::Context cl = stan::math::opencl_context.context();
  EXPECT_NE("", stan::math::opencl_context.description());
  EXPECT_EQ(256, stan::math::opencl_context.maxWorkgroupSize());

  EXPECT_EQ(21, stan::math::opencl_context.kernel_groups.size());
  EXPECT_EQ(6, stan::math::opencl_context.kernel_strings.size());
  EXPECT_EQ(1, stan::math::opencl_context.kernels.size());
  EXPECT_EQ(6, stan::math::opencl_context.compiled_kernels.size());
}
#endif
