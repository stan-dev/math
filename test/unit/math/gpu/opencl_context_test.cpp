#ifdef STAN_OPENCL
#include <stan/math/gpu/opencl_context.hpp>
#include <gtest/gtest.h>

TEST(MathGpu, setup) {
  EXPECT_NE("", stan::math::opencl_context.get_description());
  EXPECT_EQ(1024, stan::math::opencl_context.get_maximum_workgroup_size());

  EXPECT_EQ(21, stan::math::kernel_groups.size());
  EXPECT_EQ(6, stan::math::kernel_strings.size());
  EXPECT_EQ(1, stan::math::kernels.size());
  EXPECT_EQ(6, stan::math::compiled_kernels.size());
}
#endif
