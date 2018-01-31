#ifdef STAN_OPENCL
#include <stan/math/gpu/opencl_context.hpp>
#include <gtest/gtest.h>

TEST(MathGpu, setup) {
  EXPECT_EQ(0, stan::math::initialized_) << "starts uninitialized";
  EXPECT_EQ(0, stan::math::kernel_groups.size());
  EXPECT_EQ(0, stan::math::kernel_strings.size());
  EXPECT_EQ(0, stan::math::kernels.size());
  EXPECT_EQ(0, stan::math::compiled_kernels.size());


  cl::Context cl = stan::math::get_context();
  EXPECT_EQ(1, stan::math::initialized_) << "initializes after get_context()";
  EXPECT_NE("", stan::math::get_description());
  EXPECT_EQ(256, stan::math::get_maximum_workgroup_size());
  EXPECT_EQ(1, stan::math::initialized_) << "starts uninitialized";

  EXPECT_EQ(21, stan::math::kernel_groups.size());
  EXPECT_EQ(6, stan::math::kernel_strings.size());
  EXPECT_EQ(1, stan::math::kernels.size());
  EXPECT_EQ(6, stan::math::compiled_kernels.size());
}
#endif
