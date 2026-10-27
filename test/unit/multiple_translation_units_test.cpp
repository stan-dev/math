#include <gtest/gtest.h>

TEST(multiple_translation_units, compile) {
  SUCCEED() << "this test compiling indicates that compiling the math library "
            << "public umbrella headers with multiple translation units is ok.";
}

#ifdef STAN_OPENCL
double opencl_cached_kernel1();
double opencl_cached_kernel2();
void select_opencl_device2();

TEST(multiple_translation_units, opencl_cached_kernels_after_device_switch) {
  EXPECT_DOUBLE_EQ(8.0, opencl_cached_kernel1());
  EXPECT_DOUBLE_EQ(16.0, opencl_cached_kernel2());

  // Switching in another translation unit must invalidate both kernel caches.
  select_opencl_device2();

  EXPECT_DOUBLE_EQ(8.0, opencl_cached_kernel1());
  EXPECT_DOUBLE_EQ(16.0, opencl_cached_kernel2());
}
#endif
