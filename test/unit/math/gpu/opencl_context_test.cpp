#ifdef STAN_OPENCL
#include <stan/math/gpu/opencl_context.hpp>
#include <gtest/gtest.h>

// TEST(MathGpu, setup) {
//   cl::Context cl = stan::math::opencl_context.context();
//   EXPECT_NE("", stan::math::opencl_context.description());
//   EXPECT_EQ(256, stan::math::opencl_context.maxWorkgroupSize());

//   EXPECT_EQ(21, stan::math::opencl_context.kernel_groups.size());
//   EXPECT_EQ(6, stan::math::opencl_context.kernel_strings.size());
//   EXPECT_EQ(1, stan::math::opencl_context.kernels.size());
//   EXPECT_EQ(6, stan::math::opencl_context.compiled_kernels.size());
// }

TEST(opencl_context, construction) {
  stan::math::opencl_context.debug(std::cout);
}

TEST(opencl_context, platform) {
  std::vector<cl::Platform> all_platforms;
  cl::Platform::get(&all_platforms);
  std::stringstream msg;

  msg << "all_platforms: " << all_platforms.size() << std::endl;
  for (auto platform : all_platforms) {
    msg << "platform name: " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;
  }

  EXPECT_EQ(1, all_platforms.size())
      << "expecting to find one platform" << std::endl
      << msg.str();
}
TEST(opencl_context, devices) {
  cl::Platform platform = cl::Platform::get();
  std::vector<cl::Device> all_devices;
  platform.getDevices(DEVICE_FILTER, &all_devices);

  std::stringstream msg;
  msg << "all_devices: " << all_devices.size() << std::endl;
  for (auto device : all_devices) {
    msg << "- device name: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
  }

  EXPECT_EQ(1, all_devices.size())
      << "expecting to find one device" << std::endl
      << msg.str();
}
#endif
