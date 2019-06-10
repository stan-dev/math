#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl_context.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <fstream>

TEST(MathGpu, to_size_t_test) {
  cl::size_t<3> a = stan::math::opencl::to_size_t<3>({1, 2, 3});
  EXPECT_EQ(a[0], 1);
  EXPECT_EQ(a[1], 2);
  EXPECT_EQ(a[2], 3);

  EXPECT_THROW(stan::math::opencl::to_size_t<4>({1, 2, 3, 4}),
               std::domain_error);
}

TEST(MathGpu, getInfo) {
  cl::Context cl = stan::math::opencl_context.context();
  EXPECT_NE("", stan::math::opencl_context.description());
  EXPECT_NE("", stan::math::opencl_context.capabilities());
  EXPECT_GT(stan::math::opencl_context.max_thread_block_size(), 0);
}

TEST(MathGpu, context_construction) {
  cl::Context cv = stan::math::opencl_context.context();
  cl::CommandQueue cq = stan::math::opencl_context.queue();
  std::vector<cl::Device> dv = stan::math::opencl_context.device();
  std::vector<cl::Platform> pl = stan::math::opencl_context.platform();
}

TEST(opencl_context, platform) {
  std::vector<cl::Platform> all_platforms;
  cl::Platform::get(&all_platforms);
  std::stringstream msg;

  msg << "all_platforms: " << all_platforms.size() << std::endl;
  for (auto platform : all_platforms) {
    msg << "platform name: " << platform.getInfo<CL_PLATFORM_NAME>()
        << std::endl;
  }

  EXPECT_GE(all_platforms.size(), 1)
      << "expecting to find at least one platform" << std::endl
      << msg.str();
}

TEST(opencl_context, devices) {
  try {
    std::vector<cl::Platform> all_platforms;
    cl::Platform::get(&all_platforms);
    std::vector<cl::Device> all_devices;
    all_platforms[OPENCL_PLATFORM_ID].getDevices(DEVICE_FILTER, &all_devices);

    std::stringstream msg;
    msg << "all_devices: " << all_devices.size() << std::endl;
    for (auto device : all_devices) {
      msg << "- device name: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
    }

    EXPECT_GE(all_devices.size(), 1)
        << "expecting to find at least one device" << std::endl
        << msg.str();

    msg.str("");
    msg << "max_thead_block_sizes: " << std::endl;
    for (auto device : all_devices) {
      size_t thread_block_size;
      device.getInfo<size_t>(CL_DEVICE_MAX_WORK_GROUP_SIZE, &thread_block_size);
      msg << "- thread_block_size: " << thread_block_size << std::endl;
    }
  } catch (const cl::Error& e) {
    stan::math::check_opencl_error("listing_devices_test", e);
  }
}

TEST(opencl_context, compile_kernel_rawcode) {
  // build dummy kernel
  cl::Context cl = stan::math::opencl_context.context();
  std::vector<cl::Device> dv = stan::math::opencl_context.device();
  const char* dummy_kernel_src
      = "__kernel void dummy(__global const int* foo) { };";
  cl::Program::Sources source(
      1, std::make_pair(dummy_kernel_src, strlen(dummy_kernel_src)));
  cl::Program program_ = cl::Program(cl, source);
  program_.build(dv);
  cl::Kernel dummy_kernel = cl::Kernel(program_, "dummy", NULL);
}
#endif
