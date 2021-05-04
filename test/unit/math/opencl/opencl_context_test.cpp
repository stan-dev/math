#ifdef STAN_OPENCL
#include <stan/math/opencl/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <fstream>
#include <string>

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
  const std::string dummy_kernel_src
      = "__kernel void dummy(__global const int* foo) { };";
  cl::Program program_(cl, {dummy_kernel_src});
  program_.build({stan::math::opencl_context.device()});
  cl::Kernel dummy_kernel = cl::Kernel(program_, "dummy");
}

TEST(opencl_context, switch_devices_errors) {
  EXPECT_NO_THROW(stan::math::opencl_context.select_device(0, 0));
  EXPECT_THROW(stan::math::opencl_context.select_device(-1, 0),
               std::system_error);
  EXPECT_THROW(stan::math::opencl_context.select_device(0, -1),
               std::system_error);
  EXPECT_THROW(stan::math::opencl_context.select_device(99999, 0),
               std::system_error);
  EXPECT_THROW(stan::math::opencl_context.select_device(0, 99999),
               std::system_error);
}

// Checks that select_device() works for all devices found on the system. If
// there are multiple devices, this also tests that select_device() calls work
// after another device was already in use.
TEST(opencl_context, switch_devices) {
  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);
  for (int i = 0; i < platforms.size(); i++) {
    std::vector<cl::Device> devices;
    platforms[i].getDevices(DEVICE_FILTER, &devices);
    for (int j = 0; j < devices.size(); j++) {
      stan::math::opencl_context.select_device(i, j);
      Eigen::MatrixXd m(2, 2);
      m << 1, 2, 3, 4;
      Eigen::MatrixXd correct = (2 * m) * (2 * m);
      stan::math::matrix_cl<double> m_cl(m);
      m_cl = 2 * m_cl;
      m_cl = m_cl * m_cl;
      Eigen::MatrixXd res = stan::math::from_matrix_cl(m_cl);
      EXPECT_MATRIX_EQ(res, correct);
    }
  }
  stan::math::opencl_context.select_device(OPENCL_PLATFORM_ID,
                                           OPENCL_DEVICE_ID);
}

#endif
