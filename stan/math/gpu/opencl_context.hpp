#ifndef STAN_MATH_GPU_OPENCL_CONTEXT_HPP
#define STAN_MATH_GPU_OPENCL_CONTEXT_HPP
#ifdef STAN_OPENCL
#define __CL_ENABLE_EXCEPTIONS

#include <stan/math/prim/arr/err/check_opencl.hpp>
#include <stan/math/prim/scal/err/logic_error.hpp>
#include <CL/cl.hpp>
#include <cmath>
#include <fstream>
#include <map>
#include <vector>

#define DEVICE_FILTER CL_DEVICE_TYPE_GPU
#ifndef OPENCL_DEVICE
#define OPENCL_DEVICE 0
#endif

/**
 *  @file stan/math/gpu/opencl_context.hpp
 *  @brief Initialization for OpenCL:
 *    1. find OpenCL platforms and devices available
 *    2. create context
 *    3. set up job queue
 *    4. initialize kernel groups
 */

namespace stan {
namespace math {


// TODO(Rok): select some other platform/device than 0
// TODO(Rok): option to turn profiling OFF
/**
 * The <code>opencl_context</code> class represents the OpenCL context.
 *
 * See the OpenCL specification glossary for a list of terms:
 * https://www.khronos.org/registry/OpenCL/specs/opencl-1.2.pdf. The
 * context includes the set of devices available on the host, command
 * queues, manages kernels.
 *
 * This is designed so there's only one instance running on the host.
 *
 * Some design decisions that may need to be addressed later:
 * - we are assuming a single OpenCL platform. We may want to run on multiple platforms simulatenously
 * - we are assuming a single OpenCL device. We may want to run on multiple devices simulatenously
 */
class opencl_context {
  // what we need:
  // - description_? opencl platform and device
  // - oclContext_ -> context_
  // - oclQueue_ -> command_queue_: needs context + device
  // - ** oclPlatform_ -> platform_: needed to get devices and platform name; not needed
  // - oclDevice_ -> device_: should be multiple?
  // - max_workgroup_size_: ?
  // - kernels_groups: kernel function name -> which kernel group it actually lives
  // - kernels_strings: kernel group -> actual string of the kernel
  // - kernels: kernel group -> cl::Kernel
  // - compiled_kernels: kernel function name -> bool
  // - ** dummy_kernel: just make it a local
  //
  // methods:
  // - constructor: should be light, initialize things very lightly, test if OpenCL is available
  //     grabs context_, queue_, devices_, max_workgroup_size_?, and then the kernel stuff, and the description
  // - disable move, copy, and assign operators
  // - context(): returns the context
  // - command_queue(): returns the command_queue; it should be one per device
  // - get_kernel(const char* name): returns the appropriate kernel
 private:
  const char* description_;
  size_t max_workgroup_size_;
  std::string platform_name_;
  cl::Device device_;
  std::string device_name_;
  cl::Context context_;
  cl::CommandQueue command_queue_;
  /**
   * Construct the opencl_context by initializing the
   * OpenCL context, devices, command queues, and kernel
   * groups.
   *
   * @throw std::domain_error if an OpenCL error occurs
   */
  opencl_context() {
    try {
      cl::Platform platform = cl::Platform::get();
      platform_name_ = platform.getInfo<CL_PLATFORM_NAME>();

      // FIXME(dl): this logic is wrong! I have two devices on my machine. How do we choose
      //  which device to use?
      std::vector<cl::Device> all_devices;
      platform.getDevices(DEVICE_FILTER, &all_devices);
      if (all_devices.size() == 0) {
        logic_error("OpenCL Initialization", "[Device]", platform_name_,
                    "No OpenCL devices found on the selected platform: ");
      }
      device_ = all_devices[OPENCL_DEVICE];
      context_ = cl::Context(all_devices);
      command_queue_ = cl::CommandQueue(context_, device_,
         CL_QUEUE_PROFILING_ENABLE, nullptr);

      const char* dummy_kernel_src
       = "__kernel void dummy(__global const int* foo) { };";
      cl::Program::Sources source(
          1, std::make_pair(dummy_kernel_src, strlen(dummy_kernel_src)));
      cl::Program program_ = cl::Program(context_, source);
      // build dummy kernel
      try {
        program_.build(all_devices);
        cl::Kernel dummy_kernel = cl::Kernel(program_, "dummy", NULL);
      } catch (const cl::Error &e) {
        logic_error(
            "OpenCL Initialization", e.what(), e.err(),
            "\nRetrieving build log\n",
            program_.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device_).c_str());
      }
      // device info
      std::ostringstream description_message;
      device_.getInfo<size_t>(CL_DEVICE_MAX_WORK_GROUP_SIZE,
                              &max_workgroup_size_);
      device_name_ = device_.getInfo<CL_DEVICE_NAME>();
      description_message << "Device " << device_name_ <<
      " on the platform " << platform_name_;
      description_ = description_message.str().c_str();

    } catch (const cl::Error &e) {
      check_ocl_error("build", e);
    }
  }

 public:
  void debug(std::ostream& s) {
    s << "inside opencl_context" << std::endl;
    s << " - platform_name_: " << platform_name_ << std::endl;
  }
  static opencl_context& getInstance() {
      static opencl_context instance;
      return instance;
  }
  /*!
    The copy and move constructors and assign operators are
    disabled
  */
  //opencl_context(opencl_context const&) = delete;
  //opencl_context(opencl_context&&) = delete;
  //opencl_context& operator = (opencl_context const&) = delete;
  //opencl_context& operator = (opencl_context &&) = delete;

  /**
  * Returns the description of the OpenCL
  * platform and device that is used.
  *
  */
  inline const char* description() const {
     return description_;
  }

  /**
  * Returns the reference to the
  * OpenCL context. If no context was created,
  * a new context is created.
  *
  */
  inline cl::Context &context() {
   return context_;
  }
  /**
  * Returns the reference to the active
  * OpenCL command queue. If no context
  * and queue were created,
  * a new context and queue are created and
  * the reference to the new queue is returned.
  *
  */
  inline cl::CommandQueue &queue() {
   return command_queue_;
  }
  /**
  * Returns the maximum workgroup size for the
  * device in the context.
  */
  inline int maxWorkgroupSize() { return max_workgroup_size_; }
};


static opencl_context opencl_context = stan::math::opencl_context::getInstance();

}  // namespace math
}  // namespace stan

#endif
#endif
