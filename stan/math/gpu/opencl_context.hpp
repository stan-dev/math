#ifndef STAN_MATH_GPU_OPENCL_CONTEXT_HPP
#define STAN_MATH_GPU_OPENCL_CONTEXT_HPP
#ifdef STAN_OPENCL
#define __CL_ENABLE_EXCEPTIONS

#include <stan/math/prim/arr/err/check_opencl.hpp>
#include <stan/math/prim/scal/err/system_error.hpp>
#include <CL/cl.hpp>
#include <cmath>
#include <fstream>
#include <map>
#include <vector>
#include <system_error>
#include <cerrno>

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
 * - we are assuming a single OpenCL platform. We may want to run on multiple
 * platforms simulatenously
 * - we are assuming a single OpenCL device. We may want to run on multiple
 * devices simulatenously
 */
class opencl_context {
  /*
  what we need:
  - description_ opencl platform and device
  - context_: Holds the device, queue and platform info
  - command_queue_: job queue for devices in context, one queue per device
  - platform_: The platform such as NVIDIA OpenCL or AMD SDK;
  - device_: The GPU device(s)
  - max_workgroup_size_: The maximum size of a block of workers on GPU
  - kernels_groups: map for groups of kernels
  - kernels_strings: map holding kernel code for groups
  - kernels: map holding compiled kernels in each group
  - compiled_kernels: map for checking whether a kernel group is compiled
  //
  methods:
  - constructor: should be light, initialize things very lightly, test if
  OpenCL is available
      grabs context_, queue_, devices_, max_workgroup_size_?, and then the
      kernel stuff, and the description
  - disable move, copy, and assign operators
  - context(): returns the context
  - command_queue(): returns the command_queue; it should be one per device
  - get_kernel(const char* name): returns the appropriate kernel
  */
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
   * @throw std::system_error if an OpenCL error occurs
   */
  opencl_context() {
    try {
      // platform
      cl::Platform platform = cl::Platform::get();
      platform_name_ = platform.getInfo<CL_PLATFORM_NAME>();
      // device setup
      std::vector<cl::Device> all_devices;
      platform.getDevices(DEVICE_FILTER, &all_devices);
      if (all_devices.size() == 0) {
        system_error("OpenCL Initialization", "[Device]", -1,
                     "CL_DEVICE_NOT_FOUND");
      }
      device_ = all_devices[OPENCL_DEVICE];
      // context and queue
      context_ = cl::Context(all_devices);
      command_queue_ = cl::CommandQueue(context_, device_,
                                        CL_QUEUE_PROFILING_ENABLE, nullptr);
      // build dummy kernel
      const char* dummy_kernel_src
          = "__kernel void dummy(__global const int* foo) { };";
      cl::Program::Sources source(
          1, std::make_pair(dummy_kernel_src, strlen(dummy_kernel_src)));
      cl::Program program_ = cl::Program(context_, source);
      try {
        program_.build(all_devices);
        cl::Kernel dummy_kernel = cl::Kernel(program_, "dummy", NULL);
      } catch (const cl::Error& e) {
        system_error(
            "OpenCL Initialization", e.what(), e.err(),
            "\nRetrieving build log\n",
            program_.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device_).c_str());
      }
      // device info
      std::ostringstream description_message;
      device_.getInfo<size_t>(CL_DEVICE_MAX_WORK_GROUP_SIZE,
                              &max_workgroup_size_);
      device_name_ = device_.getInfo<CL_DEVICE_NAME>();
      description_message << "Device " << device_name_ << " on the platform "
                          << platform_name_;
      description_ = description_message.str().c_str();
      // setup kernel groups
      init_kernel_groups();

    } catch (const cl::Error& e) {
      check_ocl_error("build", e);
    }
  }

 public:
  typedef std::map<const char*, const char*> map_string;
  typedef std::map<const char*, cl::Kernel> map_kernel;
  typedef std::map<const char*, bool> map_bool;
  std::map<const char*, std::tuple<bool, char*, cl::Kernel>> kernels_;
  map_string kernel_groups;
  map_string kernel_strings;
  map_kernel kernels;
  map_bool compiled_kernels;
  void init_kernel_groups();
  void compile_kernel_group(const char* group);
  cl::Kernel get_kernel(const char* name);
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
  // opencl_context(opencl_context const&) = delete;
  // opencl_context(opencl_context&&) = delete;
  // opencl_context& operator = (opencl_context const&) = delete;
  // opencl_context& operator = (opencl_context &&) = delete;

  /**
   * Returns the description of the OpenCL
   * platform and device that is used.
   *
   */
  inline const char* description() const { return description_; }

  /**
   * Returns the reference to the
   * OpenCL context. If no context was created,
   * a new context is created.
   *
   */
  inline cl::Context& context() { return context_; }
  /**
   * Returns the reference to the active
   * OpenCL command queue. If no context
   * and queue were created,
   * a new context and queue are created and
   * the reference to the new queue is returned.
   *
   */
  inline cl::CommandQueue& queue() { return command_queue_; }
  /**
   * Returns the maximum workgroup size for the
   * device in the context.
   */
  inline int max_workgroup_size() { return max_workgroup_size_; }
};


/**
 * Initalizes the global std::map variables that
 * hold the OpenCL kernel sources, the groups to
 * which each kernel is assigned to and the
 * information about which kernel was already compiled.
 *
 */
inline void opencl_context::init_kernel_groups() {

  kernel_groups["dummy"] = "timing";

  /*Kernel group strings
  the dummy kernel is the only one not included in files
  so it is treated before the loop that iterates
  through  kernels to load all
  */
  kernel_strings["timing"] =
   "__kernel void dummy(__global const int* foo) { };";
  compiled_kernels["timing"] = false;
}

/**
 * Compiles all the kernels in the specified group
 *
 * @param group The kernel group name
 *
 * @throw std::system_error if there are compilation errors
 * when compiling the specified kernel group sources
 *
 */
inline void opencl_context::compile_kernel_group(const char* group) {
  cl::Context &ctx = context();
  std::vector<cl::Device> devices = ctx.getInfo<CL_CONTEXT_DEVICES>();
  const char* kernel_source = kernel_strings[group];
  cl::Program::Sources source(
      1, std::make_pair(kernel_source, strlen(kernel_source)));
  cl::Program program_ = cl::Program(ctx, source);
  try {
    char temp[100];
    int local = 32;
    int gpu_local_max = sqrt(max_workgroup_size());
    if (gpu_local_max < local)
      local = gpu_local_max;
    /*
    parameters that have special limits are for now handled here
    kernels with parameters will be compiled separately
    for now we have static parameters, so this will be OK
    */
    snprintf(temp, sizeof(temp), "-D TS=%d -D TS1=%d -D TS2=%d ",
     local, local, local);
    program_.build(devices, temp);

    cl_int err = CL_SUCCESS;
    /*
    Iterate over the kernel list
    and get all the kernels from this group
    */
    for (auto it : kernel_groups) {
      if (strcmp(group, it.second) == 0) {
        kernels[(it.first)]
            = cl::Kernel(program_, it.first, &err);
      }
    }
  } catch (const cl::Error &e) {
    system_error(
        "OpenCL Initialization", e.what(), e.err(),
        "\nRetrieving build log\n",
        program_.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device_).c_str());
  }
}

/**
 * Returns the reference to the compiled kernel.
 * If the kernel has not yet been compiled,
 * the kernel group is compiled first.
 *
 * @param name The kernel name
 *
 * @return a copy of the cl::Kernel object
 */
inline cl::Kernel opencl_context::get_kernel(const char* name) {
  // Compile the kernel group and return the kernel
  if (!compiled_kernels[kernel_groups[name]]) {
    compile_kernel_group(kernel_groups[name]);
    compiled_kernels[kernel_groups[name]] = true;
  }
  return kernels[name];
}

static opencl_context opencl_context
    = stan::math::opencl_context::getInstance();

}  // namespace math
}  // namespace stan

#endif
#endif
