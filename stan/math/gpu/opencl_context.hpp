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
 *    1. create context
 *    2. Find OpenCL platforms and devices available
 *    3. set up command queue
 *    4. initialize kernel groups
 */

namespace stan {
namespace math {

/**
 * The <code>opencl_context</code> class represents the OpenCL context.
 *
 * See the OpenCL specification glossary for a list of terms:
 * https://www.khronos.org/registry/OpenCL/specs/opencl-1.2.pdf.
 * The context includes the set of devices available on the host, command
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
 /**
  methods:
  * constructor: should be light, initialize things very lightly, test if
  * OpenCL is available grabs context_, queue_, devices_, max_workgroup_size_,
  * and then the  kernel stuff, and the description
  * disable move, copy, and assign operators
  * context(): returns the context
  */
 private:
  const char* description_; // string of platform and device info
  size_t max_workgroup_size_; // The maximum size of a block of workers on GPU
  std::string platform_name_; // The platform such as NVIDIA OpenCL or AMD SDK
  cl::Device device_; // The GPU device(s)
  std::string device_name_; // The name of the GPU
  cl::Context context_; // Manages the the device, queue, platform, memory, etc.
  cl::CommandQueue command_queue_; // job queue for device, one per device
  typedef std::map<const char*, const char*> map_string;
  typedef std::map<const char*, cl::Kernel> map_kernel;
  typedef std::map<const char*, bool> map_bool;
  // TODO(Dan): Did you put this here?
  std::map<const char*, std::tuple<bool, const char*, const char*>> kernel_source;
  void init_kernel_groups();
  void compile_kernel_group(const char* group);

  /**
   * Construct the opencl_context by initializing the
   * OpenCL context, devices, command queues, and kernel
   * groups.
   *
   * For information about
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
  map_string kernel_groups; // map for groups of kernels
  map_string kernel_strings; // map holding kernel code for each kernel
  map_kernel kernels; // map holding compiled kernels in each group
  map_bool check_compiled_kernels; //map for whether a kernel group is compiled

  cl::Kernel get_kernel(const char* name);
  // FIXME:(Steve/Dan) should probably be deleted before merge
  void debug(std::ostream& s) {
    s << "inside opencl_context" << std::endl;
    s << " * platform_name_: " << platform_name_ << std::endl;
  }
  /**
  * Initializes the OpenCL context. This is made with a static local singleton
  * design so that only one context is available.
  */
  static opencl_context& getInstance() {
    static opencl_context instance;
    return instance;
  }

  /**
   * Returns the description of the OpenCL platform and device that is used.
   * Devices will be a GPU and Platforms are a specific OpenCL implimenation
   * such as AMD SDK's or Nvidia's OpenCL implimentation.
   */
  inline const char* description() const { return description_; }

  /**
   * Returns the reference to the OpenCL context. The OpenCL context manages
   * objects such as the device, memory, command queue, program, and kernel
   * objects. For stan, there should only be one context, queue, device, and
   * program with multiple kernels.
   */
  inline cl::Context& context() { return context_; }
  /**
   * Returns the reference to the active OpenCL command queue for the device.
   * One command queue will exist per device where
   * kernels are placed on the command queue and by default executed in order.
   */
  inline cl::CommandQueue& queue() { return command_queue_; }
  /**
   * Returns the maximum workgroup size defined by CL_DEVICE_MAX_WORK_GROUP_SIZE
   * for the device in the context. This is the maximum product of work group
   * dimensions for a particular device. IE a max workgoup of 256 would allow
   * work groups of sizes (16,16), (128,2), (8, 32), etc.
   */
  inline int max_workgroup_size() { return max_workgroup_size_; }
};


/**
 * Initalizes the global std::map variables that
 * hold the OpenCL kernel sources (<code> kernel_strings </code>), the groups to
 * which kernel is assigned (<code> kernel_groups </code>), and the
 * information about which kernel was already compiled
 * (<code> check_compiled_kernels </code>). std::tuple<bool, char*, cl::Kernel>
 */
inline void opencl_context::init_kernel_groups() {
  // Messing around with how kernel_source would look
  kernel_source["dummy"] = {0, "timing",
   "__kernel void dummy(__global const int* foo) { };"};

  kernel_groups["dummy"] = "timing";
  kernel_strings["timing"] =
   "__kernel void dummy(__global const int* foo) { };";
  check_compiled_kernels["timing"] = false;
}

/**
 * Compiles all the kernels in the specified kernel group defined by
 * <code> kernel_groups </code>. The side effect of this
 *  method places all compiled kernels for a group inside of <code> kernels
 *  </code>.
 *
 * @param group_name[in] The kernel group name
 *
 * @throw std::system_error if there are compilation errors
 * when compiling the specified kernel group sources
 */
inline void opencl_context::compile_kernel_group(const char* group_name) {
  cl::Context &ctx = context();
  std::vector<cl::Device> devices = ctx.getInfo<CL_CONTEXT_DEVICES>();
  const char* kernel_source = kernel_strings[group_name];
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
    // Iterate over the kernel list and get all the kernels from this group
    for (auto kernel_group : kernel_groups) {
      if (strcmp(group_name, kernel_group.second) == 0) {
        kernels[(kernel_group.first)]
            = cl::Kernel(program_, kernel_group.first, &err);
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
 * If the kernel has not yet been compiled,
 * the kernel group is compiled first.
 *
 * @brief Passing the name of a kernel from <code> kernel_groups </code> will
 * compile all kernels in the same group as the selected kernel.
 * OpenCL kernels are compiled JIT, instead of compiling each kernel
 * individually this function will compile all kernels
 * in a predefined group. Groupings are made such that kernels commonly
 * called with on another will be compiled at the same time. For example,
 * An arithmetic group of kernels all compiled together contains the kernels
 * for <code> add() </code>, <code> subtract() </code>,
 * and <code> multiply() </code>. This function will only return the kernel
 * which was called, but when a user asks for a kernel within the group those
 * kernels will already be compiled.
 *
 * @param[in] kernel_name The kernel name
 *
 * @return a copy of the cl::Kernel object
 */
inline cl::Kernel opencl_context::get_kernel(const char* kernel_name) {
  // Compile the kernel group and return the kernel
  if (!check_compiled_kernels[kernel_groups[kernel_name]]) {
    compile_kernel_group(kernel_groups[kernel_name]);
    check_compiled_kernels[kernel_groups[kernel_name]] = true;
  }
  return kernels[kernel_name];
}

static opencl_context opencl_context
    = stan::math::opencl_context::getInstance();

}  // namespace math
}  // namespace stan

#endif
#endif
