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
  std::string platform_name_;
  cl::Device device_;
  std::string device_name_;
  cl::Context context_;
  cl::CommandQueue command_queue_;
 //  const char* description_;
 //  size_t max_workgroup_size_;
 //  typedef std::map<const char*, const char*> map_string;
 //  typedef std::map<const char*, cl::Kernel> map_kernel;
 //  typedef std::map<const char*, bool> map_bool;

 // public:
 //  std::map<const char*, std::tuple<bool, char*, cl::Kernel>> kernels_;
 //   map_string kernel_groups;
 //   map_string kernel_strings;
 //   map_kernel kernels;
 //   map_bool compiled_kernels;
 //   const char* dummy_kernel;


 //   void init_kernel_groups();
 //   void init_platforms();
 //   void init_devices();
 //   void init_context_queue();
 //   void init_program();
 //   void compile_kernel_group(const char* group);
 //   cl::Kernel get_kernel(const char* name);

 public:
  void debug(std::ostream& s) {
    s << "inside opencl_context" << std::endl;
    s << " - platform_name_: " << platform_name_ << std::endl;
  }

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
      device_ = all_devices[0];
      device_name_ = device_.getInfo<CL_DEVICE_NAME>();

      context_ = cl::Context(device_);
      command_queue_ = cl::CommandQueue(context_, CL_QUEUE_PROFILING_ENABLE, nullptr);
//   std::ostringstream message;
//   // hack to remove -Waddress, -Wnonnull-compare warnings from GCC 6
//   message << "Device " << oclDevice_.getInfo<CL_DEVICE_NAME>() <<
//    " on the platform " << oclPlatform_.getInfo<CL_PLATFORM_NAME>();
//   std::string description_ = message.str();

    } catch (const cl::Error &e) {
      check_ocl_error("build", e);
    }
  }
 //   opencl_context() {
 //     dummy_kernel =
 //       "__kernel void dummy(__global const int* foo) { };";
 //     try {
 //        init_kernel_groups();
 //        init_platforms();
 //        init_devices();
 //        init_context_queue();
 //        init_program();
 //      } catch (const cl::Error &e) {
 //        check_ocl_error("build", e);
 //      }
 //    }

 //    /*!
 //      The copy and move constructors and assign operators are
 //      disabled
 //    */
 //    opencl_context(opencl_context const&) = delete;
 //    opencl_context(opencl_context&&) = delete;
 //    opencl_context& operator = (opencl_context const&) = delete;
 //    opencl_context& operator = (opencl_context &&) = delete;

 //    /**
 //     * Returns the description of the OpenCL
 //     * platform and device that is used.
 //     *
 //     */
 //    inline const char* description() const {
 //        return description_;
 //    }

 //    /**
 //     * Returns the reference to the
 //     * OpenCL context. If no context was created,
 //     * a new context is created.
 //     *
 //     */
 //    inline cl::Context &context() {
 //      return oclContext_;
 //    }
 //    /**
 //     * Returns the reference to the active
 //     * OpenCL command queue. If no context
 //     * and queue were created,
 //     * a new context and queue are created and
 //     * the reference to the new queue is returned.
 //     *
 //     */
 //    inline cl::CommandQueue &queue() {
 //      return oclQueue_;
 //    }
 //    /**
 //     * Returns the maximum workgroup size for the
 //     * device in the context.
 //     */
 //    inline int maxWorkgroupSize() { return max_workgroup_size_; }
};
// /**
//  * Retrieves the OpenCL platforms on the system and
//  * assigns the first platform as the target platform
//  *
//  * @throw std::logic_error if no OpenCL platforms are found
//  */
// inline void opencl_context::init_platforms() {
//   cl::Platform::get(&allPlatforms);
//   if (allPlatforms.size() == 0) {
//     logic_error("OpenCL Initialization", "[Platform]", "",
//                  "No OpenCL platforms found");
//   }
//   oclPlatform_ = allPlatforms[0];
// }

// /**
//  * Retrieves the devices from the platform.
//  * and assigns the first found device
//  * as the target device.
//  *
//  * @throw std::logic_error if no OpenCL supported devices
//  * are found on the target platform
//  */
// inline void opencl_context::init_devices() {
//   oclPlatform_.getDevices(DEVICE_FILTER, &allDevices);
//   // TODO(Steve): This should throw which platform
//   if (allDevices.size() == 0) {
//     logic_error("OpenCL Initialization", "[Device]", "",
//                  "No OpenCL devices found on the selected platform.");
//   }
//   oclDevice_ = allDevices[0];
// }
// /**
//  * Initializes the OpenCL context and queue.
//  * This function also retrieves the information
//  * on the description and the maximum workgroup size
//  * of the device.
//  *
//  */
// inline void opencl_context::init_context_queue() {
//   std::ostringstream message;
//   // hack to remove -Waddress, -Wnonnull-compare warnings from GCC 6
//   message << "Device " << oclDevice_.getInfo<CL_DEVICE_NAME>() <<
//    " on the platform " << oclPlatform_.getInfo<CL_PLATFORM_NAME>();
//   std::string description_ = message.str();
//   allDevices[0].getInfo<size_t>(CL_DEVICE_MAX_WORK_GROUP_SIZE,
//                                 &max_workgroup_size_);
//   oclContext_ = cl::Context(allDevices);
//   oclQueue_ = cl::CommandQueue(oclContext_, oclDevice_,
//                                CL_QUEUE_PROFILING_ENABLE, NULL);
// }

// /**
//  * Compiles the dummy kernel that is used for
//  * profiling purposes.
//  *
//  * @throw std::logic_error if the dummy kernel has errors
//  *
//  */
// inline void opencl_context::init_program() {
//   // Compile the dummy kernel used for timing purposes
//   cl::Program::Sources source(
//       1, std::make_pair(dummy_kernel, strlen(dummy_kernel)));
//   cl::Program program_ = cl::Program(oclContext_, source);

//   try {
//     program_.build(allDevices);
//     kernels["dummy"] = cl::Kernel(program_, "dummy", NULL);
//     compiled_kernels["timing"] = true;
//   } catch (const cl::Error &e) {
//     logic_error(
//         "OpenCL Initialization", e.what(), e.err(),
//         "\nRetrieving build log\n",
//         program_.getBuildInfo<CL_PROGRAM_BUILD_LOG>(allDevices[0]).c_str());
//   }
// }

// /**
//  * Initalizes the global std::map variables that
//  * hold the OpenCL kernel sources, the groups to
//  * which each kernel is assigned to and the
//  * information about which kernel was already compiled.
//  *
//  */
// inline void opencl_context::init_kernel_groups() {
//   // To identify the kernel group
//   kernel_groups["transpose"] = "basic_matrix";
//   kernel_groups["copy"] = "basic_matrix";
//   kernel_groups["zeros"] = "basic_matrix";
//   kernel_groups["identity"] = "basic_matrix";
//   kernel_groups["copy_triangular"] = "basic_matrix";
//   kernel_groups["copy_triangular_transposed"] = "basic_matrix";
//   kernel_groups["add"] = "basic_matrix";
//   kernel_groups["subtract"] = "basic_matrix";
//   kernel_groups["copy_submatrix"] = "basic_matrix";
//   kernel_groups["scalar_mul_diagonal"] = "matrix_multiply";
//   kernel_groups["scalar_mul"] = "matrix_multiply";
//   kernel_groups["basic_multiply"] = "matrix_multiply";
//   kernel_groups["multiply_self_transposed"] = "matrix_multiply";
//   kernel_groups["lower_tri_inv_step1"] = "matrix_inverse";
//   kernel_groups["lower_tri_inv_step2"] = "matrix_inverse";
//   kernel_groups["lower_tri_inv_step3"] = "matrix_inverse";
//   kernel_groups["cholesky_block"] = "cholesky_decomposition";
//   kernel_groups["check_nan"] = "check_gpu";
//   kernel_groups["check_symmetric"] = "check_gpu";
//   kernel_groups["check_diagonal_zeros"] = "check_gpu";

//   kernel_groups["dummy"] = "timing";

//   // Kernel group strings
//   // the dummy kernel is the only one not included in files
//   // so it is treated before the loop that iterates
//   // through  kernels to load all

//   kernel_strings["timing"] = dummy_kernel;
//   kernel_strings["check_gpu"] =
//   #include <stan/math/prim/mat/fun/kern_gpu/check_gpu.cl>  // NOLINT
//       ;                                                  // NOLINT
//   kernel_strings["cholesky_decomposition"] =
//   #include <stan/math/prim/mat/fun/kern_gpu/cholesky_decomposition.cl>  // NOLINT
//       ;                                                               // NOLINT
//   kernel_strings["matrix_inverse"] =
//   #include <stan/math/prim/mat/fun/kern_gpu/matrix_inverse.cl>  // NOLINT
//       ;                                                       // NOLINT
//   kernel_strings["matrix_multiply"] =
//   #include <stan/math/prim/mat/fun/kern_gpu/matrix_multiply.cl>  // NOLINT
//       ;                                                        // NOLINT
//   kernel_strings["basic_matrix"] =
//   #include <stan/math/prim/mat/fun/kern_gpu/basic_matrix.cl>  // NOLINT
//       ;                                                     // NOLINT

//   // Check if the kernels were already compiled
//   compiled_kernels["basic_matrix"] = false;
//   compiled_kernels["matrix_multiply"] = false;
//   compiled_kernels["timing"] = false;
//   compiled_kernels["matrix_inverse"] = false;
//   compiled_kernels["cholesky_decomposition"] = false;
//   compiled_kernels["check_gpu"] = false;
// }

// /**
//  * Compiles all the kernela in the specified group
//  *
//  * @param group The kernel group name
//  *
//  * @throw std::logic_error if there are compilation errors
//  * when compiling the specified kernel group sources
//  *
//  */
// inline void opencl_context::compile_kernel_group(const char* group) {
//   cl::Context &ctx = context();
//   std::vector<cl::Device> devices = ctx.getInfo<CL_CONTEXT_DEVICES>();
//   const char* kernel_source = kernel_strings[group];
//   cl::Program::Sources source(
//       1, std::make_pair(kernel_source, strlen(kernel_source)));
//   cl::Program program_ = cl::Program(ctx, source);
//   try {
//     char temp[100];
//     int local = 32;
//     int gpu_local_max = sqrt(maxWorkgroupSize());
//     if (gpu_local_max < local)
//       local = gpu_local_max;
//     // parameters that have special limits are for now handled here
//     // kernels with paramters will be compiled separately
//     // for now we have static parameters, so this will be OK
//     snprintf(temp, sizeof(temp), "-D TS=%d -D TS1=%d -D TS2=%d ",
//      local, local, local);
//     program_.build(devices, temp);

//     cl_int err = CL_SUCCESS;
//     // Iterate over the kernel list
//     // and get all the kernels from this group
//     for (auto it : kernel_groups) {
//       if (strcmp(group, it.second) == 0) {
//         kernels[(it.first)]
//             = cl::Kernel(program_, it.first, &err);
//       }
//     }
//   } catch (const cl::Error &e) {
//     logic_error(
//         "\n OpenCL Initialization", e.what(), e.err(),
//         "\nRetrieving build log\n",
//         program_.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]).c_str());
//   }
// }
// /**
//  * Returns the reference to the compiled kernel.
//  * If the kernel has not yet been compiled,
//  * the kernel group is compiled first.
//  *
//  * @param name The kernel name
//  *
//  * @return a copy of the cl::Kernel object
//  */
// inline cl::Kernel opencl_context::get_kernel(const char* name) {
//   // Compile the kernel group and return the kernel
//   if (!compiled_kernels[kernel_groups[name]]) {
//     compile_kernel_group(kernel_groups[name]);
//     compiled_kernels[kernel_groups[name]] = true;
//   }
//   return kernels[name];
// }

static opencl_context opencl_context;



}  // namespace math
}  // namespace stan

#endif
#endif
