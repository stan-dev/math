#ifndef STAN_MATH_PRIM_MAT_FUN_OCL_HPP
#define STAN_MATH_PRIM_MAT_FUN_OCL_HPP

#define __CL_ENABLE_EXCEPTIONS

#include <stan/math/prim/arr/err/check_opencl.hpp>
#include <CL/cl.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#ifdef STAN_DEVICE_CPU
#define DEVICE_FILTER CL_DEVICE_TYPE_CPU
#else
#define DEVICE_FILTER CL_DEVICE_TYPE_GPU
#endif

/*
 *    @file stan/math/prim/mat/fun/ocl.hpp
 *    @brief Initialization for OpenCL: find platforms, devices,
 *      create context, command queue etc.
 */

namespace stan {
namespace math {

inline void compile_kernel_group(std::string group);

static int initialized_ = 0;
typedef std::map<std::string, std::string> map_string;
typedef std::map<std::string, cl::Kernel> map_kernel;
typedef std::map<std::string, bool> map_bool;
static map_string kernel_groups;
static map_string kernel_strings;
static map_kernel kernels;
static map_bool compiled_kernels;
static std::string dummy_kernel = "__kernel void dummy() { };";

/**
 * Initalizes the global std::map variables that
 * hold the OpenCL kernel sources, the groups to
 * which each kernel is assigned to and the
 * information about which kernel was already compiled.
 *
 */
inline void init_kernel_groups() {
  // To identify the kernel group
  kernel_groups["transpose"] = "basic_matrix";
  kernel_groups["copy"] = "basic_matrix";
  kernel_groups["zeros"] = "basic_matrix";
  kernel_groups["identity"] = "basic_matrix";
  kernel_groups["copy_triangular"] = "basic_matrix";
  kernel_groups["copy_triangular_transposed"] = "basic_matrix";
  kernel_groups["add"] = "basic_matrix";
  kernel_groups["subtract"] = "basic_matrix";
  kernel_groups["copy_submatrix"] = "basic_matrix";
  kernel_groups["scalar_mul_diagonal"] = "matrix_multiply";
  kernel_groups["scalar_mul"] = "matrix_multiply";
  kernel_groups["basic_multiply"] = "matrix_multiply";
  kernel_groups["multiply_self_transposed"] = "matrix_multiply";
  kernel_groups["lower_tri_inv_step1"] = "matrix_inverse";
  kernel_groups["lower_tri_inv_step2"] = "matrix_inverse";
  kernel_groups["lower_tri_inv_step3"] = "matrix_inverse";
  kernel_groups["cholesky_block"] = "cholesky_decomposition";
  kernel_groups["check_nan"] = "check_gpu";
  kernel_groups["check_symmetric"] = "check_gpu";
  kernel_groups["check_diagonal_zeros"] = "check_gpu";

  kernel_groups["dummy"] = "timing";

  // Kernel group strings
  // the dummy kernel is the only one not included in files
  // so it is treated before the loop that iterates
  // through  kernels to load all

  kernel_strings["timing"] = dummy_kernel;
  kernel_strings["check_gpu"] =
#include <stan/math/prim/mat/kern_gpu/check_gpu.cl>  // NOLINT
      ;                                              // NOLINT
  kernel_strings["cholesky_decomposition"] =
#include <stan/math/prim/mat/kern_gpu/cholesky_decomposition.cl>  // NOLINT
      ;                                                           // NOLINT
  kernel_strings["matrix_inverse"] =
#include <stan/math/prim/mat/kern_gpu/matrix_inverse.cl>  // NOLINT
      ;                                                   // NOLINT
  kernel_strings["matrix_multiply"] =
#include <stan/math/prim/mat/kern_gpu/matrix_multiply.cl>  // NOLINT
      ;                                                    // NOLINT
  kernel_strings["basic_matrix"] =
#include <stan/math/prim/mat/kern_gpu/basic_matrix.cl>  // NOLINT
      ;                                                 // NOLINT

  // Check if the kernels were already compiled
  compiled_kernels["basic_matrix"] = false;
  compiled_kernels["matrix_multiply"] = false;
  compiled_kernels["timing"] = false;
  compiled_kernels["matrix_inverse"] = false;
  compiled_kernels["cholesky_decomposition"] = false;
  compiled_kernels["check_gpu"] = false;
}

// TODO(Rok): select some other platform/device than 0
// TODO(Rok): option to turn profiling OFF
/**
 * The class that represents the OpenCL runtime
 *
 * The class is used to find the OpenCL supported
 * platforms/devices and initialize the OpenCL
 * context/command queue.
 *
 */
class ocl {
 private:
  std::string description_;
  cl::Context oclContext_;
  cl::CommandQueue oclQueue_;
  cl::Platform oclPlatform_;
  cl::Device oclDevice_;

  size_t max_workgroup_size;

  inline void init() {
    try {
      std::vector<cl::Platform> allPlatforms;
      cl::Platform::get(&allPlatforms);
      if (allPlatforms.size() == 0) {
        std::cout << " No platforms found. " << std::endl;
        exit(1);
      }
      oclPlatform_ = allPlatforms[0];
      std::vector<cl::Device> allDevices;
      oclPlatform_.getDevices(DEVICE_FILTER, &allDevices);
      if (allDevices.size() == 0) {
        std::cout << " No devices found on the selected platform." << std::endl;
        exit(1);
      }
      oclDevice_ = allDevices[0];
      description_ = "Device " + oclDevice_.getInfo<CL_DEVICE_NAME>()
                     + " on the platform "
                     + oclPlatform_.getInfo<CL_PLATFORM_NAME>();
      allDevices[0].getInfo<size_t>(CL_DEVICE_MAX_WORK_GROUP_SIZE,
                                    &max_workgroup_size);
      oclContext_ = cl::Context(allDevices);
      oclQueue_ = cl::CommandQueue(oclContext_, oclDevice_,
                                   CL_QUEUE_PROFILING_ENABLE, NULL);
      init_kernel_groups();
      // Compile the dummy kernel used for timing purposes
      cl::Program::Sources source(
          1, std::make_pair(dummy_kernel.c_str(), dummy_kernel.size()));
      cl::Program program_ = cl::Program(oclContext_, source);

      try {
        program_.build(allDevices);
        kernels["dummy"] = cl::Kernel(program_, "dummy", NULL);
        compiled_kernels["timing"] = true;
      } catch (const cl::Error &e) {
        std::cout << "Building failed, " << e.what() << "(" << e.err() << ")"
                  << "\nRetrieving build log\n"
                  << program_.getBuildInfo<CL_PROGRAM_BUILD_LOG>(allDevices[0]);
      }
    } catch (const cl::Error &e) {
      check_ocl_error("build", e);
    }
  }

 public:
  inline std::string description() const {
    if (initialized_ == 1) {
      return description_;
    }
    return "No device selected yet.";
  }

  inline cl::Context &context() {
    if (initialized_ == 0) {
      init();
      initialized_ = 1;
    }
    return oclContext_;
  }

  inline cl::CommandQueue &queue() {
    if (initialized_ == 0) {
      init();
      initialized_ = 1;
    }
    return oclQueue_;
  }
  inline int maxWorkgroupSize() { return max_workgroup_size; }
};

static ocl ocl_context_queue;

/**
 * Returns the description of the OpenCL
 * platform and device that is used.
 *
 */
inline std::string get_description() { return ocl_context_queue.description(); }

/**
 * Returns the reference to the
 * OpenCL context. If no context was created,
 * a new context is created.
 *
 */
inline cl::Context &get_context() { return ocl_context_queue.context(); }
/**
 * Returns the reference to the active
 * OpenCL command queue. If no context
 * and queue were created,
 * a new context and queue are created and
 * the reference to the new queue is returned.
 *
 */
inline cl::CommandQueue &get_queue() { return ocl_context_queue.queue(); }
/**
 * Returns the reference to the active
 * OpenCL command queue. If no context
 * and queue were created,
 * a new context and queue are created and
 * the reference to the new queue is returned.
 *
 */
inline int get_maximum_workgroup_size() {
  return ocl_context_queue.maxWorkgroupSize();
}

/**
 * Compiles all the kernel in the specified group
 *
 * @param group The kernel group name
 *
 */
inline void compile_kernel_group(std::string group) {
  cl::Context &ctx = get_context();
  std::vector<cl::Device> devices = ctx.getInfo<CL_CONTEXT_DEVICES>();
  std::string kernel_source = kernel_strings[group];
  cl::Program::Sources source(
      1, std::make_pair(kernel_source.c_str(), kernel_source.size()));
  cl::Program program_ = cl::Program(ctx, source);
  try {
    char temp[100];
    int local = 32;
    int gpu_local_max = sqrt(get_maximum_workgroup_size());
    if (gpu_local_max < local)
      local = gpu_local_max;
    // parameters that have special limits are for now handled here
    // kernels with paramters will be compiled separately
    // for now we have static parameters, so this will be OK
    snprintf(temp, sizeof(temp), "-D TS=%d -D TS1=%d -D TS2=%d ", local, local,
             local);
    program_.build(devices, temp);

    cl_int err = CL_SUCCESS;
    // Iterate over the kernel list
    // and get all the kernels from this group
    for (std::map<std::string, std::string>::iterator it
         = kernel_groups.begin();
         it != kernel_groups.end(); ++it) {
      if (group.compare((it->second).c_str()) == 0) {
        kernels[(it->first).c_str()]
            = cl::Kernel(program_, (it->first).c_str(), &err);
      }
    }
  } catch (const cl::Error &e) {
    std::cout << "Building failed, " << e.what() << "(" << e.err() << ")"
              << "\nRetrieving build log\n"
              << program_.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
  }
}

/**
 * Returns the reference to the compiled kernel.
 * If the kernel has not yet been compiled,
 * the kernel group is compiled first.
 *
 * @param name The kernel name
 *
 */
inline cl::Kernel get_kernel(std::string name) {
  // Compile the kernel group and return the kernel
  if (!compiled_kernels[kernel_groups[name]]) {
    compile_kernel_group(kernel_groups[name]);
    compiled_kernels[kernel_groups[name]] = true;
  }
  return kernels[name];
}

}  // namespace math
}  // namespace stan

#endif
