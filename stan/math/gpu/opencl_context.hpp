#ifndef STAN_MATH_GPU_OPENCL_CONTEXT_HPP
#define STAN_MATH_GPU_OPENCL_CONTEXT_HPP

#define __CL_ENABLE_EXCEPTIONS

#include <stan/math/prim/arr/err/check_opencl.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <CL/cl.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
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

namespace {
  typedef std::map<std::string, std::string> map_string;
  typedef std::map<std::string, cl::Kernel> map_kernel;
  typedef std::map<std::string, bool> map_bool;
  static map_string kernel_groups;
  static map_string kernel_strings;
  static map_kernel kernels;
  static map_bool compiled_kernels;
  static std::string dummy_kernel =
    "__kernel void dummy(__global const int* foo) { };";

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
  #include <stan/math/prim/mat/fun/kern_gpu/check_gpu.cl>  // NOLINT
        ;                                                  // NOLINT
    kernel_strings["cholesky_decomposition"] =
  #include <stan/math/prim/mat/fun/kern_gpu/cholesky_decomposition.cl>  // NOLINT
        ;                                                               // NOLINT
    kernel_strings["matrix_inverse"] =
  #include <stan/math/prim/mat/fun/kern_gpu/matrix_inverse.cl>  // NOLINT
        ;                                                       // NOLINT
    kernel_strings["matrix_multiply"] =
  #include <stan/math/prim/mat/fun/kern_gpu/matrix_multiply.cl>  // NOLINT
        ;                                                        // NOLINT
    kernel_strings["basic_matrix"] =
  #include <stan/math/prim/mat/fun/kern_gpu/basic_matrix.cl>  // NOLINT
        ;                                                     // NOLINT

    // Check if the kernels were already compiled
    compiled_kernels["basic_matrix"] = false;
    compiled_kernels["matrix_multiply"] = false;
    compiled_kernels["timing"] = false;
    compiled_kernels["matrix_inverse"] = false;
    compiled_kernels["cholesky_decomposition"] = false;
    compiled_kernels["check_gpu"] = false;
  }
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
class opencl_context {
 private:
  std::string description_;
  cl::Context oclContext_;
  cl::CommandQueue oclQueue_;
  cl::Platform oclPlatform_;
  cl::Device oclDevice_;
  std::vector<cl::Platform> allPlatforms;
  std::vector<cl::Device> allDevices;
  size_t max_workgroup_size_;


 public:
   explicit opencl_context() {
     try {
        cl::Platform::get(&allPlatforms);
        if (allPlatforms.size() == 0) {
          domain_error("OpenCL Initialization", "[Platform]", "",
                       "No OpenCL platforms found");
        }
        oclPlatform_ = allPlatforms[0];

        oclPlatform_.getDevices(DEVICE_FILTER, &allDevices);
        // TODO(Steve): This should throw which platform
        if (allDevices.size() == 0) {
          domain_error("OpenCL Initialization", "[Device]", "",
                       "No OpenCL devices found on the selected platform.");
        }
        oclDevice_ = allDevices[0];
        description_ = "Device " + oclDevice_.getInfo<CL_DEVICE_NAME>()
                       + " on the platform "
                       + oclPlatform_.getInfo<CL_PLATFORM_NAME>();
        allDevices[0].getInfo<size_t>(CL_DEVICE_MAX_WORK_GROUP_SIZE,
                                      &max_workgroup_size_);
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
          domain_error(
              "OpenCL Initialization", e.what(), e.err(),
              "\nRetrieving build log\n",
              program_.getBuildInfo<CL_PROGRAM_BUILD_LOG>(allDevices[0]).c_str());
        }
      } catch (const cl::Error &e) {
        check_ocl_error("build", e);
      }
    }
    // delete copy and move constructors and assign operators
    opencl_context(opencl_context const&) = delete;
    opencl_context(opencl_context&&) = delete;
    opencl_context& operator = (opencl_context const&) = delete;
    opencl_context& operator = (opencl_context &&) = delete;
    inline std::string description() const {
        return description_;
    }

    inline cl::Context &context() {
      return oclContext_;
    }

    inline cl::CommandQueue &queue() {
      return oclQueue_;
    }
    inline int maxWorkgroupSize() { return max_workgroup_size_; }


    /**
     * Returns the description of the OpenCL
     * platform and device that is used.
     *
     */
    inline std::string get_description() { return description(); }

    /**
     * Returns the reference to the
     * OpenCL context. If no context was created,
     * a new context is created.
     *
     */
    inline cl::Context &get_context() { return context(); }
    /**
     * Returns the reference to the active
     * OpenCL command queue. If no context
     * and queue were created,
     * a new context and queue are created and
     * the reference to the new queue is returned.
     *
     */
    inline cl::CommandQueue &get_queue() { return queue(); }
    /**
     * Returns the reference to the active
     * OpenCL command queue. If no context
     * and queue were created,
     * a new context and queue are created and
     * the reference to the new queue is returned.
     *
     */
    inline int get_maximum_workgroup_size() {
      return maxWorkgroupSize();
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
        domain_error(
            "OpenCL Initialization", e.what(), e.err(), "\nRetrieving build log\n",
            program_.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]).c_str());
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

};

static opencl_context opencl_context;


}  // namespace math
}  // namespace stan

#endif
