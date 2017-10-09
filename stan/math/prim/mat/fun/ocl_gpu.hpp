#ifndef STAN_MATH_PRIM_MAT_FUN_OCL_HPP
#define STAN_MATH_PRIM_MAT_FUN_OCL_HPP

#define __CL_ENABLE_EXCEPTIONS

#include "stan/math/prim/arr/err/check_gpu.hpp"
#include "stan/math/prim/mat/kern_gpu/basic_matrix_kernels.hpp"
#include "stan/math/prim/mat/kern_gpu/matrix_multiply_kernels.hpp"
#include "stan/math/prim/mat/kern_gpu/matrix_inverse_kernels.hpp"
#include "stan/math/prim/mat/kern_gpu/cholesky_kernels.hpp"

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.hpp>
#else
#include <CL/cl.hpp>
#endif
#include <iostream>
#include <string>
#include <map>
#include <vector>

#define DEVICE_FILTER CL_DEVICE_TYPE_GPU

/*
*    @file stanmathcl/ocl.hpp
*    @brief Initialization for OpenCL: find platforms, devices,
*      create context, command queue etc.
*/

namespace stan {
  namespace math {

    void compile_kernel_group(std::string group);

    int initialized_ = 0;
    std::map<std::string,std::string> kernel_groups;
    std::map<std::string,std::string> kernel_strings;
    std::map<std::string,cl::Kernel> kernels;
    std::map<std::string, bool> compiled_kernels;
    std::string dummy_kernel = "__kernel void dummy() { };";

    void init_kernel_groups() {
      //to identify the kernel group
      kernel_groups["transpose"] = "basic_matrix";
      kernel_groups["copy"] = "basic_matrix";
      kernel_groups["zeros"] = "basic_matrix";
      kernel_groups["identity"] = "basic_matrix";
      kernel_groups["copy_triangular"] = "basic_matrix";
      kernel_groups["copy_triangular_transposed"] = "basic_matrix";
      kernel_groups["add"] = "basic_matrix";
      kernel_groups["subtract"] = "basic_matrix";
      kernel_groups["scalar_mul_diagonal"] = "matrix_multiply";
      kernel_groups["scalar_mul"] = "matrix_multiply";
      kernel_groups["basic_multiply"] = "matrix_multiply";
      kernel_groups["lower_tri_inv_step1"] = "matrix_inverse";
      kernel_groups["lower_tri_inv_step2"] = "matrix_inverse";
      kernel_groups["lower_tri_inv_step3"] = "matrix_inverse";
      kernel_groups["cholesky_block"] = "cholesky_decomposition";
      kernel_groups["cholesky_left_update"] = "cholesky_decomposition";
      kernel_groups["cholesky_mid_update"] = "cholesky_decomposition";
      kernel_groups["cholesky_zero"] = "cholesky_decomposition";
      kernel_groups["dummy"] = "timing";
      //kernel group strings
      kernel_strings["basic_matrix"] = stan::math::kernel_sources::transpose +
       stan::math::kernel_sources::copy +
       stan::math::kernel_sources::zeros;
      kernel_strings["basic_matrix"] += stan::math::kernel_sources::identity +
       stan::math::kernel_sources::copy_triangular +
       stan::math::kernel_sources::copy_triangular_transposed;
      kernel_strings["basic_matrix"] +=
       stan::math::kernel_sources::add +
       stan::math::kernel_sources::subtract;
      kernel_strings["matrix_multiply"] =
       stan::math::kernel_sources::scalar_mul_diagonal +
       stan::math::kernel_sources::scalar_mul +
       stan::math::kernel_sources::basic_multiply;
      kernel_strings["matrix_inverse"] =
       stan::math::kernel_sources::lower_tri_inv_step1 +
       stan::math::kernel_sources::lower_tri_inv_step2_3;
      kernel_strings["cholesky_decomposition"] =
       stan::math::kernel_sources::cholesky_block +
       stan::math::kernel_sources::cholesky_left_mid_update +
       stan::math::kernel_sources::cholesky_zero;
      kernel_strings["timing"] = stan::math::dummy_kernel;

      //note if the kernels were already compiled
      compiled_kernels["basic_matrix"] = false;
      compiled_kernels["matrix_multiply"] = false;
      compiled_kernels["timing"] = false;
      compiled_kernels["matrix_inverse"] = false;
      compiled_kernels["cholesky_decomposition"] = false;

    }

    //TODO: select some other platform/device than 0
    //TODO: option to turn profiling OFF

    class ocl {
      private:
        std::string description_;
        cl::Context oclContext;
        cl::CommandQueue oclQueue;
        cl::Platform oclPlatform;
        cl::Device oclDevice;

        void init() {
          try {
            std::vector<cl::Platform> allPlatforms;
            cl::Platform::get(&allPlatforms);
            if(allPlatforms.size()  ==  0) {
              std::cout << " No platforms found. " << std::endl;
              exit(1);
            }
            oclPlatform = allPlatforms[0];
            std::vector<cl::Device> allDevices;
            oclPlatform.getDevices(DEVICE_FILTER, &allDevices);
            if(allDevices.size() == 0) {
              std::cout<<" No devices found on the selected platform." <<
               std::endl;
              exit(1);
            }
            oclDevice = allDevices[0];
            description_ = "Device " + oclDevice.getInfo<CL_DEVICE_NAME>() +
             " on the platform " + oclPlatform.getInfo<CL_PLATFORM_NAME>();
            oclContext = cl::Context(allDevices);
            oclQueue = cl::CommandQueue(oclContext, oclDevice,
             CL_QUEUE_PROFILING_ENABLE, NULL);
            init_kernel_groups();
            //compile the dummy kernel used for timing purposes
            cl::Program::Sources source(1,
             std::make_pair(dummy_kernel.c_str(), dummy_kernel.size()));
            cl::Program program_ = cl::Program(oclContext, source);

            try {
              program_.build(allDevices);
              kernels["dummy"] = cl::Kernel(program_, "dummy", NULL);
              compiled_kernels["timing"] = true;
            } catch (cl::Error& e) {
                std::cout << "Building failed, " << e.what() << "(" << e.err()
                 << ")" << "\nRetrieving build log\n" <<
                 program_.getBuildInfo<CL_PROGRAM_BUILD_LOG>(allDevices[0]);
            }
          } catch (cl::Error& e) {
            check_ocl_error(e);
          }
        }

      public:

        std::string description() const {

          if(initialized_ == 1) {
            return description_;
          }
          return "No device selected yet.";
        }

        cl::Context& context() {
          if(initialized_ == 0) {
            init();
            initialized_ = 1;
          }
          return oclContext;
        }

        cl::CommandQueue queue() {
          if(initialized_ == 0) {
            init();
            initialized_ = 1;
          }
          return oclQueue;
        }

    };

    stan::math::ocl ocl_context_queue;



    std::string get_description() {
      return ocl_context_queue.description();
    }

    cl::Context& get_context() {
      return ocl_context_queue.context();
    }

    cl::CommandQueue get_queue() {
      return ocl_context_queue.queue();
    }

    void compile_kernel_group(std::string group) {
        cl::Context ctx = stan::math::get_context();
        std::vector<cl::Device> devices = ctx.getInfo<CL_CONTEXT_DEVICES>();
        std::string kernel_source = kernel_strings[group];
        cl::Program::Sources source(1,
        std::make_pair(kernel_source.c_str(), kernel_source.size()));
        cl::Program program_ = cl::Program(ctx, source);
        try{
          program_.build(devices);

          cl_int err = CL_SUCCESS;
          //iterate over the kernel list and get all the kernels from this group
          for (std::map<std::string,std::string>::iterator
           it = kernel_groups.begin(); it!= kernel_groups.end(); ++it) {
            if(group.compare((it->second).c_str()) == 0) {
              kernels[(it->first).c_str()] = cl::Kernel(program_,
               (it->first).c_str(), &err);
            }
          }
        } catch (cl::Error& e) {
            std::cout << "Building failed, " << e.what() << "(" << e.err() <<
             ")" << "\nRetrieving build log\n" <<
              program_.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
        }



    }

    cl::Kernel& get_kernel(std::string name) {
      //check which kernel group the kernel is assigned to
      cl::Context & ctx = stan::math::get_context();
      //compile the kernel group and return the kernel
      if(!compiled_kernels[kernel_groups[name]]) {
        compile_kernel_group(kernel_groups[name]);
        compiled_kernels[kernel_groups[name]] = true;
      }
      return kernels[name];
    }

  }
}

#endif
