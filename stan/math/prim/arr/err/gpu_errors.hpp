#ifndef STAN_MATH_PRIM_ARR_ERR_GPU_ERRORS_HPP_
#define STAN_MATH_PRIM_ARR_ERR_GPU_ERRORS_HPP_

#include <iostream>
#include <string>

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.hpp>
#else
#include <CL/cl.hpp>
#endif

/** @file stan/math/prim/arr/err/gpu_errors.hpp
*    @brief errors - utility functions for error checks: 
*    OpenCL and user errors (input dimensions,etc..)
*/

namespace stan {
  namespace math {

    //OpenCL errors - this lists all OpenCL erros, inlcuding the ones from OpenCL 2.0
    void check_ocl_error(cl::Error& e) {   
	    std::string error_string;
      switch (e.err()) {
        case 0:
          error_string="CL_SUCCESS";
          break;
        case -1:
          error_string="CL_DEVICE_NOT_FOUND";
          break;
        case -2:
          error_string="CL_DEVICE_NOT_AVAILABLE";
          break;
        case -3:
          error_string="CL_COMPILER_NOT_AVAILABLE";
          break;
        case -4:
          error_string="CL_MEM_OBJECT_ALLOCATION_FAILURE";
          break;
        case -5:
          error_string="CL_OUT_OF_RESOURCES";
          break;
        case -6:
          error_string="CL_OUT_OF_HOST_MEMORY";
          break;
        case -7:
          error_string="CL_PROFILING_INFO_NOT_AVAILABLE";
          break;
        case -8:
          error_string="CL_MEM_COPY_OVERLAP";
          break;
        case -9:
          error_string="CL_IMAGE_FORMAT_MISMATCH";
          break;
        case -10:
          error_string="CL_IMAGE_FORMAT_NOT_SUPPORTED";
          break;
        case -11:
          error_string="CL_BUILD_PROGRAM_FAILURE";
          break;
        case -12:
          error_string="CL_MAP_FAILURE";
          break;
        case -13:
          error_string="CL_MISALIGNED_SUB_BUFFER_OFFSET";
          break;
        case -14:
          error_string="CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
          break;
        case -15:
          error_string="CL_COMPILE_PROGRAM_FAILURE";
          break;
        case -16:
          error_string="CL_LINKER_NOT_AVAILABLE";
          break;
        case -17:
          error_string="CL_LINK_PROGRAM_FAILURE";
          break;
        case -18:
          error_string="CL_DEVICE_PARTITION_FAILED";
          break;
        case -19:
          error_string="CL_KERNEL_ARG_INFO_NOT_AVAILABLE";
          break;
        case -30:
          error_string="CL_INVALID_VALUE";
          break;
        case -31:
          error_string="CL_INVALID_DEVICE_TYPE";
          break;
        case -32:
          error_string="CL_INVALID_PLATFORM";
          break;
        case -33:
          error_string="CL_INVALID_DEVICE";
          break;
        case -34:
          error_string="CL_INVALID_CONTEXT";
          break;
        case -35:
          error_string="CL_INVALID_QUEUE_PROPERTIES";
          break;
        case -36:
          error_string="CL_INVALID_COMMAND_QUEUE";
          break;
        case -37:
          error_string="CL_INVALID_HOST_PTR";
          break;
        case -38:
          error_string="CL_INVALID_MEM_OBJECT";
          break;
        case -39:
          error_string="CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
          break;
        case -40:
          error_string="CL_INVALID_IMAGE_SIZE";
          break;
        case -41:
          error_string="CL_INVALID_SAMPLER";
          break;
        case -42:
          error_string="CL_INVALID_BINARY";
          break;
        case -43:
          error_string="CL_INVALID_BUILD_OPTIONS";
          break;
        case -44:
          error_string="CL_INVALID_PROGRAM";
          break;
        case -45:
          error_string="CL_INVALID_PROGRAM_EXECUTABLE";
          break;
        case -46:
          error_string="CL_INVALID_KERNEL_NAME";
          break;
        case -47:
          error_string="CL_INVALID_KERNEL_DEFINITION";
          break;
        case -48:
          error_string="CL_INVALID_KERNEL";
          break;
        case -49:
          error_string="CL_INVALID_ARG_INDEX";
          break;
        case -50:
          error_string="CL_INVALID_ARG_VALUE";
          break;
        case -51:
          error_string="CL_INVALID_ARG_SIZE";
          break;
        case -52:
          error_string="CL_INVALID_KERNEL_ARGS";
          break;
        case -53:
          error_string="CL_INVALID_WORK_DIMENSION";
          break;
        case -54:
          error_string="CL_INVALID_WORK_GROUP_SIZE";
          break;
        case -55:
          error_string="CL_INVALID_WORK_ITEM_SIZE";
          break;
        case -56:
          error_string="CL_INVALID_GLOBAL_OFFSET";
          break;
        case -57:
          error_string="CL_INVALID_EVENT_WAIT_LIST";
          break;
        case -58:
          error_string="CL_INVALID_EVENT";
          break;
        case -59:
          error_string="CL_INVALID_OPERATION";
          break;
        case -60:
          error_string="CL_INVALID_GL_OBJECT";
          break;
        case -61:
          error_string="CL_INVALID_BUFFER_SIZE";
          break;
        case -62:
          error_string="CL_INVALID_MIP_LEVEL";
          break;
        case -63:
          error_string="CL_INVALID_GLOBAL_WORK_SIZE";
          break;
        case -64:
          error_string="CL_INVALID_PROPERTY";
          break;
        case -65:
          error_string="CL_INVALID_IMAGE_DESCRIPTOR";
          break;
        case -66:
          error_string="CL_INVALID_COMPILER_OPTIONS";
          break;
        case -67:
          error_string="CL_INVALID_LINKER_OPTIONS";
          break;
        case -68:
          error_string="CL_INVALID_DEVICE_PARTITION_COUNT";
          break;
        case -69:
          error_string="CL_INVALID_PIPE_SIZE";
          break;
        case -70:
          error_string="CL_INVALID_DEVICE_QUEUE";
          break;
        case -1000:
          error_string="CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
          break;
        case -1001:
          error_string="CL_PLATFORM_NOT_FOUND_KHR";
          break;
        case -1002:
          error_string="CL_INVALID_D3D10_DEVICE_KHR";
          break;
        case -1003:
          error_string="CL_INVALID_D3D10_RESOURCE_KHR";
          break;
        case -1004:
          error_string="CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
          break;
        case -1005:
          error_string="CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
          break;
        case -1006:
          error_string="CL_INVALID_D3D11_DEVICE_KHR";
          break;
        case -1007:
          error_string="CL_INVALID_D3D11_RESOURCE_KHR";
          break;
        case -1008:
          error_string="CL_D3D11_RESOURCE_ALREADY_ACQUIRED_KHR";
          break;
        case -1009:
          error_string="CL_D3D11_RESOURCE_NOT_ACQUIRED_KHR";
          break;
        case -101:
          error_string="CL_INVALID_D3D9_DEVICE_NV CL_INVALID_DX9_DEVICE_INTEL";
          break;
        case -1011:
          error_string="CL_INVALID_D3D9_RESOURCE_NV CL_INVALID_DX9_RESOURCE_INTEL";
          break;
        case -1012:
          error_string="CL_D3D9_RESOURCE_ALREADY_ACQUIRED_NV CL_DX9_RESOURCE_ALREADY_ACQUIRED_INTEL";
          break;
        case -1013:
          error_string="CL_D3D9_RESOURCE_NOT_ACQUIRED_NV CL_DX9_RESOURCE_NOT_ACQUIRED_INTEL";
          break;
        case -1092:
          error_string="CL_EGL_RESOURCE_NOT_ACQUIRED_KHR";
          break;
        case -1093:
          error_string="CL_INVALID_EGL_OBJECT_KHR";
          break;
        case -1094:
          error_string="CL_INVALID_ACCELERATOR_INTEL";
          break;
        case -1095:
          error_string="CL_INVALID_ACCELERATOR_TYPE_INTEL";
          break;
        case -1096:
          error_string="CL_INVALID_ACCELERATOR_DESCRIPTOR_INTEL";
          break;
        case -1097:
          error_string="CL_ACCELERATOR_TYPE_NOT_SUPPORTED_INTEL";
          break;
        case -1098:
          error_string="CL_INVALID_VA_API_MEDIA_ADAPTER_INTEL";
          break;
        case -1099:
          error_string="CL_INVALID_VA_API_MEDIA_SURFACE_INTEL";
          break;
        case -1100:
          error_string="CL_VA_API_MEDIA_SURFACE_ALREADY_ACQUIRED_INTEL";
          break;
        case -1101:
          error_string="CL_VA_API_MEDIA_SURFACE_NOT_ACQUIRED_INTEL";
          break;
        case -9999:
          error_string="ILLEGAL_READ_OR_WRITE_NVIDIA";
          break;
        default:
          error_string="number "+e.err();
          break;
      }
      std::cout << "The OpenCL application ended with the error: " << error_string << std::endl;  
      exit(1);
    }
    
    void app_error(std::string text) {
      std::cout << "ERROR: " << text << std::endl;
      exit(1);
    }
  }
}
#endif
