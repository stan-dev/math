#ifndef STAN_MATH_PRIM_ARR_ERR_CHECK_OPENCL_HPP_
#define STAN_MATH_PRIM_ARR_ERR_CHECK_OPENCL_HPP_

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include <iostream>
#include <string>
#include <exception>

/** @file stan/math/prim/arr/err/check_opencl.hpp
*    @brief checking OpenCL error numbers
*/

namespace stan {
  namespace math {

    /**
     * Function that throws the OpenCL exception with the
     * given explanation
     *
     * @param function the name of the function where the error occured
     * @param msg information on the OpenCL error
     * 
     * @throw std::domain_error Always.
     */
    inline void throw_openCL(const std::string& function, const char* msg) {
      throw std::domain_error(function +
        ": The OpenCL application ended with the error: " + msg);
    }
    /**
     * Throws the domain error with specifying the OpenCL error that
     * occured. It outputs the OpenCL errors that are specified
     * in OpenCL 2.0. If no matching error number is found,
     * it throws the error with the number. 
     *
     * @param function the name of the function where the error occured
     * @param e The error number
     * 
     * @throw std::domain_error Always.
     */
    inline void check_ocl_error(const std::string& function,
                         const cl::Error& e) {
      switch (e.err()) {
        case 0:
          // CL_SUCCESS - no need to throw
          return;
        case -1:
          throw_openCL(function,
            "CL_DEVICE_NOT_FOUND");
        case -2:
          throw_openCL(function,
            "CL_DEVICE_NOT_AVAILABLE");
        case -3:
          throw_openCL(function,
            "CL_COMPILER_NOT_AVAILABLE");
        case -4:
          throw_openCL(function,
            "CL_MEM_OBJECT_ALLOCATION_FAILURE");
        case -5:
          throw_openCL(function,
            "CL_OUT_OF_RESOURCES");
        case -6:
          throw_openCL(function,
            "CL_OUT_OF_HOST_MEMORY");
        case -7:
          throw_openCL(function,
            "CL_PROFILING_INFO_NOT_AVAILABLE");
        case -8:
          throw_openCL(function,
            "CL_MEM_COPY_OVERLAP");
        case -9:
          throw_openCL(function,
            "CL_IMAGE_FORMAT_MISMATCH");
        case -10:
          throw_openCL(function,
            "CL_IMAGE_FORMAT_NOT_SUPPORTED");
        case -11:
          throw_openCL(function,
            "CL_BUILD_PROGRAM_FAILURE");
        case -12:
          throw_openCL(function,
            "CL_MAP_FAILURE");
        case -13:
          throw_openCL(function,
            "CL_MISALIGNED_SUB_BUFFER_OFFSET");
        case -14:
          throw_openCL(function,
            "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST");
        case -15:
          throw_openCL(function,
            "CL_COMPILE_PROGRAM_FAILURE");
        case -16:
          throw_openCL(function,
            "CL_LINKER_NOT_AVAILABLE");
        case -17:
          throw_openCL(function,
            "CL_LINK_PROGRAM_FAILURE");
        case -18:
          throw_openCL(function,
            "CL_DEVICE_PARTITION_FAILED");
        case -19:
          throw_openCL(function,
            "CL_KERNEL_ARG_INFO_NOT_AVAILABLE");
        case -30:
          throw_openCL(function,
            "CL_INVALID_VALUE");
        case -31:
          throw_openCL(function,
            "CL_INVALID_DEVICE_TYPE");
        case -32:
          throw_openCL(function,
            "CL_INVALID_PLATFORM");
        case -33:
          throw_openCL(function,
            "CL_INVALID_DEVICE");
        case -34:
          throw_openCL(function,
            "CL_INVALID_CONTEXT");
        case -35:
          throw_openCL(function,
            "CL_INVALID_QUEUE_PROPERTIES");
        case -36:
          throw_openCL(function,
            "CL_INVALID_COMMAND_QUEUE");
        case -37:
          throw_openCL(function,
            "CL_INVALID_HOST_PTR");
        case -38:
          throw_openCL(function,
            "CL_INVALID_MEM_OBJECT");
        case -39:
          throw_openCL(function,
            "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR");
        case -40:
          throw_openCL(function,
            "CL_INVALID_IMAGE_SIZE");
        case -41:
          throw_openCL(function,
            "CL_INVALID_SAMPLER");
        case -42:
          throw_openCL(function,
            "CL_INVALID_BINARY");
        case -43:
          throw_openCL(function,
            "CL_INVALID_BUILD_OPTIONS");
        case -44:
          throw_openCL(function,
            "CL_INVALID_PROGRAM");
        case -45:
          throw_openCL(function,
            "CL_INVALID_PROGRAM_EXECUTABLE");
        case -46:
          throw_openCL(function,
            "CL_INVALID_KERNEL_NAME");
        case -47:
          throw_openCL(function,
            "CL_INVALID_KERNEL_DEFINITION");
        case -48:
          throw_openCL(function,
            "CL_INVALID_KERNEL");
        case -49:
          throw_openCL(function,
            "CL_INVALID_ARG_INDEX");
        case -50:
          throw_openCL(function,
            "CL_INVALID_ARG_VALUE");
        case -51:
          throw_openCL(function,
            "CL_INVALID_ARG_SIZE");
        case -52:
          throw_openCL(function,
            "CL_INVALID_KERNEL_ARGS");
        case -53:
          throw_openCL(function,
            "CL_INVALID_WORK_DIMENSION");
        case -54:
          throw_openCL(function,
            "CL_INVALID_WORK_GROUP_SIZE");
        case -55:
          throw_openCL(function,
            "CL_INVALID_WORK_ITEM_SIZE");
        case -56:
          throw_openCL(function,
            "CL_INVALID_GLOBAL_OFFSET");
        case -57:
          throw_openCL(function,
            "CL_INVALID_EVENT_WAIT_LIST");
        case -58:
          throw_openCL(function,
            "CL_INVALID_EVENT");
        case -59:
          throw_openCL(function,
            "CL_INVALID_OPERATION");
        case -60:
          throw_openCL(function,
            "CL_INVALID_GL_OBJECT");
        case -61:
          throw_openCL(function,
            "CL_INVALID_BUFFER_SIZE");
        case -62:
          throw_openCL(function,
            "CL_INVALID_MIP_LEVEL");
        case -63:
          throw_openCL(function,
            "CL_INVALID_GLOBAL_WORK_SIZE");
        case -64:
          throw_openCL(function,
            "CL_INVALID_PROPERTY");
        case -65:
          throw_openCL(function,
            "CL_INVALID_IMAGE_DESCRIPTOR");
        case -66:
          throw_openCL(function,
            "CL_INVALID_COMPILER_OPTIONS");
        case -67:
          throw_openCL(function,
            "CL_INVALID_LINKER_OPTIONS");
        case -68:
          throw_openCL(function,
            "CL_INVALID_DEVICE_PARTITION_COUNT");
        case -69:
          throw_openCL(function,
            "CL_INVALID_PIPE_SIZE");
        case -70:
          throw_openCL(function,
            "CL_INVALID_DEVICE_QUEUE");
        case -1000:
          throw_openCL(function,
           "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR");
        case -1001:
          throw_openCL(function,
            "CL_PLATFORM_NOT_FOUND_KHR");
        case -1002:
          throw_openCL(function,
            "CL_INVALID_D3D10_DEVICE_KHR");
        case -1003:
          throw_openCL(function,
            "CL_INVALID_D3D10_RESOURCE_KHR");
        case -1004:
          throw_openCL(function,
            "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR");
        case -1005:
          throw_openCL(function,
            "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR");
        case -1006:
          throw_openCL(function,
            "CL_INVALID_D3D11_DEVICE_KHR");
        case -1007:
          throw_openCL(function,
            "CL_INVALID_D3D11_RESOURCE_KHR");
        case -1008:
          throw_openCL(function,
            "CL_D3D11_RESOURCE_ALREADY_ACQUIRED_KHR");
        case -1009:
          throw_openCL(function,
            "CL_D3D11_RESOURCE_NOT_ACQUIRED_KHR");
        case -101:
          throw_openCL(function,
           "CL_INVALID_D3D9_DEVICE_NV ");
        case -1011:
          throw_openCL(function, "CL_INVALID_D3D9_RESOURCE_NV ");
        case -1012:
          throw_openCL(function, "CL_D3D9_RESOURCE_ALREADY_ACQUIRED_NV "
           "CL_DX9_RESOURCE_ALREADY_ACQUIRED_INTEL");
        case -1013:
          throw_openCL(function, "CL_D3D9_RESOURCE_NOT_ACQUIRED_NV "
           "CL_DX9_RESOURCE_NOT_ACQUIRED_INTEL");
        case -1092:
          throw_openCL(function,
            "CL_EGL_RESOURCE_NOT_ACQUIRED_KHR");
        case -1093:
          throw_openCL(function,
            "CL_INVALID_EGL_OBJECT_KHR");
        case -1094:
          throw_openCL(function,
            "CL_INVALID_ACCELERATOR_INTEL");
        case -1095:
          throw_openCL(function,
            "CL_INVALID_ACCELERATOR_TYPE_INTEL");
        case -1096:
          throw_openCL(function,
           "CL_INVALID_ACCELERATOR_DESCRIPTOR_INTEL");
        case -1097:
          throw_openCL(function,
            "CL_ACCELERATOR_TYPE_NOT_SUPPORTED_INTEL");
        case -1098:
          throw_openCL(function,
            "CL_INVALID_VA_API_MEDIA_ADAPTER_INTEL");
        case -1099:
          throw_openCL(function,
            "CL_INVALID_VA_API_MEDIA_SURFACE_INTEL");
        case -1100:
          throw_openCL(function,
            "CL_VA_API_MEDIA_SURFACE_ALREADY_ACQUIRED_INTEL");
        case -1101:
          throw_openCL(function,
           "CL_VA_API_MEDIA_SURFACE_NOT_ACQUIRED_INTEL");
        case -9999:
          throw_openCL(function,
            "ILLEGAL_READ_OR_WRITE_NVIDIA");
        default:
          throw_openCL(function,
            std::to_string(e.err()).c_str());
      }
    }

  }
}
#endif
