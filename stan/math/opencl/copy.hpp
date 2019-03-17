#ifndef STAN_MATH_PRIM_MAT_FUN_OPENCL_COPY_HPP
#define STAN_MATH_PRIM_MAT_FUN_OPENCL_COPY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/copy.hpp>
#include <stan/math/opencl/kernels/pack.hpp>
#include <stan/math/opencl/kernels/unpack.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <CL/cl.hpp>
#include <iostream>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Copies the source Eigen matrix to
 * the destination matrix that is stored
 * on the OpenCL device.
 *
 * @tparam T type of data in the Eigen matrix
 * @param dst destination matrix on the OpenCL device
 * @param src source Eigen matrix
 *
 * @throw <code>std::invalid_argument</code> if the
 * matrices do not have matching dimensions
 */
template <int R, int C>
void copy(matrix_cl& dst, const Eigen::Matrix<double, R, C>& src) {
  check_size_match("copy (Eigen -> (OpenCL))", "src.rows()", src.rows(),
                   "dst.rows()", dst.rows());
  check_size_match("copy (Eigen -> (OpenCL))", "src.cols()", src.cols(),
                   "dst.cols()", dst.cols());
  if (src.size() > 0) {
    cl::CommandQueue queue = opencl_context.queue();
    try {
      /**
       * Writes the contents of src to the OpenCL buffer
       * starting at the offset 0
       * CL_TRUE denotes that the call is blocking
       * We do not want to execute any further kernels
       * on the device until we are sure that the data is transferred)
       */
      queue.enqueueWriteBuffer(dst.buffer(), CL_TRUE, 0,
                               sizeof(double) * dst.size(), src.data());
    } catch (const cl::Error& e) {
      check_opencl_error("copy Eigen->(OpenCL)", e);
    }
  }
}

/**
 * Copies the source matrix that is stored
 * on the OpenCL device to the destination Eigen
 * matrix.
 *
 * @tparam T type of data in the Eigen matrix
 * @param dst destination Eigen matrix
 * @param src source matrix on the OpenCL device
 *
 * @throw <code>std::invalid_argument</code> if the
 * matrices do not have matching dimensions
 */
template <int R, int C>
void copy(Eigen::Matrix<double, R, C>& dst, const matrix_cl& src) {
  check_size_match("copy ((OpenCL) -> Eigen)", "src.rows()", src.rows(),
                   "dst.rows()", dst.rows());
  check_size_match("copy ((OpenCL) -> Eigen)", "src.cols()", src.cols(),
                   "dst.cols()", dst.cols());
  if (src.size() > 0) {
    cl::CommandQueue queue = opencl_context.queue();
    try {
      /**
       * Reads the contents of the OpenCL buffer
       * starting at the offset 0 to the Eigen
       * matrix
       * CL_TRUE denotes that the call is blocking
       * We do not want to execute any further kernels
       * on the device until we are sure that the data is transferred)
       */
      queue.enqueueReadBuffer(src.buffer(), CL_TRUE, 0,
                              sizeof(double) * dst.size(), dst.data());
    } catch (const cl::Error& e) {
      check_opencl_error("copy (OpenCL)->Eigen", e);
    }
  }
}

/**
 * Packs the flat triagnular matrix on the OpenCL device and
 * copies it to the std::vector.
 *
 * @tparam triangular_view the triangularity of the source matrix
 * @param dst the destination std::vector
 * @param src the flat triangular source matrix on the OpenCL device
 */
template <TriangularViewCL triangular_view>
inline void packed_copy(std::vector<double>& dst, const matrix_cl& src) {
  if (src.size() > 0) {
    int packed_size = src.rows() * (src.rows() + 1) / 2;
    dst.reserve(packed_size);
    cl::CommandQueue queue = opencl_context.queue();
    try {
      cl::Buffer oclBuffer_packed
          = cl::Buffer(opencl_context.context(), CL_MEM_READ_WRITE,
                       sizeof(double) * packed_size);
      stan::math::opencl_kernels::pack(cl::NDRange(src.rows(), src.rows()),
                                       oclBuffer_packed, src.buffer(),
                                       src.rows(), src.rows(), triangular_view);
      queue.enqueueReadBuffer(oclBuffer_packed, CL_TRUE, 0,
                              sizeof(double) * packed_size, dst.data());
    } catch (const cl::Error& e) {
      check_opencl_error("packed_copy (OpenCL)->std::vector", e);
    }
  }
}

/**
 * Copies and unpacks the packed triangular matrix from
 * the source std::vector to the flat matrix_cl on the OpenCL device.
 *
 * @tparam triangular_view the triangularity of the source matrix
 * @param src the packed source std::vector
 * @param dst the destination flat matrix on the OpenCL device
 * @throw <code>std::invalid_argument</code> if the
 * size of the vector does not match the expected size
 * for the packed triangular matrix
 */
template <TriangularViewCL triangular_view>
inline void packed_copy(matrix_cl& dst, const std::vector<double>& src) {
  int packed_size = dst.rows() * (dst.rows() + 1) / 2;
  check_size_match("copy (packed std::vector -> OpenCL)", "src.size()",
                   src.size(), "dst.rows() * (dst.rows() + 1) / 2",
                   packed_size);
  if (src.size() > 0) {
    cl::CommandQueue queue = opencl_context.queue();
    try {
      cl::Buffer oclBuffer_packed
          = cl::Buffer(opencl_context.context(), CL_MEM_READ_WRITE,
                       sizeof(double) * packed_size);
      queue.enqueueWriteBuffer(oclBuffer_packed, CL_TRUE, 0,
                               sizeof(double) * packed_size, src.data());
      stan::math::opencl_kernels::unpack(
          cl::NDRange(dst.rows(), dst.rows()), dst.buffer(), oclBuffer_packed,
          dst.rows(), dst.rows(), triangular_view);
    } catch (const cl::Error& e) {
      check_opencl_error("packed_copy std::vector->OpenCL", e);
    }
  }
}

/**
 * Copies the source matrix to the
 * destination matrix. Both matrices
 * are stored on the OpenCL device.
 *
 * @param dst destination matrix
 * @param src source matrix
 *
 * @throw <code>std::invalid_argument</code> if the
 * matrices do not have matching dimensions
 */
inline void copy(matrix_cl& dst, const matrix_cl& src) {
  check_size_match("copy ((OpenCL) -> (OpenCL))", "src.rows()", src.rows(),
                   "dst.rows()", dst.rows());
  check_size_match("copy ((OpenCL) -> (OpenCL))", "src.cols()", src.cols(),
                   "dst.cols()", dst.cols());
  if (src.size() > 0) {
    try {
      /**
       * Copies the contents of the src buffer to the dst buffer
       * see the matrix_cl(matrix_cl&) constructor
       *  for explanation
       */
      opencl_kernels::copy(cl::NDRange(dst.rows(), dst.cols()), src.buffer(),
                           dst.buffer(), dst.rows(), dst.cols());
    } catch (const cl::Error& e) {
      std::cout << e.err() << std::endl;
      check_opencl_error("copy (OpenCL)->(OpenCL)", e);
    }
  }
}

}  // namespace math
}  // namespace stan
#endif
#endif
