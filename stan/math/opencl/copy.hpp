#ifndef STAN_MATH_OPENCL_COPY_HPP
#define STAN_MATH_OPENCL_COPY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/cache_copy.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/copy.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <CL/cl.hpp>
#include <iostream>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Copies the source primitive type Eigen matrix to
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
inline void copy(matrix_cl& dst, const Eigen::Matrix<double, R, C>& src) {
  check_size_match("copy (Eigen -> (OpenCL))", "src.rows()", src.rows(),
                   "dst.rows()", dst.rows());
  check_size_match("copy (Eigen -> (OpenCL))", "src.cols()", src.cols(),
                   "dst.cols()", dst.cols());
  if (src.size() > 0) {
    internal::cache_copy(dst.buffer(), src);
  }
}

/**
 * Copies the source var type Eigen matrix to
 * the destination matrix that is stored
 * on the GPU.
 *
 * @tparam T type of data in the Eigen matrix
 * @param dst destination matrix on the GPU
 * @param src source Eigen matrix
 *
 * @throw <code>std::invalid_argument</code> if the
 * matrices do not have matching dimensions
 */
template <int R, int C>
inline void copy(matrix_cl& dst, const Eigen::Matrix<var, R, C>& src) {
  check_size_match("copy (Eigen -> GPU)", "src.rows()", src.rows(),
                   "dst.rows()", dst.rows());
  check_size_match("copy (Eigen -> GPU)", "src.cols()", src.cols(),
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
                               sizeof(double) * dst.size(),
                               value_of_rec(src).data());
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
inline void copy(Eigen::Matrix<double, R, C>& dst, const matrix_cl& src) {
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
      cl::CommandQueue queue = opencl_context.queue();
      /**
       * Copies the contents of the src buffer to the dst buffer
       * see the matrix_cl(matrix_cl&) constructor
       *  for explanation
       */
      cl::Event copy_event;
      queue.enqueueCopyBuffer(src.buffer(), dst.buffer(), 0, 0,
                              sizeof(double) * dst.size(), NULL, &copy_event);
      copy_event.wait();
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
