#ifndef STAN_MATH_OPENCL_COPY_HPP
#define STAN_MATH_OPENCL_COPY_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/kernels/copy.hpp>
#include <stan/math/opencl/kernels/pack.hpp>
#include <stan/math/opencl/kernels/unpack.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/err/check_triangular.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/arr/fun/vec_concat.hpp>

#include <CL/cl.hpp>
#include <algorithm>
#include <iostream>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

/**
 * Copies the source Eigen matrix to
 * the destination matrix that is stored
 * on the OpenCL device.
 *
 * @tparam R Compile time rows of the Eigen matrix
 * @tparam C Compile time columns of the Eigen matrix
 * @param src source Eigen matrix
 * @return matrix_cl with a copy of the data in the source matrix
 */
template <typename T, int R, int C, typename = enable_if_arithmetic<T>>
inline matrix_cl<T> to_matrix_cl(const Eigen::Matrix<T, R, C>& src) {
  matrix_cl<T> dst(src.rows(), src.cols());
  if (src.size() == 0) {
    return dst;
  }
  try {
    /**
     * Writes the contents of src to the OpenCL buffer
     * starting at the offset 0
     * CL_FALSE denotes that the call is non-blocking
     * This means that future kernels need to know about the copy_event
     * So that they do not execute until this transfer finishes.
     */
    cl::Event copy_event;
    const cl::CommandQueue& queue = opencl_context.queue();
    queue.enqueueWriteBuffer(dst.buffer(), CL_FALSE, 0, sizeof(T) * dst.size(),
                             src.data(), &dst.write_events(), &copy_event);
    dst.add_write_event(copy_event);
  } catch (const cl::Error& e) {
    check_opencl_error("copy Eigen->(OpenCL)", e);
  }
  return dst;
}

/**
 * Copies the source matrix that is stored
 * on the OpenCL device to the destination Eigen
 * matrix.
 *
 * @param src source matrix on the OpenCL device
 * @return Eigen matrix with a copy of the data in the source matrix
 */
template <typename T, typename = enable_if_arithmetic<T>>
inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> from_matrix_cl(
    const matrix_cl<T>& src) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dst(src.rows(), src.cols());
  if (src.size() == 0) {
    return dst;
  }
  try {
    /**
     * Reads the contents of the OpenCL buffer
     * starting at the offset 0 to the Eigen
     * matrix
     * CL_TRUE denotes that the call is blocking
     * We do not want to pass data back to the CPU until all of the jobs
     * called on the source matrix are finished.
     */
    cl::Event copy_event;
    const cl::CommandQueue queue = opencl_context.queue();
    queue.enqueueReadBuffer(src.buffer(), CL_FALSE, 0,
                            sizeof(double) * dst.size(), dst.data(),
                            &src.write_events(), &copy_event);
    copy_event.wait();
    src.clear_write_events();
  } catch (const cl::Error& e) {
    check_opencl_error("copy (OpenCL)->Eigen", e);
  }
  return dst;
}

/**
 * Packs the flat triagnular matrix on the OpenCL device and
 * copies it to the std::vector.
 *
 * @param src the flat triangular source matrix on the OpenCL device
 * @return the packed std::vector
 * @throw <code>std::invalid_argument</code> if the matrix is not triangular
 */
template <typename T, typename = enable_if_arithmetic<T>>
inline std::vector<T> packed_copy(const matrix_cl<T>& src) {
  check_triangular("packed_copy", "src", src);
  const int packed_size = src.rows() * (src.rows() + 1) / 2;
  std::vector<T> dst(packed_size);
  if (dst.size() == 0) {
    return dst;
  }
  try {
    const cl::CommandQueue queue = opencl_context.queue();
    matrix_cl<T> packed(packed_size, 1);
    stan::math::opencl_kernels::pack(cl::NDRange(src.rows(), src.rows()),
                                     packed, src, src.rows(), src.rows(),
                                     src.view());
    const std::vector<cl::Event> mat_events
        = vec_concat(packed.read_write_events(), src.write_events());
    cl::Event copy_event;
    queue.enqueueReadBuffer(packed.buffer(), CL_FALSE, 0,
                            sizeof(T) * packed_size, dst.data(), &mat_events,
                            &copy_event);
    copy_event.wait();
    src.clear_write_events();
  } catch (const cl::Error& e) {
    check_opencl_error("packed_copy (OpenCL->std::vector)", e);
  }
  return dst;
}

/**
 * Copies the packed triangular matrix from
 * the source std::vector to an OpenCL buffer and
 * unpacks it to a flat matrix on the OpenCL device.
 *
 * @tparam matrix_view the triangularity of the source matrix
 * @param src the packed source std::vector
 * @param rows the number of rows in the flat matrix
 * @return the destination flat matrix on the OpenCL device
 * @throw <code>std::invalid_argument</code> if the
 * size of the vector does not match the expected size
 * for the packed triangular matrix
 */
template <matrix_cl_view matrix_view, typename T,
          typename = enable_if_arithmetic<T>>
inline matrix_cl<T> packed_copy(const std::vector<T>& src, int rows) {
  const int packed_size = rows * (rows + 1) / 2;
  check_size_match("copy (packed std::vector -> OpenCL)", "src.size()",
                   src.size(), "rows * (rows + 1) / 2", packed_size);
  matrix_cl<T> dst(rows, rows, matrix_view);
  if (dst.size() == 0) {
    return dst;
  }
  try {
    matrix_cl<T> packed(packed_size, 1);
    cl::Event packed_event;
    const cl::CommandQueue queue = opencl_context.queue();
    queue.enqueueWriteBuffer(packed.buffer(), CL_FALSE, 0,
                             sizeof(T) * packed_size, src.data(), NULL,
                             &packed_event);
    packed.add_write_event(packed_event);
    stan::math::opencl_kernels::unpack(cl::NDRange(dst.rows(), dst.rows()), dst,
                                       packed, dst.rows(), dst.rows(),
                                       matrix_view);
  } catch (const cl::Error& e) {
    check_opencl_error("packed_copy (std::vector->OpenCL)", e);
  }
  return dst;
}

/**
 * Copies the source matrix to the
 * destination matrix. Both matrices
 * are stored on the OpenCL device.
 *
 * @param src source matrix
 * @return matrix_cl with copies of values in the source matrix
 * @throw <code>std::invalid_argument</code> if the
 * matrices do not have matching dimensions
 */
template <typename T, typename = enable_if_arithmetic<T>>
inline matrix_cl<T> copy_cl(const matrix_cl<T>& src) {
  matrix_cl<T> dst(src.rows(), src.cols(), src.view());
  if (src.size() == 0) {
    return dst;
  }
  try {
    /**
     * Copies the contents of the src buffer to the dst buffer
     * see the matrix_cl(matrix_cl&) constructor
     *  for explanation
     */
    cl::CommandQueue queue = opencl_context.queue();
    const std::vector<cl::Event> mat_events
        = vec_concat(dst.read_write_events(), src.write_events());
    cl::Event copy_event;
    queue.enqueueCopyBuffer(src.buffer(), dst.buffer(), 0, 0,
                            sizeof(T) * src.size(), &mat_events, &copy_event);
    dst.add_write_event(copy_event);
    src.add_read_event(copy_event);
  } catch (const cl::Error& e) {
    check_opencl_error("copy_cl (OpenCL)->(OpenCL)", e);
  }
  return dst;
}

/**
 * Copy A 1 by 1 source matrix from the Device to  the host.
 * @tparam T An arithmetic type to pass the value from the OpenCL matrix to.
 * @param src A 1x1 matrix on the device.
 * @return dst Arithmetic to receive the matrix_cl value.
 */
template <typename T, typename = enable_if_arithmetic<T>>
inline T from_matrix_cl_error_code(const matrix_cl<T>& src) {
  T dst;
  check_size_match("copy_error_code ((OpenCL) -> (OpenCL))", "src.rows()",
                   src.rows(), "dst.rows()", 1);
  check_size_match("copy_error_code ((OpenCL) -> (OpenCL))", "src.cols()",
                   src.cols(), "dst.cols()", 1);
  try {
    cl::Event copy_event;
    const cl::CommandQueue queue = opencl_context.queue();
    queue.enqueueReadBuffer(src.buffer(), CL_FALSE, 0, sizeof(T), &dst,
                            &src.write_events(), &copy_event);
    copy_event.wait();
    src.clear_write_events();
  } catch (const cl::Error& e) {
    check_opencl_error("copy_error_code (OpenCL)->(OpenCL)", e);
  }
  return dst;
}

/**
 * Copy an arithmetic type to the device.
 * @tparam T An arithmetic type to pass the value from the OpenCL matrix to.
 * @param src Arithmetic to receive the matrix_cl value.
 * @return A 1x1 matrix on the device.
 */
template <typename T, typename = enable_if_arithmetic<T>>
inline matrix_cl<T> to_matrix_cl(const T& src) {
  matrix_cl<T> dst(1, 1);
  check_size_match("to_matrix_cl ((OpenCL) -> (OpenCL))", "src.rows()",
                   dst.rows(), "dst.rows()", 1);
  check_size_match("to_matrix_cl ((OpenCL) -> (OpenCL))", "src.cols()",
                   dst.cols(), "dst.cols()", 1);
  try {
    cl::Event copy_event;
    const cl::CommandQueue queue = opencl_context.queue();
    queue.enqueueWriteBuffer(dst.buffer(), CL_FALSE, 0, sizeof(T), &src,
                             &dst.write_events(), &copy_event);
    dst.add_write_event(copy_event);
  } catch (const cl::Error& e) {
    check_opencl_error("to_matrix_cl (OpenCL)->(OpenCL)", e);
  }
  return dst;
}

}  // namespace math
}  // namespace stan
#endif
#endif
