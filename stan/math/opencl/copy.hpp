#ifndef STAN_MATH_OPENCL_COPY_HPP
#define STAN_MATH_OPENCL_COPY_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/vec_concat.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/value_type.hpp>
#include <stan/math/opencl/kernels/copy.hpp>
#include <stan/math/opencl/kernels/pack.hpp>
#include <stan/math/opencl/kernels/unpack.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/err/check_triangular.hpp>

#include <CL/cl2.hpp>
#include <algorithm>
#include <iostream>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

/** \ingroup opencl
 * Copies the source Eigen matrix, `std::vector` or scalar to the destination
 * matrix that is stored on the OpenCL device. The function also accepts
 * `matrix_cl`s in which case it just returns the argument. If a lvalue matrix
 * is passed to this function the caller must make sure that the matrix does not
 * go out of scope before copying is complete.
 *
 * That means `.wait()` must be called on the event associated on copying or
 * any other event that requires completion of this event. This can be done by
 * calling `.wait_for_write_events()` or `.wait_for_read_write_events()` on
 * returned matrix or any matrix that is calculated from that one.
 *
 * @param src source Eigen matrix
 * @return matrix_cl with a copy of the data in the source matrix
 */
template <typename T, require_vt_arithmetic<T>* = nullptr>
inline matrix_cl<value_type_t<T>> to_matrix_cl(T&& src) {
  return matrix_cl<value_type_t<T>>(std::forward<T>(src));
}

/** \ingroup opencl
 * Copies the source matrix that is stored on the OpenCL device to the
 * destination Eigen matrix.
 *
 * @param src source matrix on the OpenCL device
 * @return Eigen matrix with a copy of the data in the source matrix
 */
template <int R = Eigen::Dynamic, int C = Eigen::Dynamic, typename T,
          typename = require_arithmetic_t<T>>
inline Eigen::Matrix<T, R, C> from_matrix_cl(const matrix_cl<T>& src) {
  Eigen::Matrix<T, R, C> dst(src.rows(), src.cols());
  if (src.size() == 0) {
    return dst;
  }
  try {
    /** \ingroup opencl
     * Reads the contents of the OpenCL buffer
     * starting at the offset 0 to the Eigen
     * matrix
     * CL_TRUE denotes that the call is blocking
     * We do not want to pass data back to the CPU until all of the jobs
     * called on the source matrix are finished.
     */
    cl::Event copy_event;
    const cl::CommandQueue queue = opencl_context.queue();
    queue.enqueueReadBuffer(src.buffer(), opencl_context.in_order(), 0,
                            sizeof(T) * dst.size(), dst.data(),
                            &src.write_events(), &copy_event);
    copy_event.wait();
    src.clear_write_events();
  } catch (const cl::Error& e) {
    check_opencl_error("copy (OpenCL)->Eigen", e);
  }
  return dst;
}

/** \ingroup opencl
 * Copies result of a kernel generator expression to the
 * destination Eigen matrix.
 *
 * @tparam R rows type of the destination
 * @tparam C cols type of the destination
 * @tparam T type of expression
 * @param src source expression
 * @return Eigen matrix with a copy of the data in the source matrix
 */
template <int R = Eigen::Dynamic, int C = Eigen::Dynamic, typename T,
          require_all_kernel_expressions_t<T>* = nullptr,
          require_not_matrix_cl_t<T>* = nullptr>
inline Eigen::Matrix<value_type_t<T>, R, C> from_matrix_cl(const T& src) {
  return from_matrix_cl<R, C>(src.eval());
}

/** \ingroup opencl
 * Packs the flat triangular matrix on the OpenCL device and
 * copies it to the std::vector.
 *
 * @param src the flat triangular source matrix on the OpenCL device
 * @return the packed std::vector
 * @throw <code>std::invalid_argument</code> if the matrix is not triangular
 */
template <typename T, typename = require_arithmetic_t<T>>
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
    queue.enqueueReadBuffer(packed.buffer(), opencl_context.in_order(), 0,
                            sizeof(T) * packed_size, dst.data(), &mat_events,
                            &copy_event);
    copy_event.wait();
    src.clear_write_events();
  } catch (const cl::Error& e) {
    check_opencl_error("packed_copy (OpenCL->std::vector)", e);
  }
  return dst;
}

/** \ingroup opencl
 * Copies the packed triangular matrix from the source std::vector to an OpenCL
 * buffer and unpacks it to a flat matrix on the OpenCL device. If a lvalue is
 * passed to this constructor the caller must make sure that it does not go out
 * of scope before copying is complete.
 *
 * That means `.wait()` must be called on the event associated on copying or
 * any other event that requires completion of this event. This can be done by
 * calling `.wait_for_write_events()` or `.wait_for_read_write_events()` on
 * returned matrix or any matrix that is calculated from that one.
 *
 * @tparam matrix_view the triangularity of the source matrix
 * @param src the packed source std::vector
 * @param rows the number of rows in the flat matrix
 * @return the destination flat matrix on the OpenCL device
 * @throw <code>std::invalid_argument</code> if the
 * size of the vector does not match the expected size
 * for the packed triangular matrix
 */
template <matrix_cl_view matrix_view, typename Vec,
          typename Vec_scalar = scalar_type_t<Vec>,
          require_vector_vt<std::is_arithmetic, Vec>* = nullptr>
inline matrix_cl<Vec_scalar> packed_copy(Vec&& src, int rows) {
  const int packed_size = rows * (rows + 1) / 2;
  check_size_match("copy (packed std::vector -> OpenCL)", "src.size()",
                   src.size(), "rows * (rows + 1) / 2", packed_size);
  matrix_cl<Vec_scalar> dst(rows, rows, matrix_view);
  if (dst.size() == 0) {
    return dst;
  }
  try {
    matrix_cl<Vec_scalar> packed(packed_size, 1);
    cl::Event packed_event;
    const cl::CommandQueue queue = opencl_context.queue();
    queue.enqueueWriteBuffer(
        packed.buffer(),
        opencl_context.in_order() || std::is_rvalue_reference<Vec&&>::value, 0,
        sizeof(Vec_scalar) * packed_size, src.data(), nullptr, &packed_event);
    packed.add_write_event(packed_event);
    stan::math::opencl_kernels::unpack(cl::NDRange(dst.rows(), dst.rows()), dst,
                                       packed, dst.rows(), dst.rows(),
                                       matrix_view);
  } catch (const cl::Error& e) {
    check_opencl_error("packed_copy (std::vector->OpenCL)", e);
  }
  return dst;
}

/** \ingroup opencl
 * Copies the source matrix to the
 * destination matrix. Both matrices
 * are stored on the OpenCL device.
 *
 * @tparam T An arithmetic type to pass the value from the OpenCL matrix to.
 * @param src the source matrix
 * @return matrix_cl with copies of values in the source matrix
 * @throw <code>std::invalid_argument</code> if the
 * matrices do not have matching dimensions
 */
template <typename T, typename = require_arithmetic_t<T>>
inline matrix_cl<T> copy_cl(const matrix_cl<T>& src) {
  return matrix_cl<T>(src);
}

/** \ingroup opencl
 * Copy A 1 by 1 source matrix from the Device to  the host.
 * @tparam T An arithmetic type to pass the value from the OpenCL matrix to.
 * @param src A 1x1 matrix on the device.
 * @return dst Arithmetic to receive the matrix_cl value.
 */
template <typename T, typename = require_arithmetic_t<T>>
inline T from_matrix_cl_error_code(const matrix_cl<T>& src) {
  T dst;
  check_size_match("copy_error_code ((OpenCL) -> (OpenCL))", "src.rows()",
                   src.rows(), "dst.rows()", 1);
  check_size_match("copy_error_code ((OpenCL) -> (OpenCL))", "src.cols()",
                   src.cols(), "dst.cols()", 1);
  try {
    cl::Event copy_event;
    const cl::CommandQueue queue = opencl_context.queue();
    queue.enqueueReadBuffer(src.buffer(), opencl_context.in_order(), 0,
                            sizeof(T), &dst, &src.write_events(), &copy_event);
    copy_event.wait();
    src.clear_write_events();
  } catch (const cl::Error& e) {
    check_opencl_error("copy_error_code (OpenCL)->(OpenCL)", e);
  }
  return dst;
}

}  // namespace math
}  // namespace stan
#endif
#endif
