#ifndef STAN_MATH_OPENCL_COPY_HPP
#define STAN_MATH_OPENCL_COPY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/value_type.hpp>
#include <stan/math/opencl/scalar_type.hpp>
#include <stan/math/opencl/kernels/pack.hpp>
#include <stan/math/opencl/kernels/unpack.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/err/check_triangular.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/vec_concat.hpp>

#include <CL/opencl.hpp>
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
template <typename T, require_st_arithmetic<T>* = nullptr>
inline matrix_cl<scalar_type_t<T>> to_matrix_cl(T&& src) {
  return matrix_cl<scalar_type_t<T>>(std::forward<T>(src));
}

/** \ingroup opencl
 * Copies the source matrix that is stored on the OpenCL device to the
 * destination Eigen matrix.
 *
 * @tparam T_ret destination type
 * @tparam T scalar in the source matrix
 * @param src source matrix on the OpenCL device
 * @return Eigen matrix with a copy of the data in the source matrix
 */
template <typename T_ret, typename T, require_eigen_t<T_ret>* = nullptr,
          require_matrix_cl_t<T>* = nullptr,
          require_st_same<T_ret, T>* = nullptr>
inline auto from_matrix_cl(const T& src) {
  using T_val = value_type_t<T>;
  using T_ret_col_major
      = Eigen::Matrix<scalar_type_t<T_ret>, T_ret::RowsAtCompileTime,
                      T_ret::ColsAtCompileTime>;
  T_ret_col_major dst(src.rows(), src.cols());
  if (src.size() == 0) {
    return dst;
  }
  if ((src.view() == matrix_cl_view::Lower
       || src.view() == matrix_cl_view::Upper)
      && src.rows() == src.cols()) {
    using T_not_bool
        = std::conditional_t<std::is_same<T_val, bool>::value, char, T_val>;
    std::vector<T_not_bool> packed = packed_copy(src);

    size_t pos = 0;
    if (src.view() == matrix_cl_view::Lower) {
      for (int j = 0; j < src.cols(); ++j) {
        for (int k = 0; k < j; ++k) {
          dst.coeffRef(k, j) = 0;
        }
        for (int i = j; i < src.cols(); ++i) {
          dst.coeffRef(i, j) = packed[pos++];
        }
      }
    } else {
      for (int j = 0; j < src.cols(); ++j) {
        for (int i = 0; i <= j; ++i) {
          dst.coeffRef(i, j) = packed[pos++];
        }
        for (int k = j + 1; k < src.cols(); ++k) {
          dst.coeffRef(k, j) = 0;
        }
      }
    }
  } else {
    try {
      cl::Event copy_event;
      const cl::CommandQueue queue = opencl_context.queue();
      queue.enqueueReadBuffer(src.buffer(), opencl_context.in_order(), 0,
                              sizeof(T_val) * dst.size(), dst.data(),
                              &src.write_events(), &copy_event);
      copy_event.wait();
      src.clear_write_events();
    } catch (const cl::Error& e) {
      check_opencl_error("copy (OpenCL)->Eigen", e);
    }
    if (!contains_nonzero(src.view(), matrix_cl_view::Lower)) {
      dst.template triangularView<Eigen::StrictlyLower>()
          = T_ret_col_major::Zero(dst.rows(), dst.cols());
    }
    if (!contains_nonzero(src.view(), matrix_cl_view::Upper)) {
      dst.template triangularView<Eigen::StrictlyUpper>()
          = T_ret_col_major::Zero(dst.rows(), dst.cols());
    }
  }
  return dst;
}

/** \ingroup opencl
 * Copies result of a kernel generator expression to the specified
 * destination type.
 *
 * @tparam T_ret destination type
 * @tparam T type of expression
 * @param src source expression
 * @return Eigen matrix with a copy of the data in the source matrix
 */
template <typename T_ret, typename T,
          require_all_kernel_expressions_t<T>* = nullptr,
          require_not_matrix_cl_t<T>* = nullptr>
inline auto from_matrix_cl(const T& src) {
  return from_matrix_cl<T_ret>(src.eval());
}

/** \ingroup opencl
 * Copy A 1 by 1 source matrix from the Device to the host.
 * @tparam T An arithmetic type to pass the value from the OpenCL matrix to.
 * @param src A 1x1 matrix on the device.
 * @return dst Arithmetic to receive the matrix_cl value.
 */
template <typename T_dst, typename T, require_arithmetic_t<T>* = nullptr,
          require_same_t<T_dst, T>* = nullptr>
inline T_dst from_matrix_cl(const matrix_cl<T>& src) {
  T dst;
  check_size_match("from_matrix_cl<scalar>", "src.rows()", src.rows(),
                   "dst.rows()", 1);
  check_size_match("from_matrix_cl<scalar>", "src.cols()", src.cols(),
                   "dst.cols()", 1);
  try {
    cl::Event copy_event;
    const cl::CommandQueue queue = opencl_context.queue();
    queue.enqueueReadBuffer(src.buffer(), opencl_context.in_order(), 0,
                            sizeof(T), &dst, &src.write_events(), &copy_event);
    copy_event.wait();
    src.clear_write_events();
  } catch (const cl::Error& e) {
    check_opencl_error("from_matrix_cl<scalar>", e);
  }
  return dst;
}

/** \ingroup opencl
 * Copies the source matrix that is stored on the OpenCL device to the
 * destination `std::vector`.
 *
 * @tparam T_dst destination type
 * @tparam T scalar in the source matrix
 * @param src source matrix on the OpenCL device
 * @return `std::vector` with a copy of the data in the source matrix
 */
template <typename T_dst, typename T,
          require_std_vector_vt<std::is_arithmetic, T_dst>* = nullptr,
          require_all_st_same<T_dst, T>* = nullptr>
inline T_dst from_matrix_cl(const matrix_cl<T>& src) {
  check_size_match("from_matrix_cl<std::vector>", "src.cols()", src.cols(),
                   "dst.cols()", 1);
  T_dst dst(src.rows());
  if (src.rows() == 0) {
    return dst;
  }
  try {
    cl::Event copy_event;
    const cl::CommandQueue queue = opencl_context.queue();
    queue.enqueueReadBuffer(src.buffer(), opencl_context.in_order(), 0,
                            sizeof(T) * src.rows(), dst.data(),
                            &src.write_events(), &copy_event);
    copy_event.wait();
    src.clear_write_events();
  } catch (const cl::Error& e) {
    check_opencl_error("from_matrix_cl<std::vector>", e);
  }
  return dst;
}

/** \ingroup opencl
 * Copies the source matrix that is stored on the OpenCL device to the
 * destination `std::vector` containing Eigen vectors.
 *
 * @tparam T_dst destination type
 * @tparam T scalar in the source matrix
 * @param src source matrix on the OpenCL device
 * @return `std::vector` containing Eigen vectors with a copy of the data in the
 * source matrix
 */
template <typename T_dst, typename T,
          require_std_vector_vt<is_eigen_vector, T_dst>* = nullptr,
          require_all_st_same<T_dst, T>* = nullptr>
inline T_dst from_matrix_cl(const matrix_cl<T>& src) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp = from_matrix_cl(src);
  T_dst dst;
  dst.reserve(src.cols());
  for (int i = 0; i < src.cols(); i++) {
    dst.emplace_back(tmp.col(i));
  }
  return dst;
}

/** \ingroup opencl
 * Copies the source kernel generator expression or matrix that is stored on the
 * OpenCL device to the destination Eigen matrix.
 *
 * @tparam T source type
 * @param src expression or source matrix on the OpenCL device
 * @return Eigen matrix with a copy of the data in the source matrix
 */
template <typename T, require_all_kernel_expressions_t<T>* = nullptr>
auto from_matrix_cl(const T& src) {
  return from_matrix_cl<
      Eigen::Matrix<scalar_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>>(src);
}

/** \ingroup opencl
 * Packs the square flat triangular matrix on the OpenCL device and
 * copies it to the std::vector.
 *
 * @param src the flat triangular source matrix on the OpenCL device
 * @return the packed std::vector
 * @throw <code>std::invalid_argument</code> if the matrix is not triangular
 */
template <typename T, require_matrix_cl_t<T>* = nullptr>
inline auto packed_copy(const T& src) {
  check_triangular("packed_copy", "src", src);
  const int packed_size = src.rows() * (src.rows() + 1) / 2;
  using T_val = value_type_t<T>;
  using T_not_bool
      = std::conditional_t<std::is_same<T_val, bool>::value, char, T_val>;
  std::vector<T_not_bool> dst(packed_size);
  if (dst.size() == 0) {
    return dst;
  }
  try {
    const cl::CommandQueue queue = opencl_context.queue();
    matrix_cl<T_val> packed(packed_size, 1);
    stan::math::opencl_kernels::pack(cl::NDRange(src.rows(), src.rows()),
                                     packed, src, src.rows(), src.rows(),
                                     src.view());
    const std::vector<cl::Event> mat_events
        = vec_concat(packed.read_write_events(), src.write_events());
    cl::Event copy_event;
    queue.enqueueReadBuffer(packed.buffer(), opencl_context.in_order(), 0,
                            sizeof(T_val) * packed_size, dst.data(),
                            &mat_events, &copy_event);
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
template <typename T, require_matrix_cl_t<T>* = nullptr>
inline plain_type_t<T> copy_cl(const T& src) {
  return plain_type_t<T>(src);
}

}  // namespace math
}  // namespace stan
#endif
#endif
