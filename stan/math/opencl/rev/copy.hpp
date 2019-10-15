#ifndef STAN_MATH_OPENCL_REV_COPY_HPP
#define STAN_MATH_OPENCL_REV_COPY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/rev/matrix_cl.hpp>
#include <stan/math/opencl/kernels/copy.hpp>
#include <stan/math/opencl/kernels/pack.hpp>
#include <stan/math/opencl/kernels/unpack.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/arr/fun/vec_concat.hpp>

#include <CL/cl.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <type_traits>

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
  template <int R, int C>
  inline matrix_cl<var> to_matrix_cl(const Eigen::Matrix<var, R, C>& src) {
    return {src};
  }


  /**
   * Packs the flat triagnular matrix on the OpenCL device and
   * copies it to the std::vector.
   *
   * @tparam triangular_view the triangularity of the source matrix
   * @param src the flat triangular source matrix on the OpenCL device
   * @return the packed std::vector
   */
  template <matrix_cl_view triangular_view>
  inline std::vector<double> packed_copy(const matrix_cl<var>& src) {
    const int packed_size = src.rows() * (src.rows() + 1) / 2;
    std::vector<double> dst(packed_size);
    if (dst.size() == 0) {
      return dst;
    }
    try {
      const cl::CommandQueue queue = opencl_context.queue();
      matrix_cl<double> packed(packed_size, 1, src.view());
      stan::math::opencl_kernels::pack(cl::NDRange(src.rows(), src.rows()),
                                       packed, src.val(), src.rows(), src.rows(),
                                       triangular_view);
      const std::vector<cl::Event> mat_events
          = vec_concat(packed.read_write_events(), src.adj().write_events());
      cl::Event copy_event;
      queue.enqueueReadBuffer(packed.buffer(), CL_FALSE, 0,
                              sizeof(double) * packed_size, dst.data(),
                              &mat_events, &copy_event);
      copy_event.wait();
      src.adj().clear_write_events();
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
* @tparam triangular_view the triangularity of the source matrix
* @param src the packed source std::vector
* @param rows the number of rows in the flat matrix
* @return the destination flat matrix on the OpenCL device
* @throw <code>std::invalid_argument</code> if the
* size of the vector does not match the expected size
* for the packed triangular matrix
*/
template <matrix_cl_view triangular_view>
inline matrix_cl<var> packed_copy(vari** src, int rows) try {
 const int packed_size = rows * (rows + 1) / 2;
   // TODO(Steve): idk how to check size match here
   //check_size_match("copy (packed std::vector -> OpenCL)", "src.size()",
   //                 rows, "rows * (rows + 1) / 2", packed_size);
   matrix_cl<var> dst(rows, rows);
   if (dst.size() == 0) {
     return dst;
   }
   matrix_cl<var> packed(src, rows, 1, triangular_view);
   cl::Event packed_event;
   stan::math::opencl_kernels::unpack(cl::NDRange(dst.rows(), dst.rows()), dst.val(),
                                      packed.val(), dst.rows(), dst.rows(),
                                      triangular_view);
  stan::math::opencl_kernels::unpack(cl::NDRange(dst.rows(), dst.rows()), dst.adj(),
                                     packed.adj(), dst.rows(), dst.rows(),
                                     triangular_view);
   return dst;
 } catch (const cl::Error& e) {
   check_opencl_error("packed_copy (std::vector->OpenCL)", e);
   matrix_cl<var> dst(rows, rows);
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
  inline matrix_cl<var> copy_cl(const matrix_cl<var>& src) {
    matrix_cl<var> dst(src);
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
inline Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> from_matrix_cl(
    const matrix_cl<var>& src) try {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dst(src.rows(),
                                                            src.cols());
  if (src.size() == 0) {
   return dst;
  }
  dst = from_matrix_cl(src.adj());
  return dst;
} catch (const cl::Error& e) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dst(src.rows(),
                                                            src.cols());
  check_opencl_error("copy (OpenCL)->Eigen", e);
  return dst;
}

}
}

#endif
#endif
