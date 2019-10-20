#ifndef STAN_MATH_OPENCL_MULTIPLY_HPP
#define STAN_MATH_OPENCL_MULTIPLY_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/kernels/scalar_mul.hpp>
#include <stan/math/opencl/kernels/matrix_multiply.hpp>
#include <stan/math/opencl/kernels/add.hpp>
#include <stan/math/opencl/sub_block.hpp>
#include <stan/math/opencl/zeros.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <algorithm>

namespace stan {
namespace math {
namespace opencl {
/**
 * Computes the product of the specified matrices with the option
 * of specifying the triangularity of either input matrices.
 *
 * Computes the matrix multiplication C[M, K] = A[M, N] x B[N, K]
 *
 * @param A first matrix
 * @param B second matrix
 * @tparam partial_view_A specifies whether the matrix A is a
 *  lower/upper triangular or a rectangular matrix
 * @tparam partial_view_B specifies whether the matrix B is a
 *  lower/upper triangular or a rectangular matrix
 * @return the product of the first and second matrix
 *
 * @throw <code>std::invalid_argument</code> if the
 *   number of columns in A and rows in B do not match
 */

template <typename T1, typename T2, typename = require_all_arithmetic_t<T1, T2>>
inline matrix_cl<return_type_t<T1, T2>> multiply(const matrix_cl<T1>& A,
                                                 const matrix_cl<T2>& B) {
  check_size_match("multiply ((OpenCL))", "A.cols()", A.cols(), "B.rows()",
                   B.rows());
  matrix_cl<return_type_t<T1, T2>> temp(A.rows(), B.cols(),
                                        either(A.view(), B.view()));
  if (A.size() == 0 || B.size() == 0) {
    temp.zeros();
    return temp;
  }
  if (A.rows() == 1) {
    const int local_size
        = opencl_kernels::row_vector_matrix_multiply.make_functor.get_opts().at(
            "LOCAL_SIZE_");
    try {
      opencl_kernels::row_vector_matrix_multiply(
          cl::NDRange(temp.cols() * local_size), cl::NDRange(local_size), A, B,
          temp, B.rows(), B.cols(), A.view(), B.view());
    } catch (cl::Error& e) {
      check_opencl_error("row_vector - matrix multiply", e);
    }
    return temp;
  }
  if (B.cols() == 1) {
    try {
      opencl_kernels::matrix_vector_multiply(cl::NDRange(temp.rows()), A, B,
                                             temp, A.rows(), A.cols(), A.view(),
                                             B.view());
    } catch (cl::Error& e) {
      check_opencl_error("matrix - vector multiply", e);
    }
    return temp;
  }
  int local = opencl_kernels::matrix_multiply.make_functor.get_opts().at(
      "THREAD_BLOCK_SIZE");
  const int Mpad = ((A.rows() + local - 1) / local) * local;
  const int Npad = ((B.cols() + local - 1) / local) * local;
  const int wpt = opencl_kernels::matrix_multiply.make_functor.get_opts().at(
      "WORK_PER_THREAD");
  const int wgs = Mpad / local * Npad / local;
  const int split = std::min(
      A.cols() / local,
      (opencl_context.tuning_opts().multiply_wgs_per_compute_unit
           * static_cast<int>(opencl_context.device()[0]
                                  .getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>())
       + wgs - 1)
          / wgs);
  try {
    if (split <= 1) {
      opencl_kernels::matrix_multiply(
          cl::NDRange(Mpad, Npad / wpt), cl::NDRange(local, local / wpt), A, B,
          temp, A.rows(), B.cols(), B.rows(), A.view(), B.view());
    } else {
      matrix_cl<return_type_t<T1, T2>> tempSplit(A.rows(), B.cols() * split);
      opencl_kernels::matrix_multiply(cl::NDRange(Mpad, Npad / wpt, split),
                                      cl::NDRange(local, local / wpt, 1), A, B,
                                      tempSplit, A.rows(), B.cols(), B.rows(),
                                      A.view(), B.view());
      opencl_kernels::add_batch(cl::NDRange(A.rows(), B.cols()), temp,
                                tempSplit, A.rows(), B.cols(), split);
    }
  } catch (cl::Error& e) {
    check_opencl_error("multiply", e);
  }
  return temp;
}
}  // namespace opencl
}  // namespace math
}  // namespace stan
#endif
#endif
