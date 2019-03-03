#ifndef STAN_MATH_OPENCL_MULTIPLY_TRANSPOSE_HPP
#define STAN_MATH_OPENCL_MULTIPLY_TRANSPOSE_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/multiply_transpose.hpp>
#include <stan/math/opencl/err/check_square.hpp>
#include <Eigen/Dense>

namespace stan {
namespace math {
/**
 * Computes the product of a square OpenCL matrix with its transpose.
 *
 * Computes the matrix multiplication C = A x A^T
 *
 * @param A input matrix
 * @return the product of the input matrix and its transpose
 *
 */
inline matrix_cl multiply_transpose(const matrix_cl& A) {
  matrix_cl temp(A.rows(), A.rows());
  if (A.size() == 0)
    return temp;
  // padding the matrices so the dimensions are divisible with local
  // improves performance becasuse we can omit if statements in the
  // multiply kernel
  int local = opencl_kernels::multiply_transpose.make_functor.get_opts().at(
      "THREAD_BLOCK_SIZE");
  int Mpad = ((A.rows() + local - 1) / local) * local;
  int Npad = ((A.cols() + local - 1) / local) * local;
  matrix_cl tempPad(Mpad, Mpad);
  matrix_cl Apad(Mpad, Npad);
  opencl_kernels::zeros(cl::NDRange(Mpad, Npad), Apad.buffer(), Mpad, Npad,
                        TriangularViewCL::Entire);
  Apad.sub_block(A, 0, 0, 0, 0, A.rows(), A.cols());
  int wpt = opencl_kernels::multiply_transpose.make_functor.get_opts().at(
      "WORK_PER_THREAD");
  try {
    opencl_kernels::multiply_transpose(
        cl::NDRange(Mpad, Mpad / wpt), cl::NDRange(local, local / wpt),
        Apad.buffer(), tempPad.buffer(), Apad.rows(), Apad.cols());
  } catch (cl::Error& e) {
    check_opencl_error("multiply self transpose", e);
  }
  temp.sub_block(tempPad, 0, 0, 0, 0, temp.rows(), temp.cols());
  return temp;
}
}  // namespace math
}  // namespace stan

#endif
#endif
