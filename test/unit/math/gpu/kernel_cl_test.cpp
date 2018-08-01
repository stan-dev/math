#ifdef STAN_OPENCL

#include <stan/math/prim/mat.hpp>
#include <CL/cl.hpp>
#include <stan/math/gpu/opencl_context.hpp>
#include <stan/math/gpu/kernel_cl.hpp>
#include <stan/math/gpu/copy.hpp>
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/gpu/transpose.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathGpu, make_kernel) {
  stan::math::vector_d v0;
  stan::math::row_vector_d rv0;
  stan::math::matrix_d m0;
  v0 << 1, 2, 3;
  rv0 << 1, 2, 3;
  m0 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::math::matrix_gpu v00(v0);
  stan::math::matrix_gpu rv00(rv0);
  stan::math::matrix_gpu m00(m0);

  stan::math::matrix_gpu v00_dst(v0.cols(), v0.rows());
  stan::math::matrix_gpu rv00_dst(rv0.cols(), rv0.rows());
  stan::math::matrix_gpu m00_dst(m0.cols(), m0.rows());
  cl::CommandQueue cmdQueue = stan::math::opencl_context.queue();
  stan::math::kernel_cl kernel("tranpose");
//  kernel.set_args(v00_dst.buffer(), v00.buffer(), v00.rows(), v00.cols());
//  cmdQueue.enqueueNDRangeKernel(kernel.compiled_, cl::NullRange,
//                                cl::NDRange(v00.rows(), v00.cols()),
//                                cl::NullRange, NULL, NULL);
}
#endif
