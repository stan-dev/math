#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/gpu/opencl_context.hpp>
#include <stan/math/gpu/matrix_gpu.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

TEST(MathMatrixGPU, matrix_gpu_copy) {
  stan::math::vector_d d1;
  stan::math::vector_d d1_a;
  stan::math::vector_d d1_b;  
  stan::math::matrix_d d2;
  stan::math::matrix_d d2_a;
  stan::math::matrix_d d2_b;
  stan::math::matrix_d d0;
  d1.resize(3);
  d1_a.resize(3);
  d1_b.resize(3);
  d1 << 1, 2, 3;
  d2.resize(2, 3);
  d2_a.resize(2, 3);
  d2_b.resize(2, 3);
  d2 << 1, 2, 3, 4, 5, 6;
  // vector
  stan::math::matrix_gpu d11(3, 1);
  stan::math::matrix_gpu d111(3, 1);
  EXPECT_NO_THROW(stan::math::copy(d11, d1));
  EXPECT_NO_THROW(stan::math::copy(d111, d11));
  EXPECT_NO_THROW(stan::math::copy(d1_a, d11));
  EXPECT_NO_THROW(stan::math::copy(d1_b, d111));  
  EXPECT_EQ(1, d1_a(0));
  EXPECT_EQ(2, d1_a(1));
  EXPECT_EQ(3, d1_a(2));
  EXPECT_EQ(1, d1_b(0));
  EXPECT_EQ(2, d1_b(1));
  EXPECT_EQ(3, d1_b(2));
  // matrix
  stan::math::matrix_gpu d00;
  stan::math::matrix_gpu d000;
  stan::math::matrix_gpu d22(2, 3);
  stan::math::matrix_gpu d222(2, 3);
  EXPECT_NO_THROW(stan::math::copy(d22, d2));
  EXPECT_NO_THROW(stan::math::copy(d222, d22));
  EXPECT_NO_THROW(stan::math::copy(d2_a, d22));
  EXPECT_NO_THROW(stan::math::copy(d2_b, d222));  
  EXPECT_EQ(1, d2_a(0, 0));
  EXPECT_EQ(2, d2_a(0, 1));
  EXPECT_EQ(3, d2_a(0, 2));
  EXPECT_EQ(4, d2_a(1, 0));
  EXPECT_EQ(5, d2_a(1, 1));
  EXPECT_EQ(6, d2_a(1, 2));
  EXPECT_EQ(1, d2_b(0, 0));
  EXPECT_EQ(2, d2_b(0, 1));
  EXPECT_EQ(3, d2_b(0, 2));
  EXPECT_EQ(4, d2_b(1, 0));
  EXPECT_EQ(5, d2_b(1, 1));
  EXPECT_EQ(6, d2_b(1, 2));
  // zero sized copy
  EXPECT_NO_THROW(stan::math::copy(d00, d0));
  EXPECT_NO_THROW(stan::math::copy(d0, d00));
  EXPECT_NO_THROW(stan::math::copy(d000, d00));
}

TEST(MathMatrixGPU, barebone_buffer_copy) {
  // a barebone OpenCL example of copying
  // a vector of doubles to the GPU and back
  size_t size = 512;
  std::vector<double> cpu_buffer(size);
  for (unsigned int i = 0; i < size; i++) {
    cpu_buffer[i] = i * 1.0;
  }
  std::vector<double> cpu_dst_buffer(size);
  // retrieve the command queue
  cl::CommandQueue queue = stan::math::opencl_context.queue();
  // retrieve the context
  cl::Context& ctx = stan::math::opencl_context.context();
  // create the gpu buffer of the same size
  cl::Buffer gpu_buffer
      = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * size);

  // write the cpu_buffer to the GPU (gpu_buffer)
  queue.enqueueWriteBuffer(gpu_buffer, CL_TRUE, 0, sizeof(double) * size,
                           &cpu_buffer[0]);

  // write the gpu buffer back to the cpu_dst_buffer
  queue.enqueueReadBuffer(gpu_buffer, CL_TRUE, 0, sizeof(double) * size,
                          &cpu_dst_buffer[0]);

  for (unsigned int i = 0; i < size; i++) {
    EXPECT_EQ(i * 1.0, cpu_dst_buffer[i]);
  }
}

#endif
