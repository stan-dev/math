#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <gtest/gtest.h>
#include <CL/cl.hpp>
#include <algorithm>
#include <vector>

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
  cl::Event write_event;
  queue.enqueueWriteBuffer(gpu_buffer, CL_TRUE, 0, sizeof(double) * size,
                           &cpu_buffer[0], NULL, &write_event);
  write_event.wait();
  cl::Event read_event;
  // write the gpu buffer back to the cpu_dst_buffer
  queue.enqueueReadBuffer(gpu_buffer, CL_TRUE, 0, sizeof(double) * size,
                          &cpu_dst_buffer[0], NULL, &read_event);
  read_event.wait();
  for (unsigned int i = 0; i < size; i++) {
    EXPECT_EQ(i * 1.0, cpu_dst_buffer[i]);
  }
}

TEST(MathMatrixGPU, matrix_cl_vector_copy) {
  stan::math::vector_d d1_cpu;
  stan::math::vector_d d1_a_cpu;
  stan::math::vector_d d1_b_cpu;
  d1_cpu.resize(3);
  d1_a_cpu.resize(3);
  d1_b_cpu.resize(3);
  d1_cpu << 1, 2, 3;
  // vector
  stan::math::matrix_cl d11_cl(3, 1);
  stan::math::matrix_cl d111_cl(3, 1);
  EXPECT_NO_THROW(stan::math::copy(d11_cl, d1_cpu));
  EXPECT_NO_THROW(stan::math::copy(d111_cl, d11_cl));
  EXPECT_NO_THROW(stan::math::copy(d1_a_cpu, d11_cl));
  EXPECT_NO_THROW(stan::math::copy(d1_b_cpu, d111_cl));
  EXPECT_EQ(1, d1_a_cpu(0));
  EXPECT_EQ(2, d1_a_cpu(1));
  EXPECT_EQ(3, d1_a_cpu(2));
  EXPECT_EQ(1, d1_b_cpu(0));
  EXPECT_EQ(2, d1_b_cpu(1));
  EXPECT_EQ(3, d1_b_cpu(2));
}

TEST(MathMatrixCL, matrix_cl_matrix_copy) {
  stan::math::matrix_d d2_cpu;
  stan::math::matrix_d d2_a_cpu;
  stan::math::matrix_d d2_b_cpu;
  stan::math::matrix_d d0_cpu;
  d2_cpu.resize(2, 3);
  d2_a_cpu.resize(2, 3);
  d2_b_cpu.resize(2, 3);
  d2_cpu << 1, 2, 3, 4, 5, 6;
  // matrix
  stan::math::matrix_cl d00_cl;
  stan::math::matrix_cl d000_cl;
  stan::math::matrix_cl d22_cl(2, 3);
  stan::math::matrix_cl d222_cl(2, 3);
  EXPECT_NO_THROW(stan::math::copy(d22_cl, d2_cpu));
  EXPECT_NO_THROW(stan::math::copy(d222_cl, d22_cl));
  EXPECT_NO_THROW(stan::math::copy(d2_a_cpu, d22_cl));
  EXPECT_NO_THROW(stan::math::copy(d2_b_cpu, d222_cl));
  EXPECT_EQ(1, d2_a_cpu(0, 0));
  EXPECT_EQ(2, d2_a_cpu(0, 1));
  EXPECT_EQ(3, d2_a_cpu(0, 2));
  EXPECT_EQ(4, d2_a_cpu(1, 0));
  EXPECT_EQ(5, d2_a_cpu(1, 1));
  EXPECT_EQ(6, d2_a_cpu(1, 2));
  EXPECT_EQ(1, d2_b_cpu(0, 0));
  EXPECT_EQ(2, d2_b_cpu(0, 1));
  EXPECT_EQ(3, d2_b_cpu(0, 2));
  EXPECT_EQ(4, d2_b_cpu(1, 0));
  EXPECT_EQ(5, d2_b_cpu(1, 1));
  EXPECT_EQ(6, d2_b_cpu(1, 2));
  // zero sized copy
  EXPECT_NO_THROW(stan::math::copy(d00_cl, d0_cpu));
  EXPECT_NO_THROW(stan::math::copy(d0_cpu, d00_cl));
  EXPECT_NO_THROW(stan::math::copy(d000_cl, d00_cl));
}

TEST(MathMatrixCL, matrix_cl_matrix_copy_arithmetic) {
  int test_val = 5;
  // Use this for successful copy
  stan::math::matrix_cl d22_cl(1, 1);
  // Use this for failure copy
  stan::math::matrix_cl d222_cl(2, 3);
  EXPECT_NO_THROW(stan::math::copy(d22_cl, test_val));
  EXPECT_NO_THROW(stan::math::copy(test_val, d22_cl));
  EXPECT_THROW(copy(d222_cl, test_val), std::invalid_argument);
  EXPECT_THROW(copy(test_val, d222_cl), std::invalid_argument);
}

#endif
