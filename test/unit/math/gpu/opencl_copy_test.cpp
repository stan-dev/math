#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/gpu/opencl_context.hpp>
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/rev/mat.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/gpu/copy.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

void test_cache_speed() {
  auto m = stan::math::matrix_d::Random(100, 100).eval();
  stan::math::matrix_gpu d11(100, 100);
    stan::math::matrix_gpu d12(100, 100);
  std::chrono::steady_clock::time_point first_begin = std::chrono::steady_clock::now();
  stan::math::copy(d11, m);
  std::chrono::steady_clock::time_point first_end = std::chrono::steady_clock::now();
  size_t first_pass = std::chrono::duration_cast<std::chrono::nanoseconds>(first_end - first_begin).count();
  std::chrono::steady_clock::time_point second_begin = std::chrono::steady_clock::now();
  stan::math::copy(d12, m);
  std::chrono::steady_clock::time_point second_end = std::chrono::steady_clock::now();
  size_t second_pass = std::chrono::duration_cast<std::chrono::nanoseconds>(second_end - second_begin).count();
  ASSERT_GT(first_pass, second_pass);
  ASSERT_FALSE(m.opencl_buffer_() == NULL);
}

TEST(MathMatrixGPU, matrix_gpu_copy_cache) {
  test_cache_speed();
}

TEST(MathMatrixGPU, matrix_gpu_var_copy) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> m(5, 5);
  double pos_ = 1.1;
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      m(i, j) = pos_++;

  stan::math::matrix_gpu d1_gpu(5, 5);
  stan::math::copy(d1_gpu, m);
  EXPECT_TRUE(m.opencl_buffer_() == NULL);
  stan::math::matrix_d d1_cpu_return(5, 5);
  stan::math::copy(d1_cpu_return, d1_gpu);
  EXPECT_EQ(1.1, d1_cpu_return(0, 0));
  EXPECT_EQ(6.1, d1_cpu_return(1, 0));
  EXPECT_EQ(11.1, d1_cpu_return(2, 0));
  EXPECT_EQ(16.1, d1_cpu_return(3, 0));
  EXPECT_EQ(21.1, d1_cpu_return(4, 0));
}

TEST(MathMatrixGPU, matrix_gpu_copy) {
  stan::math::vector_d d1_cpu;
  stan::math::vector_d d1_a_cpu;
  stan::math::vector_d d1_b_cpu;
  stan::math::matrix_d d2_cpu;
  stan::math::matrix_d d2_a_cpu;
  stan::math::matrix_d d2_b_cpu;
  stan::math::matrix_d d0_cpu;
  d1_cpu.resize(3);
  d1_a_cpu.resize(3);
  d1_b_cpu.resize(3);
  d1_cpu << 1, 2, 3;
  d2_cpu.resize(2, 3);
  d2_a_cpu.resize(2, 3);
  d2_b_cpu.resize(2, 3);
  d2_cpu << 1, 2, 3, 4, 5, 6;
  // vector
  stan::math::matrix_gpu d1_gpu(3, 1);
  stan::math::matrix_gpu d11_gpu(3, 1);
  EXPECT_NO_THROW(stan::math::copy(d1_gpu, d1_cpu));
  EXPECT_NO_THROW(stan::math::copy(d11_gpu, d1_gpu));
  EXPECT_NO_THROW(stan::math::copy(d1_a_cpu, d1_gpu));
  EXPECT_EQ(1, d1_a_cpu(0));
  EXPECT_EQ(2, d1_a_cpu(1));
  EXPECT_EQ(3, d1_a_cpu(2));
  EXPECT_NO_THROW(stan::math::copy(d1_b_cpu, d11_gpu));
  EXPECT_EQ(1, d1_b_cpu(0));
  EXPECT_EQ(2, d1_b_cpu(1));
  EXPECT_EQ(3, d1_b_cpu(2));
  // matrix
  stan::math::matrix_gpu d00_gpu;
  stan::math::matrix_gpu d000_gpu;
  stan::math::matrix_gpu d22_gpu(2, 3);
  stan::math::matrix_gpu d222_gpu(2, 3);
  EXPECT_NO_THROW(stan::math::copy(d22_gpu, d2_cpu));
  EXPECT_NO_THROW(stan::math::copy(d222_gpu, d22_gpu));
  EXPECT_NO_THROW(stan::math::copy(d2_a_cpu, d22_gpu));
  EXPECT_NO_THROW(stan::math::copy(d2_b_cpu, d222_gpu));
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
  EXPECT_NO_THROW(stan::math::copy(d00_gpu, d0_cpu));
  EXPECT_NO_THROW(stan::math::copy(d0_cpu, d00_gpu));
  EXPECT_NO_THROW(stan::math::copy(d000_gpu, d00_gpu));
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
