#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/rev/mat.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/gpu/opencl_context.hpp>
#include <stan/math/gpu/matrix_gpu.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

void test_cache_speed() {
  auto m = stan::math::matrix_d::Random(100, 100).eval();
  std::chrono::steady_clock::time_point begin1
      = std::chrono::steady_clock::now();
  stan::math::matrix_gpu d33(m);
  std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
  size_t first_pass
      = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1)
            .count();
  std::chrono::steady_clock::time_point begin
      = std::chrono::steady_clock::now();
  stan::math::matrix_gpu d33b(m);
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  size_t second_pass
      = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin)
            .count();
  ASSERT_GT(first_pass, second_pass);
  ASSERT_FALSE(m.opencl_buffer_() == NULL);
}

TEST(MathMatrixGPU, matrix_gpu_cache) {
  for (int i = 0; i <= 10; i++) {
    test_cache_speed();
  }
}

TEST(MathMatrixGPU, matrix_gpu_primitive_creation) {
  stan::math::vector_d d1;
  stan::math::matrix_d d2;
  stan::math::matrix_d d3;

  d1.resize(3);
  d2.resize(2, 3);
  EXPECT_NO_THROW(stan::math::matrix_gpu A(1, 1));
  EXPECT_NO_THROW(stan::math::matrix_gpu d11(d1));
  EXPECT_NO_THROW(stan::math::matrix_gpu d22(d2));
  stan::math::matrix_gpu d33(d3);
}

TEST(MathMatrixGPU, matrix_gpu_var_creation) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> m(5, 5);
  double pos_ = 1.1;
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      m(i, j) = pos_++;

  EXPECT_NO_THROW(stan::math::matrix_gpu d11(m));
  ASSERT_TRUE(m.opencl_buffer_() == NULL);
}

#endif
