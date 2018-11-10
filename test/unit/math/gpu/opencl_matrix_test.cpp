#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/rev/mat.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/gpu/opencl_context.hpp>
#include <stan/math/gpu/matrix_gpu.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

void testy_mctest() {
  double pos = 1.1;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> m(5000, 5000);
  for (int i = 0; i < 5000; ++i)
    for (int j = 0; j < 5000; ++j)
      m(i, j) = stan::math::var(pos++);
  std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();
  stan::math::matrix_gpu d33(m);
  std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
  size_t current1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count();
    std::cout << "First Copy Time: " << current1 << "\n";
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    stan::math::matrix_gpu d33b(m);
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  size_t current = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
  std::cout << "Second Copy Time:" << current << "\n";
}

TEST(MathMatrixGPU, matrix_gpu_creation) {
  std::cout << "Var test:  \n";

  for (int i = 0; i <= 10; i++) {
    testy_mctest();
  }
  std::cout << "Normie test:\n";
  for (int i = 0; i <= 10; i++) {
    auto m = stan::math::matrix_d::Random(5000, 5000).eval();
    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();
    stan::math::matrix_gpu d33(m);
    std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    size_t current1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count();
      std::cout << "First Copy Time: " << current1 << "\n";
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
      stan::math::matrix_gpu d33b(m);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    size_t current = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
    std::cout << "Second Copy Time:  " << current << "\n";
  }
}

#endif
