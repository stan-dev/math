#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/gpu/multiply.hpp>
#include <stan/math/gpu/copy.hpp>
#include <stan/math/gpu/lower_tri_inverse.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <algorithm>

void populate_lower_triangular_inverse_input(stan::math::matrix_d& A) {
  boost::random::mt19937 rng;
  for (int i = 0; i < A.rows(); i++) {
    for (int j = 0; j < i; j++) {
      A(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    A(i, i) = 10000.0;
    for (int j = i + 1; j < A.rows(); j++) {
      A(i, j) = 0.0;
    }
  }
}

bool are_timings_steady(std::vector<size_t>& t) {
  if (t.size() < 20) {
    return false;
  }
  size_t temp_old = 0;
  size_t temp = 0;
  for (size_t i = 0; i < (t.size() - 1); i++) {
    temp_old += t[i];
  }
  temp = temp_old + t[t.size() - 1];
  temp_old /= (t.size() - 1);
  temp /= t.size();

  double diff = abs(static_cast<double>(temp_old) - static_cast<double>(temp));
  double diff_pct = diff / static_cast<double>(temp_old);
  if (diff_pct > 0.001) {
    return false;
  } else {
    std::cout << temp << std::endl;
    return true;
  }
}

TEST(MathMatrixGPU, lower_tri_inverse) {
  std::vector<int> sizes
      = {32,   48,   64,   80,   96,   112,  128,  160,  192,  224,  256,
         288,  320,  352,  384,  416,  448,  480,  512,  576,  640,  704,
         768,  832,  896,  960,  1024, 1152, 1280, 1408, 1536, 1664, 1798,
         1920, 2048, 2304, 2560, 2816, 3072, 3328, 3584, 3840, 4096, 4608,
         5120, 5632, 6144, 6656, 7168, 7680, 8192, 9216, 10240};
  for (auto s : sizes) {
    std::cout << s << "\t";
    auto m1 = stan::math::matrix_d(s, s);
    auto m5 = stan::math::matrix_d(s, s);
    populate_lower_triangular_inverse_input(m1);
    std::vector<size_t> temp_cpu;
    while (!are_timings_steady(temp_cpu)) {
      std::chrono::steady_clock::time_point begin
          = std::chrono::steady_clock::now();
      auto m1_cpu = stan::math::mdivide_left_tri<Eigen::Lower>(m1).eval();
      std::chrono::steady_clock::time_point end
          = std::chrono::steady_clock::now();
      size_t current
          = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin)
                .count();
      temp_cpu.push_back(current);
    }
  }
}
#endif