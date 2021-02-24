#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto append_array_functor = [](const auto& a, const auto& b) {
  return stan::math::append_array(a, b);
};

TEST(OpenCLAppendArray, prim_rev_values_small) {
  std::vector<double> a{-2.2, -0.8, 0.5, 1, 1.5, 3, 3.4, 4};
  std::vector<double> b{1, 2, 3, 4, 5};
  stan::math::test::compare_cpu_opencl_prim_rev(append_array_functor, a, b);
}

TEST(OpenCLAppendArray, prim_rev_size_0) {
  int N = 0;

  std::vector<double> a{};
  std::vector<double> b{};
  stan::math::test::compare_cpu_opencl_prim_rev(append_array_functor, a, b);
}

TEST(OpenCLAppendArray, prim_rev_values_large) {
  int N = 71;
  int M = 87;

  std::vector<double> a;
  std::vector<double> b;
  for (int i = 0; i < N; i++) {
    a.push_back(Eigen::VectorXd::Random(1)[0]);
  }
  for (int i = 0; i < M; i++) {
    b.push_back(Eigen::VectorXd::Random(1)[0]);
  }
  stan::math::test::compare_cpu_opencl_prim_rev(append_array_functor, a, b);
}

#endif
