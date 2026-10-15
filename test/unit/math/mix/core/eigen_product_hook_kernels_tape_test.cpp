#include <stan/math/mix.hpp>
#include <gtest/gtest.h>

// Eigen's own kernels (triangular products and solves, Cholesky, LU) form
// coefficient-wise and scalar products internally, which Stan's product
// hooks evaluate. These bounds are the number of reverse-pass (chain) varis
// each kernel recorded before the hooks existed, for n = 8; the hooks must
// not increase them.

namespace {
using stan::math::fvar;
using stan::math::var;

template <typename T>
T kernel_scalar(double x) {
  return T(x);
}
template <>
fvar<var> kernel_scalar<fvar<var>>(double x) {
  return fvar<var>(var(x), var(0.5));
}

struct kernel_chain_counts {
  std::size_t tri_solve;
  std::size_t tri_mul;
  std::size_t llt;
  std::size_t lu;
};

template <typename T>
kernel_chain_counts measure_kernels() {
  const int n = 8;
  Eigen::Matrix<T, -1, -1> a(n, n);
  Eigen::Matrix<T, -1, 1> b(n);
  for (int j = 0; j < n; ++j) {
    b(j) = kernel_scalar<T>(0.1 * j + 0.3);
    for (int i = 0; i < n; ++i) {
      a(i, j) = kernel_scalar<T>((i == j ? n : 0) + 0.01 * (i + 1) * (j + 2));
    }
  }
  const auto& chain_stack = stan::math::ChainableStack::instance_->var_stack_;
  kernel_chain_counts counts;
  {
    const std::size_t start = chain_stack.size();
    Eigen::Matrix<T, -1, 1> x
        = a.template triangularView<Eigen::Lower>().solve(b);
    counts.tri_solve = chain_stack.size() - start;
  }
  {
    const std::size_t start = chain_stack.size();
    Eigen::Matrix<T, -1, 1> x = a.template triangularView<Eigen::Lower>() * b;
    counts.tri_mul = chain_stack.size() - start;
  }
  {
    const std::size_t start = chain_stack.size();
    Eigen::Matrix<T, -1, 1> x = a.llt().solve(b);
    counts.llt = chain_stack.size() - start;
  }
  {
    const std::size_t start = chain_stack.size();
    Eigen::Matrix<T, -1, 1> x = a.partialPivLu().solve(b);
    counts.lu = chain_stack.size() - start;
  }
  stan::math::recover_memory();
  return counts;
}
}  // namespace

TEST(mathMixEigenProductHook, varKernelChainCountsDoNotGrow) {
  const auto counts = measure_kernels<var>();
  EXPECT_LE(counts.tri_solve, 64);
  EXPECT_LE(counts.tri_mul, 82);
  EXPECT_LE(counts.llt, 340);
  EXPECT_LE(counts.lu, 484);
}

TEST(mathMixEigenProductHook, fvarVarKernelChainCountsDoNotGrow) {
  const auto counts = measure_kernels<fvar<var>>();
  EXPECT_LE(counts.tri_solve, 216);
  EXPECT_LE(counts.tri_mul, 256);
  EXPECT_LE(counts.llt, 1414);
  EXPECT_LE(counts.lu, 1504);
}
