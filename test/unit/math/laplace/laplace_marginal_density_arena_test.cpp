#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
#include <vector>

namespace laplace_marginal_density_arena_test {

// Highest arena use reached at the likelihood calls of one gradient of a
// Poisson GP with the given solver, with every allocation in one arena block
size_t peak_arena_bytes(int solver, int n) {
  auto& memalloc = stan::math::ChainableStack::instance_->memalloc_;
  memalloc.alloc(16 << 20);
  stan::math::recover_memory();
  memalloc.alloc(1 << 20);  // skip the small first block
  size_t peak = 0;
  auto ll_fun = [&](auto&& theta, auto&& y, std::ostream*) {
    peak = std::max(peak, memalloc.approx_bytes_used());
    return stan::math::poisson_log_lpmf(y, theta);
  };
  std::vector<double> x(n);
  std::iota(x.begin(), x.end(), 0.0);
  std::vector<int> y(n, 2);
  Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(n);
  stan::math::var alpha = 1.2, rho = 0.7;
  stan::math::laplace_marginal_tol(
      ll_fun, std::forward_as_tuple(y), 1,
      stan::math::test::sqr_exp_kernel_functor{},
      std::forward_as_tuple(x, alpha, rho),
      std::make_tuple(theta_0, 1e-8, 200, solver, 0, false), nullptr);
  stan::math::recover_memory();
  return peak;
}

}  // namespace laplace_marginal_density_arena_test

TEST(laplace_marginal_density, solvers_2_and_3_allocate_no_unused_buffers) {
  using laplace_marginal_density_arena_test::peak_arena_bytes;
  const int n = 50;
  const size_t matrix = 8 * n * n;  // bytes of one n x n matrix
  size_t peak_1 = peak_arena_bytes(1, n);
  // Beyond solver 1's n x n arena matrices, solver 2 needs none and solver 3
  // one (A); a pre-sized buffer that is then reassigned adds a whole one.
  EXPECT_LT(peak_arena_bytes(2, n), peak_1 + matrix / 2);
  EXPECT_LT(peak_arena_bytes(3, n), peak_1 + matrix * 3 / 2);
}
