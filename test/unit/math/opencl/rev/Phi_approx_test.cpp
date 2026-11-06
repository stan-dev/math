#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/math/mix/fun/Phi_approx_tail_refs.hpp>
#include <vector>

auto Phi_approx_functor
    = [](const auto& a) { return stan::math::Phi_approx(a); };

TEST(OpenCL_Phi_approx, prim_rev_values_small) {
  Eigen::VectorXd a(8);
  a << -0.22, -0.8, 0.5, 1, 0.15, 0.3, 0.34, -0.04;
  stan::math::test::compare_cpu_opencl_prim_rev(Phi_approx_functor, a);
}

/**
 * Device derivative in the tails against fixed references, which
 * `compare_cpu_opencl_prim_rev` cannot check: its relative tolerance floors
 * the denominator at 1. References in Phi_approx_tail_refs.hpp.
 */
TEST(OpenCL_Phi_approx, rev_tail_derivative_against_references) {
  const Eigen::VectorXd x = stan::math::to_vector(Phi_approx_tail_refs::args());
  const std::vector<double> expected = Phi_approx_tail_refs::derivs();

  stan::math::var_value<stan::math::matrix_cl<double>> x_cl(
      stan::math::to_matrix_cl(x));
  stan::math::sum(stan::math::Phi_approx(x_cl)).grad();

  const Eigen::VectorXd adj = stan::math::from_matrix_cl(x_cl.adj());
  for (int i = 0; i < x.size(); ++i) {
    Phi_approx_tail_refs::expect_rel(expected[i], adj[i], x[i]);
  }
  stan::math::recover_memory();
}

TEST(OpenCL_Phi_approx, prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(Phi_approx_functor, a);
}

TEST(OpenCL_Phi_approx, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(Phi_approx_functor, a);
}

#endif
