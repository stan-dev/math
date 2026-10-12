#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

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
 * the denominator at 1. References by mpmath at 60 digits.
 */
TEST(OpenCL_Phi_approx, rev_tail_derivative_against_references) {
  const int N = 8;
  Eigen::VectorXd x(N);
  x << -40, -10, -7, 1, 7, 8, 10, 40;
  // at x = -40 and 40 the true value, about 4e-1987, underflows to 0
  Eigen::VectorXd expected(N);
  expected << 0.0, 5.9589784505245725e-37, 5.1340418352583332e-15,
      0.24152728992165639, 5.1340418352583332e-15, 8.7097694481792379e-21,
      5.9589784505245725e-37, 0.0;

  stan::math::var_value<stan::math::matrix_cl<double>> x_cl(
      stan::math::to_matrix_cl(x));
  stan::math::sum(stan::math::Phi_approx(x_cl)).grad();

  const Eigen::VectorXd adj = stan::math::from_matrix_cl(x_cl.adj());
  for (int i = 0; i < N; ++i) {
    EXPECT_NEAR(
        expected[i], adj[i],
        1e-13 * std::fmax(expected[i], std::numeric_limits<double>::min()))
        << "x = " << x[i];
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
