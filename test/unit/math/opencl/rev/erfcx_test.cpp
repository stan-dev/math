#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

auto erfcx_functor = [](const auto& a) { return stan::math::erfcx(a); };

TEST(OpenCLerfcx, prim_rev_values_small) {
  Eigen::VectorXd a(7);
  a << -2.6, -2, -1, -0.2, 1, 1.3, 2.6;
  stan::math::test::compare_cpu_opencl_prim_rev(erfcx_functor, a);
}

/**
 * The tail, `x >= 4`, where the value uses the Cody rational and the
 * derivative uses its correction factor directly.
 *
 * Every other case here stays below 2.6, so before this the tail branch was
 * never executed on the device, for the value or for the derivative.
 */
TEST(OpenCLerfcx, prim_rev_values_tail) {
  Eigen::VectorXd a(8);
  a << 3.9, 4.0, 4.1, 6.0, 10.0, 100.0, 1000.0, 10000.0;
  stan::math::test::compare_cpu_opencl_prim_rev(erfcx_functor, a);
}

/**
 * Device derivative in the tail against fixed references.
 *
 * `compare_cpu_opencl_prim_rev` cannot check this. It uses
 * `EXPECT_NEAR_REL`, whose relative tolerance floors its denominator at 1,
 * so for a derivative of magnitude 5.6e-09 even a 2.4e-08 relative error is
 * a 1.4e-16 absolute difference and passes. Verified: the differential test
 * above passes against both the old and the new device derivative.
 *
 * References are `2 * x * erfcx(x) - 2 / sqrt(pi)` at 60 significant
 * digits, rounded to double. They are the same values the CPU suite uses in
 * `test/unit/math/mix/fun/erfcx_derivative_test.cpp`.
 */
TEST(OpenCLerfcx, rev_tail_derivative_against_references) {
  const int N = 6;
  Eigen::VectorXd x(N);
  x << 4.0, 6.0, 10.0, 100.0, 1000.0, 10000.0;
  Eigen::VectorXd expected(N);
  expected << -0.032383506095021455, -0.015060353489052321,
      -0.0055593122190608565, -5.6410497625993184e-05, -5.641887372654967e-07,
      -5.641895750849127e-09;

  stan::math::var_value<stan::math::matrix_cl<double>> x_cl(
      stan::math::to_matrix_cl(x));
  auto y = stan::math::erfcx(x_cl);
  stan::math::var total = stan::math::sum(y);
  total.grad();

  const Eigen::VectorXd adj = stan::math::from_matrix_cl(x_cl.adj());
  for (int i = 0; i < N; ++i) {
    EXPECT_LT(std::fabs(adj[i] / expected[i] - 1.0), 1e-12)
        << "device derivative lost precision at x = " << x[i];
  }
  stan::math::set_zero_all_adjoints();
}

TEST(OpenCLerfcx, prim_rev_size_0) {
  int N = 0;

  Eigen::MatrixXd a(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(erfcx_functor, a);
}

TEST(OpenCLerfcx, prim_rev_values_large) {
  int N = 71;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(erfcx_functor, a);
}

#endif
