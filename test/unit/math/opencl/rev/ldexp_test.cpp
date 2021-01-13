#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

auto ldexp_functor
    = [](const auto& a, const auto& b) { return stan::math::ldexp(a, b); };

TEST(OpenCLPrim, ldexp_small_zero) {
  stan::math::matrix_d d1(3, 3);
  d1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> d2(3, 3);
  d2 << 0, 8, 7, 6, 5, 4, 3, 2, 1;
  stan::math::test::compare_cpu_opencl_prim_rev(ldexp_functor, d1, d2);

  stan::math::matrix_d d0(0, 0);
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> i0(0, 0);
  stan::math::test::compare_cpu_opencl_prim_rev(ldexp_functor, d0, i0);
}

TEST(OpenCLPrim, ldexp_rev_exceptions) {
  stan::math::matrix_d md1(2, 2);
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> md2(3, 3);

  stan::math::var_value<stan::math::matrix_cl<double>> md11
      = stan::math::to_matrix_cl(md1);
  stan::math::matrix_cl<int> md22 = stan::math::to_matrix_cl(md2);
  EXPECT_THROW(stan::math::ldexp(md11, md22), std::invalid_argument);
}

TEST(OpenCLPrim, ldexp_big) {
  int N = 71;
  stan::math::matrix_d m = 10 * Eigen::MatrixXd::Random(N, N);
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> b = m.cast<int>();

  stan::math::matrix_d a = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(ldexp_functor, a, b);
}

TEST(OpenCLPrim, ldexp_scalar_exponent_big) {
  int N = 71;
  stan::math::matrix_d m = 10 * Eigen::MatrixXd::Random(N, N);
  int b = 2;

  stan::math::matrix_d a = Eigen::MatrixXd::Random(N, N);
  stan::math::test::compare_cpu_opencl_prim_rev(ldexp_functor, a, b);
}

TEST(OpenCLPrim, ldexp_scalar_significand_big) {
  int N = 71;
  stan::math::matrix_d m = 10 * Eigen::MatrixXd::Random(N, N);
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> b = m.cast<int>();

  double a = 1.1;
  stan::math::test::compare_cpu_opencl_prim_rev(ldexp_functor, a, b);
}
#endif
