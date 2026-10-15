#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <stan/math/opencl/kernels/device_functions/digamma.hpp>
#include <stan/math/opencl/kernels/device_functions/digamma_diff.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/expect_near_rel.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

static const std::string test_digamma_diff_kernel_code
    = STRINGIFY(__kernel void test(__global double *C, __global double *A,
                                   __global double *B) {
        const int i = get_global_id(0);
        C[i] = digamma_diff(A[i], B[i]);
      });

const stan::math::opencl_kernels::kernel_cl<
    stan::math::opencl_kernels::out_buffer,
    stan::math::opencl_kernels::in_buffer,
    stan::math::opencl_kernels::in_buffer>
    digamma_diff_kernel(
        "test", {stan::math::opencl_kernels::digamma_device_function,
                 stan::math::opencl_kernels::digamma_diff_device_function,
                 test_digamma_diff_kernel_code});

namespace {
Eigen::VectorXd digamma_diff_cl(const Eigen::VectorXd &x,
                                const Eigen::VectorXd &d) {
  stan::math::matrix_cl<double> x_cl(x);
  stan::math::matrix_cl<double> d_cl(d);
  stan::math::matrix_cl<double> res_cl(x.size(), 1);
  digamma_diff_kernel(cl::NDRange(x.size()), res_cl, x_cl, d_cl);
  return stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);
}
}  // namespace

TEST(MathMatrixCL, digamma_diff) {
  Eigen::VectorXd x = Eigen::VectorXd::Random(1000).array() * 15 + 15.01;
  Eigen::VectorXd d = Eigen::VectorXd::Random(1000).array() * 50 + 50;
  Eigen::VectorXd res = digamma_diff_cl(x, d);
  stan::test::expect_near_rel("digamma_diff (OpenCL)", res,
                              stan::math::digamma_diff(x, d));
}

TEST(MathMatrixCL, digamma_diff_large_shapes_match_cpu) {
  // log-uniform x in [1e-3, 1e16] and d in [1e-3, 1e8]; the CPU function
  // has a relative error of a few ulp
  const int n = 4096;
  Eigen::VectorXd x
      = (Eigen::VectorXd::Random(n).array() * 9.5 + 6.5) * std::log(10.0);
  Eigen::VectorXd d
      = (Eigen::VectorXd::Random(n).array() * 5.5 + 2.5) * std::log(10.0);
  x = x.array().exp();
  d = d.array().exp();
  Eigen::VectorXd res = digamma_diff_cl(x, d);
  Eigen::VectorXd cpu = stan::math::digamma_diff(x, d);
  double max_rel = 0;
  for (int i = 0; i < n; ++i) {
    const double rel = std::fabs(res(i) - cpu(i)) / std::fabs(cpu(i));
    max_rel = std::max(max_rel, rel);
    EXPECT_NEAR(res(i), cpu(i), 1e-14 * std::fabs(cpu(i)))
        << "x = " << x(i) << ", d = " << d(i);
  }
  std::cout << "digamma_diff OpenCL vs CPU, max relative difference " << max_rel
            << " (" << max_rel / std::numeric_limits<double>::epsilon()
            << " eps)" << std::endl;
}

TEST(MathMatrixCL, digamma_diff_reference_values) {
  // psi(x + d) - psi(x) from mp.digamma at 120 digits for the exact double
  // arguments, checked against the same computation at 80 digits
  struct TestValue {
    double x;
    double d;
    double value;
  };
  const std::vector<TestValue> values = {
      {0x1.56e1fc2f8f359p-997, 0x1.0000000000000p+0, 9.9999999999999997e+299},
      {0x1.0624dd2f1a9fcp-10, 0x1.0000000000000p+0, 9.9999999999999998e+2},
      {0x1.0000000000000p-1, 0x0.0p+0, 0.0},
      {0x1.0000000000000p-1, 0x1.b7cdfd9d7bdbbp-34, 4.9348021997032397e-10},
      {0x1.0000000000000p-1, 0x1.8000000000000p+1, 3.0666666666666667},
      {0x1.4000000000000p+1, 0x1.0000000000000p-2, 1.1574438433018941e-1},
      {0x1.d333333333333p+2, 0x1.c800000000000p+5, 2.2379430906763218},
      {0x1.3ff7ced916873p+3, 0x1.0000000000000p+0, 1.000100010001e-1},
      {0x1.4000000000000p+3, 0x1.0000000000000p+0, 1.0e-1},
      {0x1.9000000000000p+3, 0x1.e848000000000p+19, 1.1330326906617404e+1},
      {0x1.f400000000000p+9, 0x1.0000000000000p+0, 1.0e-3},
      {0x1.7d78400000000p+26, 0x1.8000000000000p+1, 2.9999999700000005e-8},
      {0x1.d1a94a2000000p+39, 0x1.c800000000000p+5, 5.6999999998404e-11},
      {0x1.c6bf526340000p+49, 0x1.8000000000000p+1, 2.999999999999997e-15},
      {0x1.c6bf526340000p+49, 0x1.c6bf526340000p+49, 6.9314718055994556e-1},
      {0x1.550f7dca70000p+51, 0x1.d400000000000p+6, 3.8999999999999246e-14},
      {0x1.5af1d78b58c40p+66, 0x1.0000000000000p+0, 1.0e-20},
      {0x1.0000000000000p+0, 0x1.7e43c8800759cp+996, 6.9135274356311524e+2},
      {0x1.7e43c8800759cp+996, 0x1.7e43c8800759cp+996, 6.9314718055994531e-1},
  };
  const int n = values.size();
  Eigen::VectorXd x(n);
  Eigen::VectorXd d(n);
  for (int i = 0; i < n; ++i) {
    x(i) = values[i].x;
    d(i) = values[i].d;
  }
  Eigen::VectorXd res = digamma_diff_cl(x, d);
  double max_rel = 0;
  for (int i = 0; i < n; ++i) {
    const double ref = values[i].value;
    if (ref == 0) {
      EXPECT_EQ(res(i), 0.0) << "x = " << x(i) << ", d = " << d(i);
      continue;
    }
    const double rel = std::fabs(res(i) - ref) / std::fabs(ref);
    max_rel = std::max(max_rel, rel);
    EXPECT_NEAR(res(i), ref, 1e-14 * std::fabs(ref))
        << "x = " << x(i) << ", d = " << d(i);
  }
  std::cout << "digamma_diff OpenCL vs mpmath, max relative error " << max_rel
            << " (" << max_rel / std::numeric_limits<double>::epsilon()
            << " eps)" << std::endl;
}

TEST(MathMatrixCL, digamma_diff_edge_cases) {
  const double inf = std::numeric_limits<double>::infinity();
  Eigen::VectorXd x(7);
  x << NAN, 1.5, 0.0, -1.0, 1.5, inf, 1.5;
  Eigen::VectorXd d(7);
  d << 1.0, NAN, 1.0, 1.0, -1.0, 3.0, inf;
  Eigen::VectorXd res = digamma_diff_cl(x, d);
  for (int i = 0; i < 5; ++i) {
    EXPECT_TRUE(std::isnan(res(i))) << "x = " << x(i) << ", d = " << d(i);
  }
  EXPECT_EQ(res(5), 0.0);
  EXPECT_EQ(res(6), inf);
}

#endif
