#include <stan/math/opencl/double_d.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <Eigen/Core>

#define EXPECT_NORMALIZED(a) \
  EXPECT_LT(std::abs(a.low), \
            std::abs(a.high) * std::numeric_limits<double>::epsilon());

TEST(double_d, add_dd_dd_test) {
  using stan::math::internal::add_dd_dd;
  using stan::math::internal::double_d;
  double eps = std::numeric_limits<double>::epsilon();
  double h_eps = eps * 0.5;
  double_d a{1.0, h_eps};
  double_d c{h_eps, 0.0};
  double_d d{-1.0, -h_eps + eps * h_eps};

  // simple
  double_d res = add_dd_dd(a, a);
  EXPECT_NORMALIZED(res);
  EXPECT_EQ(res.high, 2.0);
  EXPECT_EQ(res.low, eps);

  // carry
  res = add_dd_dd(a, c);
  EXPECT_NORMALIZED(res);
  EXPECT_EQ(res.high, 1.0 + eps);
  EXPECT_EQ(res.low, 0.0);

  // cancelation
  res = add_dd_dd(a, d);
  EXPECT_NORMALIZED(res);
  EXPECT_EQ(res.high, eps * h_eps);
  EXPECT_EQ(res.low, 0.0);
}

TEST(double_d, mul_dd_dd_test) {
  using stan::math::internal::double_d;
  using stan::math::internal::mul_dd_dd;
  double eps = std::numeric_limits<double>::epsilon();
  double h_eps = eps * 0.5;
  double_d a{1.0, h_eps * 0.001};
  double_d b{h_eps, h_eps * h_eps * h_eps * h_eps};
  double_d c{1.0, h_eps};
  double_d d{1.0, -h_eps};

  // simple
  double_d res = mul_dd_dd(a, b);
  EXPECT_NORMALIZED(res);
  EXPECT_EQ(res.high, h_eps);
  EXPECT_EQ(res.low, h_eps * h_eps * 0.001);

  // carry
  res = mul_dd_dd(c, c);
  EXPECT_NORMALIZED(res);
  EXPECT_EQ(res.high, 1.0 + eps);
  EXPECT_EQ(res.low, 0.0);

  // cancelation
  res = mul_dd_dd(c, d);
  EXPECT_NORMALIZED(res);
  EXPECT_EQ(res.high, 1.0);
  EXPECT_EQ(res.low, 0.0);
}

TEST(double_d, div_dd_dd_test) {
  using stan::math::internal::div_dd_dd;
  using stan::math::internal::double_d;
  double eps = std::numeric_limits<double>::epsilon();
  double h_eps = eps * 0.5;
  double_d a{1.0, h_eps};
  double_d b{1.0, h_eps * (1 - eps)};
  double_d c{0.0, 0.0};
  double_d d{1.0, -h_eps};

  // simple
  double_d res = div_dd_dd(a, a);
  EXPECT_NORMALIZED(res);
  EXPECT_EQ(res.high, 1.0);
  EXPECT_EQ(res.low, 0.0);

  res = div_dd_dd(a, b);
  EXPECT_NORMALIZED(res);
  EXPECT_EQ(res.high, 1.0);
  EXPECT_EQ(res.low, h_eps * eps);

  // div by zero
  res = div_dd_dd(a, c);
  EXPECT_EQ(res.high, std::numeric_limits<double>::infinity());
  EXPECT_EQ(res.low, 0);
}

// this is the best test we have, but it relies on the nonstandard gcc extension
// quadmath, so it is not run by default
#ifdef GCC_QUADMATH_TEST
#include <quadmath.h>

#define EXPECT_DD_F128_EQ(dd, f128) \
  tmp = (dd);                       \
  EXPECT_NEAR((__float128)tmp.high + tmp.low, (f128), 1e-30 * 1e30)

TEST(double_d, all) {
  using stan::math::internal::double_d;
  for (int i = 0; i < 10000000; i++) {
    double_d tmp;
    double high = Eigen::MatrixXd::Random(0, 0).coeff(0, 0) * 1e30;
    double low = Eigen::MatrixXd::Random(0, 0).coeff(0, 0) * 1e-20 * 1e30;
    double_d dd = high;
    dd.low = low;
    __float128 f128 = high;
    f128 += low;
    EXPECT_DD_F128_EQ(dd, f128);

    double high2 = Eigen::MatrixXd::Random(0, 0).coeff(0, 0);
    double low2 = Eigen::MatrixXd::Random(0, 0).coeff(0, 0) * 1e-20;
    double_d dd2 = high2;
    dd2.low = low2;
    __float128 f1282 = high2;
    f1282 += low2;
    EXPECT_DD_F128_EQ(dd2, f1282);

    EXPECT_DD_F128_EQ(dd + dd2, f128 + f1282);
    EXPECT_DD_F128_EQ(dd - dd2, f128 - f1282);
    EXPECT_DD_F128_EQ(dd * dd2, f128 * f1282);
    EXPECT_DD_F128_EQ(dd / dd2, f128 / f1282);

    __float128 f1283 = high;
    __float128 f1284 = high2;
    EXPECT_DD_F128_EQ(stan::math::internal::mul_d_d(high, high2),
                      f1283 * f1284);
  }
}
#endif

#ifdef STAN_OPENCL

#include <stan/math/opencl/prim.hpp>
static const std::string double_d_test_kernel_code
    = STRINGIFY(__kernel void double_d_test_kernel(__global double_d *C,
                                                   const __global double_d *A,
                                                   const __global double *B) {
        const int i = get_global_id(0);
        C[i] = mul_dd_d(A[i], B[i]);
      });

const stan::math::opencl_kernels::kernel_cl<
    stan::math::opencl_kernels::out_buffer,
    stan::math::opencl_kernels::in_buffer,
    stan::math::opencl_kernels::in_buffer>
    double_d_test_kernel("double_d_test_kernel",
                         {stan::math::internal::double_d_src,
                          double_d_test_kernel_code});

TEST(double_d, opencl) {
  using stan::math::internal::double_d;
  using VectorXdd = Eigen::Matrix<double_d, -1, 1>;
  int n = 10;
  VectorXdd a(n);
  Eigen::VectorXd b = Eigen::VectorXd::Random(n);
  for (int i = 0; i < n; i++) {
    a[i].high = Eigen::VectorXd::Random(1)[0];
    a[i].low = Eigen::VectorXd::Random(1)[0] * 1e-17;
  }
  stan::math::matrix_cl<double_d> a_cl(a);
  stan::math::matrix_cl<double> b_cl(b);
  stan::math::matrix_cl<double_d> c_cl(n, 1);
  double_d_test_kernel(cl::NDRange(n), c_cl, a_cl, b_cl);
  VectorXdd c = stan::math::from_matrix_cl(c_cl);
  for (int i = 0; i < n; i++) {
    double_d correct = a[i] * b[i];
    EXPECT_EQ(c[i].high, correct.high);
    EXPECT_NEAR(c[i].low, correct.low, 1e-30);
  }
}

#endif
