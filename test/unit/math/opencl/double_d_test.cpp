
#include <stan/math/opencl/double_d.hpp>
//#include <test/unit/math/expect_near_rel.hpp>
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
  EXPECT_TRUE(res.high == 2.0);
  EXPECT_TRUE(res.low == eps);

  // carry
  res = add_dd_dd(a, c);
  EXPECT_NORMALIZED(res);
  EXPECT_TRUE(res.high == 1.0 + eps);
  EXPECT_TRUE(res.low == 0.0);

  // cancelation
  res = add_dd_dd(a, d);
  EXPECT_NORMALIZED(res);
  EXPECT_TRUE(res.high == eps * h_eps);
  EXPECT_TRUE(res.low == 0.0);
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
  EXPECT_TRUE(res.high == h_eps);
  EXPECT_TRUE(res.low == h_eps * h_eps * 0.001);

  // carry
  res = mul_dd_dd(c, c);
  EXPECT_NORMALIZED(res);
  EXPECT_TRUE(res.high == 1.0 + eps);
  EXPECT_TRUE(res.low == 0.0);

  // cancelation
  res = mul_dd_dd(c, d);
  EXPECT_NORMALIZED(res);
  EXPECT_TRUE(res.high == 1.0);
  EXPECT_TRUE(res.low == 0.0);
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
  EXPECT_TRUE(res.high == 1.0);
  EXPECT_TRUE(res.low == 0.0);

  res = div_dd_dd(a, b);
  EXPECT_NORMALIZED(res);
  EXPECT_TRUE(res.high == 1.0);
  EXPECT_TRUE(res.low == h_eps * eps);

  // div by zero
  res = div_dd_dd(a, c);
  EXPECT_TRUE(res.high == std::numeric_limits<double>::infinity());
  EXPECT_TRUE(res.low == 0);
}

#include <quadmath.h>

#define EXPECT_DD_F128_EQ(dd, f128) \
  tmp = (dd);              \
  EXPECT_NEAR((__float128)tmp.high + tmp.low, (f128), 1e-30*1e30)

TEST(double_d, all) {
  using stan::math::internal::double_d;
  for (int i = 0; i < 10000000; i++) {
    double_d tmp;
    double high = Eigen::MatrixXd::Random(0, 0).coeff(0, 0)*1e30;
    double low = Eigen::MatrixXd::Random(0, 0).coeff(0, 0) * 1e-20*1e30;
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
    EXPECT_DD_F128_EQ(stan::math::internal::mul_d_d(high, high2), f1283 * f1284);
  }
}




#ifdef STAN_OPENCL

#endif
