
#include <stan/math/opencl/double_d.hpp>
#include <gtest/gtest.h>
#include <limits>

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
  double_d b{h_eps, h_eps * h_eps * h_eps* h_eps};
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
  using stan::math::internal::double_d;
  using stan::math::internal::div_dd_dd;
  double eps = std::numeric_limits<double>::epsilon();
  double h_eps = eps * 0.5;
  double_d a{1.0, h_eps};
  double_d b{1.0, h_eps*(1-eps)};
  double_d c{1.0, h_eps};
  double_d d{1.0, -h_eps};

  // simple
  double_d res = div_dd_dd(a, a);
  EXPECT_NORMALIZED(res);
  EXPECT_TRUE(res.high == 1.0);
  EXPECT_TRUE(res.low == 0.0);

  res = div_dd_dd(a, b);
  EXPECT_NORMALIZED(res);
  EXPECT_TRUE(res.high == 1.0);
  EXPECT_TRUE(res.low == h_eps*eps);
}

#ifdef STAN_OPENCL

#endif
