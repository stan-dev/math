#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <limits>

TEST(AgradRev, log_a) {
  AVAR a(5.0);
  AVAR f = log(a);
  EXPECT_FLOAT_EQ(log(5.0), f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(1.0 / 5.0, g[0]);
}

TEST(AgradRev, log_inf) {
  AVAR a = std::numeric_limits<double>::infinity();
  EXPECT_TRUE(boost::math::isinf(log(a)));
}

TEST(AgradRev, log_0) {
  AVAR a(0.0);
  EXPECT_TRUE(boost::math::isinf(log(a)) && (log(a) < 0.0));
}

TEST(AgradRev, log_neg) {
  AVAR a(0.0 - stan::math::EPSILON);
  EXPECT_TRUE(std::isnan(log(a)));
}

struct log_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    return log(arg1);
  }
};

TEST(AgradRev, log_NaN) {
  log_fun log_;
  test_nan(log_, false, true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a(5.0);
  test::check_varis_on_stack(stan::math::log(a));
}

TEST(AgradRev, complex) {
  stan::math::var x = stan::math::pi();
  std::complex<stan::math::var> z(x, 2.0 / 3);
  EXPECT_TRUE(log(conj(z)) == conj(log(z)));

  double h = 1e-8;
  z = std::complex<stan::math::var>(x, h);
  auto f = log(z);

  AVEC v = createAVEC(real(z));
  VEC g;
  real(f).grad(v, g);
  EXPECT_FLOAT_EQ(g[0], imag(f).val() / h);
}
