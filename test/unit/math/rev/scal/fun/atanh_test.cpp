#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev,atanh) {
  AVAR a = 0.3;
  AVAR f = atanh(a);
  EXPECT_FLOAT_EQ(atanh(0.3), f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x,g);
  EXPECT_FLOAT_EQ(1.0/(1.0 - 0.3 * 0.3), g[0]);
}

TEST(AgradRev,atanh_1) {
  double inf = std::numeric_limits<double>::infinity();
  AVAR a = 1;
  AVAR f = atanh(a);
  EXPECT_FLOAT_EQ(inf, f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x,g);
  EXPECT_FLOAT_EQ(1.0/(1.0 - 1.0 * 1.0), g[0]);
}

TEST(AgradRev,atanh_neg_1) {
  double inf = std::numeric_limits<double>::infinity();
  AVAR a = -1;
  AVAR f = atanh(a);
  EXPECT_FLOAT_EQ(-inf, f.val());

  AVEC x = createAVEC(a);
  VEC g;
  f.grad(x,g);
  EXPECT_FLOAT_EQ(1.0/(1.0 - (-1.0 * -1.0)), g[0]);
}

TEST(AgradRev,atanh_out_of_bounds) {
  using stan::math::atanh;
  EXPECT_THROW(atanh(AVAR(-2)), std::domain_error);
  EXPECT_THROW(atanh(AVAR(1001.2)), std::domain_error);
}

struct atanh_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return atanh(arg1);
  }
};

TEST(AgradRev,atanh_NaN) {
  atanh_fun atanh_;
  test_nan(atanh_,false,true);
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 0.3;
  test::check_varis_on_stack(stan::math::atanh(a));
}
