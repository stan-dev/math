#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev, modified_bessel_first_kind_int_var) {
  int a(1);
  AVAR b(4.0);
  AVAR f = stan::math::modified_bessel_first_kind(a, b);
  EXPECT_FLOAT_EQ(9.75946515370444990947519256731268090, f.val());

  AVEC x = createAVEC(a, b);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0, g[0]);
  EXPECT_FLOAT_EQ(8.862055663710218018987472041388932, g[1]);

  a = -1;
  b = -3.0;
  f = stan::math::modified_bessel_first_kind(a, b);

  EXPECT_FLOAT_EQ(-3.95337021740260939647863574058058,
                  f.val());

  x = createAVEC(a, b);
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0, g[0]);
  EXPECT_FLOAT_EQ(3.5630025133974876201183569658, g[1]);
}

struct modified_bessel_first_kind_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return modified_bessel_first_kind(1, arg1);
  }
};

TEST(AgradRev, modified_bessel_first_kind_NaN) {
  modified_bessel_first_kind_fun modified_bessel_first_kind_;
  test_nan(modified_bessel_first_kind_, true, false);
}

TEST(AgradRev, check_varis_on_stack) {
  int a(1);
  AVAR b(4.0);

  test::check_varis_on_stack(stan::math::modified_bessel_first_kind(a, b));
}
