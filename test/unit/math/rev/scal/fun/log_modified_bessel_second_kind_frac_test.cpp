#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

TEST(AgradRev, log_modified_bessel_second_kind_frac_base) {
  std::ostringstream err_out;

  double f1 = stan::math::log_modified_bessel_second_kind_frac(0.92, 1.73, err_out);
  EXPECT_FLOAT_EQ(-1.6402036, f1);

  double f2 = stan::math::log_modified_bessel_second_kind_frac(1.0, 2.0, err_out);
  EXPECT_FLOAT_EQ(-1.9670713, f2);
}

TEST(AgradRev, log_modified_bessel_second_kind_frac_int_var) {
  std::ostringstream err_out;

  AVAR a(1);
  AVAR b(2);
  AVAR f = stan::math::log_modified_bessel_second_kind_frac(a, b, err_out);
  EXPECT_FLOAT_EQ(-1.9670713, f.val());

  AVEC x = createAVEC(a, b);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0.40715387938, g[0]);
  EXPECT_FLOAT_EQ(-1.314307758763, g[1]);
}

TEST(AgradRev, log_modified_bessel_second_kind_frac_int_v) {
  std::ostringstream err_out;

  AVAR a(3);
  AVAR b(6.378);
  AVAR f = stan::math::log_modified_bessel_second_kind_frac(a, b, err_out);
  EXPECT_FLOAT_EQ(-6.4474869, f.val());

  AVEC x = createAVEC(a, b);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0.42755982, g[0]);
  EXPECT_FLOAT_EQ(-1.1689174, g[1]);

}

TEST(AgradRev, log_modified_bessel_second_kind_frac_int_z) {
  std::ostringstream err_out;

  AVAR a(0.11);
  AVAR b(8.0);
  AVAR f = stan::math::log_modified_bessel_second_kind_frac(a, b, err_out);
  EXPECT_FLOAT_EQ(-8.82797, f.val());

  AVEC x = createAVEC(a, b);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0.012987931, g[0]);
  EXPECT_FLOAT_EQ(-1.0608374, g[1]);

}


TEST(AgradRev, log_modified_bessel_second_kind_frac_general) {
  std::ostringstream err_out;

  AVAR a(0.92);
  AVAR b(1.73);
  AVAR f = stan::math::log_modified_bessel_second_kind_frac(a, b, err_out);
  EXPECT_FLOAT_EQ(-1.6402, f.val());

  AVEC x = createAVEC(a, b);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0.423179, g[0]);
  EXPECT_FLOAT_EQ(-1.35463, g[1]);
}

TEST(AgradRev, log_modified_bessel_second_kind_frac_general_2) {
  std::ostringstream err_out;

  AVAR a(0.93);
  AVAR b(1.71);
  AVAR f = stan::math::log_modified_bessel_second_kind_frac(a, b, err_out);
  EXPECT_FLOAT_EQ(-1.6402, f.val());

  AVEC x = createAVEC(a, b);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0.423179, g[0]);
  EXPECT_FLOAT_EQ(-1.35463, g[1]);
}
