#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <limits>

TEST(AgradRev, ldexp) {
  AVAR a = 1.3;
  AVAR f = stan::math::ldexp(a, 5);
  EXPECT_FLOAT_EQ(stan::math::ldexp(1.3, 5), f.val());

  AVEC x = createAVEC(a);
  VEC grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(stan::math::exp2(5), grad_f[0]);

  a = std::numeric_limits<AVAR>::infinity();
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(),
                  stan::math::ldexp(a, 5).val());
}

TEST(AgradRev, check_varis_on_stack) {
  AVAR a = 1.3;
  test::check_varis_on_stack(stan::math::ldexp(a, 5));
}
