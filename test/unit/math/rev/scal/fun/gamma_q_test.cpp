#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <boost/math/special_functions/gamma.hpp>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>

// leaving regression tests here for value
// main tests moved to mix/scal/fun

TEST(AgradRevGammaQ, infLoopInVersion2_0_1_var_var) {
  AVAR a = 8.01006;
  AVAR b = 2.47579e+215;
  AVEC x = createAVEC(a, b);

  AVAR f = gamma_q(a, b);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0, g[0]);
}
TEST(AgradRevGammaQ, infLoopInVersion2_0_1_var_double) {
  AVAR a = 8.01006;
  double b = 2.47579e+215;
  AVEC x = createAVEC(a);

  AVAR f = gamma_q(a, b);
  VEC g;
  f.grad(x, g);
  EXPECT_FLOAT_EQ(0, g[0]);
}
