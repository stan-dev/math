#include <gtest/gtest.h>
#include <stan/math/prim/scal/fun/F32.hpp>

TEST(MathPrimScalFun, F32_converges_by_z) { // converge
  EXPECT_NEAR(2.5,
              stan::math::F32(1.0, 1.0, 1.0, 1.0, 1.0, 0.6, 1e-10),
              1e-8);
}

TEST(MathPrimScalFun, F32_polynomial) { // terminate by zero numerator, no sign-flip
  EXPECT_NEAR(11.28855722705942,
              stan::math::F32(1.0, 31.0, -27.0, 19.0, -41.0, .99999, 1e-10),
              1e-8);
}

TEST(MathPrimScalFun, F32_short_polynomial) { // terminate by zero numerator, single-step, no sign-flip
  EXPECT_NEAR(-0.08000000000000007,
              stan::math::F32(1.0, 12.0, -1.0, 10.0, 1.0, .9),
              1e-8);
}

TEST(MathPrimScalFun, F32_short_polynomial_undef) { // at pole, should throw 
  EXPECT_THROW(
              stan::math::F32(1.0, 12.0, -1.0, 10.0, -1.0, 1.0),
              std::domain_error);
}

TEST(MathPrimScalFun, F32_sign_flip_numerator) { // converge, single sign flip via numerator
  EXPECT_NEAR(0.96935324630667443905,
              stan::math::F32(1.0, -.5, 2.0, 10.0, 1.0, 0.3, 1e-10),
              1e-8);
}

TEST(MathPrimScalFun, F32_diverge_by_z) {
  // This should throw (Mathematica claims the answer is -10 but... ?
  EXPECT_THROW(
              stan::math::F32(1.0, 12.0, 1.0, 10.0, 1.0, 1.1),
              std::domain_error);
}

TEST(MathPrimScalFun, F32_double_sign_flip) { // convergence, double sign flip 
  EXPECT_NEAR(1.03711889198028226149,
              stan::math::F32(1.0, -.5, -2.5, 10.0, 1.0, 0.3, 1e-10),
              1e-8);
  EXPECT_NEAR(1.06593846110441323674,
              stan::math::F32(1.0, -.5, -4.5, 10.0, 1.0, 0.3, 1e-10),
              1e-8);
}



