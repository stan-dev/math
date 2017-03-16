#include <gtest/gtest.h>
#include <stan/math/prim/scal/fun/F32.hpp>

TEST(MathPrimScalFun, F32) {
  // This should terminate by convergence
  EXPECT_NEAR(2.5,
              stan::math::F32(1.0, 1.0, 1.0, 1.0, 1.0, 0.6),
              1e-6);

  // This should terminate by zero numerator, no sign-flip
  EXPECT_NEAR(11.28915378492300834453857665243661995978358572684,
              stan::math::F32(1.0, 31.0, -27.0, 19.0, -41.0, 1.0),
              1e-6);

  // This should terminate by zero numerator, single-step, no sign-flip
  EXPECT_NEAR(-0.2,
              stan::math::F32(1.0, 12.0, -1.0, 10.0, 1.0, 1.0),
              1e-6);

  // This should terminate by NaN increment, single-step, no sign-flip
  EXPECT_NEAR(2.2,
              stan::math::F32(1.0, 12.0, -1.0, 10.0, -1.0, 1.0),
              1e-6);

  // This should terminate by convergence, single sign flip via numerator
  EXPECT_NEAR(0.9693532463066744390558635404413378640,
              stan::math::F32(1.0, -.5, 2.0, 10.0, 1.0, 0.3, 1e-8),
              1e-6);

  // This should throw (Mathematica claims the answer is -10 but... ?
  EXPECT_THROW(
              stan::math::F32(1.0, 12.0, 1.0, 10.0, 1.0, 1.1),
              std::domain_error);

  // This should terminate by convergence, double sign flip via
  // numerator
  EXPECT_NEAR(1.0371188919802822614906491960448352702,
              stan::math::F32(1.0, -.5, -2.5, 10.0, 1.0, 0.3, 1e-10),
              1e-6);

  // This should terminate by convergence, double sign flip via
  // numerator
  EXPECT_NEAR(1.0659384611044132367442134708091552077,
              stan::math::F32(1.0, -.5, -4.5, 10.0, 1.0, 0.3, 1e-10),
              1e-6);

  // This should terminate by convergence, double sign flip via
  // numerator
  EXPECT_NEAR(0.5310817155578250445980684217963963819041890167,
              stan::math::F32(-6.5, -.5, -4.5, 10.0, 1.0, 0.3, 1e-10),
              1e-6);
}


