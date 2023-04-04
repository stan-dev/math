#include <stan/math/prim/functor.hpp>
#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/prob.hpp>
#include <stan/math/rev.hpp>

#include <gtest/gtest.h>
#include <vector>

TEST(mathPrimScalProbWienerFullScal, valid) {
  using stan::math::INFTY;
  using stan::math::wiener_full_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_NO_THROW(wiener_full_lpdf(rt, a, t0, w, v, sv, sw, st0));
  rt = 5;
  a = 1;
  v = 1;
  w = 0.5;
  t0 = 0.0;
  sv = 0.0;
  sw = 0.0;
  st0 = 0.0;
  EXPECT_NO_THROW(wiener_full_lpdf(rt, a, t0, w, v, sv, sw, st0));
}

// rt
TEST(mathPrimScalProbWienerFullScal, invalid_rt) {
  using stan::math::INFTY;
  using stan::math::wiener_full_lpdf;
  double a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_full_lpdf(0, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(-1, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(INFTY, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(-INFTY, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(NAN, a, t0, w, v, sv, sw, st0),
               std::domain_error);
}

// a
TEST(mathPrimScalProbWienerFullScal, invalid_a) {
  using stan::math::INFTY;
  using stan::math::wiener_full_lpdf;
  double rt = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_full_lpdf(rt, 0, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, -1, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, INFTY, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, -INFTY, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, NAN, t0, w, v, sv, sw, st0),
               std::domain_error);
}

// v
TEST(mathPrimScalProbWienerFullScal, invalid_v) {
  using stan::math::INFTY;
  using stan::math::wiener_full_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, INFTY, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, -INFTY, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, NAN, sv, sw, st0),
               std::domain_error);
}

// w
TEST(mathPrimScalProbWienerFullScal, invalid_w) {
  using stan::math::INFTY;
  using stan::math::wiener_full_lpdf;
  double rt = 1, a = 1, v = -1, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, -0.1, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, 0, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, 1, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, 1.1, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, INFTY, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, -INFTY, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, NAN, v, sv, sw, st0),
               std::domain_error);
}

// t0
TEST(mathPrimScalProbWienerFullScal, invalid_t0) {
  using stan::math::INFTY;
  using stan::math::wiener_full_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_full_lpdf(rt, a, 2, w, v, sv, sw, st0),
               std::domain_error);  // rt must be greater than t0
  EXPECT_THROW(wiener_full_lpdf(rt, a, -1, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, INFTY, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, -INFTY, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, NAN, w, v, sv, sw, st0),
               std::domain_error);
}

// sv
TEST(mathPrimScalProbWienerFullScal, invalid_sv) {
  using stan::math::INFTY;
  using stan::math::wiener_full_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, v, -1, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, v, INFTY, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, v, -INFTY, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, v, NAN, sw, st0),
               std::domain_error);
}

// sw
TEST(mathPrimScalProbWienerFullScal, invalid_sw) {
  using stan::math::INFTY;
  using stan::math::wiener_full_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, v, sv, -1, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, 0.8, v, sv, 0.5, st0),
               std::domain_error);  // sw must be smaller than 2*(1-w)
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, 0.3, v, sv, 0.7, st0),
               std::domain_error);  // sw must be smaller than 2*w
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, v, sv, INFTY, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, v, sv, -INFTY, st0),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, v, sv, NAN, st0),
               std::domain_error);
}

// st0
TEST(mathPrimScalProbWienerFullScal, invalid_st0) {
  using stan::math::INFTY;
  using stan::math::wiener_full_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2;
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, v, sv, sw, -1),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, v, sv, sw, INFTY),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, v, sv, sw, -INFTY),
               std::domain_error);
  EXPECT_THROW(wiener_full_lpdf(rt, a, t0, w, v, sv, sw, NAN),
               std::domain_error);
}


