#include <stan/math/prim/functor.hpp>
#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/prob.hpp>
#include <stan/math/rev.hpp>

#include <gtest/gtest.h>
#include <vector>

TEST(mathPrimScalProbWienerFullLccdfScal, valid) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_NO_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, st0));
  rt = 5;
  a = 1;
  v = 1;
  w = 0.5;
  t0 = 0.0;
  sv = 0.0;
  sw = 0.0;
  st0 = 0.0;
  EXPECT_NO_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, st0));
}

// rt
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_rt) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(0, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(-1, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(INFTY, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(-INFTY, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(NAN, a, t0, w, v, sv, sw, st0),
               std::domain_error);
}

// a
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_a) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, 0, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, -1, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, INFTY, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, -INFTY, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, NAN, t0, w, v, sv, sw, st0),
               std::domain_error);
}

// v
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_v) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, INFTY, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, -INFTY, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, NAN, sv, sw, st0),
               std::domain_error);
}

// w
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_w) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, -0.1, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 0, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 1, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 1.1, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, INFTY, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, -INFTY, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, NAN, v, sv, sw, st0),
               std::domain_error);
}

// t0
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_t0) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, 2, w, v, sv, sw, st0),
               std::domain_error);  // rt must be greater than t0
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, -1, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, INFTY, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, -INFTY, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, NAN, w, v, sv, sw, st0),
               std::domain_error);
}

// sv
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_sv) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, -1, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, INFTY, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, -INFTY, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, NAN, sw, st0),
               std::domain_error);
}

// sw
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_sw) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, -1, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 0.8, v, sv, 0.5, st0),
               std::domain_error);  // sw must be smaller than 2*(1-w)
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 0.3, v, sv, 0.7, st0),
               std::domain_error);  // sw must be smaller than 2*w
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, INFTY, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, -INFTY, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, NAN, st0),
               std::domain_error);
}

// st0
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_st0) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, -1),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, INFTY),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, -INFTY),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, NAN),
               std::domain_error);
}

TEST(mathPrimScalProbWienerFullLccdfPrecScal, valid) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_NO_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, st0, 1e-4));
  rt = 5;
  a = 1;
  v = 1;
  w = 0.5;
  t0 = 0.0;
  sv = 0.0;
  sw = 0.0;
  st0 = 0.0;
  EXPECT_NO_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, st0, 1e-4));
}

// rt
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_rt) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(0, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(-1, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(INFTY, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(-INFTY, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(NAN, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// a
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_a) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, 0, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, -1, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, INFTY, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, -INFTY, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, NAN, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// v
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_v) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, INFTY, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, -INFTY, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, NAN, sv, sw, st0, 1e-4),
               std::domain_error);
}

// w
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_w) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, -0.1, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 0, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 1, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 1.1, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, INFTY, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, -INFTY, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, NAN, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// t0
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_t0) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, 2, w, v, sv, sw, st0, 1e-4),
               std::domain_error);  // rt must be greater than t0
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, -1, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, INFTY, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, -INFTY, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, NAN, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// sv
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_sv) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, -1, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, INFTY, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, -INFTY, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, NAN, sw, st0, 1e-4),
               std::domain_error);
}

// sw
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_sw) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, -1, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 0.8, v, sv, 0.5, st0, 1e-4),
               std::domain_error);  // sw must be smaller than 2*(1-w)
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 0.3, v, sv, 0.7, st0, 1e-4),
               std::domain_error);  // sw must be smaller than 2*w
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, INFTY, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, -INFTY, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, NAN, st0, 1e-4),
               std::domain_error);
}

// st0
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_st0) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, -1, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, INFTY, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, -INFTY, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, NAN, 1e-4),
               std::domain_error);
}
