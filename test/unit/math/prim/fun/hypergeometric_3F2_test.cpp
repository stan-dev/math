#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

// converge
TEST(MathPrimScalFun, F32_converges_by_z) {
  using stan::math::hypergeometric_3F2;
  using stan::math::to_row_vector;
  using stan::math::to_vector;
  std::vector<double> a = {1.0, 1.0, 1.0};
  std::vector<double> b = {1.0, 1.0};
  double z = 0.6;

  EXPECT_NEAR(2.5, hypergeometric_3F2(a, b, z), 1e-8);
  EXPECT_NEAR(2.5, hypergeometric_3F2(to_vector(a), to_vector(b), z), 1e-8);
  EXPECT_NEAR(2.5, hypergeometric_3F2(to_row_vector(a), to_row_vector(b), z),
              1e-8);
}
// terminate by zero numerator, no sign-flip
TEST(MathPrimScalFun, F32_polynomial) {
  EXPECT_NEAR(
      11.28855722705942,
      stan::math::hypergeometric_3F2({1.0, 31.0, -27.0}, {19.0, -41.0}, .99999),
      1e-8);
}
// terminate by zero numerator, single-step, no sign-flip
TEST(MathPrimScalFun, F32_short_polynomial) {
  EXPECT_NEAR(
      -0.08000000000000007,
      stan::math::hypergeometric_3F2({1.0, 12.0, -1.0}, {10.0, 1.0}, .9), 1e-8);
}
// terminate by zero numerator at k = 1; the denominator (b2)_k is zero only
// from k = 2 on, so the polynomial is defined
TEST(MathPrimScalFun, F32_short_polynomial_denominator_equal_numerator) {
  EXPECT_NEAR(
      2.2, stan::math::hypergeometric_3F2({1.0, 12.0, -1.0}, {10.0, -1.0}, 1.0),
      1e-14);
}
// at pole, should throw
TEST(MathPrimScalFun, F32_short_polynomial_undef) {
  EXPECT_THROW(
      stan::math::hypergeometric_3F2({1.0, 12.0, -2.0}, {10.0, -1.0}, 1.0),
      std::domain_error);
}
// converge, single sign flip via numerator
TEST(MathPrimScalFun, F32_sign_flip_numerator) {
  EXPECT_NEAR(0.96935324630667443905,
              stan::math::hypergeometric_3F2({1.0, -.5, 2.0}, {10.0, 1.0}, 0.3),
              1e-8);
}

TEST(MathPrimScalFun, F32_diverge_by_z) {
  // This should throw (Mathematica claims the answer is -10 but... ?
  EXPECT_THROW(
      stan::math::hypergeometric_3F2({1.0, 12.0, 1.0}, {10.0, 1.0}, 1.1),
      std::domain_error);
}
// convergence, double sign flip
TEST(MathPrimScalFun, F32_double_sign_flip) {
  EXPECT_NEAR(
      1.03711889198028226149,
      stan::math::hypergeometric_3F2({1.0, -.5, -2.5}, {10.0, 1.0}, 0.3), 1e-8);
  EXPECT_NEAR(
      1.06593846110441323674,
      stan::math::hypergeometric_3F2({1.0, -.5, -4.5}, {10.0, 1.0}, 0.3), 1e-8);
}

// The tests below use z = 1 and sum(b) <= sum(a). The reference values are
// the exact finite sums.

// terminate by zero numerator at k = 6, with large a2 and b2
TEST(MathPrimScalFun, F32_polynomial_large_parameters) {
  EXPECT_NEAR(17.62967394122564,
              stan::math::hypergeometric_3F2({1.0, 302.30333970, -6.0},
                                             {5.0, -151.75306521}, 1.0),
              1e-12);
}
// a numerator parameter equal to zero: only the first term is not zero
TEST(MathPrimScalFun, F32_zero_numerator) {
  EXPECT_EQ(1.0, stan::math::hypergeometric_3F2({1.0, 302.30333970, 0.0},
                                                {5.0, -151.75306521}, 1.0));
  EXPECT_EQ(1.0,
            stan::math::hypergeometric_3F2({1.0, 101.0, 0.0}, {2.0, 0.5}, 1.0));
}
// terminate by zero numerator at k = 6; b2 = a3, so the denominator (b2)_k
// is zero only from k = 7 on
TEST(MathPrimScalFun, F32_polynomial_denominator_equal_numerator) {
  EXPECT_NEAR(
      14.118912760416667,
      stan::math::hypergeometric_3F2({1.0, 6.5, -6.0}, {5.0, -6.0}, 1.0),
      1e-12);
}
// terminate by zero numerator at k = 6, with sum(b) == sum(a)
TEST(MathPrimScalFun, F32_polynomial_equal_parameter_sums) {
  EXPECT_NEAR(
      9.5967841682127396,
      stan::math::hypergeometric_3F2({1.0, 4.5, -6.0}, {5.0, -5.5}, 1.0),
      1e-12);
}
// a1 = -1 ends the series at k = 1, before (b1)_k is zero from k = 3 on; the
// larger non-positive integer a2 = -5 does not matter
TEST(MathPrimScalFun, F32_polynomial_smallest_numerator_ends) {
  EXPECT_NEAR(
      -0.25,
      stan::math::hypergeometric_3F2({-1.0, -5.0, 1.0}, {-2.0, 1.0}, 0.5),
      1e-15);
}
// the sign of z enters every term of the sum
TEST(MathPrimScalFun, F32_infsum_negative_z) {
  Eigen::VectorXd a(3);
  a << 1.0, 2.0, -3.0;
  Eigen::VectorXd b(2);
  b << 4.0, 5.0;
  EXPECT_NEAR(204.0 / 175.0,
              stan::math::internal::hypergeometric_3F2_infsum(a, b, -0.5),
              1e-15);
}
// the sum throws when the steps run out before the end of the series
TEST(MathPrimScalFun, F32_infsum_max_steps) {
  Eigen::VectorXd a(3);
  a << 1.0, 4.5, -6.0;
  Eigen::VectorXd b(2);
  b << 5.0, -5.5;
  EXPECT_THROW(
      stan::math::internal::hypergeometric_3F2_infsum(a, b, 1.0, 1e-6, 3),
      std::domain_error);
}
// terminate by zero numerator at k = 116; the terms fall below 1e-6 before
// the end of the series, so the sum must not stop at an absolute tolerance
TEST(MathPrimScalFun, F32_polynomial_small_terms) {
  EXPECT_NEAR(
      3.1748878200738662,
      stan::math::hypergeometric_3F2({1.0, 1.1, -116.0}, {2.0, -125.0}, 1.0),
      1e-13);
}
// all terms are positive; the term after 1 is about 2e-21, and the sum is
// about 6e36, so the sum must not stop at a small term
TEST(MathPrimScalFun, F32_polynomial_small_term_before_large_terms) {
  const double F = 5.7668140904936036e+36;
  EXPECT_NEAR(F,
              stan::math::hypergeometric_3F2({1e-23, 300.0, -100.0},
                                             {1.0, -150.5}, 1.0),
              1e-13 * F);
}
