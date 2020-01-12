#include <test/unit/math/test_ad.hpp>
#include <limits>

auto f = [](const auto& a, const auto& b, const auto& z) {
  typedef stan::return_type_t<decltype(a), decltype(b), decltype(z)> T;
  T digamma_a = stan::math::digamma(a);
  T digamma_ab = stan::math::digamma(a + b);
  T aT = a;
  T bT = b;
  T zT = z;
  return stan::math::inc_beta_dda(aT, bT, zT, digamma_a, digamma_ab);
};

TEST(mathMixScalFun, inc_beta_dda) {
  // replicate pre-existing test cases
  stan::test::expect_ad(f, 0.6, 0.3, 0.5);
  stan::test::expect_ad(f, 0.9, 0.9, 0.9);
  stan::test::expect_ad(f, 1.0, 1.0, 0.4);
}

TEST(mathMixScalFun, inc_beta_dda_integer) {
  stan::test::expect_ad(f, 3, 2, 0.2);
  stan::test::expect_ad(f, 3, 2.0, 0.2);
  stan::test::expect_ad(f, 3.0, 2, 0.2);
  stan::test::expect_ad(f, 3.0, 2.0, 0.2);
}

TEST(mathMixScalFun, inc_beta_dda_a_boundary) {
  const double inf = stan::math::INFTY;
  stan::test::expect_ad(f, inf, 0.5, 0.5);
}

TEST(mathMixScalFun, inc_beta_dda_b_boundary) {
  const double inf = stan::math::INFTY;
  stan::test::expect_ad(f, 0.5, inf, 0.5);
}

TEST(mathMixScalFun, inc_beta_dda_a_b_boundary) {
  stan::test::expect_ad(f, 0.0, 0.0, 0.5);
}
