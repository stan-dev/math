#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/constraint/lub_constrain_helpers.hpp>
#include <test/unit/math/mix/fun/inv_logit_tail_refs.hpp>

TEST(mathMixMatFun, lub_constrain_scalars) {
  double x1 = 0.7;
  double x2 = -38.1;
  double lb = -2.0;
  double ub = 3.5;
  lub_constrain_tests::expect(x1, lb, ub);
  lub_constrain_tests::expect(x2, lb, ub);
  lub_constrain_tests::expect(x1, lb, lb);
  lub_constrain_tests::expect(x2, lb, lb);
  // ub inf
  auto ub_inf = stan::math::INFTY;
  lub_constrain_tests::expect(x1, lb, ub_inf);
  lub_constrain_tests::expect(x2, lb, ub_inf);

  // lb inf
  auto lb_inf = stan::math::NEGATIVE_INFTY;
  lub_constrain_tests::expect(x1, lb_inf, ub);
  lub_constrain_tests::expect(x2, lb_inf, ub);

  // both inf
  lub_constrain_tests::expect(x1, lb_inf, ub_inf);
  lub_constrain_tests::expect(x2, lb_inf, ub_inf);
}

TEST(mathMixMatFun, lub_constrain_fwd_derivative_tails) {
  using inv_logit_tail_refs::ref;
  // derivatives of -2 + 5 inv_logit(x)
  inv_logit_tail_refs::expect_fwd_derivatives(
      [](const auto& x) { return stan::math::lub_constrain(x, -2.0, 3.0); },
      [](const ref& r) { return 5 * r.d1; },
      [](const ref& r) { return 5 * r.d2; });
}
