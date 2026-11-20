#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/inv_logit_tail_refs.hpp>

TEST(mathMixMatFun, log1mInvLogit) {
  auto f = [](const auto& x1) { return stan::math::log1m_inv_logit(x1); };
  stan::test::expect_common_nonzero_unary_vectorized<
      stan::test::ScalarSupport::Real>(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -1, -0.5, -0.2, 0.5, 1, 1.3,
                                      3, 5);
}

TEST(mathMixMatFun, log1minvlogit_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_args;
  auto f = [](const auto& x1) {
    using stan::math::log1m_inv_logit;
    return log1m_inv_logit(x1);
  };
  std::vector<double> com_args = common_args();
  std::vector<double> args{-2.6, -2, -1, -0.5, -0.2, 0.5, 1, 1.3, 3, 5};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}

TEST(mathMixMatFun, log1mInvLogitFwdDerivativeTails) {
  using inv_logit_tail_refs::ref;
  // d/dx = -inv_logit(x), d^2/dx^2 = -inv_logit'(x)
  inv_logit_tail_refs::expect_fwd_derivatives(
      [](const auto& x) { return stan::math::log1m_inv_logit(x); },
      [](const ref& r) { return -r.s; }, [](const ref& r) { return -r.d1; });
}
