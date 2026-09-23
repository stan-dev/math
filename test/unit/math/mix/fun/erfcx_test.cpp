#include <test/unit/math/test_ad.hpp>

// Finite arguments only: erfcx(-inf) is +inf and erfcx grows like
// 2 * exp(x * x) to the left, so the common-argument set (which includes the
// infinities) is not finite-differenceable here.
TEST(mathMixMatFun, erfcx) {
  auto f = [](const auto& x1) { return stan::math::erfcx(x1); };
  stan::test::expect_unary_vectorized(f, -2.6, -2, -1, -0.2, 0, 0.5, 1, 1.3,
                                      2.6, 3.9, 4.1, 12.0);
}

TEST(mathMixMatFun, erfcxmatvar) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  auto f = [](const auto& x1) { return stan::math::erfcx(x1); };
  std::vector<double> args{-2.6, -2, -1, -0.2, 0, 0.5, 1, 1.3, 2.6, 3.9, 4.1};
  Eigen::VectorXd A(args.size());
  for (int i = 0; i < A.size(); ++i) {
    A(i) = args[i];
  }
  expect_ad_vector_matvar(f, A);
}
