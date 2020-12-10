#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, cbrt) {
  auto f = [](const auto& x1) { return stan::math::cbrt(x1); };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, 1, 1.3, 3);
}

TEST(mathMixMatFun, cbrt_varmat) {
  auto f = [](const auto& x1) {
    using stan::math::cbrt;
    return cbrt(x1);
  };
  auto com_args = stan::test::internal::common_nonzero_args();
  std::vector<double> extra_args{-2.6, -2, 1, 1.3, 3};
  Eigen::VectorXd A(com_args.size() + extra_args.size());
  int i = 0;
  for (double x : com_args) {
    A(i) = x;
    ++i;
  }
  for (double x : extra_args) {
    A(i) = x;
    ++i;
  }
  stan::test::expect_ad_matvar(f, A);
}
