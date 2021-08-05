#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, erfc) {
  auto f = [](const auto& x1) { return stan::math::erfc(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -1, -0.2, 1, 1.3, 2.6);
}

TEST(mathMixMatFun, erfcmatvar) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_args;
  auto f = [](const auto& x1) { return stan::math::erfc(x1); };
  std::vector<double> com_args = common_args();
  std::vector<double> args{-2.6, -2, -1, -0.2, 1, 1.3, 2.6};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
