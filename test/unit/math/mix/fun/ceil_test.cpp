#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, ceil) {
  auto f = [](const auto& x1) { return stan::math::ceil(x1); };
  // can't autodiff ceil through integers
  for (auto x : stan::test::internal::common_nonzero_args())
    stan::test::expect_unary_vectorized(f, x);
  stan::test::expect_unary_vectorized(f, -2.6, -2.1, -0.2, 1.1, 1.51, 3.1);
}

TEST(mathMixMatFun, ceilmatvar) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_args;
  auto f = [](const auto& x1) { return stan::math::ceil(x1); };
  std::vector<double> com_args = common_args();
  std::vector<double> args{-2.6, -2.1, -0.2, 1.1, 1.51, 3.1};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
