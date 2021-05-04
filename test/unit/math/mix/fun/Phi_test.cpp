#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, Phi) {
  auto f = [](const auto& x1) { return stan::math::Phi(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -27.5, 27.5);
  for (double x = -37.5; x <= 10; x += 0.5)
    stan::test::expect_unary_vectorized(x);
}

TEST(mathMixMatFun, Phi_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) {
    using stan::math::Phi;
    return Phi(x1);
  };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{-27.5, -0.5, 0.0, 1.1, 27.5};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
