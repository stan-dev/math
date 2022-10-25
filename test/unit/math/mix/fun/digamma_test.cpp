#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, digamma) {
  auto f = [](const auto& x1) { return stan::math::digamma(x1); };
  stan::test::expect_common_nonzero_unary_vectorized<
      stan::test::ScalarSupport::Real>(f);
  stan::test::expect_unary_vectorized(f, -25, -10.2, -1.2, -1, 2.3, 5.7);
}

TEST(mathMixMatFun, digamma_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) {
    using stan::math::digamma;
    return digamma(x1);
  };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{-25, -10.2, -1.2, -1, 2.3, 5.7};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
