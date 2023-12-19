#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, logInvLogit) {
  auto f = [](const auto& x1) { return stan::math::log_inv_logit(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -1.0, -0.5, 0.5, 1.3, 5);
}

TEST(mathMixMatFun, logInvLogitVarMat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) {
    using stan::math::acos;
    return stan::math::log_inv_logit(x1);
  };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{
      -2.2, -0.8, 0.5, 1 + std::numeric_limits<double>::epsilon(),
      1.5,  3,    3.4, 4};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
