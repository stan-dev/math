#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, PhiApprox) {
  auto f = [](const auto& x1) { return stan::math::Phi_approx(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -3.0, 1, 1.3, 3);
}

TEST(mathMixMatFun, PhiApprox_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) {
    using stan::math::Phi_approx;
    return Phi_approx(x1);
  };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{-3.0, 1, 1.3, 3};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
