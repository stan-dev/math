#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, invPhi) {
  auto f = [](const auto& x1) { return stan::math::inv_Phi(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, 0.02425, 0.97575);  // breakpoints
  stan::test::expect_unary_vectorized(f, -100.25, -2, 0.01, 0.1, 0.98, 0.5,
                                      2.0);
}

TEST(mathMixMatFun, invPhi_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_args;
  auto f = [](const auto& x1) {
    using stan::math::inv_Phi;
    return inv_Phi(x1);
  };
  std::vector<double> com_args = common_args();
  std::vector<double> args{0.02425, 0.97575, -2, 0.1, 0.5, 2.0};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
