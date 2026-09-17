#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, inv) {
  auto f = [](const auto& x1) { return stan::math::inv(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -0.2, 1.3, 3);
}

TEST(mathMixMatFun, inv_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_args;
  auto f = [](const auto& x1) {
    using stan::math::inv;
    return inv(x1);
  };
  std::vector<double> com_args = common_args();
  std::vector<double> args{-2.6, -2, -0.2, 1.3, 3};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}

TEST(mathMixMatFun, inv_extreme_derivative) {
  using namespace stan::math;
  nested_rev_autodiff nested;
  var x = 1e160;
  inv(x).grad();
  EXPECT_NEAR(-1e-320, x.adj(), 1e-323);
  EXPECT_NEAR(-1e-320, inv(fvar<double>(1e160, 1)).d_, 1e-323);
  var_value<Eigen::VectorXd> xv(Eigen::VectorXd::Constant(1, 1e160));
  sum(inv(xv)).grad();
  EXPECT_NEAR(-1e-320, xv.adj()[0], 1e-323);
}
