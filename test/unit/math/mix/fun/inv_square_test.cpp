#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, invSquare) {
  auto f = [](const auto& x1) { return stan::math::inv_square(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -0.2, 1, 1.3, 3);
}

TEST(mathMixMatFun, invsquare_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_args;
  auto f = [](const auto& x1) {
    using stan::math::inv_square;
    return inv_square(x1);
  };
  std::vector<double> com_args = common_args();
  std::vector<double> args{-2.6, -2, -0.2, 1, 1.3, 3};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}

TEST(mathMixMatFun, inv_square_extreme_value_and_derivative) {
  using namespace stan::math;
  EXPECT_NEAR(1e-320, inv_square(1e160), 1e-323);
  EXPECT_NEAR(1e-320, inv_square(Eigen::VectorXd::Constant(1, 1e160))[0],
              1e-323);
  nested_rev_autodiff nested;
  var x = 1e104;
  inv_square(x).grad();
  EXPECT_NEAR(-2e-312, x.adj(), 1e-323);
  EXPECT_NEAR(-2e-312, inv_square(fvar<double>(1e104, 1)).d_, 1e-323);
  var_value<Eigen::VectorXd> xv(Eigen::VectorXd::Constant(1, 1e104));
  sum(inv_square(xv)).grad();
  EXPECT_NEAR(-2e-312, xv.adj()[0], 1e-323);
}
