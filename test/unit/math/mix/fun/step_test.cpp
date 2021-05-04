#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST(mathMixScalFun, step) {
  auto f = [](const auto& x) { return stan::math::step(x); };

  stan::test::expect_common_prim([](auto x) { return x < 0.0 ? 0 : 1; }, f);

  stan::test::expect_ad(f, std::numeric_limits<double>::quiet_NaN());
}

TEST(mathMixMatFun, step_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) {
    using stan::math::asin;
    return asin(x1);
  };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{-2.6, -2, -0.2, 0.7, 1.3, 3.4, 5};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
