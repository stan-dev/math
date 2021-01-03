#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST(mathMixScalFun, step) {
  auto f = [](const auto& x) { return stan::math::step(x); };
  stan::test::expect_ad(f, -18765.3);
  stan::test::expect_ad(f, -1.0);
  stan::test::expect_ad(f, -0.001);
  stan::test::expect_ad(f, 0.001);
  stan::test::expect_ad(f, 1.0);
  stan::test::expect_ad(f, 3.5);
  stan::test::expect_ad(f, 2713.5);

  stan::test::expect_ad(f, std::numeric_limits<double>::quiet_NaN());
}

TEST(mathMixMatFun, step_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) {
    using stan::math::step;
    return step(x1);
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
