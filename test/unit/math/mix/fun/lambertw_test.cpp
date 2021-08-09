#include <test/unit/math/test_ad.hpp>
#include <stan/math/mix.hpp>
#include <gtest/gtest.h>

TEST(mathMixMatFun, lambert_w0) {
  auto f = [](const auto& x1) {
    using stan::math::lambert_w0;
    return lambert_w0(x1);
  };
  stan::test::expect_unary_vectorized(f, -0.3, -0.1, 0.0, 1, 10, 20);

  // Test bounds
  stan::test::expect_all_throw(f, -0.38);

  stan::math::recover_memory();
}

TEST(mathMixMatFun, lambert_w0_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_args;
  auto f = [](const auto& x1) {
    using stan::math::lambert_w0;
    return lambert_w0(x1);
  };
  std::vector<double> com_args = common_args();
  std::vector<double> args{-0.3, -0.1, 0.0, 1, 10, 20, -0.38};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}

TEST(mathMixMatFun, lambert_wm1) {
  auto f = [](const auto& x1) {
    using stan::math::lambert_wm1;
    return lambert_wm1(x1);
  };
  stan::test::expect_unary_vectorized(f, -0.35, -0.3, -0.1, -0.01);

  // Test bounds
  stan::test::expect_all_throw(f, -0.38);
  stan::test::expect_all_throw(f, 0.001);

  stan::math::recover_memory();
}

TEST(mathMixMatFun, lambert_wm1_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_args;
  auto f = [](const auto& x1) {
    using stan::math::lambert_wm1;
    return lambert_wm1(x1);
  };
  std::vector<double> com_args = common_args();
  std::vector<double> args{-0.35, -0.3, -0.1, -0.01, -0.38, 0.001};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
